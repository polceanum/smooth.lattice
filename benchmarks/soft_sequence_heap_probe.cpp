#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

// Auditable soft-sequence-heap probe after Brodal's sequence-heap
// presentation of Chazelle/Kaplan-style soft heaps. This is intended as a
// correctness/corruption-semantics bridge, not yet as the optimized selector
// primitive used by the published soft-heap X+Y algorithms.

class SoftSequenceHeap {
public:
    struct Extracted {
        std::uint64_t value = 0;
        double real_key = 0.0;
        double current_key = 0.0;
        bool corrupt = false;
        std::vector<std::uint64_t> newly_corrupt;
    };

    explicit SoftSequenceHeap(double epsilon)
        : epsilon_(epsilon), r0_(rank_threshold(epsilon)) {
        if (!(epsilon > 0.0 && epsilon <= 0.5)) {
            throw std::invalid_argument("epsilon must be in (0, 0.5]");
        }
    }

    void insert(double key, std::uint64_t value) {
        int idx = (int)items_.size();
        items_.push_back(Item{key, value});
        sequences_.insert(sequences_.begin(), Sequence{0, {idx}, 0});
        inserted_++;

        while (sequences_.size() >= 2 && sequences_[0].rank == sequences_[1].rank) {
            Sequence merged = merge_sequences(sequences_[0], sequences_[1]);
            sequences_.erase(sequences_.begin(), sequences_.begin() + 2);
            sequences_.insert(sequences_.begin(), std::move(merged));
        }
    }

    bool empty() const {
        return live_count() == 0;
    }

    std::size_t inserted_count() const {
        return inserted_;
    }

    std::size_t live_count() const {
        std::size_t n = 0;
        for (const auto& item : items_) {
            if (item.live) n++;
        }
        return n;
    }

    std::size_t corrupt_count() const {
        std::size_t n = 0;
        for (const auto& item : items_) {
            if (item.live && item.c_owner >= 0 && item.w_owner < 0) n++;
        }
        return n;
    }

    std::size_t corrupt_bound() const {
        return (std::size_t)std::floor(epsilon_ * (double)inserted_ + 1e-12);
    }

    Extracted extract_min() {
        while (true) {
            int si = min_sequence_index();
            if (si < 0) throw std::runtime_error("extract_min from empty heap");
            Sequence& seq = sequences_[(std::size_t)si];
            int head = seq.items[seq.head];
            Item& h = items_[(std::size_t)head];

            if (!h.cset.empty()) {
                int x = h.cset.back();
                h.cset.pop_back();
                Item& ix = items_[(std::size_t)x];
                if (!ix.live || ix.c_owner != head || ix.w_owner >= 0) {
                    throw std::runtime_error("invalid corrupted-set head during extract");
                }
                ix.live = false;
                ix.c_owner = -1;
                return Extracted{ix.value, ix.key, h.key, true, {}};
            }

            std::vector<std::uint64_t> newly;
            newly.reserve(h.wset.size());
            for (int x : h.wset) {
                Item& ix = items_[(std::size_t)x];
                if (!ix.live || ix.w_owner != head || ix.c_owner < 0) {
                    throw std::runtime_error("invalid witness-set member during extract");
                }
                ix.w_owner = -1;
                newly.push_back(ix.value);
            }
            h.wset.clear();

            h.live = false;
            seq.head++;
            if (seq.head >= seq.items.size()) {
                sequences_.erase(sequences_.begin() + si);
            }
            return Extracted{h.value, h.key, h.key, false, std::move(newly)};
        }
    }

    bool validate(std::string* why = nullptr) const {
        auto fail = [&](const std::string& message) {
            if (why) *why = message;
            return false;
        };

        for (std::size_t i = 1; i < sequences_.size(); i++) {
            if (sequences_[i - 1].rank >= sequences_[i].rank) {
                return fail("sequence ranks are not strictly increasing");
            }
        }

        std::vector<int> in_sequence(items_.size(), 0);
        std::vector<int> in_cset(items_.size(), 0);
        std::vector<int> in_wset(items_.size(), 0);

        for (const auto& seq : sequences_) {
            if (seq.head >= seq.items.size()) return fail("empty sequence retained");
            double last = -std::numeric_limits<double>::infinity();
            for (std::size_t p = seq.head; p < seq.items.size(); p++) {
                int x = seq.items[p];
                if (x < 0 || (std::size_t)x >= items_.size()) return fail("bad sequence item index");
                const Item& ix = items_[(std::size_t)x];
                if (!ix.live) return fail("dead item retained in sequence");
                if (ix.c_owner >= 0 || ix.w_owner >= 0) return fail("sequence item has C/W owner");
                if (ix.key + 1e-15 < last) return fail("sequence is not sorted");
                last = ix.key;
                in_sequence[(std::size_t)x]++;
            }
        }

        for (std::size_t owner = 0; owner < items_.size(); owner++) {
            const Item& io = items_[owner];
            for (int x : io.cset) {
                if (x < 0 || (std::size_t)x >= items_.size()) return fail("bad C index");
                const Item& ix = items_[(std::size_t)x];
                if (!io.live || !ix.live) return fail("dead item in C relation");
                if (ix.c_owner != (int)owner) return fail("C owner mismatch");
                if (ix.key > io.key + 1e-15) return fail("C order violation");
                in_cset[(std::size_t)x]++;
            }
            for (int x : io.wset) {
                if (x < 0 || (std::size_t)x >= items_.size()) return fail("bad W index");
                const Item& ix = items_[(std::size_t)x];
                if (!io.live || !ix.live) return fail("dead item in W relation");
                if (ix.w_owner != (int)owner) return fail("W owner mismatch");
                if (ix.c_owner < 0) return fail("witnessed item has no C owner");
                if (ix.key + 1e-15 < io.key) return fail("W order violation");
                in_wset[(std::size_t)x]++;
            }
        }

        for (std::size_t i = 0; i < items_.size(); i++) {
            const Item& item = items_[i];
            if (!item.live) {
                if (in_sequence[i] || in_cset[i] || in_wset[i]) return fail("dead item is still referenced");
                continue;
            }
            if (in_sequence[i] + in_cset[i] != 1) return fail("live item is not in exactly one primary place");
            if (in_wset[i] > 1) return fail("item has multiple witnesses");
            if (item.w_owner >= 0 && in_wset[i] != 1) return fail("w_owner without W membership");
            if (item.w_owner < 0 && in_wset[i] != 0) return fail("W membership without w_owner");
        }

        if (corrupt_count() > corrupt_bound()) {
            std::ostringstream os;
            os << "corrupt count " << corrupt_count() << " exceeds bound " << corrupt_bound();
            return fail(os.str());
        }

        return true;
    }

private:
    struct Item {
        double key = 0.0;
        std::uint64_t value = 0;
        bool live = true;
        int c_owner = -1;
        int w_owner = -1;
        std::vector<int> cset;
        std::vector<int> wset;
    };

    struct Sequence {
        int rank = 0;
        std::vector<int> items;
        std::size_t head = 0;
    };

    static int rank_threshold(double epsilon) {
        return (int)std::ceil(std::log2(1.0 / epsilon));
    }

    int min_sequence_index() const {
        int best = -1;
        double best_key = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < sequences_.size(); i++) {
            const Sequence& seq = sequences_[i];
            if (seq.head >= seq.items.size()) continue;
            double key = items_[(std::size_t)seq.items[seq.head]].key;
            if (key < best_key) {
                best_key = key;
                best = (int)i;
            }
        }
        return best;
    }

    void add_to_c(int owner, int x) {
        items_[(std::size_t)owner].cset.push_back(x);
        items_[(std::size_t)x].c_owner = owner;
    }

    void add_to_w(int owner, int x) {
        items_[(std::size_t)owner].wset.push_back(x);
        items_[(std::size_t)x].w_owner = owner;
    }

    void prune_item(int pred, int x, int succ) {
        add_to_c(succ, x);
        std::vector<int> old_c;
        old_c.swap(items_[(std::size_t)x].cset);
        for (int y : old_c) add_to_c(succ, y);

        add_to_w(pred, x);
        std::vector<int> old_w;
        old_w.swap(items_[(std::size_t)x].wset);
        for (int y : old_w) add_to_w(pred, y);
    }

    std::vector<int> reduce_sequence(const std::vector<int>& input) {
        std::vector<int> out;
        out.reserve((input.size() + 1) / 2 + 1);
        for (std::size_t i = 0; i < input.size(); i++) {
            bool prune = (i % 2 == 1) && (i + 1 < input.size());
            if (prune) {
                prune_item(input[i - 1], input[i], input[i + 1]);
            } else {
                out.push_back(input[i]);
            }
        }
        return out;
    }

    Sequence merge_sequences(const Sequence& a, const Sequence& b) {
        std::vector<int> merged;
        merged.reserve((a.items.size() - a.head) + (b.items.size() - b.head));
        std::size_t i = a.head, j = b.head;
        while (i < a.items.size() || j < b.items.size()) {
            if (j >= b.items.size() ||
                (i < a.items.size() && items_[(std::size_t)a.items[i]].key <= items_[(std::size_t)b.items[j]].key)) {
                merged.push_back(a.items[i++]);
            } else {
                merged.push_back(b.items[j++]);
            }
        }

        int rank = a.rank + 1;
        if (rank > r0_ && ((rank - r0_) % 2 == 0)) {
            merged = reduce_sequence(merged);
        }
        return Sequence{rank, std::move(merged), 0};
    }

    double epsilon_;
    int r0_;
    std::size_t inserted_ = 0;
    std::vector<Item> items_;
    std::vector<Sequence> sequences_;
};

static double permuted_key(std::uint64_t i) {
    std::uint64_t x = i * 2862933555777941757ULL + 3037000493ULL;
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    return (double)(x % 1000000007ULL) + 1e-9 * (double)i;
}

static bool run_validation(std::size_t n, double epsilon, const std::string& pattern, bool verbose) {
    SoftSequenceHeap heap(epsilon);
    for (std::size_t i = 0; i < n; i++) {
        double key = 0.0;
        if (pattern == "increasing") key = (double)i;
        else if (pattern == "decreasing") key = (double)(n - i);
        else key = permuted_key(i);
        heap.insert(key, i);
        std::string why;
        if (!heap.validate(&why)) {
            std::cerr << "validation failed after insert i=" << i << " pattern=" << pattern << ": " << why << "\n";
            return false;
        }
    }

    double last_current = -std::numeric_limits<double>::infinity();
    std::size_t extracted = 0;
    std::size_t corrupt_returns = 0;
    std::size_t reported_corruptions = 0;
    while (!heap.empty()) {
        auto e = heap.extract_min();
        if (e.current_key + 1e-15 < last_current) {
            std::cerr << "current-key order violation pattern=" << pattern << "\n";
            return false;
        }
        if (e.real_key > e.current_key + 1e-15) {
            std::cerr << "corrupt key below real key pattern=" << pattern << "\n";
            return false;
        }
        last_current = e.current_key;
        corrupt_returns += e.corrupt ? 1 : 0;
        reported_corruptions += e.newly_corrupt.size();
        extracted++;
        std::string why;
        if (!heap.validate(&why)) {
            std::cerr << "validation failed after extract count=" << extracted << " pattern=" << pattern << ": " << why << "\n";
            return false;
        }
    }

    if (extracted != n) {
        std::cerr << "extracted " << extracted << " items from " << n << " inserts\n";
        return false;
    }
    if (verbose) {
        std::cout << "soft_sequence_heap_validate pattern=" << pattern
                  << " n=" << n
                  << " epsilon=" << epsilon
                  << " extracted=" << extracted
                  << " corrupt_returns=" << corrupt_returns
                  << " reported_corruptions=" << reported_corruptions
                  << "\n";
    }
    return true;
}

static void run_bench(std::size_t n, double epsilon) {
    std::vector<double> keys;
    keys.reserve(n);
    for (std::size_t i = 0; i < n; i++) keys.push_back(permuted_key(i));

    auto t0 = std::chrono::high_resolution_clock::now();
    SoftSequenceHeap sh(epsilon);
    for (std::size_t i = 0; i < n; i++) sh.insert(keys[i], i);
    std::size_t corrupt_returns = 0;
    while (!sh.empty()) {
        corrupt_returns += sh.extract_min().corrupt ? 1 : 0;
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    struct Node {
        double key;
        std::size_t id;
        bool operator>(const Node& other) const {
            return key > other.key || (key == other.key && id > other.id);
        }
    };
    auto h0 = std::chrono::high_resolution_clock::now();
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
    for (std::size_t i = 0; i < n; i++) pq.push(Node{keys[i], i});
    while (!pq.empty()) pq.pop();
    auto h1 = std::chrono::high_resolution_clock::now();

    std::cout << std::fixed << std::setprecision(6)
              << "soft_sequence_heap_bench n=" << n
              << " epsilon=" << epsilon
              << " soft_sec=" << std::chrono::duration<double>(t1 - t0).count()
              << " binary_heap_sec=" << std::chrono::duration<double>(h1 - h0).count()
              << " corrupt_returns=" << corrupt_returns
              << "\n";
}

int main(int argc, char** argv) {
    try {
        std::string mode = argc >= 2 ? argv[1] : "validate";
        std::size_t n = argc >= 3 ? (std::size_t)std::stoull(argv[2]) : 2048;
        double epsilon = argc >= 4 ? std::stod(argv[3]) : 0.25;
        if (mode == "validate") {
            bool ok = true;
            for (const std::string pattern : {"increasing", "decreasing", "permuted"}) {
                ok = run_validation(n, epsilon, pattern, true) && ok;
            }
            return ok ? 0 : 1;
        }
        if (mode == "bench") {
            run_bench(n, epsilon);
            return 0;
        }
        throw std::invalid_argument("usage: soft_sequence_heap_probe [validate|bench] [n] [epsilon]");
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
