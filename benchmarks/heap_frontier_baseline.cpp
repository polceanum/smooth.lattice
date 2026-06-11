#include <boost/multiprecision/cpp_int.hpp>

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
#include <vector>

using boost::multiprecision::cpp_int;
using namespace std;

struct Trace {
    uint32_t parent;
    uint16_t dim;
};

struct HeapNode {
    double log_value;
    uint32_t trace_id;
    uint16_t last_dim;
};

struct HeapGreater {
    bool operator()(const HeapNode& a, const HeapNode& b) const {
        if (a.log_value != b.log_value) return a.log_value > b.log_value;
        return a.trace_id > b.trace_id;
    }
};

vector<int> parse_primes(const string& raw) {
    vector<int> primes;
    stringstream ss(raw);
    string item;
    while (getline(ss, item, ',')) {
        if (!item.empty()) primes.push_back(stoi(item));
    }
    if (primes.empty()) throw invalid_argument("empty prime list");
    return primes;
}

vector<uint64_t> reconstruct_exps(const vector<Trace>& traces, uint32_t trace_id, size_t k) {
    vector<uint64_t> exps(k, 0);
    while (trace_id != 0) {
        const Trace& trace = traces[trace_id];
        exps[trace.dim]++;
        trace_id = trace.parent;
    }
    return exps;
}

cpp_int value_from_exps(const vector<int>& primes, const vector<uint64_t>& exps) {
    cpp_int out = 1;
    for (size_t i = 0; i < primes.size(); ++i) {
        for (uint64_t j = 0; j < exps[i]; ++j) out *= primes[i];
    }
    return out;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "usage: heap_frontier_baseline primes_csv N\n";
        return 2;
    }
    const string primes_raw = argv[1];
    const vector<int> primes = parse_primes(primes_raw);
    const uint64_t n = stoull(argv[2]);
    if (n == 0) {
        cerr << "N is 1-based\n";
        return 2;
    }
    if (primes.size() > numeric_limits<uint16_t>::max()) {
        cerr << "too many primes for trace representation\n";
        return 2;
    }

    vector<double> logs;
    logs.reserve(primes.size());
    for (int p : primes) logs.push_back(log(static_cast<double>(p)));

    vector<Trace> traces;
    traces.reserve(static_cast<size_t>(min<uint64_t>(n * primes.size() + 1, 10000000ULL)));
    traces.push_back({0, numeric_limits<uint16_t>::max()});

    priority_queue<HeapNode, vector<HeapNode>, HeapGreater> heap;
    heap.push({0.0, 0, 0});

    uint64_t pops = 0;
    uint64_t pushes = 1;
    size_t max_heap = heap.size();
    HeapNode current{0.0, 0, 0};

    auto t0 = chrono::high_resolution_clock::now();
    while (pops < n) {
        if (heap.empty()) {
            cerr << "frontier unexpectedly empty\n";
            return 1;
        }
        current = heap.top();
        heap.pop();
        ++pops;
        if (pops == n) break;

        for (uint16_t dim = current.last_dim; dim < primes.size(); ++dim) {
            if (traces.size() >= numeric_limits<uint32_t>::max()) {
                cerr << "trace id overflow\n";
                return 1;
            }
            traces.push_back({current.trace_id, dim});
            const uint32_t child_id = static_cast<uint32_t>(traces.size() - 1);
            heap.push({current.log_value + logs[dim], child_id, dim});
            ++pushes;
        }
        if (heap.size() > max_heap) max_heap = heap.size();
    }
    auto t1 = chrono::high_resolution_clock::now();

    const vector<uint64_t> exps = reconstruct_exps(traces, current.trace_id, primes.size());
    const cpp_int exact_value = value_from_exps(primes, exps);
    const size_t digits = exact_value.convert_to<string>().size();
    const double seconds = chrono::duration<double>(t1 - t0).count();

    cout.setf(ios::fixed);
    cout << setprecision(6);
    cout << "heap_frontier"
         << " P=" << primes_raw
         << " k=" << primes.size()
         << " N=" << n
         << " seconds=" << seconds
         << " last_log=" << current.log_value
         << " heap_pops=" << pops
         << " heap_pushes=" << pushes
         << " max_heap=" << max_heap
         << " trace_nodes=" << traces.size()
         << " digits=" << digits
         << " exps=[";
    for (size_t i = 0; i < exps.size(); ++i) {
        if (i) cout << ",";
        cout << exps[i];
    }
    cout << "]\n";
    return 0;
}
