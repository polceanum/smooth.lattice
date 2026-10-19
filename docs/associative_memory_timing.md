# Timing-adaptive retrieval: a candidate conversion theorem

This is a working proof record, not a contribution manuscript. The active
memory/access/retention goal is unchanged. Novelty, superiority over the
strongest applicable method, and publication readiness are **not established**.
No new neural training or correspondence is authorized by this record.

The preceding [research record](interleaving_research.md) identifies existing
lattice interpolation and Gram samplers as stronger baselines. This pass asks
a different, explicit question: can a canonical, always-correct index retain
its resource guarantees when queries and key edits depend on earlier *times*,
for an arbitrarily long interaction? A query supplied as a raw vector is not
an already supplied table address. Timings, answers, and memory-access traces
are different observation models.

## Model and the proposed claim

Use a sequential finite-word RAM. Private randomness and internal memory are
not visible; the client observes the completed operation's answer and number
of RAM steps. Padding steps are charged. This is a mathematical timing model,
not a claim about constant-time software on an actual CPU or protection against
cache, power, concurrent-access, or memory-inspection attacks.

The logical dataset has a deterministic reference representation of U words,
with initial construction/loading cost I.
It supports canonical updates/answers, a worst-case D-step input commit/queue/output
operation, and a deterministic F-step fallback query. All input and payload
widths, address bits and reference-state costs belong in these quantities.
The fixed word width must hold addresses and the phase/budget counters;
arithmetic and fresh random-bit generation are charged at that width. Wider
inputs require multiword costs. No counter or stored precision grows with Q.
Canonical means the answer is the same function of the current logical state
and input for every private seed. Merely returning an arbitrary valid neighbor
does not satisfy this assumption. Reads must leave the index and its update
randomness unchanged; any fresh query coins use a separate stream.

Take R,T,M_L as positive integer upper bounds. The base randomized index is
Las Vegas: whenever its operation completes it
is correct. It has expected build/reset cost at most R. For every fixed legal
operation history, independently of its private seed, each complete operation
has expected cost at most T, uniformly over allowed starting datasets. These
are per-operation expectations, not just a total amortized bound. Copies have
self-contained states: a delayed operation must not dereference shared key
storage that a later logical edit can overwrite. Queue records carry all
required old/new data; any persistent snapshots belong in the space and
operation costs. For a phase
of at most L fixed operations, let M_L bound the expected maximum live index
and workspace over the *entire* build and phase, including interrupted updates.
An expected size bound at one instant is not a bound on M_L.

The candidate dynamic conversion has, for integer L,k>=1:

Every returned answer is canonical and correct. Rebuild attempts terminate
almost surely because their acceptance probability is uniformly positive.
Thus every operation in an unbounded countable interaction terminates and is
correct almost surely. Exact reference updates introduce no overwrite drift;
this does not compress that reference state or establish a semantic encoder.

```text
space O(U + k*M_L + L*w + k) words,
initial expected setup <= I + 8*k*R,
E[total cost of Q operations]
  <= I + 8*k*R + Q*[D + k*(G+B) + 2*(8*k*R+F)*(p+1/L)],

G = ceil(8*T*ceil(log2(2*L))), B = ceil(16*T), p=(3/4)^k.
```

Here w is the maximum word length of an update record. The queue is shared;
each copy has a cursor and a resumable operation checkpoint. Constants in
simulation, cap checks, padding and reference maintenance must be included in
R,T,D and the displayed budgets. At most L records are retained. This is an
expected *total* bound for every deterministic finite Q with one fixed set of
algorithm parameters, not a worst-case bound for each individual operation,
and not a bound for arbitrary adversarial stopping times with E[Q] substituted.
The algorithm does not need Q in advance.

For L=max(1,ceil((R+F)/T)) and
k=4*ceil(log2(2*L)), so p<=1/(2*L), the bound is

```text
I + O(k*R) initial setup + O(Q*[D + T*log^2(2+(R+F)/T)]),
space O(U + M_L*log(2*L) + L*w).
```

This is a generic conversion under the hypotheses, not a new storage
compression result. A competitive application still has to supply useful
R,T,M_L,U,D,F,w and a canonical base index with all relevant bit costs.

## Safe rebuilding: its duration must not reveal the new seed

For each copy, try independent builds in slots of exactly ceil(4*R) steps.
Reject an attempt if it fails to finish or would exceed 16*M_L live words;
pad successful and failed slots alike. Under the original seed distribution,
Markov's inequality and a union bound give acceptance probability at least
1-1/4-1/16>1/2. Thus the expected cost of an accepted build is at most 8*R.
Cleanup/abandonment costs are part of R; old private state is inaccessible in
this observation model. This is not a cryptographic secure-erasure claim.

If A is acceptance, the installed seed has law mu(.|A). For every slot count j,

```text
Pr(number of slots=j, installed seed in E)
  = Pr(not A)^(j-1) * mu(E intersect A).
```

Consequently the installed seed is independent of the observed number of
slots, conditional on the public logical dataset. Different copies remain
independent even after their build times are observed. An unpadded successful
build could reveal the new seed and invalidate the argument immediately.

Acceptance has probability at least 1/2, so every nonnegative cost/space
expectation conditional on acceptance is at most twice its original bound.
The same reasoning applies to a future fixed history chosen after the build:
it must not depend on the installed seed. No empirical timing estimate replaces
these hypotheses.

## An epoch-count lemma without a false conditional-independence step

Suppose each fresh epoch, conditional on the history before it starts, has a
first exceptional operation tau satisfying Pr(tau<=s)<=s*a. Queries need not
remain independent of the seed after observing ordinary replies: survival
itself changes its posterior. Let N_Q count exceptions followed by fresh
epochs in Q operations. Then

```text
E[N_Q] <= 2*a*Q.
```

Here is a proof that also handles a truncated last epoch. Fix a remaining
horizon h and write z=Pr(tau<=h), so z<=min(1,a*h). For 0<=t<=h,
Pr(tau>t)>=1-min(a*t,z). Therefore

```text
E[min(tau,h)] >= integral_0^(z/a) (1-a*t) dt
                 = (z-z^2/2)/a >= z/(2*a).
```

The zero-a case is immediate. Induct on h, conditioning on the first exception
and then using the same uniform claim for the new epoch and its new history:

```text
E[N_h] <= z + 2*a*E[(h-tau)*1_{tau<=h}]
       = 2*a*h + z - 2*a*E[min(tau,h)] <= 2*a*h.
```

This argument needs only a first-exception CDF bound at each fresh epoch. It
does not assume a per-query conditional hazard <=a. That stronger statement is
false: with one hidden bad query among m, m-1 ordinary replies can identify
the last query as bad with certainty. The factor two is asymptotically sharp
under the CDF hypothesis: fresh uniform epoch lengths in {1,...,m} have a=1/m
and long-run exception rate 2/(m+1). This sharpness is for the abstract lemma,
not an architectural or ANN lower bound.

## Static read-only special case

For a static read-only index with expected build peak M, use build caps 4*R
steps and 4*M words. Acceptance is again at least 1/2. Run k independent copies
for B>=2*T steps per query; return the canonical answer if any finishes, and
pad the ordinary response to k*B steps. If none finishes, use the reference
fallback and discard/rebuild every copy. Charge any output read separately
if it is not already in T and F. Query workspace touched in B steps is at
most B words per copy; peak storage is O(U+k*(M+B)), not automatically O(U+k*M).

Within a fresh epoch, fix the client's coins and the hypothetical transcript
on which every operation is ordinary. Canonical replies and padded timings
make this query path independent of the installed seeds. At each position
on this path, Markov's inequality and independence of the copies give an
unconditional all-copy timeout probability at most p=(2*T/B)^k. Union bound
gives Pr(first timeout<=s)<=s*p. The epoch lemma yields

```text
E[total work] <= I + 8*k*R + Q*[k*B + 2*(8*k*R+F)*(2*T/B)^k].
```

For k=1, balancing terms gives O(sqrt(T*(R+F))) expected amortized query work,
up to the displayed constants and initial setup. A simple explicit choice
for a logarithmic number of copies is

```text
A=max(2,ceil((R+F)/T)), k=4*ceil(log2(2*A)), B=ceil(4*T).
```

It gives O(T*log(2+A)) expected amortized work and O(log(2+A)) copies:
2^k >= (2*A)^4 >= 8*k*A, so (8*k*R+F)*2^(-k)<=T.
For integer A>=2, 2*A<=2^A implies k<=4*A; hence
(2*A)^4>=32*A^2>=8*k*A. The dynamic choice similarly uses
(3/4)^4<1/2 to ensure p<=1/(2*L), using only integer parameter choices.
All guarantees are independent of Q. This does not make a rare fallback
operation cheap in worst-case time, or prove a speedup over another index.

## Dynamic proof: queues, cap checks, and replacement

At the start of an epoch, build k independent copies of the current reference
dataset by the safe procedure above. On every operation, first commit its
logical change in the deterministic reference representation and append its
self-contained update record, if any. Each live copy receives G charged simulation steps to
process its pending updates in order. A copy which would exceed 16*M_L words
is disabled until the next epoch. Every input record and checkpoint counts
toward storage. No partially updated copy answers a query.

For a query, give each copy that has caught up B query steps. Return the
canonical answer from any successful copy. If no copy finishes correctly
within its budget, answer using the deterministic reference fallback and
replace all copies. Pad ordinary operation times to D+k*(G+B), including
unused work slots. Start a new epoch after L operations even if there has
been no fallback. This bounds the queue and avoids a horizon-dependent G.
Rebuild after an exceptional response before processing the next input;
the client may see all rebuild timings. Their independence was proved above.

For analysis, fix a fresh epoch's hypothetical all-ordinary path. Its updates
and queries are independent of the new seeds. Consider a single copy's full
uninterrupted execution of these updates. Each update cost has expectation
at most 2*T under the accepted-seed law; no independence between costs is
needed. The backlog is the counter

```text
C_t = max(0,C_(t-1)+X_t-G), C_0=0.
```

Lemma 2.3 of Bernstein--Forster--Henzinger (primary source below) gives
Pr(C_t>0)<=1/2 for histories of length at most L with the chosen G. Query
steps do not consume or modify update randomness. Skipped reads consequently
do not alter this counterfactual update history.

For this fixed history, conditional-on-acceptance expected maximum live
space is <=2*M_L. Thus the probability of ever disabling this copy during
the epoch is <=1/8. At any fixed query position its full, caught-up query
cost exceeds B with probability <=1/8. A union bound shows that a copy is
alive, caught up and fast with probability at least
1-1/2-1/8-1/8=1/4. These events are independent *between copies* on the fixed
path, so all-copy query failure probability is <=p=(3/4)^k.

An early replacement within s<L operations therefore has probability <=s*p.
For s>=L the scheduled replacement gives probability one, which is <=s/L.
The first replacement CDF is bounded by s*(p+1/L). Apply the epoch lemma;
each replacement costs at most F+8*k*R in conditional expectation (charging
F even for scheduled replacements is conservative). This proves the stated
expected total-work bound. The argument never conditions the fresh-copy
expectation directly on survival of the current epoch.

The phase logarithm is not removable from this *generic first-moment counter*
argument alone. For L=2^J and per-step mean bound one, take a work budget
g=J/3 (take J divisible by three for integer budgets). Partition times into dyadic blocks
{L-2^(j+1)+1,...,L-2^j}, of sizes 2^j. On a disjoint outcome (j,t), put
a single work spike 3*g*2^j at t, with probability 1/(3*g*2^j).
These outcomes have total probability J/(3*g)=1; each time has expected
work at most one, yet every outcome leaves positive backlog at L.
Thus a budget of order log L can be necessary even for an oblivious history.
This is an elementary counter witness, not an index/architecture lower bound
or a newly established contribution. Additional structure of an actual index
could improve the budget.

## Primary-source boundary and strongest comparisons

* [Andoni--Haris--Kelman--Onak, arXiv:2601.00272v1](https://arxiv.org/html/2601.00272v1):
  Definitions 3.1/4.1, Theorem 4.2, Claim 4.3 and Sections 4/7 inspected.
  The formal transcript lists query/answer pairs, with oblivious key updates.
  Its timing-release construction protects annulus-selection statistics;
  do not silently identify it with exposure of every physical operation time.
  The fair-ANN bound has a density-ratio term and logarithmic horizon factors;
  its other robust constructions and for-all guarantee are also relevant.
  Under separated keys the density ratio is small, so comparing only against
  its sqrt(Q)-space construction would be misleading. No defect in its stated
  guarantees, or superiority over its whole tradeoff, is established here.
* [Ahle, full version dated June 27, 2018](https://thomasahle.com/papers/lasvegas.pdf),
  Definition 1, Theorem 1, Lemma 1 and Section 3 construction/cost discussion
  inspected; the statement on PDF page 3 was visually checked;
  [versioned preprint](https://arxiv.org/abs/1704.02054v3).
  It gives covering filters and Las Vegas Hamming near-neighbor queries.
  Its query and construction costs are expected. The covering property,
  rather than ordinary Monte Carlo LSH, is needed to reject false empty
  answers. A complete dynamic M_L/word-cost instantiation remains to audit;
  the preprint abstract alone is not that audit.
* [Wei, Optimal Las Vegas Approximate Near Neighbors](https://arxiv.org/abs/1807.07527)
  is another, stronger applicable Las Vegas baseline. The primary abstract
  and Corollary 3.11/construction discussion were inspected. The entire
  data-dependent proof and a dynamic finite-precision instantiation have not
  been audited here. No numerical exponent comparison is claimed.
* [Bernstein--Forster--Henzinger, A Deamortization Approach](https://eprints.cs.univie.ac.at/7159/1/A_Deamortization_Approach.pdf),
  Theorem 2.1, Lemma 2.3 and their queue/counter proof inspected in the
  57-page author version; PDF pages 6--7 visually checked. Independent copies,
  bounded work slices and the
  logarithmic backlog budget are already explicit. Its graph applications
  assume oblivious updates. Our dynamic derivation depends on that known
  lemma; it is not a new backlog/deamortization theorem.
* [Bethea--Reiter, Data Structures with Unpredictable Timing](https://reitermk.github.io/papers/2009/ESORICS.pdf),
  Sections 1/2/5 and the stated timing-adversary model inspected. It already
  studies timing-adaptive access and per-operation rerandomization of sets.
  Its skip-list mechanism differs from the candidate conversion.
* [Braun--Jana--Boneh, Robust and Efficient Elimination of Cache and Timing Side Channels](https://arxiv.org/abs/1506.00189),
  Section IV-A, including overtaken padding thresholds and secret replacement,
  inspected. Padding and replacing a secret after an overrun are existing
  ideas. Their use here is not by itself a new result, and this abstract RAM
  proof does not supply their hardware protections.
* [Beimel et al., Dynamic Algorithms Against an Adaptive Adversary](https://arxiv.org/abs/2111.03980):
  the later audit below reads Theorems 3.1/3.7 and the blinking-adversary
  definition beyond the abstract. Restarting fixed-length phases to obtain
  horizon-independent resource parameters is already explicit prior art.

An unsuccessful keyword search is not evidence that the conversion is new.
Further exact comparisons are needed with adaptive runtime conversions,
bounded leakage/rekeying, and dynamic robust near-neighbor indexes, including
the concurrent search-problem work cited by Andoni et al. The concrete claim
to audit is the *complete, horizon-independent timing-adaptive conversion*,
including changing canonical bindings and the phase space bound. The component
Markov, independence, renewal and queue calculations alone are insufficient
paper contributions.

## What an associative-memory application would and would not give

For raw Hamming keys with pairwise distance >2*c*r, at most one key lies within
c*r of any query. An always-correct (r,c*r)-ANN answer, verified at exact
distance <=r, therefore implements the same canonical exact-radius lookup
for every raw query: its current payload, or bottom if none qualifies.
This is an elementary reduction from a known ANN guarantee. Exact radius
comparisons, full raw keys, all independent payloads, indexes, construction,
and overwrite costs must be charged. Pairwise separation is an input promise;
checking it or learning such an encoder is not free. If the key set changes,
the promise must continue to hold.

There is no two-latent projection in this specification: the index receives
the full raw vector. Payload overwrites can live in an exact current table.
But all n*v payload coordinates remain stored. This supplies no smooth
attention function, training advantage, general semantic encoder, or automatic
compression. Also, using cr approximately d/2 with this separation would
leave room for very few keys; do not combine a dense-code retention promise
with a favorable near-neighbor exponent from an incompatible geometric regime.

A less restrictive canonical target is the smallest-ID point within r in a
bounded-density dataset. That would require complete candidate enumeration,
not an ANN routine that stops at its first valid point. The resulting
reporting cost, key-edit cost and phase maximum space must be proved before
using this as an application. No such performance theorem is claimed here.

**Decision.** Preserve the conversion proof as a candidate, rather than train
another model. The next gate is an audited dynamic instantiation and a genuine
advantage for the same observation/workload model. If the complete guarantee
is already a published conversion, or its extra canonical/separation promises
erase its advantage, classify it as a baseline and move on. The original goal
remains active and unachieved.

The bounded controls are reproducible with:

```bash
python3 -B scripts/run_timing_restart_checks.py
```

They check optimal adaptive timeout policies in finite hidden-seed games,
the sharp renewal example, setup-count independence, and the need to retain
the phase-length logarithm in a first-moment-only backlog argument. They are
not an ANN implementation, a performance benchmark, a formal proof assistant,
or evidence of novelty.

## Concrete instantiation: exact local attention with synchronous key edits

The next pass replaces the unspecified index parameters with a complete
construction. This is a **candidate guarantee**, not an established new
algorithm or paper result. In particular, Pagh's original covering-hash
paper already combines query timeouts, independent copies and a brute-force
fallback. Neither that combination nor the covering family is our invention.
The precise remaining priority question concerns the *whole* dynamic,
timing-adaptive, all-prefix resource guarantee stated here.

Let n>=2 be a fixed slot capacity, 1<=r<=t<=d, and 1<=s<=n. A present slot i
stores an arbitrary raw key x_i in {0,1}^d and an independent v-vector of
fixed-grid values V_i=z_i/2^S, |z_i[j]|<=2^S. The promise is that, at every
operation boundary, every radius-t Hamming ball contains at most s present
slots. Multiplicities count as separate slots. This promise is not verified
for free; the theorem provides no learned encoder that enforces it.

For an arbitrary raw query q, define N_r(q)={i:dist(q,x_i)<=r}. The requested
*complete* reader is

```text
a_i(q) = r+1-dist(q,x_i),
F(q) = sum_{i in N_r(q)} a_i(q)*V_i / sum_{i in N_r(q)} a_i(q),
```

or bottom if N_r(q) is empty. All nearby slots contribute. Returning just an
arbitrary approximate neighbor does not implement this function. If keys have
pairwise distance >2*r, a query within r of a stored key returns its current
value exactly. Without that extra separation, the specified answer is the
local mixture, not a promise to retrieve one arbitrarily designated neighbor.

Use a fixed word width w large enough for addresses, r-bit parameters, and
S+O(log(s*(r+1)+1)) bit accumulators. For example a sufficiently large
constant times log(n*H*(d+v+1)+2)+S+log(s*(r+1)+2) suffices for the chosen
logarithmic number of copies, where
H=2^(r+1)-1 and mu=ceil(d/w). Set ell=ceil(log2(n+1)). All storage bounds
below are words of this explicitly charged width; multiply by w for bits.
The bounds use ordinary multiword comparisons/bit operations and conservatively
charge O(d) for a Hamming-distance calculation. They do not assume a free
whole-vector comparison, free popcount, or unlimited-precision arithmetic.

### Known covering family and exact base reader

Draw independent uniform (r+1)-bit labels m_j for j in [d]. For every nonzero
u in F_2^(r+1), use mask A_u[j]=<m_j,u> over F_2 and projection
h_u(x)=x AND A_u. The H tables are indexed by u; repeated masks may occur
and are charged with multiplicity. For every seed, every difference support
of size <=r has a nonzero orthogonal u. Thus every radius-r neighbor collides
in some table, simultaneously for all pairs of raw keys/queries.

For every *fixed* pair at Hamming distance D, every nonzero u independently
across coordinates has collision probability 2^(-D). Consequently

```text
E[number of colliding table positions] = H*2^(-D).
```

These are Pagh's established covering and collision calculations. The equality
counts table positions with multiplicity, rather than distinct masks.

Each table is a deterministic balanced search tree of IDs, ordered by
(h_u(x_i),i). Store the full raw keys once in a shared canonical table;
tree comparisons obtain projected chunks from these keys and the current
mask, costing O(mu) per comparison. There is no stored copy of the d-bit
projection at every tree node. An exact equal-range query costs
O((mu+1)*ell) plus candidate visits. Enumerate all H tables, deduplicate IDs
in query scratch, verify full distances, and return the sorted list of all
IDs/distances within r. Complete enumeration certifies an empty answer too.
An ordinary Monte Carlo LSH index cannot certify absence merely by failing
to find a collision.

Generate masks in reflected Gray-code order from the r+1 packed bit planes.
Only O((r+1)*mu) seed/basis words and a current mu-word mask are needed;
there are O(H*mu) XOR/control operations, with O(H) total Gray-code control
work by its recursive traversal. All trees together contain at most n*H
nodes of O(1) metadata words each. Query scratch has at most n IDs/marks;
initialize its marks once and clear only touched entries after completion
or abortion. Clearing costs at most a constant times work already performed.
Scratch is restored before another copy is queried.

For a fixed legal history and a fixed next query, at most s keys are within
t. Their collision visits are at most H*s. Every other key has integer
distance >=t+1, so expected far visits are at most H*n*2^(-t-1). The uniform
base index-query expectation is therefore bounded by a sufficiently
large constant C times

```text
T0 = C*H*[(mu+1)*ell + (d+ell)*(s+n*2^(-t-1))].
```

Take an integer upper bound, rounding upward. This includes deduplication,
distance checks, mask generation, sorted output IDs, and scratch cleanup.
The expectation is over the private covering seed. Tree shapes may depend
on the fixed update history, but their worst-case height bound does not.

After one copy produces the exact list, form F in a common output stage in
A0=O(s*v) word operations. Return the integer numerator vector and the
common denominator 2^S*sum_i a_i. Their bit bounds were included in w.
Exact sums make the answer independent of enumeration order/private seeds;
floating accumulation with seed-dependent rounding would not establish this.
Only this final stage accesses payloads: at most s*v value words, including
arbitrary independent current values, are used.

### Key edits have bounded cost, so no update queues are needed

Updates are supplied a slot ID; queries are supplied the full raw q, without
a target ID. An insertion, deletion or key replacement has deterministic
worst-case index cost

```text
W = O(H*(mu+1)*ell)
```

per copy, by changing one entry in every relevant table. For shared raw-key
storage, ordering matters: remove the old ID from *every copy's* trees before
changing its raw key, then insert it into every copy using the new key.
There are no queries during this synchronous operation. Mutating the shared
key first would invalidate the comparisons used to find the old entries.
No paused copy observes overwritten key storage, and no hidden snapshots are
needed. Empty slots and duplicate keys use explicit occupancy/ID metadata.

Payload overwrite alone writes O(v) current-table words and never changes
the covering indexes. A key move with a simultaneous payload replacement
pays both costs. All updates preserve the stated density promise; checking
that promise is separate work. If writes instead name a unique existing raw
key, an ordinary deterministic balanced dictionary resolves its ID in
O((mu+1)*ell) time with O(n) metadata words, which must be added. In particular,
a raw-key-named payload overwrite costs O((mu+1)*ell+v), not O(v).
Duplicate keys require the explicit slot semantics used above. Ordinary
key-edit timings are padded to the deterministic bound, and
value-edit timings to O(v), so neither reveals a private index seed.

A complete index build/reset has a deterministic padded upper bound

```text
R = O(n*W + H*mu + d*(r+1)),
M = O(n*H + (r+1)*mu + n) private index/workspace words,
F0 = O(n*d) deterministic fallback for all IDs/distances.
```

The d*(r+1) term charges generation of fresh unbiased seed bits. The
unbounded-horizon theorem consumes fresh independent bits at every rebuild;
a single finite pseudorandom-generator seed is not covered. Raw-key,
payload, occupancy and common scratch storage is O(n*(mu+v)+n) words;
initial loading/initialization I is paid. There are no build retries or
expected-space substitutions in this instantiation. M is a worst-case bound
throughout arbitrary key histories of capacity n. Since H>=r+1 and s>=1,
R+F0=O(n*T0). The reference fallback computes the same complete list, and
the common output stage computes the same F on either path.

### Dynamic timing-adaptive bound

Use k independent indexes. A query runs each in turn for at most 2*T0 body
steps; allow another 2*T0 steps for restoring touched scratch. On completion
use its exact list, and pad the index portion of an ordinary response to
4*k*T0. If all copies time out, use the reference fallback and rebuild every
index before accepting the next operation. Padding includes unused copies;
no parallel hardware or uncharged waiting is assumed. All raw reference
data and payloads survive a reset.

On an epoch's hypothetical all-ordinary transcript, canonical outputs and
fixed operation-type timing budgets make the entire key/value/query history
independent of its fresh seeds. Markov and independence between copies give
per-position all-copy timeout probability <=p=2^(-k). A union bound gives
Pr(first timeout within a queries)<=a*p, despite arbitrarily many intervening
synchronous key/value edits. Independence is between families, not masks
within a family. The earlier epoch lemma is applied in *query*
counts, not by pretending surviving seeds remain independent of the transcript.

For deterministic prescribed numbers Q queries, J key edits and P payload
overwrites (in an adaptive order), the expected total work is at most

```text
I + k*R + O(k*W*J + v*P)
  + Q*[4*k*T0 + A0 + 2*p*(k*R+F0)].
```

Choose A=max(2,ceil((R+F0)/T0)) and k=4*ceil(log2(2*A)). Then k=O(log(n+1))
and p*(k*R+F0)<=T0. Thus queries cost O(T0*log(n+1)+s*v) in expectation
after initial setup, key edits O(W*log(n+1)), and independent value overwrites
O(v). The total persistent/workspace bound is

```text
O(n*(mu+v) + log(n+1)*[n*H+(r+1)*mu]) words.
```

For t=ceil(log2 n), r=floor(t/c)>=1 with fixed c>1, H=O(n^(1/c)). For
bounded s, this is sublinear index-query work in n, up to logarithms and
the explicitly charged full-input dimension, with a superlinear index.
More generally require n^(1/c)*(s+1)*log(n+1)=o(n) for the displayed query
index work to beat a full O(n*d) scan in the usual d>=log n regime. The
output term s*v remains. This is a specified locality/density workload,
not unrestricted softmax attention or a semantic learning guarantee.

### A bound for all prefixes, including a client-chosen stopping point

The synchronous instantiation also admits a stronger cost statement than
the earlier general conversion's fixed-horizon expected bound. It uses
classical exponential-supermartingale reasoning, not a new probability tool.
Let m=1/p=2^k, z=1-1/m, and beta=z*(1-z^m). A positive integer epoch query
length tau with Pr(tau<=a)<=a/m stochastically dominates a uniform integer
in {1,...,m}; an epoch that never times out has tau=infinity. Therefore

```text
E[z^tau | history before epoch] <= (1/m)*sum_{a=1}^m z^a = beta.
```

If S_j is the cumulative number of queries at timeout j, then
L_j=z^(S_j)*beta^(-j) is a nonnegative supermartingale starting at one.
Stop it at its first crossing of 1/delta or at a deterministic epoch count.
Its stopped expectation is <=1, so crossing probability is <=delta.
Letting that count increase proves the same bound for any crossing. Thus,
with probability at least 1-delta, *simultaneously for every operation prefix*
with Q completed queries, its timeout count N satisfies

```text
N*log(1/beta) <= Q*log(1/z) + log(1/delta).
```

For m>=2, (1-1/m)^m>=1/4, so beta<=3/4;
log(1/z)<=2/m and log(4/3)>=1/4. A convenient looser bound is
N<=8*p*Q+4*log(1/delta). Combining with the deterministic cost per normal
operation, fallback and rebuild proves the simultaneous all-prefix bound

```text
work <= I+k*R + O(k*W*J+v*P)
        + Q*(4*k*T0+A0+8*p*(k*R+F0))
        + 4*(k*R+F0)*log(1/delta).
```

Here Q,J,P are the realized counts in each prefix, so the client can choose
when to stop based on earlier answers/timings. The algorithm, storage, word
precision and number of copies do not depend on a horizon or delta. The
confidence statement concerns compute cost: all answers are correct for
every seed. It does not promise a small worst-case latency for an individual
fallback, or turn the earlier general expected bound into an arbitrary
stopping-time theorem without this separate argument.

### Prior-art audit and current decision

* [Pagh's SODA 2016 author paper](https://www.rasmuspagh.net/papers/coveringLSH.pdf),
  Section 3, Lemma 3.2, Theorem 3.1 and Section 3.1 were inspected; PDF pages
  5--6 visually checked. It already supplies the covering family, collision
  expectation, all-neighbor extension, timeouts, independent copies and
  brute-force fallback. No exponent or covering-family novelty is claimed.
  Its 2018 journal version and the faster fcLSH implementation must still be
  compared before declaring a priority boundary for the complete guarantee.
* [Beimel et al.](https://arxiv.org/abs/2111.03980), Theorems 3.1/3.7, their
  algorithms/proofs and the blinking-adversary definition were inspected;
  PDF page 11 visually checked. The generic reduction concerns approximate
  scalar estimates and uses a private median. Fixed-length phase restart is
  already explicit. Blinking adversaries receive outputs only after refresh;
  this differs from the ordinary-answer-and-timing transcript used here.
  These differences alone are not proof that our conversion is new.
* [Feng et al., ICML 2025 proceedings](https://proceedings.mlr.press/v267/feng25c.html),
  Definitions/Assumptions 1.1--1.3, Theorems 1.5/1.6 and Appendix C's
  Theorem C.1 were inspected; PDF page 4 visually checked. The search
  guarantee returns a valid point/subset, not necessarily the complete exact
  neighborhood needed by F. The inspected update theorem deletes returned
  points. Do not silently replace it by an arbitrary-key-replacement theorem,
  or claim its scalar-regression guarantees solve the present workload.
  Its phase restart observations also prevent novelty claims based solely on
  having parameters independent of an indefinitely long stream.

The direct base-index construction and complete cost proof now remove the
previous unspecified M_L/dynamic-index hypothesis for this workload. The
all-prefix bound removes a horizon from the cost confidence statement.
However, the construction is built from established indexing and probability
ingredients. **A substantive new advantage over the strongest applicable
published guarantee remains unestablished.** In particular, neither a comparison
only with sqrt(Q)-copy search, nor the absence of an identical keyword result,
is a defensible novelty/impact argument. The next gate is the complete
same-workload comparison with later covering/range indexes and timing-adaptive
conversions. No neural experiment or contribution manuscript follows yet.

## Spectral concentration and constant-family robustification

This next iteration strengthens our preceding construction. Its mathematical
claim is stated below with a proof; **priority, a substantive advantage over
the strongest applicable index, and publication readiness remain open**.
The improvement over our own preceding bound is not evidence of such an
advantage. The published fast hash factorization is explicitly retained.

### The base family has more than a first-moment guarantee

Let x_1,...,x_a be the occupied keys, a<=n, with slot multiplicities retained.
Define the positive semidefinite matrix

```text
K_ij = 2^(-dist(x_i,x_j)/2), lambda = ||K||_2.
```

For no occupied keys take lambda=0. Otherwise lambda>=1. Fix any raw binary
query q, let S_i be its difference support from x_i, p_i=2^(-|S_i|),
mu=sum_i p_i, and v_i=sqrt(p_i). Use the indexed covering family from above,
H=2^(r+1)-1, with independent uniform coordinate labels m_j in F_2^(r+1).
Let X_u be the number of colliding slots in table u and C=sum_{u!=0} X_u,
counting repeated table/slot visits. Then, exactly,

```text
E C = H*mu,
Var C = H*(v^T K v - mu^2),
mu <= lambda,
Var C <= H*lambda^2/4.
```

For the proof, each nonzero u gives independent fair mask bits, so E X_u=mu.
For two slots i,j, a simultaneous collision requires the mask to vanish on
S_i union S_j. Since

```text
|S_i union S_j| = (|S_i|+|S_j|+dist(x_i,x_j))/2,
```

we have E X_u^2=v^T K v. Distinct nonzero u,u' are linearly independent over
F_2. Their mask bits are uniform independent pairs in every coordinate,
and the coordinate labels are independent. Thus the entire two masks,
and hence X_u and X_u', are independent. This is pairwise independence;
the whole family is not mutually independent. For example, the masks for
u,u',u+u' always have bitwise XOR zero. Adding pairwise covariances proves
the exact variance identity, even when a realized seed repeats some masks.

To bound mu, use unit vectors phi(0)=(1,0), phi(1)=(1/sqrt(2),1/sqrt(2)),
and their d-fold tensor products Phi(x). Their inner products equal
2^(-dist(x,y)/2). If B has rows Phi(x_i), then K=B B^T and

```text
mu = ||B Phi(q)||^2 <= ||B||^2 = lambda,
v^T K v <= lambda*mu,
lambda*mu-mu^2 <= lambda^2/4.
```

This exponentially large feature space is only a proof device: it is never
computed or stored. The tensor Gram argument is standard linear algebra,
not a new tensor-factorization algorithm. For any supplied upper bound
A>=max(1,lambda), Chebyshev now gives

```text
Pr(C > 2*H*A) <= 1/(4*H).
```

The variance inequality is sharp: take eight identical one-bit keys, q=1,
r=1. Then lambda=8, mu=4, H=3 and Var C=48=H*lambda^2/4.
It would be false to omit the spectral factor or assume independence of
all H tables. This is an unconditional finite-space statement for every
fixed dataset/query; applying it to a posterior seed law after observing
times would still be invalid.

The bound allows a spectral promise instead of a hard distance partition.
A simpler sufficient promise, directly comparable to the bounded-overlap
assumption in Feng et al., is that each key has at most Delta other keys
at distance <=2*t. With t=ceil(log2 n), each row of K sums to at most

```text
Delta+1+n*2^(-t) <= Delta+2.
```

The spectral norm is at most the maximum row sum. Thus A=Delta+2 works,
and every radius-t ball has at most Delta+1 occupied slots: any two keys
in that ball are at distance <=2*t. This is a stronger promise than just
the previous global ball-density bound, and is charged as a hypothesis.
The spectral lemma itself does not require that overlap promise. A bound
s on the radius-r neighborhood size is a separate output/aggregation cost;
without an additional promise one can use s<=min(n,2^r*A).

The spectral promise is strictly broader. For n=2B keys, let
t=ceil(log2 n), d=B*t+1. Key x_(i,z) has ones in block i of t coordinates,
zeros in the other blocks, and final bit z in {0,1}. Between blocks,
distances are 2t or 2t+1; within a block they are 1. Its Gram matrix is
the tensor product of the B-by-B matrix with diagonal 1/off-diagonal
2^(-t), and [[1,1/sqrt(2)],[1/sqrt(2),1]]. Consequently

```text
lambda=(1+(B-1)*2^(-t))*(1+1/sqrt(2)) < 3.
```

Yet each key has B other keys within distance 2t, and the radius-t ball
about zero contains B keys. For 1<=r<t a query can meet only one block's
two keys, so s<=2 and reads can be genuine weighted mixtures. Thus the
new condition is not a restatement of constant ball density/overlap.
This example has d=Theta(n log n) and charges all those raw bits; its known
block structure also permits specialized indexes. It witnesses different
promises, not a separation from all indexing algorithms or a learned-model
benefit.

### Retain the published fast factorization, with finite fingerprints

[Pham--Pagh's fcLSH](https://arxiv.org/abs/1602.02620v2), Section 3,
Algorithm 2, Equation (5) and Lemma 3, already compute the family's
universal fingerprints in O(d+H*log(H+1)) arithmetic operations using
a fast Hadamard transform. Those sections and PDF pages 6--7 were
inspected, with page 7 visually checked. Their modulus is prime. We use
the same factorization below with a finite dyadic modulus; this arithmetic
adaptation is not proposed as a separate novelty claim.
We always use independent coordinate labels, including when d<=H+1;
the paper's optional column-permutation variant is not being substituted
for this distribution.

Set b=ceil(log2(4*n*H^2)). Choose independent uniform b-bit coefficients
c_j. The fingerprint in table u is

```text
g_u(x) = sum_j c_j*x_j*a_u[j] mod 2^b.
```

For any two distinct binary projections, a differing coordinate has
coefficient +1 or -1, which is invertible modulo 2^b. Conditioning on
the other coefficients proves collision probability exactly 2^(-b).
This needs no prime search, rejection sampling or unlimited precision.
Conditional on all covering masks, a union bound over at most n*H
unequal query/slot projections gives

```text
Pr(any additional fingerprint collision) <= n*H/2^b <= 1/(4*H).
```

The coefficients are shared across tables. Do not infer independence of
compressed-table counts from the earlier independence of the masks.
Instead combine this union bound with the uncompressed variance bound:
the total fingerprint candidate count is at most 2*H*A except with
probability at most 1/(2*H). Fingerprint collisions affect work only;
full Hamming verification and deduplication still return the exact
complete neighborhood for every seed.

The transform can be done with only b+1 working bits. Aggregate
T_l=sum_{j:m_j=l} c_j*x_j modulo 2^(b+1), take its Hadamard transform
HT at that modulus, and put S=sum_l T_l modulo the same modulus. Then

```text
g_u(x) = ((S-(HT)_u) mod 2^(b+1))/2,  u!=0.
```

The representative in parentheses is even. The integer identity before
reduction is twice the desired sum, so reducing before dividing retains
the required low b bits. Using only b bits during the transform is wrong:
for a one-bit vector and coefficient 1 it can turn fingerprint 1 into 0.
The checker compares this modular transform with independently computed
direct projected dot products and enumerates extra fingerprint collisions.

### Concrete dynamic costs and the number of families

Use deterministic balanced trees of (fingerprint,slot ID), with one-word
comparisons, instead of projected full-key comparisons. Store raw keys and
payloads once. Cache fingerprints in the tree entries and maintain
slot/table handles, charging their O(n*H) words. Before changing a
raw key, remove its old entries in every family; then insert its new
fingerprints. No asynchronous query runs during a key edit. Pure payload
overwrite by supplied slot ID costs O(v); insertion also pays for its
v payload words. An update named by a raw key instead of an ID additionally
pays for a deterministic exact-key dictionary, as in the preceding section.

Here w is word width, mu_w=ceil(d/w), ell=ceil(log2(n+1)), and D_h is the
worst-case cost of verifying one full Hamming distance. A conservative
ordinary word RAM can take D_h=O(d); O(mu_w) requires explicitly including
word population count as an allowed primitive. Fix w large enough for
addresses, b+1 bits, labels, and the exact payload numerators/denominators
S+O(log(s*(r+1)+1)) bits. Otherwise pay multiword arithmetic costs.
The quantities below include fixed simulation/padding constants and can
be rounded up to integers:

```text
W_f = O(d+H*log(H+1)+H*ell),             worst-case per-family key edit,
T_f = O(d+H*log(H+1)+H*ell+H*A*(D_h+1)), fast complete candidate search,
R_f = O(n*W_f+d*(b+r+1)+H),             worst-case padded family rebuild,
F_0 = O(n*D_h),                         exact fallback neighbor scan,
common aggregation = O(s*v).
```

R_f explicitly charges every fresh random bit, including d*b coefficient
bits and d*(r+1) label bits, and clearing abandoned state. All transforms
operate on bounded words. Labels and coefficients take O(d) words per
family; no full H*d mask array is stored. Tree entries, raw arrays,
marks, touched-ID lists, transform buffers and answer lists are included:

```text
space = O(n*(mu_w+v)+k*(n*H+d)) words, plus word width w bits per word;
initial work = I+k*R_f, with I charging raw key/value initialization.
```

Range iteration starts with O(ell) tree search and O(1) amortized work
per entry. A reusable n-entry mark array and touched list deduplicate slots
without clearing all n entries at each query. On abort, clear only touched
marks; reserve this cleanup in the query budget. A complete search stores
distances/IDs, then forms the exact commutative weighted sums in common
storage. The answer is the same numerator vector and denominator for all
seeds, with no seed-dependent ordering or random valid-neighbor output.
All successful reads have the same padded index time, key edits are padded
to their worst-case budget, and rebuilds are padded to R_f.

Choose k independent complete families, pad their trials, and scan/reset
all families only when every trial times out. The previous canonical
no-timeout-branch and first-timeout-CDF proof applies with

```text
p=(2*H)^(-k),  Pr(first timeout within a queries)<=a*p.
```

This covers key edits depending on earlier completed-operation times:
along the all-normal branch, the entire logical history is independent of
the fresh seed, and fixed-seed updates produce exactly the index for its
current keys. There is no conditional-hazard independence claim.

An explicit non-circular choice is

```text
L_c=max(2,ceil((R_f+F_0)/T_f)),
l = smallest positive integer with (2*H)^l >= 2*L_c,
k=4*l.
```

For H>=3, l<=ceil(log2(2*L_c))<=L_c. Hence
(2*H)^k >= (2*L_c)^4 >= k*L_c and p*(k*R_f+F_0)<=T_f.
The preceding all-prefix proof uses only integer m=1/p, so it applies to
m=(2*H)^k, even though m need not be a power of two. With probability
at least 1-delta, simultaneously at every completed-operation prefix,

```text
work <= I+k*R_f + O(k*W_f*J+v*P)
        + Q*O(k*T_f+s*v+8*p*(k*R_f+F_0))
        + 4*(k*R_f+F_0)*log(1/delta).
```

As before, Q,J,P count queries, key edits and fixed-width vector writes in
that prefix; v*P charges those P writes. Insertions
pay their vector-write term as well. Exact answers hold for every seed;
only this compute bound is probabilistic. True fresh random bits, private
internals and the abstract timing model are assumptions; a fixed finite
PRNG or actual-CPU constant-time guarantee is not supplied.

For fixed c>1, t=ceil(log2 n), r=floor(t/c)>=1, H=Theta_c(n^(1/c)),
and D_h=O(d), we have R_f+F_0=O(n*T_f). Therefore k=O_c(1), rather than
the O(log n) families in our previous first-moment conversion. This does
not say k is constant for fixed r as n grows. Under bounded overlap,
A=Delta+2 and s=Delta+1; the resulting bounds become

```text
query O_c(d+H*log n+H*(Delta+2)*(d+1)+(Delta+1)*v),
key edit O_c(d+H*log n), vector overwrite O(v),
space O_c(n*(mu_w+v)+n*H+d) words.
```

These are expected amortized query costs at deterministic horizons, with
the separate all-prefix confidence statement above. An individual timeout
still costs a complete scan and rebuild. If d>=log n and
H*(Delta+1)=o(n), the index-search bound is sublinear relative to scanning
all n raw d-bit keys in the conservative verification model. The v payload
cost and superlinear index storage are still explicit. No free retained
value compression, attention approximation, learned geometric promise,
or general transformer speedup follows.

### What this changes, and what the priority gate still needs

| Source/guarantee | What has been checked | Boundary for this candidate |
|---|---|---|
| [Pagh's author-linked preprint v3](https://arxiv.org/abs/1507.03225v3), 16 pages, January 2016 | Targeted construction/collision and representation passages in Sections 3/6, Lemma 3.3 and Theorem 3.5; covering, first moment and faster representations. Earlier SODA timeout/copy/fallback ingredients are recorded above | These are prior art. The author's publication list links this preprint to the 2018 journal article; the journal body was not separately compared. No exact journal/preprint identity is asserted. |
| [fcLSH v2](https://arxiv.org/abs/1602.02620v2), 14 pages | Section 3/Algorithm 2/Lemma 3/Eq. (5), plus all-neighbor/output-size discussion | The O(d+H log H) factorization and complete-neighbor retrieval are prior art. Our stronger candidate is the spectral second moment and its constant-family timing-adaptive consequence, not those algorithms. |
| [Andoni--Haris--Kelman--Onak v1](https://arxiv.org/html/2601.00272v1), Theorems 1.2/1.5 and Definitions 1.1/3.1 | Fair/density method already uses O(n^(1+rho) log(nQ)) stated space and O(d*(n^rho+density ratio)*log(nQ)) stated query time; for-all Hamming uses the stronger exponent rho=1/(2c-1). Formal game and update paragraph inspected | Do not compare only against its sqrt(Q) constructions or ignore its better exponent. Its returned point, oblivious-key-update and formal answer transcript guarantees differ; a complete exact-operator timing/adaptive-key tradeoff comparison remains open. Space units/raw storage and word costs must be reconciled before comparing formulas. |
| [Feng et al. 2025](https://proceedings.mlr.press/v267/feng25c.html), Assumption 1.3 and Appendix C | Bounded overlap and adaptive search are existing research topics; the inspected update theorem deletes returned points | The overlap hypothesis is not our discovery. A different output/update specification is not itself an impact claim. |
| [Ahle](https://thomasahle.com/papers/lasvegas.pdf) and [Wei](https://arxiv.org/abs/1807.07527) | Established Las Vegas near-neighbor guarantees, including stronger data-dependent exponents | Complete exact-range, dynamic-key, finite-word and transcript costs must be matched. This record does not establish exponent superiority or a general lower bound against them. |

The spectral identity/inequality, modular arithmetic and repetition proof
have now been derived explicitly and checked on exact finite controls.
The potentially useful advance is a robustification whose overhead matches
the underlying fast covering index up to a c-dependent constant in the
stated regime. **We have not established that this guarantee is new or
that it beats the strongest same-workload published construction.**
In particular, the identities use elementary covariance/frame bounds;
neither unsuccessful searches nor a logarithm removed from our own earlier
construction makes a paper. The next decision requires auditing optimal
Las Vegas range-search constructions and timing-adaptive transformations
against this full statement, then either identifying a substantive gap
or classifying the whole result as a baseline. No training or claimed-
contribution manuscript is warranted on the current evidence.

Reproduce the extended certificates and unchanged core checks with:

```sh
python3 -B scripts/run_timing_restart_checks.py
./build.sh
python3 tests/test_smoke.py
```

The new controls enumerate 25,152 covering-label families, verify 100 exact
mean/variance/tail identities, and compare 430,896 modular-transform hash
values with direct projected sums. All 94,480 complete weighted reads agree
with direct Hamming scans, including 58,410 reads with extra
fingerprint collisions. The report also retains sharpness/extra-bit/false-
independence negative controls and 582 integer copy-choice checks. These are
finite supporting certificates, not a production index or empirical speedup;
the uniform assertion rests on the proof above. Historical certificate blocks
and the other research artifacts are preserved.
