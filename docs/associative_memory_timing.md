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

## XOR compression: direct buckets without losing the tail bound

The next candidate removes the balanced-tree lookup and key-edit logarithm
from the preceding construction. It changes the fingerprint algebra, rather
than claiming a faster implementation of exactly the same integer hashes.
The proof below is uniform; its priority and advantage against all applicable
published methods remain unestablished. No additional neural training follows.

### Linear-time evaluation of the entire family

Write h=r+1, H=2^h-1 and choose a power of two M=2^b. A coordinate label
m_j is uniform in F_2^h as before. Independently choose each coefficient c_j
uniformly from F_2^b. Replace integer sums by bitwise XOR:

```text
g_u(x) = XOR_{j: x_j=1 and <m_j,u>=1} c_j,  u!=0.
```

This is a standard random linear fingerprint of each selected binary
projection. Unequal projections collide with probability 1/M. A true
projection collision is preserved for every coefficient choice, so covering
and exact full-distance verification still give every inner neighbor.
XOR tabulation is established prior art; see
[Patrascu--Thorup](https://arxiv.org/abs/1011.5200v2). We do not import
its concentration results as if the H correlated projected keys were
independent inputs to independent tables.

All H fingerprints can be evaluated in O(d+H) word operations, with O(H)
scratch and O(d) seed words. First form T_l=XOR_{j:m_j=l,x_j=1} c_j in
O(d+H) steps. Let P be its prefix-XOR array. For each coordinate t of a
label, compute

```text
G_t = XOR_{l: bit_t(l)=1} T_l
    = XOR_{a=2^t, 3*2^t, ... < 2^h} (P[a+2^t] XOR P[a]).
```

There are sum_t 2^(h-t-1)=H such interval queries in total. Finally
g_u(x)=XOR_{t:bit_t(u)=1} G_t. Generate the outputs by repeated doubling:
start with g_0=0; for t=0,...,h-1, set g_(2^t+a)=g_a XOR G_t for
0<=a<2^t. This uses exactly H output XORs. No per-output scan of all
h basis words, integer Hadamard transform, division, or large mask array
is needed. The RAM algorithm scans packed input words to extract their
bits; the small Python certificate uses integers and is not a timing
implementation of that input scan.

The factorization has h basis words and the original raw keys remain
available for exact checking. It is an indexing fingerprint, not a
two-scalar learned representation or a promise that arbitrary labels can
be retained in h words without the stored values and index.

### Shared coefficients: the missing covariance is explicit

Direct buckets require fingerprints of only about log2 n bits. The earlier
union bound forbidding *any* extra collision would require larger buckets
and waste space. Instead allow extra candidates and analyze their count.
They are correlated; pairwise independence of the masks is insufficient.

Fix a dataset of a<=n occupied slots and q. Retain p_i, mu, K, lambda and
v from the preceding proof, and put

```text
Q_ij=2^(-|S_i union S_j|), L_ij=2^(-dist(x_i,x_j)),
V=sum_ij Q_ij-mu^2=v^T K v-mu^2,
J=sum_ij L_ij, theta=1/M, beta=1-theta, alpha=theta*beta.
```

Let X_u be the compressed bucket's number of slots and C_xor=sum_u X_u,
including repeated visits in different tables. The exact identities are

```text
E C_xor = H*(beta*mu+a/M),
Cov(X_u,X_w)=alpha*V,                         u!=w,
Var C_xor = H*beta*(1+(H-3)/M)*V + H*alpha*(J-mu^2).
```

Here is a derivation that keeps the shared randomness. Conditional on all
coefficients, put p_i(c)=2^(-rank{c_j:j in S_i}) and mu_c=sum_i p_i(c).
The masks in distinct tables are independent conditional on c, giving
Cov(X_u,X_w)=Var_c(mu_c). Express p_i(c) as the probability that the
linear coefficient map annihilates a uniform binary subset of S_i. For
two such subsets, distinct nonzero vectors have independent uniform
images; equal nonzero vectors have the same uniform image. Therefore

```text
E_c p_i(c) = beta*p_i+theta,
Cov_c(p_i(c),p_j(c)) = alpha*(Q_ij-p_i*p_j),
Var_c(mu_c)=alpha*V.
```

For one table's shared mask, its two masked difference vectors are both
zero with probability Q_ij, equal with probability L_ij, and individually
zero with probabilities p_i,p_j. Splitting into both-zero, one-zero,
equal-nonzero and distinct-nonzero cases gives

```text
E[I_i I_j] = (1-3*theta+2*theta^2)*Q_ij
            + alpha*(p_i+p_j+L_ij)+theta^2.
```

Summing, subtracting the squared mean and adding the H*(H-1) cross-table
covariances proves the displayed total variance. For a one-bit key 0,
query 1, r=1, M=4, the cross-table covariance is 3/64: the actual total
variance is 63/64, whereas incorrectly independent compressed tables
would give 45/64. Thus the covariance term is necessary.

Suppose M>=max(n,H), H>=3 and A>=max(1,lambda). We already have
mu<=A and V<=A^2/4. Since 0<=K_ij<=1, entrywise L_ij<=K_ij, whence

```text
J <= 1^T K 1 <= a*lambda <= n*A.
```

The variance identity now gives

```text
E C_xor <= H*(A+1),
Var C_xor <= H*(A^2/2+A) <= H*(A+1)^2/2,
Pr(C_xor > 2*H*(A+1)) <= 1/(2*H).
```

All inequalities hold with multiplicities and for every fixed q. They
bound the actual compressed candidate visits, including false positives;
no fingerprint-collision-free event or unjustified conditioning on an
observed surviving seed is used. Exceptional seeds still exist: sixteen
simplex-code keys in dimension 15 have lambda=31/16, but all-zero
coefficients give 48 visits for r=1, above the A=2 budget of 18. The
fallback/rebuild mechanism remains necessary for the runtime theorem.

### A complete direct-address index and its cost

Take H<=n and M=2^ceil(log2 n), so n<=M<2n. Store an H-by-M array of
bucket heads and H preallocated nodes per slot. Each node has an ID and
previous/next handles; the slot/table cache records its current bucket.
Push and remove a node in O(1) worst-case operations, even when every
fingerprint collides. Initialization and resets clear H*M heads and
n*H nodes, and charge those operations. No second-level randomized
dictionary, hidden lookup oracle, or tree comparison remains.

To edit a key, remove its H old nodes using their handles, evaluate the
new H hashes in O(d+H), and insert its H new nodes. Remove old entries in
all private families before changing the shared key array; no query
interleaves with that synchronous update. Raw-key-named updates still pay
the separate exact dictionary cost. Pure value writes by slot ID still
cost O(v) and leave all index metadata unchanged.

Read the H bucket heads, visit candidates, deduplicate IDs and verify full
Hamming distances. The existing touched-mark cleanup and common exact
weighted aggregation apply. With the same finite-word and verification
model as above, the per-family bounds are

```text
W_x = O(d+H),                         worst-case key edit,
T_x = O(d+H+H*(A+1)*(D_h+1)),         fast complete candidate search,
R_x = O(n*(d+H)+H*M+d*(b+r+1)),       worst-case padded rebuild,
F_0 = O(n*D_h),                      exact fallback scan,
common aggregation = O(s*v).

space = O(n*(ceil(d/w)+v)+k*(n*H+d)) words,
initial work = I+k*R_x.
```

The word width w holds addresses, b-bit coefficients, label bits, bounded
operation counters and exact output numerators/denominators, as specified
earlier. Fresh coefficient/label bits are all charged in R_x. Lookup and
unlinking use only ordinary integer indexing and pointer operations; XOR
does not assume free large-integer arithmetic. Population count, when
used for D_h=O(ceil(d/w)), remains an explicit optional word operation.
With ordinary bit counting one can conservatively take D_h=O(d).

Use the preceding independent-family padding/reset scheme with
p=(2H)^(-k), the explicit copy choice based on (R_x+F_0)/T_x, and the
same all-prefix bound, replacing R_f,T_f,W_f by R_x,T_x,W_x. The fixed
seed's updated bucket index equals the index built on the current keys,
so the canonical all-normal-branch argument also covers keys chosen using
earlier times. Space is worst-case bounded at every point of the history.
Exact answers hold for every seed; runtime remains expected amortized,
with the simultaneous all-prefix confidence bound and occasional expensive
scan/reset. A single operation's query latency is not worst-case fast.

For fixed c>1, r=floor(ceil(log2 n)/c), H=Theta_c(n^(1/c))<=n for all
sufficiently large n, R_x+F_0=O(n*T_x). Thus k=O_c(1), with

```text
key edit O_c(d+H), vector write O(v),
query O_c(d+H+H*(A+1)*(D_h+1)+s*v),
space O_c(n*(ceil(d/w)+v)+n*H+d) words.
```

The supplied A and s promises still must hold at every operation boundary;
the algorithm does not verify them for free. The full raw-key/value cost,
superlinear n*H metadata and possible fingerprint-noise work remain.
This is a construction/proof, not an implemented production memory layer
or a practical transformer claim.

### Stronger-baseline audit and the remaining publication decision

The previous goal turn was **progress**: it saved and verified the spectral
tail proof. This turn's audit rules out several weak novelty arguments:

* [Ahle's full version](https://thomasahle.com/papers/lasvegas.pdf), Section 2,
  PDF pages 9--10, explicitly stores all filter memberships, states close-
  pair covering and discusses insertions/deletions and dynamization.
  All-neighbor enumeration from such a covering is a routine extension;
  exact vector payload sums do not themselves create a new search theorem.
  Page 10 was visually checked. A complete enumeration bound must include
  every visited membership, rather than reuse the first-hit query bound.
* [Christiani's filtering framework](https://arxiv.org/abs/1605.02687),
  Theorem 1.1, Appendix A's correctness calculation, and Appendix F were
  inspected. Appendix F already supplies dynamization and fixed-horizon
  logarithmic repetition; its stated word RAM uses Theta(log n)-bit words
  and allows negligible sampling/rounding error. PDF page 15 was visually
  checked. It does not supply our zero-error, all-time timing transcript
  theorem merely by saying the adversary does not know its randomness.
  Conversely, these differences are not evidence of novelty by themselves.
* [Wei v1](https://arxiv.org/abs/1807.07527v1), Appendix A.3's final proof
  paragraph, explicitly replaces representative leaves with full point
  lists so search can continue after a false positive; PDF page 21 was
  visually checked. Thus the earlier representative-leaf objection is
  already addressed. The stronger data-dependent exponent is retained
  as a serious comparator, with exact-range enumeration, adaptive key
  changes, preprocessing and finite precision still to match fully.
* [Schmitt et al., VLDB 2023](https://www.vldb.org/pvldb/vol16/p2686-schmitt.pdf),
  Section 7.3, was inspected. It reuses fcLSH's Hadamard computation and
  adds join-specific filtering of singleton lists. This targeted inspection
  is not a full audit of every later covering-signature algorithm.

There is now a precise comparison target. The published fcLSH signature
algorithm takes O(d+H log H) operations; this XOR signature construction
takes O(d+H) with a different universal fingerprint and the displayed
compressed-candidate bound. The direct buckets also replace our previous
balanced-tree key-edit bound O(d+H log n) by O(d+H), at the same asymptotic
index space. This is a better bound for our construction, not yet a claim
against the best full-workload method.

In particular, packed raw-mask evaluation is already O(H*ceil(d/w)) in a
suitable word model and can be competitive; when d fits in one word there
is no general logarithmic improvement to claim. With d=Theta(log^2 n),
w=Theta(log n), and H=n^rho much larger than d, XOR signature evaluation
has O(H) work versus the published Hadamard evaluation's O(H log n) and
straight packed projection's O(H log n). That is a concrete primitive-level
comparison, not a lower bound against every signature algorithm. Full read
time may still be dominated by candidate distance checks and payload work.

Random linear XOR hashing, binary rank factorization, prefix sums, intrusive
bucket lists and ordinary covering are established ingredients. An honest
paper would need to establish that their complete combination and shared-
coefficient spectral tail give a new, useful tradeoff against the strongest
applicable methods. The next audit must check later covering implementations
and optimal Las Vegas range/filter data structures on that full statement.
The objective, strongest-baseline gate and publication standard are unchanged.

The same three reproduction commands above run the added controls. They
compare 2,051,936 XOR fingerprint values with direct selected-coordinate
sums, enumerate 41,296 joint label/coefficient seed choices in the stated
finite cases, and verify 24 exact mean/variance/covariance identities.
The intrusive-list oracle checks 768 dynamic states, 512 key/occupancy
edits, 256 independent payload writes and 32,256 budgeted reads, including
abort cleanup and nontrivial mixtures. It does not implement the complete
timing wrapper or measure hardware performance. All earlier certificate
blocks, manuscripts, core kernels and learned-recall results are preserved.

## Stronger comparator audit and a compact sparse-key regime

The preceding goal turn was **progress**: it saved the XOR/direct-bucket
proof and independent checks. The full original publication goal remains
unachieved. This iteration addresses two concrete gaps: the strongest
data-dependent comparator and whether the spectral promise has a useful
regime in moderate dimension. The probability and coding arguments below
are supporting lemmas built from classical tools, not separate novelty
claims or grounds to declare the goal complete.

### What the optimal hashing paper already supplies

The [ALRW full author version](https://www.mit.edu/~andoni/papers/tradeoffs-full.pdf)
is 63 pages, dated March 21, 2017. Theorem 1.1 and the paragraph after the
Hamming/ell_1 tradeoff explicitly include insertion/deletion support using
dynamization, with logarithmic losses in space/query time and update
exponent rho_u. Thus a claim that this comparator is inherently static
would be wrong. The nominal Hamming tradeoff is

```text
c*sqrt(rho_q)+(c-1)*sqrt(rho_u) >= sqrt(2*c-1),
rho_q=rho_u=1/(2*c-1) in the balanced case.
```

Sections 4.2--4.4 and Lemma 4.13 were inspected in detail. Its proof
separately counts visited nodes and tested level-K point occurrences,
including multiplicities. So it is also wrong to dismiss this bound as
counting only tree navigation. However, exact range enumeration still
needs care: the query may stop at an acceptable outer neighbor; some
earlier base cases store just one point. As recorded above, Wei already
changes those representatives into complete lists to continue past false
positives. Full lists remove a correctness obstacle but do not by themselves
prove that scanning them has the first-answer running time.

Here is the accounting a complete comparison must supply. Separate the
root query's points at distance at most the outer threshold from all other
points. If each point is stored in at most B_p expected leaf occurrences,
all visited near occurrences contribute at most sum_near B_p. If expected
far occurrences are at most B_far and navigation costs B_nav, a complete
verified range read costs

```text
O(B_nav + D_h*(B_far + sum_near B_p) + s*v).
```

No independence of occurrences is needed for this elementary upper bound.
Claim 4.12 gives the relevant expected replication bound; Lemma 4.13 gives
the navigation and candidate recurrences. A fully justified extension
would have to handle all of the following before plugging in their optimal
exponents:

* Each visited coarse base case must either certify all its points are
  within the root outer threshold or charge the other scanned points to
  B_far. A valid-neighbor representative alone is not such a certificate.
  For a local ball B(o,R), testing ||q-o||<=R+r_1 first is safe for inner
  search; if r_1+2R<=r_2, all its points are then within local r_2.
  Radial rounding and threshold propagation must connect this local
  assertion back to the original Hamming distance.
* The one-sided dimensional reductions used by Wei can introduce false
  positives. A first-hit bound does not automatically bound all replicated
  false positives in a full read. Their expected multiplicity and the
  precision/rounding costs belong in B_far and B_p.
* The dynamization, build/rebuild space, private operation times, accepted
  seed law and canonical full outputs need matching. Supporting updates
  nominally does not by itself prove the present timing-adaptive theorem;
  differing formal guarantees are also not evidence of novelty.

For orientation only, if these steps yielded B_nav+B_far=n^(rho_q+o(1))
and B_p=n^(rho_u+o(1)), the complete cost would be
O(D_h*[n^(rho_q+o(1))+k_outer*n^(rho_u+o(1))]+s*v), where k_outer counts
root outer neighbors. This is a **conditional accounting consequence**, not
an established full-range/finite-word theorem of ALRW or Wei. Under bounded
k_outer its balanced exponent would be polynomially better than 1/c;
for c=2, 1/3 replaces 1/2. We therefore retain this as a serious comparator,
not just the published sqrt(Q)-copy robust search alternatives. The remaining
full-workload proof is explicit rather than an inferred nonexistence result.

### A uniform sparse-shell feasibility theorem

Take an even integer t>=16, n=2^t and d=16*t^2. Independently choose each
key's support S_i uniformly among all t-subsets of [d], initially allowing
repetitions. Define K_ij=2^(-dist(x_i,x_j)/2) as before. With probability at
least

```text
1 - delta_pair(t) - n^(-2),
delta_pair(t)=n^2*(3/(8*t))^(t/2),
```

all distinct slots have distance greater than t and ||K||_2<=3. This is a
single dataset event and hence supports the previous bounds for **every
raw query**, including adaptively chosen queries. No query-net union bound
or random-query assumption is used. For t=16, delta_pair=6561/2^24;
it decreases by more than a factor of three when t increases by two.

For a proof, conditional on any one support S_i, its intersections
I_j=|S_i intersect S_j| with the other keys are independent hypergeometric
variables. Their factorial moments satisfy, for 0<=l<=t,

```text
E binom(I_j,l) = binom(t,l)^2/binom(d,l)
              <= binom(t,l)/(16*t)^l.
```

Expanding z^I=(1+(z-1))^I and summing gives

```text
E[2^I] <= (1+1/(16*t))^t <= 16/15 < 3/2,
E[4^I] <= (1+3/(16*t))^t <= 16/13 < 2.
```

The last bounds follow by bounding binomial coefficients by t^l/l! and
then the exponential series by the corresponding geometric series.
With l=t/2, the same factorial moment and binom(t,l)<=(e*t/l)^l give
Pr(I>=t/2)<=(3/(8*t))^(t/2), since e<3. A union bound over pairs yields
delta_pair. Off this event every pair has distance 2t-2I>t, and duplicates
are excluded as well.

It remains to control row sums without assuming the entire matrix has
independent entries. Put Y_j=2^(I_j) if I_j<t/2 and zero otherwise.
Conditional on S_i, the Y_j are independent, bounded by sqrt(n), have
mean below 3/2 and second moment below 2. The usual bounded-variable
exponential-moment calculation gives

```text
Pr(sum_{j!=i} Y_j > 2*n | S_i)
  <= exp(-n/(16+(4/3)*sqrt(n))).
```

For completeness, a centered variable Z with |Z|<=L obeys
E exp(aZ)<=exp(a^2*E[Z^2]/(2*(1-aL/3))) for 0<a<3/L, by expanding the
series and using k!>=2*3^(k-2). Multiplying these bounds for independent
variables and choosing the exponential parameter proves Bernstein's bound
with deviation n/2, variance sum at most 2n, and L=sqrt(n), giving the
displayed exponent. This is an application of a classical concentration
argument, not a new concentration inequality.

Union bounding over n rows costs at most n^(-2): the exponent
3n/(48+4sqrt(n)) is at least 3t for even t>=16, and e>2. Its value at t=16
satisfies this, and when t increases by two the exponent at least doubles.
On the complement of both bad events, K_ij=2^(I_j)/n off the diagonal,
so every row sums to at most 1+2n/n=3. The maximum-row-sum bound completes
the spectral proof. Rows themselves were never assumed independent.

This is substantially more compact than the earlier disjoint-block witness:
d=Theta(log^2 n), while the radius-t ball centered at zero contains **all n
keys**, and every key has n-1 others within distance 2t. Thus both the global
outer-density and bounded-overlap promises can fail maximally while this
spectral promise holds with high probability. The zero query has mu=1
exactly. It also admits a trivial weight-based pruning rule, so this query
or dataset is not a lower-bound witness against every other index.

### The dimension order is necessary within this model

For any n=2^t binary keys of common weight t in dimension d, repetitions
allowed, let f_j be the number using coordinate j. Then

```text
(1/n^2)*sum_ij |S_i intersect S_j| = sum_j f_j^2/n^2 >= t^2/d.
```

The Rayleigh quotient of K at the all-ones vector, followed by Jensen,
therefore gives

```text
||K||_2 >= (1/n)*sum_ij K_ij
        = (1/n^2)*sum_ij 2^(|S_i intersect S_j|)
        >= 2^(t^2/d).
```

Consequently ||K||_2<=A for a fixed A>1 requires
d>=t^2/log2(A). The order Theta(log^2 n) is thus both sufficient and
necessary for bounded norm in this particular n=2^t, weight-t model.
This is an elementary Rayleigh/Cauchy--Schwarz/Jensen consequence. It is
not a dimension lower bound for arbitrary codes, encoders, attention
architectures or associative memories.

### Retention, mixtures, probability and the full cost

For the unpaired keys on the good event, any r<=t/2 Hamming ball contains
at most one key, so arbitrary independent fixed-grid payloads are recalled
exactly throughout each key's radius-r ball, with no overwrite interference.
This is noisy raw-vector lookup, not a supplied target ID.

There is also a two-neighbor version with the specified exact mixture
semantics. Replace each center x_i by (x_i,0) and (x_i,1). The key matrix
becomes K tensor [[1,1/sqrt(2)],[1/sqrt(2),1]], whose norm is below 6.
The total capacity is N=2n and dimension d+1=Theta(log^2 N). For
1<=r<=t/2, a query meets at most one center's pair, so s<=2. At the key
(x_i,0) the weights of the pair are r+1 and r, yielding a genuine weighted
blend of independently replaceable vector payloads. Do not call that blend
the individual value at the first key. Each raw vector is still retained
explicitly. All N keys lie within distance ceil(log2 N)=t+1 of zero.

For fixed c>=2 and r=floor(log2(N)/c), the preceding direct-bucket theorem
therefore has A=6, s=2, H=Theta_c(N^(1/c)), and in its stated word model

```text
space O_c(N*(ceil(d/w)+v)+N*H+d) words,
worst-case key edit O_c(d+H), vector write O(v),
expected amortized read O_c(d+H*(D_h+1)+v),
initial work I+O_c(N*(d+H)) including fresh index bits and cleared heads.
```

The already proved all-prefix cost confidence statement applies conditional
on the supplied geometry. Geometry failure and index-runtime failure are
different probability spaces; their probabilities may be union bounded.
The dataset theorem treats raw keys as input. Generating random supports,
loading all raw bits/values, or verifying a proposed pool is not free.
An unbounded sequence of fresh random replacements is not certified by
this one-draw probability bound: future key edits must still preserve the
stated geometry, or pay for an admission/verification rule. No learned
encoder or semantic accuracy result is supplied.

The compact feasibility/dimension theorem closes a model-regime gap. It
does **not** establish novelty of the XOR/direct-bucket theorem, superiority
to optimal exact-range structures, or a publication-ready ML result. The
next decision remains a full same-workload prior-art comparison; supporting
coding and probability lemmas do not replace that gate.

## Timing attacks with bounded kernel norm and identical answers

### A closer published reduction and the exact attack boundary

The [July 2016 version of Ahle--Aumuller--Pagh](https://arxiv.org/abs/1605.02673v2),
Appendix G (PDF pages 20--21, printed pages 19--20), already discusses
Indyk's repeated-neighbor-search/deletion reduction and a stronger
list-of-points range algorithm. It scans selected lists and removes each
reported inner neighbor from every list, obtaining the stated expected
O(n^rho_q+t*n^rho_u) bound and using the optimal Euclidean hashing
tradeoff. The latter display was visually checked. Its Definition 1 allows
each inner point to be reported with constant probability. It does not
give the present all-seed exact weighted output with private-operation-time
feedback. Nevertheless, range enumeration and deletion are existing
techniques; a different output definition is insufficient for novelty.
The full finite-word, timing-adaptive optimal-Las-Vegas comparison remains
open in this record. This appendix substantially narrows the comparison;
it is not evidence that a stronger complete comparator does not exist.

[Kapralov--Makarov--Sohler, March 2024 version](https://arxiv.org/html/2402.09707v2),
Section 2, explicitly raises possible adaptive runtime inflation of
Las Vegas LSH and says investigation of such adversaries is outside its
scope. Its proved attack concerns false negatives of Monte Carlo LSH.
The later [Andoni--Haris--Kelman--Onak preprint](https://arxiv.org/html/2601.00272v1),
Section 2, describes Las Vegas runtime vulnerability while citing that
paper. Thus the vulnerability itself is already in the literature. The
earlier paper is not a proof of the compact, constant-norm construction
below. These two source statements must not be collapsed into a claim
that we discovered the issue or that the earlier paper proves this theorem.

### A uniform obstruction for permanent covering families

Fix an integer s>=1. Consider s independent copies of the earlier covering
family, whose coordinate labels are uniform in F_2^(r+1). Keep these
families permanently installed, without timing padding, refreshes, or a
timeout/fallback rule. A complete family reader visits all selected bucket
entries and verifies the original keys. It can use raw projections or the
previous XOR compression. Selecting the fastest of the s complete readers
or running them in a fixed fair race does not evade the result. This is a
restriction on this reader, not on all geometric indexes: extra geometric
certificates, different indexes and algorithms are outside the obstruction.

For every even k>=16 there is a fixed dataset and a fixed list of raw
queries, chosen independently of the installed labels, with

```text
t=(8*s+1)*k, n=2^t, m=2^k, L=8*k, G=2^(s*L)=n/m,
r=k/2, H=2^(r+1)-1, c=2*(8*s+1), so r=log2(n)/c,
D=16*(2^s-1)*L^2, d=8*t+D=O_s(log^2 n).
```

The n keys are distinct and have a common Hamming weight, as do all G
queries, with the same weight as the keys. Their kernel matrix satisfies
||K||_2<3. Every listed query has no radius-r neighbor, so every canonical
answer is exactly bottom, for arbitrary independently stored payloads.
Nevertheless, with probability at least 15/31 over the installed covering
labels, at least one listed query forces every complete family reader to
visit at least m distinct keys. This remains true for every choice of the
XOR compression coefficients. Thus m>(H^2)/4 visits are needed on that
query even though the expected candidate count of each fixed query in
each compressed family is at most 4H. The bad query depends on the labels;
the dataset and the finite query list do not.

An adversary can submit all G listed queries, observe their completed RAM
times and repeat a query of maximum observed time. Deterministic read-only
queries under the same installed seed reproduce that time. On the event
above each subsequent query costs at least m word operations, while all
answers remain bottom. The adversary neither inspects memory nor receives
a varying answer. Its offline dataset/query-list construction and G probes
are explicit costs. This is not an efficient polylogarithmic-probe attack.
This initial attack also assumes no result memoization. Caching the probed
answers defeats its repetition step; the following extension removes that
particular limitation.

### Proof of the geometry and existence assertions

First, for every integer a>=16 there is a constant-weight binary code of
length 8a, weight 4a, minimum distance at least 2a and at least 2^a+1
words. A greedy packing proof suffices. There are at least
2^(8a)/(8a+1) weight-4a words. A forbidden Hamming ball of radius 2a-1
has at most

```text
V <= sum_{j=0}^{2a} binom(8a,j)
  <= 3^(2a)*(1+1/3)^(8a) = (65536/729)^a
```

words, including words of other weights. The second inequality follows
because (1/3)^j >= (1/3)^(2a) in this sum. The greedy algorithm therefore
selects at least (729/256)^a/(8a+1) words. This bound exceeds 2^a+1 at
a=16 by an exact rational comparison. Its ratio to 2^a is increasing for
a>=16, since (729/512)*(8a+1)/(8a+9)>1, while 1+2^(-a) decreases.
This proves the uniform claim without a numerical entropy approximation.

Use that code with a=k for local words z_0,z_1,...,z_m. Use it with a=t
for at least G distinct coarse identifiers y_1,...,y_G. These identifiers
have weight 4t and pairwise distance at least 2t. For each cluster j choose
an L-subset B_j of [D], and embed the local words in its coordinates,
using increasing coordinate order. Define

```text
x_(j,i) = (y_j, embed_(B_j)(z_i)), i=1,...,m,
q_j     = (y_j, embed_(B_j)(z_0)).
```

Every word has weight 4t+4k. Within a cluster the distinct words are at
distance at least 2k. Between clusters the coarse identifier alone supplies
distance at least 2t. In particular every q_j is outside every inner ball.
For any choice of the B_j, the kernel's maximum row sum is bounded by

```text
1 + (m-1)*2^(-k) + (n-m)*2^(-t) < 3.
```

The norm bound follows from symmetry and the row-sum bound, with no
random-geometric assumption. The arbitrary payloads need no special
encoding and remain stored once as before. The attack concerns reads;
there are no hidden deletion or key-admission operations.

To choose the B_j, temporarily sample them independently and uniformly
among L-subsets. In each installed family, consider its table indexed by
one fixed nonzero vector u. For each coordinate its mask bit is fair,
independently across coordinates and families. Let E_j be the event that
all s mask vectors vanish on B_j. Then Pr(E_j)=2^(-sL)=1/G. On E_j, all
m keys in cluster j collide with q_j in that table of every family. This
is equality of the raw projections and hence survives all XOR compression.

For two independent L-subsets let I be their intersection size. For
0<=ell<=L, the exact factorial moment is

```text
E binom(I,ell) = binom(L,ell)^2/binom(D,ell)
               <= binom(L,ell)*(L/D)^ell.
```

Expand 2^(sI)=(1+(2^s-1))^I. The chosen D gives

```text
E[2^(sI)] <= (1+1/(16L))^L <= 16/15,
Pr(E_j intersect E_j') <= (16/15)*G^(-2), j!=j'.
```

For the last inequality in the first line, bound each binomial coefficient
by L^ell and sum the geometric series. Writing W=sum_j 1_(E_j) gives
E W=1 and E W^2<31/15. Cauchy--Schwarz applied to W*1_(W>0) yields
Pr(W>0)>=15/31. This probability currently averages over both the subsets
and the private masks. Consequently at least one fixed choice of all B_j
has mask-event probability at least 15/31. This is an existence proof of a
seed-independent dataset, not an efficient deterministic codebook/subset
construction. Once such a dataset is fixed the preceding geometry remains
valid for every mask and query.

### Fresh queries: answer memoization does not repair the reader

There is a stronger version for the XOR/direct-bucket reader with its
specified random XOR coefficients, and a fixed
fair race that stops when one of the s families finishes. Permit a cache
of previously answered queries. On a cache miss, the algorithm still uses
the specified complete reader, rather than a geometric certificate or a
precomputed region of answers. Cache lookups must have the usual bounded
word cost. All adversarial queries in this extension are distinct, so
memoization of previously answered inputs never supplies an answer.

Append an auxiliary block of U=D coordinates, equal in every key and base
query, with half its bits zero and half one. Let d'=8t+D+U. Pairwise key
distances, the norm bound, the collapse events E_j and their second-moment
proof are unchanged. The keys and base queries still have one common
weight. All query modifications below are confined to this auxiliary
block, so they cannot reduce distance to any key: these coordinates
initially agree with every key.

For the fixed implementation choose B_s>=1 such that, for each family f,
the completed fair-race RAM time, including cache lookup, is at most

```text
B_s*d'*(H+C_f(q)+1),
```

where C_f(q) is that family's full compressed candidate count. Such a
constant exists for a fixed s under the preceding word-cost reader and
fixed round-robin simulation. It is an implementation constant, not a
measured performance parameter. Let mu=U/2^(s+1), and choose even k large
enough that

```text
mu/2 >= 2t+1,
epsilon = 8*B_s*U*d'*H/m + 8/mu <= 1/8.
```

For each fixed s and B_s these inequalities hold for all sufficiently
large k: mu grows quadratically, while U*d'*H/m is a polynomial in k times
2^(-k/2). This is a uniform asymptotic extension, not a claim that the
original k=16 endpoint meets these extra inequalities.

Conditional on E_j, call an auxiliary coordinate good if its distinguished
mask bit is zero in every family. The good coordinates in each initial-bit
half have independent binomial counts of mean mu. Chebyshev and a union
bound show that both counts are at least mu/2 except with probability
8/mu. The adversary need not observe these masks. It tests the U single-bit
modifications of q_j, and calls a coordinate good when the completed time
is at least m/2. A truly good coordinate preserves all m collisions in
every family and necessarily passes this timing test.

For a coordinate that is not good, consider the first family whose
distinguished mask bit there is one. Its distinguished table now has no
raw collisions with its own cluster. For every other nonzero table vector,
the mask bits on a local difference remain fair and independent after
conditioning on E_j, giving own-cluster raw collision probability at most
2^(-2k). Other clusters have at least 2t unconstrained coarse differences,
giving raw collision probability at most 2^(-2t) in every table. The XOR
coefficients are independent of all these mask events, so compression
adds at most 1/M collision probability per point/table. With M>=n,

```text
E[C_f(q_j with that bit flipped) | E_j, that family's bit is one]
  <= H*n/M + (H-1)*m*2^(-2k) + H*n*2^(-2t)
  <= H*(1+1/m+1/n) <= 2H.
```

Selecting the first broken family adds no condition on its other local or
coarse mask bits or on its XOR coefficients. The time bound therefore gives
expected fair-race time at most 4*B_s*d'*H. Markov's inequality bounds a
false timing classification by 8*B_s*d'*H/m. Union bound over U tests,
and add the two small-half events. Conditional on each E_j, all
classifications are correct and both halves are large with probability
at least 1-epsilon. No independence between these tests is assumed.
The tests themselves are fixed single-bit inputs; the argument never
conditions a fixed-input expectation on surviving earlier timing tests.

As in the preceding existence proof, average over subsets and masks. Since
sum_j Pr(E_j)=1, a union bound on the bad classifications/counts loses at
most epsilon, even though the cluster events are dependent. Hence there
exists a fixed subset assignment for which, with probability at least
15/31-epsilon>1/3 over labels and XOR coefficients, some collapsed cluster
has correctly classified large
halves. All auxiliary labels and compression coefficients remain private.

For each cluster in order the adversary queries its base and U single-bit
tests. If there are enough purported good coordinates in both halves,
pair 2t+1 coordinates from the zero half with 2t+1 from the one half.
It tries the distinct nonempty subsets of these pairs, flipping both bits
of each chosen pair. These queries preserve the common Hamming weight and
are distinct from the base and single-bit tests. Abandon this cluster on
the first completed time below m/2. Count every trial of time at least
m/2, and stop upon collecting n^2 such trials in total. A truly good
cluster supplies at least 2^(2t+1)-1>=n^2 distinct trials, all of time at
least m. If an incorrectly classified cluster passes many tests, those
tests already contribute the required expensive queries; it need not be
identified as incorrect in advance.

Different coarse identifiers ensure that queries from different clusters
are distinct. There are at most G*(U+1) base/single-bit tests, at most G
failed trials and n^2 successful trials. Thus at most

```text
Q_* = n^2 + G*(U+2) = (1+o_s(1))*n^2
```

distinct queries are used. With probability greater than 1/3 the total
work is at least n^2*m/2, giving expectation at least n^2*m/6. All answers
are bottom. An exact cache of prior answers cannot eliminate these
cache-miss reads. The classification queries can have weight one above or
below the common weight; every generated subset query has exactly the
keys' weight. Thus a simple total-weight rejection does not remove the
expensive generated reads. This still excludes more powerful geometric
certificates and precomputation of unasked queries, and does not claim an
attack on every possible indexing architecture.

### What advantage follows, and what does not

For fixed s and Q=n^2, the probing phase G=n/m is negligible compared with
Q. The unwrapped reader has expected total work at least
(15/31)*(Q-G)*m under the described adversary. By contrast, the previously
proved protected construction, applied to the same fixed dataset, has
expected total work

```text
O_s(n*(d+H)+Q*(d+H*(D_h+1)))
```

on these empty-output queries, including its initialization, refreshes,
fresh random bits, direct heads and candidate checks. The O_s number of
protected copies is chosen by its theorem; it need not equal the s
permanent families attacked here. Raw keys and independently stored
values, if loaded, add the previously specified common input cost; choose
v=O(1) here or keep that cost explicit. For D_h=O(d) the ratio is at least
Omega_s(n^(1/(2*(8*s+1)))/log^2 n). This is a polynomial comparison
against the specified permanent-family reader, not a comparison with the
strongest geometric index. The existing all-prefix upper bound supplies
the corresponding runtime-confidence statement for the protected reader.

The s=1 case already gives d=Theta(log^2 n), ||K||_2<3 and a slow read
with Omega(H^2) visits. Allowing every fixed s shows that simply installing
a small fixed bank of permanent families is not a universal cure as the
approximation parameter grows. The geometry changes with s and c=16s+2;
this is not an impossibility for an arbitrary constant number depending
on c, nor a matching lower bound on our protected copy count.
For the fresh-query extension replace d by d' in the protected bound and
use Q_* instead of Q. Its same polynomial comparison now survives exact
memoization of previously answered queries. Both constructions retain the
same asymptotic storage order for fixed s; the protected number of families
can have a different constant.

This supplies a concrete obstruction within the promised geometric regime
and rules out an answer-only explanation of robustness. Its ingredients
are elementary coding, collision and second-moment arguments. Priority
for this exact combination is unresolved, and the adversary uses many
probes. It neither establishes a new general indexing lower bound nor
finishes the original impact/publication goal. The decisive next gate is
still the optimal complete-range/timing comparison, not neural training
or a manuscript claiming a result beyond the evidence.

## Complete-range baseline in an arithmetic model

The previous iteration was **progress**: the compact timing obstruction,
including its fresh-query/cache-resistant extension, was proved and checked.
It compares specified permanent readers. This iteration addresses the
stronger complete-range comparator instead of extending that restricted
lower bound. The argument below is a derived baseline consequence of
published methods and elementary accounting, not a new contribution claim.

### Statement and model boundary

Fix c>1, binary raw keys in dimension d=log^(O(1)) n, and integer radius
1<=r<=d. Define k_out(q)=#{i:dist(x_i,q)<=c*r} and s(q) as the number of
inner neighbors. For fixed nonnegative rho_q,rho_u satisfying

```text
c*sqrt(rho_q)+(c-1)*sqrt(rho_u) >= sqrt(2*c-1),
```

the published Las Vegas filtering/recursion methods yield a **static-key,
pointwise expected arithmetic-model baseline** for the same complete exact
weighted operator. Its geometric preprocessing and space are
n^(1+rho_u+o(1)); its query work is

```text
n^(rho_q+o(1)) + k_out(q)*n^(rho_u+o(1)) + O(s(q)*v).
```

Polylogarithmic factors, original-distance checks and projection work are
absorbed in n^o(1). Raw keys and values are stored once in addition to the
geometric index: O(n*ceil(d/w)+n*v) finite words. Payloads may change at
cost O(v), without rebuilding any geometric index. Keys are fixed in this
statement. No two-coordinate query representation or supplied target ID
is assumed. Every inner neighbor is retained, original Hamming distances
are checked, and the numerator/denominator of the specified weighted mean
are accumulated exactly. A word width covering the S-bit payloads and
log(n*(r+1)) guard bits suffices for this aggregation.

**Geometry here uses the exact-real arithmetic model of the underlying
construction**, including rotations, square roots and filter comparisons.
Geometric cells are not asserted to fit finite words. Finite payload bits
are counted separately. This is not a finite-word implementation theorem,
an adaptive-key theorem or a new fastest-index claim. In particular, it
cannot be compared literally with our finite-word storage bound without
compiling that geometric representation. It does remove the complete-output
objection to this strong baseline in its native cost model.

### Published components used, with their roles kept separate

The [ALRW full version](https://www.mit.edu/~andoni/papers/tradeoffs-full.pdf),
Figure 5 and Section 4.2, supplies the radial recursion and its thresholds;
Claim 4.12 bounds per-point leaf replication; Lemma 4.13 separately bounds
navigation and point occurrences. Figure 5 (PDF page 29) was visually
checked. These are the components used below, rather than a black-box
assumption that first-answer time also bounds enumeration.

[Wei v1](https://arxiv.org/abs/1807.07527v1), Lemma 18, Section 6.2.1 and
Appendix A, supplies one-sided orthogonal reductions, Las Vegas spherical
filters and the complete-list modification at representative leaves.
Section A.3 calibrates the far-filter ratio with rho_q, while the final
right-hand side in Proposition 25 prints rho_u. Our asymmetric calculation
uses the displayed ratio G(cr)/G(r) and its proof, not that inconsistent
right-hand side. The balanced comparison below has rho_q=rho_u and does
not depend on this notational discrepancy. The previously audited
range-reporting Appendix G is a closely related enumeration precedent.

### Full traversal before dimension reduction

Replace every representative terminal node by its full point list, and
traverse every selected branch rather than stop after the first acceptable
point. This retains all inner pairs because the Las Vegas filters cover
every such pair. Keep raw slot IDs, mark and later clear reported IDs to
deduplicate, verify the original distance and then fetch the current
payload. Scratch cleanup is linear in touched IDs; no unbounded timestamp
counter is needed. A failed radius test never suppresses another point.

For a fixed root query classify a point as outer if its root distance is
at most r_2, and far otherwise. Each outer point is scanned at most as
often as it is stored. The per-point replication induction still bounds
its expected full-list occurrences by n^(rho_u+o(1)): at an early terminal
node a processed point now contributes one copy, just as at a level-K
terminal; cluster duplication and filtering steps are unchanged. Summing
over the outer points gives k_out*n^(rho_u+o(1)).

For the other points, the radial map has an important **outer pullback**
property. If rounded radii are R_1,R_2>0 and their rounding displacement is
at most h=delta*r_1 per point, radial projection gives

```text
projected_distance^2
  = (R_1/R_2)*(rounded_distance^2-(R_1-R_2)^2).
```

The constructor uses this formula with outer threshold r_2-2h and inner
threshold r_1+2h. Thus projected distance at most the projected outer
threshold implies rounded distance at most r_2-2h, and original distance
at most r_2 by the triangle inequality. Original inner distance at most
r_1 implies projected distance at most the projected inner threshold.
In eligible annuli |R_1-R_2|<=r_1+2h; the parameter separation r_2>r_1+4h
keeps the outer radicand positive. The published fixed-c, small-delta
regime meets this condition.

Consequently a root-far point stays outside the local outer threshold
along every query path that also contains that point. At a coarse
spherical leaf, r_2>=2R certifies every point in that leaf as local outer.
At a coarse ball leaf require the safe query test ||q-o||<=R+r_1 before
scanning; r_1+2R<=r_2 then gives the same certificate. Root-far points
cannot contribute at these leaves. This explicit ball test may prune
extra work and preserves all inner pairs. A zero-radius data child is a
single geometric location; test that location directly. A query at the
ball center has known distance to each annulus, so retain only annuli
whose distance interval can meet the inner radius. In these degenerate
cases, any scanned points satisfy the same local outer certificate; a
singular radial division is unnecessary.

The navigation recurrence A in Lemma 4.13 already counts all selected
branches. Its far-point recurrence B can likewise count only root-far
occurrences, normalized by the total node population. When r_star=r_2,
these points are local far by the pullback just proved, so the calibrated
G(r_2)/G(r_1) bound applies. For a pseudorandom remainder, the small-cap
density bound controls all points, hence also this subset. The cluster
term and the terminal bound B(u,K)<=1 remain upper bounds. No early
return is used in these recurrences. They give expected navigation and
root-far occurrences n^(rho_q+o(1)). This proves the displayed
full-traversal accounting on the preprocessed low-dimensional sphere.

### One-sided projections without an unjustified independence step

For completeness, the following construction handles our polylogarithmic
binary dimension with one orthogonal reduction, rather than copying a
first-hit false-positive bound. Choose e=(log log n)^(-1/2), a constant
integer J>rho_u+6, delta_proj=n^(-J), and block dimension
B=Theta(e^(-2)*log(1/delta_proj))=Theta(log n*log log n). Pad the dimension
to a multiple of B and split a random rotation into a blocks P_i. Apply
g_i=sqrt(a)*P_i. If the original dimension is below B, pad to B and use
the single identity block instead.

For every vector z, (1/a)*sum_i ||g_i z||^2=||z||^2, so each inner pair
has a projected inner pair in at least one block for every rotation.
For each fixed nonzero z and block, Lemma 18 bounds distortion exceeding
e by delta_proj. A root-far original pair has squared distance greater
than c*r. Its projected distance exceeds (1-e)*sqrt(c*r) unless that
block underestimates the pair. Union over n points and a blocks bounds
the probability of **any** such underestimation for a fixed query by
n*a*delta_proj. The bad event is allowed to correlate with the entire
subsequent index.

Each scaled projected point or query has norm at most U=sqrt(a*d), where
d now includes zero padding. Lift a point z to
(z,sqrt(R^2-||z||^2)) with R^2=U^2*(1+1/e). This map is noncontracting,
and for z,z' in that ball its squared distance expands by at most 1+e.
Indeed the added coordinate difference is
(||z||^2-||z'||^2)/(sqrt(R^2-||z||^2)+sqrt(R^2-||z'||^2)); Cauchy--Schwarz
bounds its square by e*||z-z'||^2. Its radius normalized by sqrt(r) is
polylogarithmic. The low-dimensional spherical construction therefore
applies with inner radius sqrt((1+e)*r), outer radius (1-e)*sqrt(c*r),
and squared ratio c*(1-e)^2/(1+e)=c-o(1). This changes the desired
exponents only by o(1). Lifting cannot introduce additional contractions
of far projected pairs. The dimension and normalized radius remain in
the required regime: replacing the original fixed polylogarithmic starting
radius by another fixed power of log n changes only the constant in the
O(log log n) radius-depth calculation. If a=1 was the identity block, there is no projection
underestimation event at all.

Build all block indexes in an attempt. Impose construction-peak geometric
space cap M equal to a constant multiple of a uniform expected peak bound.
Such a bound has the same exponent: allocated cells are charged to
construction work, whose expected exponent matches index space. Raw input
and values are retained separately. Also cap total work at a constant
multiple of the uniform expected preprocessing bound.
Abort and independently retry the entire attempt on overflow. Markov
and a union bound give acceptance probability alpha>=1/2. All completed
attempts still cover every inner pair; a space/time cap never certifies
an incomplete index. Attempts may be padded to fixed slots as in the
earlier safe-build lemma, so the retry count does not expose the accepted
seed. Expected construction work changes by only a constant factor.

On an accepted attempt, all false-point occurrences together are bounded
by M. Therefore their expected contribution under the accepted law is at
most

```text
M*Pr(any projected underestimation)/alpha
  <= M*n*a*delta_proj/alpha.
```

This is valid even with arbitrarily correlated projection error and
replication. For sufficiently large n, M<=n^(rho_u+2) and a<=n, so the
chosen J makes it negligible. The extra logarithm in 1/delta_proj changes
B only by a constant factor. The choice of J precedes the random attempt;
there is no circular selection based on its realized space.

For good original-far pairs use the spherical far-occurrence bound,
conditional on the rotation. For original outer points use the uniform
per-point storage bound, also conditional on the rotation. Restricting
these nonnegative quantities to accepted attempts and dividing by alpha
costs at most a constant. Projection work and the a block indexes cost
polylogarithmic factors. This proves the stated complete-query baseline
without multiplying two dependent expectations.

### Timing protection, changing values, and the resulting decision

Suppose k_out(q)<=K for every allowed raw query, for a supplied bound K.
The earlier generic static timing conversion applies in the same arithmetic
model, with full-data fallback and fixed query budget based on
T=n^(rho_q+o(1))+K*n^(rho_u+o(1)). Construction plus fallback is polynomial
in n, so O(log n) independent copies suffice. It gives horizon-independent
storage and expected amortized geometry work O(T*log n); the already
proved all-prefix confidence bound applies as well. Initialization includes
all copies and raw input loading. This application has no new conversion
lemma. Geometry and unfinished-query cleanup are padded; exact payload
aggregation after the geometric result is complete costs O(s(q)*v).
Only its true canonical neighbor set determines that operation count.
Pure payload writes cost O(v) and do not alter any filter, so they preserve
this static-key statement with one shared current payload array.

In particular for K=O(1), the balanced choice rho_q=rho_u=1/(2c-1) gives
smaller read and geometric-space exponents than 1/c, even after the
logarithmic timing-conversion overhead. **In this arithmetic model the
bounded-outer-density/static-key case is not an advantage for our current
construction.** Exact full outputs and the ability to change values do
not rescue that claim. The finite-word implementation still needs a
separate bit-cost proof; that remaining accounting issue alone is not a
reason to claim a new mathematical frontier.

Key deletion can physically unlink full-list memberships using stored
handles, at expected n^(rho_u+o(1)) cost for a fixed externally specified
history. The published dynamization is also relevant to insertions, but
this static proof is not a burst-cost/deamortization or adaptive-key
theorem. A complete dynamic comparison must include rebuilding slots,
phase-max space and the observation of every key-edit time. Likewise the
spectral promise alone permits unbounded outer density; this baseline's
k_out term has not been eliminated under that promise. Those are the
remaining precise research questions. No learned experiment or new paper
claim is justified by the baseline consequence proved here.

## A finite-word signature baseline for the compact sparse examples

The previous iteration was **progress**: it closed the complete-range
comparison in the published arithmetic model. The compact spectral examples
still have large outer balls, but that alone does not make them hard for
known search methods. This iteration gives a deterministic finite-word
comparator that is asymptotically stronger on a high-probability subfamily
of those same examples, including their two-payload mixture extension.
It is an elementary signature-dictionary consequence, not a new result
beyond prior art or completion of the publication goal.

### Prior art that changes the comparison

[Ahle--Aumuller--Pagh](https://arxiv.org/abs/1605.02673v2), Equation (1) and
Theorem 3, already optimize expected query work using the full collision
profile, rather than just a coarse outer-neighbor count. A large outer
ball is therefore insufficient evidence of an advantage even against
LSH-based methods.

[Arasu--Ganti--Kaushik, VLDB 2006](https://www.vldb.org/conf/2006/p918-arasu.pdf),
Sections 3 and 7, gives exact signature filtering followed by verification.
Before prefix pruning, its unit-weight minimal-subset enumeration produces
all subsets of the chosen size. Figure 8 and the accompanying WTENUM text
on PDF page 7 were visually inspected in the
[author-hosted copy](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/vldb06.pdf).
Our baseline below uses that elementary unpruned subset family. We do not
attribute its dynamic/promise-specific cost analysis to the paper.

[Christiani--Pagh, Set Similarity Search Beyond MinHash](https://arxiv.org/abs/1612.07710)
is another important comparator: its Chosen Path construction exploits set
intersection, and its model is a logarithmic-word RAM. Theorem 1.1 assumes
fixed similarity thresholds. We do not substitute a threshold tending to
zero into that statement without accounting for its parameter dependence.
The direct deterministic proof below needs no such substitution. These
sources establish that the signature idea and the benefit of sparse-set
geometry are existing research, not discoveries of this project.

### General finite-parameter dictionary lemma

Let capacity be n, raw dimension d>=2, and every occupied key have exactly
t ones. Choose an integer 1<=ell<=t. Require distinct active supports to
have intersection less than ell. For a supplied integer radius
0<=r<=t-ell, define

```text
F = binom(t,ell), G = binom(r+ell,ell) <= F,
b = ceil(log2 d), B = ell*b.
```

Store every ell-subset of each support, in increasing coordinate order,
as an exact dictionary signature. Each signature belongs to at most one
active key, by the intersection requirement. A signature is encoded by
ell fields of b bits; leading zero fields cause no ambiguity because the
field count and total length B are fixed. No probabilistic fingerprint
or unchecked hash collision is used.

A simple binary trie of depth B suffices. Preallocate at most
1+n*F*B nodes, with constant-word child/slot metadata and a bounded free
list. Lookup, insertion and deletion take O(B) word operations; deletion
prunes the empty path. No randomized dictionary or amortized allocator
is needed. A compressed trie could save space, but is unnecessary for
the stated exclusion of the covering construction.

For a raw query q let w_q be its number of ones. If |w_q-t|>r, the exact
answer is empty. Otherwise put u=floor((w_q-t+r)/2), select the first
u+ell one-coordinates in a fixed order, enumerate their ell-subsets,
look them up and verify each new candidate in the original Hamming metric.
If a key x is within radius r, then

```text
|supp(x) intersect supp(q)| = (t+w_q-dist(x,q))/2
                          >= t-r >= ell.
```

A true key omits at most u of the query's one-coordinates, so the selected
prefix contains at least ell shared coordinates and a common signature
is looked up. Also 0<=u<=r and u+ell<=w_q, making the selection valid
and giving at most G lookups. Distinct active keys have distance at least
2*(t-ell+1)>2r, so at most one can actually be an inner neighbor.
After an original-distance verification succeeds, stopping is correct
for the complete operator. False signature candidates are verified and
discarded. Candidate IDs are marked and explicitly cleared on every
return; no increasing query counter is needed.

All independently chosen payloads remain in a separate current-value
array. Query signatures never encode a supplied target ID: the full raw
query determines its family of coordinate subsets. The support condition
is a declared restriction on keys, not a free learned encoder.

### Admission, replacement, deletion and total cost

For a proposed key, first check its weight. Then look up all F signatures.
A hit belonging to another slot proves overlap at least ell and rejects
the proposal. Conversely, overlap at least ell supplies a common
signature, so this test is complete. When replacing a slot, hits on that
same old slot are allowed during this precheck. Do not mutate the index
until the whole precheck succeeds. Then remove the old signatures and
insert the new ones. Rejection leaves keys, payloads and trie unchanged.
The preallocated node bound suffices because the old key is removed
before the new one is installed. Operations are sequential and atomic.

With v finite words per payload, the explicit bounds are

```text
space O(n*(ceil(d/w)+v+1) + n*F*B + d) words,
initial work O(n*d+n*F*B+n*v),
worst-case key admission/replacement/deletion O(d+F*B),
payload-only write O(v),
worst-case complete read O(d+G*(B+d)+s*v), s<=1.
```

An accompanying new payload costs its additional O(v) to load/write.
The d factors allow straightforward bit scans, without assuming a
population-count instruction. Enumerating a signature costs O(ell),
absorbed by B. Address, coordinate and fixed-grid aggregation guard bits
are included in the word width. Input loading, trie initialization,
free-list metadata and scratch marks are charged. Finite word widths
depend on capacity and payload precision, not on the operation horizon.
Every bound is deterministic, so arbitrary timing-adaptive operation
histories do not require copies, refreshes or failure probabilities.
Retention is exact until a slot is explicitly deleted or replaced;
inadmissible keys are rejected rather than silently weakening the promise.

The mixture version of the previous compact witness is also covered.
For each center keep the keys (x,0) and (x,1), with independent current
payloads and optional occupancy flags. Decode the base query prefix to
its at-most-one radius-r center. For that center check the two suffix
distances and form exactly the specified weighted sum/denominator. This
permits two genuine neighbors and independently overwritten values; at
(x,0) their weights are r+1 and r when r>=1. Raw center bits and its
signatures are stored once, and the two payloads are both charged. A
center-key replacement changes that structured pair together. This
extension does not assert support for arbitrary independently moved
paired keys that destroy the center-code restriction.

### The same random-shell regime admits these signatures

Take even t>=256, n=2^t, d=16*t^2, and the same independent uniform
weight-t supports as before. Put

```text
a=floor(log2 t), ell=ceil(4*t/a).
```

Then a>=8 and ell<=t/2, so all the earlier radii r=floor(t/c), c>=2,
are eligible. With probability at least 1-1/n every pair intersects in
fewer than ell positions, including exclusion of duplicate keys.
To prove this, for a fixed pair let I be its hypergeometric intersection.
The factorial moment identity gives

```text
Pr(I>=ell) <= E binom(I,ell)
  = binom(t,ell)^2/binom(d,ell)
  <= binom(t,ell)/(16*t)^ell <= (3/(16*ell))^ell.
```

The final inequality uses the usual binomial bound and e<3. Union bound
over fewer than n^2 pairs gives failure at most
n^2*(3/(16*ell))^ell. Here is a uniform elementary bound, not a numerical
fit: a<=16*2^(a/4) for a>=8, by induction, since
(a+1)/a<=9/8<2^(1/4). As t>=2^a and ell>=4t/a, this implies
16*ell/3 >= (4/3)*2^(3a/4). Raising to ell therefore gives

```text
(16*ell/3)^ell >= 2^(3t),
n^2*(3/(16*ell))^ell <= 2^(-t)=1/n.
```

For even t, combine this event with the previous spectral event. Both
the subset restriction and ||K||_2<=3 then hold with probability at least
1-1/n-delta_pair(t)-n^(-2). No extra query-independence assumption is
needed. The same paired construction has norm below 6, as before.

Finally ell<=5t/a, b<=2a+6 and B<=14t. For every r<=t-ell,

```text
G <= F <= (3*t/ell)^ell <= ((3/4)*a)^(5*t/a),
log2(F)/t <= (5/a)*log2((3/4)*a) -> 0.
```

Thus these inputs have n^(1+o(1))+O(n*v) finite-word space,
n^(1+o(1))+O(n*v) initialization, n^o(1)+O(v) worst-case admitted key
replacement, and n^o(1)+O(s*v) complete reads. The raw dimension is still
Theta(log^2 n), raw queries remain unrestricted binary vectors, and the
legal-key condition is checked within the stated update budget. In the
mixture extension charge two payloads per center and s<=2. No latent
rank restriction is imposed on those payloads.

These asymptotic bounds converge slowly. They do not establish a practical
speedup at the learned pilot's sizes, and the written finite-parameter
formulas, not n^o(1) alone, determine actual cost. No new training or
hardware benchmark is supplied.

### Consequence for the current research direction

The compact random-shell witness, even with exact mixtures, therefore
does not establish a new access/storage/retention frontier for our
covering block. On a high-probability subfamily of these same inputs a
standard signature dictionary has stronger asymptotic bounds in a
finite-word model, including key admission and overwrites. This removes
the earlier precision caveat as a defense of that particular example.
It does not dominate every bounded-spectral-norm dataset: bounded norm
allows a few nearly coincident keys whose overlaps violate the signature
restriction. Nor does admission here certify the spectral norm; the
dictionary simply does not need that norm bound.

The remaining task must address the full spectral-only class or a
different substantive workload, with comparisons that include exact set
filters and dictionaries. The existing variance identity, XOR batch
factorization and timing obstruction remain valid records. Their
ingredients or a favorable synthetic outer-density example do not
establish the original publication goal.
