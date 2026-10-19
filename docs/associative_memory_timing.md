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
storage and expected amortized geometry work O(T*log n). An all-prefix
confidence statement additionally needs the random-rebuild argument in
"Random rebuilds require a compound cost bound" below: the earlier
deterministically bounded rebuild argument alone does not justify it.
Initialization includes all copies and raw input loading. The compound
argument supplies this missing supporting lemma; it does not establish a
new conversion beyond prior art. Geometry and unfinished-query cleanup
are padded; exact payload aggregation after the geometric result is
complete costs O(s(q)*v).
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

## High overlap and full rank still admit classical product decoding

The previous iteration was **progress**: a finite-word signature baseline
removed the compact sparse examples from the novelty case. The next test
must genuinely violate that dictionary's restriction, rather than add
constant coordinates that an index can simply remove. The family below
does so, has a uniformly well-conditioned full-rank kernel, and nevertheless
admits cheap exact product decoding. This is another baseline/stress-test
consequence of established ideas, not a new architecture or paper claim.

### Product indexing is already a strong comparator

[Babenko--Lempitsky, The Inverted Multi-Index](https://www.robots.ox.ac.uk/~vilem/cvpr2012.pdf),
Section 3, describes ordered Cartesian-cell traversal by a multi-sequence
algorithm. Its discussion on PDF page 5 explicitly considers more than
two dimension groups; that page was visually checked. Thus increasing
the number of product factors is not itself a new indexing principle.
Their database application also uses approximate quantization. We use the
exact product-codebook part, not an approximation guarantee for arbitrary
raw database vectors.

[Lample et al., Large Memory Layers with Product Keys](https://arxiv.org/html/1907.05242v2),
Sections 3.1--3.2, already combines exact product-key access with a full
independent value table. Its two subqueries are vectors, not the two scalar
latents of the earlier failed prototype. The main research note retains
that distinction. The construction below is a fixed binary codebook with
many factors; it supplies no trained-model comparison or new optimization
method. Its specific kernel conditioning and radius-count calculations
are derived here as diagnostics, rather than attributed to those papers.

### A spectral family that defeats the unique-subset hypothesis

For any integer p>=2 let

```text
n=2^p, L=ceil(log2(4p)), d=2pL, theta=2^(-L).
```

Encode a logical binary word a in p blocks. Its j-th block is
(a_j repeated L times, (1-a_j) repeated L times). Write x(a) for the
resulting raw binary key. Every key has weight t=pL, no coordinate is
constant across the full codebook, and

```text
dist(x(a),x(b))=2L*dist(a,b),
K_ab=2^(-dist(x(a),x(b))/2)=theta^(dist(a,b)),
K = [[1,theta],[theta,1]] tensor ... tensor [[1,theta],[theta,1]].
```

The eigenvalues are (1+theta)^(p-h)*(1-theta)^h, h=0,...,p, with the
usual binomial multiplicities. Since p*theta<=1/4, Bernoulli's inequality
and a binomial/geometric-series bound give

```text
3/4 <= (1-theta)^p <= lambda_min(K)
    <= lambda_max(K) = (1+theta)^p <= 1/(1-p*theta) <= 4/3.
```

In particular K has rank n and condition number at most 16/9. Any active
subset inherits these eigenvalue bounds by the Rayleigh quotient; no
per-update spectral computation is necessary within this fixed codebook.
This K is the exponential kernel used in our promise, not automatically
the matrix of truncated weighted reads considered below.

Two logical words differing in one bit have support intersection t-L.
Any unique-subset dictionary from the preceding section must consequently
choose ell>t-L. Its radius condition r<=t-ell then forces r<L. Thus for
r>=L there is **no** ell satisfying that lemma's two hypotheses. The graph
joining distance-2L keys is the connected p-dimensional hypercube, with
n vertices. Merging all nearest pairs therefore does not give bounded-size
components. These observations rule out two simple explanations of this
example; they are not lower bounds against other dictionaries or clustering
algorithms. The product structure remains explicit.

### Complete raw-query decoding by a finite-word heap

For a completely arbitrary raw query q, count its ones u_j,v_j in the two
halves of block j. The two possible block distances are

```text
d_j(0)=L+u_j-v_j, d_j(1)=L-u_j+v_j,
m_j=min(d_j(0),d_j(1)), Delta_j=|d_j(0)-d_j(1)|=2*(L-m_j).
```

Choose a favored bit b_j realizing m_j, breaking ties consistently, and
put M=sum_j m_j. For every codeword a,

```text
dist(q,x(a)) = M + sum_{j:a_j!=b_j} Delta_j.
```

If M>r the exact answer is empty. Otherwise sort the p nonnegative
penalties and enumerate their subset sums up to r-M. A simple heap has
one root, the empty subset. A subset whose greatest chosen sorted index
is j has children obtained by adjoining i>j. Each subset has exactly one
parent, and a child's cost never decreases. Push a child only if its
distance is at most r; stop the child loop once the sorted penalty is
too large. No visited-set hash or growing timestamp is needed.

Every codeword within radius r appears exactly once: all its ancestors
also satisfy the budget. Conversely no pushed codeword lies outside the
radius. Let Z(q,r) be their number in the entire fixed codebook, including
currently unoccupied addresses. There are Z-1 accepted child insertions
and at most one unsuccessful boundary check per emitted node. A binary
heap therefore takes O(Z*log(Z+1)) word operations after sorting. Logical
addresses are p-bit words; adjoining a coordinate is one XOR. The sorting
cost O(p log p) is absorbed by d. Complete exact geometric work is

```text
O(d + Z(q,r)*log(Z(q,r)+1)).
```

For occupied addresses fetch their independent current payloads and add
the exact weight r+1-distance to the numerator/denominator. If S(q,r)
addresses are occupied, this costs an additional O(S(q,r)*v). Empty
addresses still count in Z; no active-output-sensitive bound is claimed
for them. All raw coordinates are processed to determine the block scores.
The input is not a supplied address and has no learned low-dimensional
encoder hidden outside the accounting.

### A uniform subpolynomial bound at the earlier radius scale

Assume M<=r. A block with Delta_j<L has m_j>L/2, so there are at most
floor(2r/L) such blocks. Flipping a bit in any other block costs at least
L, so a valid subset flips at most floor(r/L) of those other blocks.
Consequently, with binomial sums truncated at p if necessary,

```text
Z(q,r) <= min(n, 2^floor(2r/L) * sum_{i=0}^{floor(r/L)} binom(p,i)).
```

For every fixed c>1 and all r<=p/c this is n^o(1), uniformly over all raw
queries. For J=floor(r/L)>=1, the usual generating-function bound gives
sum_{i<=J} binom(p,i)<=(e*p/J)^J; here J<p/2. The function
x*log(e*p/x) increases for x<=p, so replacing J by p/(cL), and using
e<3, yields

```text
log2(Z(q,r))/p <= (2+log2(3cL))/(cL) -> 0.
```

When J=0 the same upper bound follows directly. This is an analytical
uniform bound, not a sampled-query estimate. The code has
d=Theta(log n*log log n), a particularly strong kernel bound, large
overlaps, and n independent addresses, yet an existing product-decoding
principle gives subpolynomial complete geometric reads. The finite count
Z and payload work, rather than this asymptotic notation alone, determine
practical costs.

### Admission, overwrites, raw storage and horizon

Use a declared capacity n for the fixed codebook. Keep an owner array
indexed by its p-bit address, an address per occupied user slot, the raw
key bits, and the independent value array. To admit a raw key, check that
every block has exactly one of the two allowed patterns and decode its
address. Reject a duplicate belonging to another slot before mutation.
On success unlink the old address, install the new owner and store the
raw input. Deletion and payload replacement are explicit operations.
Any active subset is automatically spectrally admissible by the preceding
principal-submatrix argument.

The complete finite-word bounds are

```text
space O(n*(ceil(d/w)+v+1)+d) words, including a capacity-n heap,
initial work O(n*d+n*v), including generation/loading of all supplied keys,
worst-case key admission/replacement/deletion O(d), plus O(v) for new values,
payload-only write O(v),
complete read O(d+Z*log(Z+1)+S*v).
```

Preallocate the heap and reset its length; do not clear n entries for
each query. All address, distance, ownership and fixed-grid accumulator
guard bits fit the declared word width and have no horizon-dependent
counter. Changing a logical key can move any slot to any unoccupied valid
codeword. Rejected keys leave state unchanged. Changing the codebook
patterns, training a semantic encoder or learning the block layout is
not included in this bound and is not assumed free. The guarantee is
restricted to this explicit codebook, not all spectrally admissible keys.

### Independent retention and the limits of the conclusion

Independent stored values are not merely formal parameters. For every
0<=r<=p, put h=max(0,r-L+1). From x(a) flip one physical bit in each of
h distinct blocks to obtain q_a. Then dist(q_a,x(a))=h, while for b!=a

```text
dist(q_a,x(b)) >= h+2*(L-1).
```

Every query within eta=min(r,L-2) further bit flips of q_a has x(a) as
its only possible radius-r codeword. If a is occupied it returns that
slot's exact current value; if it is absent the answer is empty. For
r<L-1 this follows from h=0; otherwise the target is at most r-1 away
and every other codeword is at least r+1 away after the extra flips.
In particular the complete read family contains an identity submatrix
on n addresses, so its dependence on arbitrary values has rank n too.
The pure-recall basin has radius O(log log n), not Theta(log n).

At the stored key itself, r>=2L includes neighboring logical addresses
and yields a genuine mixture of independent values. Do not call that
mixture the individually stored value. The explicit q_a basins are what
justify the independent-retention statement. All claims remain exact
through arbitrary admitted updates; no implicit forgetting is used.

This example separates the spectral condition from the previous
unique-subset restriction without relying on common-coordinate padding
or a rank-two value bottleneck. Nevertheless a classical product index
already has stronger asymptotic geometric/storage costs than our covering
bound here. Full rank, good conditioning and high overlap therefore do
not supply the missing novelty argument. The original goal remains
unachieved. A useful new result must address genuinely more general
interactions or prove a substantive new tradeoff against these product
methods, with all raw-input, value, admission and precision costs included.

## Random rebuilds require a compound cost bound

The preceding comparisons remove several tempting example-based novelty
claims. This iteration instead repairs a gap in the general cost argument.
The earlier concrete covering index has a deterministically bounded rebuild
and its all-prefix proof remains valid. The complete-range arithmetic
baseline uses a random number of capped construction attempts. Its expected
rebuild cost cannot be substituted for a deterministic reset bound in that
proof. Even a single geometric build count G with success probability 1/2
has mean 2 but Pr(G>=m)=2^(-(m-1)) for every positive integer m. Bounding
the number of rebuilds does not bound the work of those rebuilds.

The following supporting lemma supplies the missing argument. Its tools
are classical exponential supermartingales and Ville's inequality; see
[Howard--Ramdas--McAuliffe--Sekhon, Time-uniform Chernoff bounds via
nonnegative supermartingales](https://arxiv.org/abs/1808.03204v7), particularly
Section 2.3, Lemma 1. No new probability inequality or strongest-baseline
advantage is claimed.

### A compound epoch lemma, allowing dependent durations and costs

An epoch starts just after the initial installation or a reset. Let tau_j
be its number of operations through the next reset, a positive integer or
infinity. Write T_0=0 and T_j=sum_{i=1}^j tau_i. Let W_0 be initial randomized build
work and W_j the fallback and rebuild work at reset j. A reset following
the last requested operation is included in that operation's cost. Ordinary
work through Q completed operations is at most I+A*Q, where I is the
deterministic raw-input/reference initialization cost. Thus

```text
Cost_Q <= I+A*Q+W_0+sum_{j=1}^{N_Q} W_j,  N_Q=max{j:T_j<=Q}.
```

The boundary history consists of the logical state and the observations
available to the client, including earlier reset times and build durations.
It does not reveal the newly installed private seed. The proof needs two
different conditional statements. For fixed eta,theta>0 and K,K_0,beta>0:

```text
E[exp(eta*W_0)] <= K_0;
E[exp(-theta*tau_j) | preceding boundary history] <= beta;
E[exp(eta*W_j) | complete triggering history before the new build] <= K;
K*beta <= 1.
```

The third bound is uniform in the triggering operation, dataset, tau_j
and any old private state. It is stronger than an unconditional bound on
W_j and allows W_j and tau_j to be dependent. If an epoch never ends,
define its next process value as zero. For finite epochs, conditioning first
on the triggering history and then on the preceding boundary gives

```text
E[exp(eta*W_j-theta*tau_j) | preceding boundary] <= K*beta <= 1.
```

Using only the two marginal moments would be false. For a uniform choice
between (tau,W)=(1,2) and (2,0), take exp(eta)=2 and exp(-theta)=1/2.
The marginal product is (5/2)*(3/8)=15/16<1, but the actual joint moment
is 9/8>1. This is why the cost bound must survive conditioning on the
trigger, rather than merely hold at the preceding boundary.

Consequently

```text
M_j = exp(eta*(W_0+sum_{i=1}^j W_i)-theta*T_j)/K_0
```

is a nonnegative supermartingale with E[M_0]<=1 in the boundary filtration.
For 0<delta<1, apply Ville's inequality (or stop at its first crossing of 1/delta and a
deterministic epoch bound, then let that bound increase). With probability
at least 1-delta, simultaneously at every finite boundary j,

```text
W_0+sum_{i=1}^j W_i <= (theta*T_j+log(K_0/delta))/eta.
```

Since T_{N_Q}<=Q, the same event gives, simultaneously for all Q>=0,

```text
Cost_Q <= I+A*Q+(theta*Q+log(K_0/delta))/eta.
```

This includes Q chosen from previous answers and times, without replacing
an expected cost by a pathwise cost. No independence between epochs or
between their durations and rebuild costs is asserted. Correctness is still
for every seed; delta concerns total work. A long single-operation stall
remains possible.

### Explicit constants for capped Las Vegas construction attempts

Assume the established fresh-epoch bound
Pr(tau_j<=s | preceding boundary history)<=a*s for every integer s>=1,
with a>0. This is an epoch CDF condition, not a per-operation hazard bound
conditional on surviving earlier queries. The all-ordinary counterfactual
path and canonical padded outputs give this condition in the conversions
above. The installed seed must remain independent of the observed number
of rejected build slots. Fresh IID trials from the same snapshot, with
successful and failed attempts padded to the same slot length, supply that
property. Merely bounding their acceptance probabilities would not supply
it.

Let S>=1 be an integer upper bound on each padded construction attempt,
including rejected-attempt cleanup and fresh-randomness generation. Each
of k>=1 copies accepts a fresh attempt with probability at least 1/2.
Let Gamma_j be the sum of their trial counts, and let F>=0 be a deterministic
integer bound for fallback and the remaining reset overhead. Require

```text
W_0 <= S*Gamma_0,   W_j <= F+S*Gamma_j.
C=F+2*k*S,   gamma=1+1/(16*C),   z=1/(1+4*a).
```

Trial counts can depend on the current dataset. Uniform conditional success
probability >=1/2 gives domination by a geometric count G on {1,2,...}.
Writing v=gamma^S<2, its exact generating function gives
E[gamma^(S*G)]=v/(2-v). Sequential conditioning over copies therefore bounds
each reset's conditional exponential moment by
gamma^F*(v/(2-v))^k. To see that this is <=16/15 using explicit constants,
the geometric series and Bernoulli inequality give

```text
gamma^S <= 1/(1-S/(16*C));
v/(2-v) <= 1/(1-S/(8*C));
(1-S/(8*C))^k >= 1-k*S/(8*C);
gamma^F <= 1/(1-F/(16*C));
(1-F/(16*C))*(1-k*S/(8*C)) >= 15/16.
```

All denominators are positive since C=F+2*k*S. The same 16/15 bound covers
initial construction, which has no fallback F. It holds after conditioning
on the entire preceding trigger; hence it can be combined with tau_j even
when the operation that caused the reset also affects construction cost.

For the epoch transform put h=floor(1/(2*a)), possibly zero. With
q=Pr(tau_j<=h | boundary)<=1/2 and z^infinity=0,

```text
E[z^tau_j | boundary] <= q+(1-q)*z^(h+1)
                         <= 1/2+(1/2)*z^(h+1) < 2/3.
```

Indeed (1+4*a)^(h+1)>=1+4*a*(h+1)>3. This proof works for every a>0;
neither 1/a nor the expected epoch duration has to be an integer. The
compound lemma applies with eta=log(gamma), theta=log(1+4*a),
K=K_0=16/15 and beta=2/3, because K*beta=32/45<1. Thus with probability
at least 1-delta, for every prefix Q,

```text
Cost_Q <= I+A*Q
          + (Q*log(1+4*a)+log(16/(15*delta)))/log(1+1/(16*C))
       <= I+Q*(A+68*a*C)+17*C*log(16/(15*delta)).
```

The second line uses 1/log(1+1/(16*C))<=16*C+1<=17*C and
log(1+4*a)<=4*a. For a=0, all epochs have infinite duration almost surely;
the initial-build moment bound alone gives the same last line with a=0.
An additional output-dependent canonical cost, such as s(q)*v for exact
payload aggregation, is added at its realized sum over the prefix.

### Consequence, limitations, and the next gate

For the generic dynamic conversion take a=p+1/L and
A=D+k*(G+B), retaining its self-contained queued updates, maximum phase
space cap, and safe-build assumptions. For the generic static conversion
take a=(2*T/B)^k and A=k*B. S is the actual padded attempt bound; it is
not the expected total reset cost R. When all attempt cleanup is included
in S=O(R), C=O(F+k*R). Choosing the previous logarithmic number of copies
makes a*C small enough for the stated amortized bounds, now with an
all-prefix confidence term O((F+k*R)*log(1/delta)). Initialization is charged
by this random-cost bound too, rather than inserting its expectation as
a deterministic additive term. No storage slot, precision parameter or
algorithmic counter depends on Q or delta; the unbounded sums above are
analysis variables only. Bounded arenas are reused on each attempt.

This closes the confidence gap for the complete-range arithmetic baseline
under its already stated caps. It strengthens that comparator, so it does
not create a separation for our covering memory. It also does not compile
the optimal filters into a finite-word implementation or establish an
adaptive-key guarantee for them. The abstract timing observation is the
completed RAM-step count, not memory access leakage or hardware time.

The exact controls cover geometric generating functions, admissible epoch
CDFs including infinite epochs, dependent duration/cost laws, and finite
adaptive timing games. They are verification controls for the written
uniform proof, not a replacement for it. The general martingale framework,
padding, rerandomization and queue bounds already have prior art. A novel
generic compiler tradeoff would require a separate comparison with the
specific timing constructions of Bethea--Reiter and Braun--Jana--Boneh,
the deamortization theorem of Bernstein--Forster--Henzinger, and adaptive
search guarantees discussed earlier. No such separation has been proved.
The original goal and publication gate remain open. Further work should
resolve that generic comparison or the full nonproduct spectral workload,
rather than counting more favorable special examples as progress toward
a new architecture.

## Near-identity keys do not make all exact mixture queries cheap

The preceding turn repaired random-rebuild accounting; it did not establish
a new access advantage. A broader claim would be that the spectral promise
itself, perhaps with excellent conditioning and robust individual recall,
makes every rich memory query cheap. The following explicit reduction rules
out a fixed-power query saving for that unrestricted claim, conditional on
the OMv conjecture. It concerns the actual scalar weighted-mean output,
not the cost of materializing an entire neighbor list.

The source of the hardness is existing mathematics:
[Henzinger--Krinninger--Nanongkai--Saranurak, version 1, Definition 2.6,
Theorem 2.7 and Corollary 2.8](https://arxiv.org/html/1511.06773v1#S2.SS2)
give the polynomial-preprocessing OuMv/uMv consequences of the OMv
conjecture. The construction below transfers that established obstruction
to our spectral and retention promises. This is a diagnostic corollary;
neither a new conditional-complexity conjecture nor publication novelty is
claimed. The earlier bilinear barrier in the main research record did not
verify these additional geometric promises.

### Fixed binary keys, including every tag and input bit

Choose m>=2, n=m^2, integers a>=1 and b>=1, and put

```text
t=ceil(log2(n)),  L=t+b,  w=t*L,
d=4*a*m+2*w,  r_0=2*a*(m-1)+w,  r=r_0+a-1=d/2-a-1.
```

Here w denotes tag weight, not RAM word width. Let e_i be the i-th unit
vector of length m. Repetition z^[a] repeats each bit a times. Assign slot
(i,j), with 0<=i,j<m, the t-bit binary label i*m+j. Replace each label bit
z by (z^L,(1-z)^L); the resulting tag T_ij has length 2*w, weight w, and
distinct tags have distance at least 2*L. Define the explicit raw key

```text
x_ij = (e_i^[a], 0^(a*m), e_j^[a], 0^(a*m), T_ij).
```

All n keys have weight 2*a+w. The tag does not depend on the stored value.
All stored values V_ij can vary independently, including arbitrary Boolean
matrices. No low-rank factorization of V, hidden learned encoder, free
query address or compressed value representation is assumed.

For distinct slots, dist(x_ij,x_kl)>=2*L, and therefore their kernel entry
2^(-dist/2)<=2^(-L). Each off-diagonal row sum is at most
(n-1)*2^(-L)<2^(-b). The symmetric matrix K-I consequently has operator
norm <2^(-b), by the maximum absolute row-sum bound. Thus every eigenvalue
of K, and of every nonempty active principal submatrix, lies in

```text
[1-2^(-b), 1+2^(-b)],
condition number <= (1+2^(-b))/(1-2^(-b)).
```

This gives rank n and a kernel arbitrarily close to the identity. Taking
b=O(log n) still uses only O(log^2 n) tag bits per key. Neither estimating
the norm numerically nor charging a global spectral test at each value
overwrite is needed: every subset of this fixed codebook satisfies the
promise by construction.

### The exact operator contains Boolean rectangle queries

For arbitrary Boolean u,v in {0,1}^m, form and supply the full raw query

```text
q(u,v)=(u^[a], (1-u)^[a], v^[a], (1-v)^[a], 0^(2*w)).
```

Its distance to a key is exactly

```text
dist(q(u,v),x_ij)=2*a*m+2*a+w-2*a*(u_i+v_j).
```

When u_i=v_j=1 this is r_0; all other keys are at least r_0+2*a away.
For our fixed radius r, the weight max(0,r+1-dist) is therefore exactly a
on the selected rectangle and zero elsewhere. If both vectors are nonempty,
the complete scalar read is

```text
F(q(u,v),V) = (sum_{i,j} u_i*V_ij*v_j)/(|u|*|v|).
```

If either vector is empty the canonical output is empty and the corresponding
Boolean product is zero. For Boolean V and nonempty vectors, the output is
zero exactly when the Boolean product u^T*V*v is zero; otherwise it is at
least 1/n. This remains true even if the algorithm preaggregates values or
never enumerates any matching slots. Exact rational output uses numerator
and denominator of O(log n) bits. Alternatively, additive error strictly
below 1/(2*n), with comparison at 1/(2*n), suffices. A fixed binary grid
with O(log n) fractional bits provides that accuracy; no infinite-precision
real arithmetic is involved.

There is also a stronger occupancy version. Occupy exactly those codewords
with M_ij=1 and give every occupied slot value 1. For a nonempty selected
rectangle the read is 1 if u^T*M*v is true and empty otherwise. Every active
principal kernel retains the same near-identity bounds. Thus deciding the
canonical empty/nonempty flag already has the conditional obstruction,
without inverse-polynomial numerical accuracy or a required neighbor list.
If the interface encodes an empty read numerically as zero, the gap is
0 versus 1, so additive error <1/2 suffices. For an interface with an empty
sentinel, this statement requires preserving that sentinel; an approximation
allowed to misreport emptiness is a different guarantee. The matrix-dependent
occupancy is installed during polynomial preprocessing; the reduction needs
no online key edits at all.

Take a to be a fixed constant >=2 and b=O(log n). Then d=O(sqrt(n)+log^2 n),
and constructing each raw query costs O(d). Materializing the raw key table
costs O(n*d), which is polynomial and is included in preprocessing, together
with loading every value. Word width O(log n) suffices for polynomial
addresses, distances and exact accumulators. The algorithm may compress the
fixed keys or use any other polynomial-size auxiliary index; that freedom
does not avoid the reduction.

Suppose a reader for this full promise had polynomial preprocessing and
worst-case read cost O(d+n^(1-epsilon)) for any fixed epsilon>0. Applying it
to this family and forming the raw queries would give uMv time
O(m+m^(2-2*epsilon)+log^2 m), a fixed-power improvement on m^2. This
contradicts the cited OMv consequence. For an always-correct randomized
reader with a uniform expected total online cost, take m query pairs:
the same substitution gives expected total
O(m^2+m^(3-2*epsilon)+m*log^2 m). Capping the entire online computation at
three times its uniform expected bound gives a bounded-time OuMv algorithm
with error probability at most 1/3 by Markov's inequality, again contradicting
the conjecture. This also rules out that uniform amortized claim. The result
also allows expected-polynomial preprocessing: cap preprocessing and online
work separately at six times their respective uniform expected bounds and
union-bound the two failure events by 1/3. The online horizon here is the
fixed m query pairs, so no expected-stopping-time substitution is involved.
The conclusion is conditional and excludes a fixed-power improvement, not
logarithmic savings or all sublinear functions. All one-time setup belongs
to preprocessing; the expected online bound covers the complete sequence
after that setup.

### Robust individual retention coexists with that obstruction

For the pure prototype q(e_i,e_j), the target distance is r_0 and every
other slot is at least r_0+2*a away. Every additional raw-bit perturbation
of at most a-1 bits leaves the target at distance <=r and all others at
distance >=r+2. The read is exactly its current independent value, or empty
if the slot is inactive. This holds through arbitrary latest-value overwrites
and deletions. It proves full dependence on all n values separately from
the rank of K. It is a basin around a specified query prototype; it does
not claim q(e_i,e_j) equals the stored key.

A classical decoder can certify these basins in O(d): for each coordinate
of u and v, combine its a direct and a complemented repetitions, then
majority-decode. At most a-1 errors cannot change any decoded bit. Require
one-hot decoded vectors and explicitly verify the full raw distance to their
prototype is <=a-1. The check includes tag bits. If it passes, ordinary
lookup returns the exact current value; otherwise the decoder makes no
claim and a complete reader can fall back. Thus fast robust individual
recall is compatible with the conditional difficulty of rich mixture reads.
Taking a=log^O(1) n permits a growing noise basin with only polylogarithmic
extra raw-input work; the same fixed-power hardness consequence survives.
This encoding/decoder is classical redundancy and lookup, not a new memory
architecture advantage.

### The radius and accuracy boundaries are essential

Here r=d/2-a-1=Theta(a*sqrt(n)+log^2 n). This is outside the earlier
r=floor(log2(n)/c) regime. The reduction does not refute that covering
guarantee, give a lower bound in its small-radius regime, or establish that
its n^(1/c) exponent is optimal. It prevents dropping the radius restriction
and replacing the missing proof by an appeal to conditioning or retention.

Accuracy and occupancy must be distinguished. The all-occupied payload
version needs inverse-polynomial accuracy; the occupancy version needs the
canonical empty/nonempty semantics, or the constant numerical gap described
above. For all-occupied, unperturbed rectangle prototypes and values in [-1,1],
sampling independent uniform selected coordinates I,J and averaging V_IJ
gives the normalized answer to constant additive accuracy with
O(epsilon^(-2)*log(1/delta)) samples, after
O(d) validation/sampler construction. Values reside in an ordinary table
and overwrite in constant time per scalar. This is the existing bounded-
variable sampling baseline already recorded in the main research note;
the payload reduction does not rule it out. With an arbitrary occupied
subset, the selected rectangle may have only one occupied slot. Uniform
rectangle sampling then has no comparable constant-sample guarantee for
emptiness. The all-occupied baseline is not a full-raw-query or
noise-basin sampling theorem, nor a simultaneous infinite-horizon error
guarantee.

The exact checker verifies the raw geometry, rational kernel row bounds,
complete weighted means, Boolean decisions after finite-grid rounding,
occupied-subset empty decisions and principal kernel row bounds,
independent values and exhaustive small noise basins. It cannot verify the
OMv conjecture. The next positive theorem still needs a quantified advantage
over applicable methods in an explicit radius/accuracy regime. This barrier
is evidence against a broad shortcut, not completion of the original goal
or a contribution manuscript.

## Exact mixture aggregation from low-intersection keys

The preceding conditional barrier does not rule out a useful structural
restriction on the keys. This iteration gives a candidate deterministic
tradeoff for that route. It extends the earlier unique-subset dictionary
from isolated-neighbor queries to the complete scalar/vector aggregate,
including queries matching many independently valued keys. The uniform
algorithm and accounting below are proved; **novelty and an advantage over
the strongest applicable aggregation methods remain unestablished**.

The initial distance-profile audit did not close the full spectral class.
Ahle's Hamming filters bound complete candidate work through intersections
of Hamming balls, while the data-dependent ALRW/Wei recursion also changes
centers and radii. The root spectral bound cannot simply be substituted for
each local distance profile. This is not evidence that those methods fail.
The restriction here is explicitly stronger than a bounded spectral norm.

### Model and finite parameters

There are at most n active distinct binary keys x_i in {0,1}^d. Each has
exact weight w, and any two active keys share fewer than k one-coordinates,
where 1<=k<=w. Values are arbitrary independently overwritable signed
fixed-grid vectors of v words, with S input bits per scalar coordinate.
For a fixed integer radius r>=0 the reader returns the exact weighted mean
with weight max(0,r+1-dist(z,x_i)), or the canonical empty answer when the
weight sum is zero. The query z is supplied as all d raw bits. Its support
and popcount are computed and charged. Changing the encoder or the key
class is not free.

Define

```text
F = sum_{j=0}^k binom(w,j),
Q_T = sum_{j=0}^k binom(T,j),  T=weight(z).
```

Take binomial coefficients with j greater than the upper argument as zero.
A necessary condition for any neighbor is |T-w|<=r; otherwise return empty
after reading the raw query. In particular only T<=w+r needs the following
algorithm. It works for every r, with the displayed T-dependent cost; the
subpolynomial consequence later assumes r=O(w).

### The maintained subsets and the exact identity

For each coordinate subset A of size less than k that occurs in an active
key, maintain its support count c_A and vector sum H_A of current values
over keys containing A. The empty subset holds the global count and sum.
Maintain a separate dictionary for size-k subsets, mapping each to its
unique containing slot. Uniqueness follows from the intersection hypothesis.
Counts, not whether a signed vector sum vanishes, determine whether an
entry is live. Do not discard a zero sum supported by active slots.

For this query put

```text
g(s)=max(0,r+1-w-T+2*s),  0<=s<=w,
b_j=sum_{h=0}^j (-1)^(j-h)*binom(j,h)*g(h),  0<=j<k,
P(s)=sum_{j=0}^{k-1} b_j*binom(s,j).
```

Newton's finite-difference formula gives P(s)=g(s) for integer
0<=s<k. The polynomial is allowed to have negative coefficients; signed
exact arithmetic is essential. For s_i=|support(z) intersect support(x_i)|,
the denominator and numerator are exactly

```text
den = sum_{j<k} b_j*sum_{A subset support(z), |A|=j} c_A
      + sum_{i:s_i>=k} (g(s_i)-P(s_i)),
num = sum_{j<k} b_j*sum_{A subset support(z), |A|=j} H_A
      + sum_{i:s_i>=k} (g(s_i)-P(s_i))*V_i.
```

Indeed, a fixed key contributes binom(s_i,j) to its j-th subset sum.
Summing first over subsets therefore gives P(s_i) for that key, which
already equals the required g(s_i) when s_i<k. The correction restores
the required coefficient for every other key. This is an identity over
the current values, not a stochastic approximation or a rank restriction.
The identity works for any integer score table g(0),...,g(w). The Hamming
hat additionally gives the population gate T<=w+r; no subpolynomial claim
for arbitrary dense softmax inputs follows from that identity alone. Any
other score table must also pay its evaluation and coefficient bit costs.

Enumerate all subsets of the query support of sizes 0 through k. Below
order k fetch the maintained aggregates, treating absent entries as zero.
At order k fetch unique owners and deduplicate them. Every key with
s_i>=k occurs because some size-k subset of the shared support was queried;
at most binom(T,k) candidate occurrences are generated. Read each distinct
candidate's raw key to compute s_i and add its correction once. An n-slot
mark array with an explicitly cleared touched list avoids horizon-dependent
timestamps. Return num/den when den>0 and empty otherwise. Although
intermediate counts can be signed, the identity proves the final den is
nonnegative. No all-neighbor list is materialized.

For a pure low-order query, exponentially many keys can contribute through
the cached sums with no candidate visits. For a high-overlap query the
candidate correction preserves dependence on individual values. Both parts
are needed: the lower-order feature state alone need not have rank n.

### Admission, updates, storage, and precision

Exact tuple dictionaries can use the earlier preallocated binary trie,
with one root per subset order. A j-tuple of sorted coordinate IDs has
j*ceil(log2(d+1)) bits and is looked up in O(k*log(d+1)) worst-case steps.
There are at most n*F live subset entries. Explicit raw keys, values, trie
nodes, support counts, unique owners, per-slot subset handles, free lists,
and query scratch space all fit in

```text
O(n*(ceil(d/W)+v+1) + n*F*(k*log(d+1)+v+1) + d+w*k) words.
```

Here W is RAM word width. Preallocate arenas for that capacity and recycle
freed entries. Cached handles require O(n*F) additional words already in
the bound. A scalar/vector overwrite changes every lower-order subset sum
for that key by its exact old/new difference, at worst O(F*(v+1)) steps
using those handles; no other independent value is altered. Deletion removes
the old contribution, decrements counts, frees unsupported entries and
unlinks size-k owners.

For an inserted or replacement raw key, first check its exact weight and
all its size-k subset owners. Ignore only the replacing slot's own entries.
Another owner exists iff another key has intersection at least k, so this
is a complete admission test. Reject before changing state. On acceptance,
remove the old slot then install all new subsets and the new value. Deletion
and admitted key replacement cost

```text
O(d + F*(k*log(d+1)+v+1)).
```

The d term conservatively charges scanning every raw input bit. Counts and
the active-key table prevent stale sums after deletes, reuse or replacement.
All operations are deterministic; there is no probabilistic horizon or
forgetting assumption. Initial loading costs O(n) times the admitted-edit
bound, plus zeroing the preallocated storage and storing all raw keys.

Precompute binom(s,j) for 0<=s<=w, 0<=j<k and the query's coefficients
and P(s) table in O(w*k+k^2) steps. The conservative complete read bound is

```text
O(d+w*k+k^2 + Q_T*(k*log(d+1)+d+v+1)).
```

It includes all subset lookups, raw candidate distances, vector arithmetic,
deduplication and scratch cleanup. Fetching one aggregate reads v words;
this cost is not hidden behind a constant-time vector operation.

Put G=r+w+T+2. Then |b_j|<=2^j*G and
|P(s)|<=G*2^k*sum_{j<k}binom(w,j). Summing at most Q_T dictionary
terms each of count <=n bounds all intermediate integer magnitudes using
S+O(log(n+1)+log G+k+log F+log Q_T) bits per scalar. Addresses must
also cover the complete raw-key/payload arrays, arenas and scratch space;
O(log(n+1)+log(d+1)+log(v+1)+log(F+1)+log(k+1)) bits suffice.
Choose W to cover these bounds and charge multiword payloads through v.
The final denominator is at most
n*(r+1); storing the exact shared denominator and signed numerators does
not require floating point. Input precision and these guard bits are fixed
before execution, using T_max=min(d,w+r) in Q_T, and do not grow with the
number of overwrites.

### A uniform asymptotic regime and an explicit family

If w=Theta(log n), d=log^O(1) n, k=O(w/log w), and r=O(w), then

```text
F,Q_T = exp(O(w*log(log(w)+2)/log(w))) = n^o(1).
```

Use the elementary binomial-sum bound with k<=w/2; bounded small parameter
cases are handled directly. Storage and preprocessing are
n^(1+o(1))*(v+1), and every admitted key edit, value overwrite and complete
read costs n^o(1)*(v+1) worst-case. All d raw query coordinates are charged.
This is an asymptotic statement with potentially large subset tables, not
a practical latency or training result. For bit-operation costs the same
subpolynomial factors apply at polylogarithmic scalar precision; larger
precision must also pay the corresponding exact-integer arithmetic costs.
The saved finite parameter examples include cases where Q_T exceeds n.

Independent retention holds throughout the generic admitted class, not
only for the codebook below. Its minimum pair distance is at least
2*(w-k+1). With r=w and k<=w/2, every raw-bit perturbation of at most
w-2*k+1 bits around any active key has that key as its sole neighbor.
The target is inside radius w and every other key is at least w+1 away.
The exact read is therefore its current arbitrary value. This invariant
survives every admitted key edit and every value overwrite. Meanwhile the
all-zero query at r=w has all active keys at weight one and returns their
full mean through the empty-subset sum, with no owner candidates at all.
The same minimum distance also gives, in every admitted state,
||K-I||_2<=(n-1)*2^(-(w-k+1)). Thus choosing w>=4*log2 n and k=o(log n)
preserves a near-identity kernel under the generic admission rule itself;
the algorithm does not need to compute that norm to accept an edit.

The class is nonempty and has substantial query binding. Let p be prime,
take every polynomial f over F_p with degree <k, and encode it by the
one-hot evaluation positions (a,f(a)), a in F_p. There are n=p^k keys,
w=p, d=p^2, and distinct keys intersect in at most k-1 coordinates. For
k=floor(p/(4*ceil(log2 p))) and sufficiently large p, log n=Theta(p),
so the regime above applies. This is a classical Reed--Solomon codebook;
admission can use either its polynomial constraints or the generic subset
test. It is not a newly discovered prime-specific construction.

Its spectral norm can also be calculated exactly. Translation by a codeword
makes every kernel row sum identical. Expanding 2^(number of zero
evaluations) over subsets and using independent evaluation constraints gives

```text
||K||_2 = 1+2^(-p)*sum_{j=0}^{k-1} binom(p,j)*(p^(k-j)-1)
         < 1+3*n*2^(-p) <= 1+3*n^(-3).
```

The last inequality uses p>=4*log2 n in the selected regime. The symmetric
off-diagonal row-sum bound gives the matching lower eigenvalue bound
1-3*n*2^(-p); every active subset inherits these bounds. No large Gram
matrix needs to be formed.

For r=p, the stored key itself returns its independent value, and every
further raw-bit perturbation of at most p-2*k+1 bits still has that key
as its only possible inner key. This follows from the minimum binary
distance 2*(p-k+1), and holds for arbitrary raw perturbations, not only
valid one-hot symbol queries. More generally the pure basin has radius
min(r,2*(p-k+1)-r-1) whenever that number is nonnegative.

Nontrivial blends are also explicit. For an odd p with p+1>4*(k-1), take
two codewords f and f+X. They agree at zero. Construct a one-hot query by
using each word on half of the nonzero evaluation positions and their
common value at zero. At r=p both have weight two; any other degree-<k
polynomial agrees at at most 2*(k-1) queried positions and is outside.
The exact read is their current two-value mean when both are active.
For any fixed ell dividing p-1, the analogous balanced query over
f+a*X for ell distinct a has exactly those ell neighbors if
1+(p-1)/ell>ell*(k-1), at radius 2*(p-1)*(1-1/ell).
Thus the architecture is not globally restricted to two values or two
latent scalars. Pure reads on all n stored keys give an identity submatrix
on their independent payloads.

### Prior-art boundary and current decision

Subset signatures and admission are already present in the earlier exact
dictionary comparator, based on
[Arasu--Ganti--Kaushik, Sections 3 and 7](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/vldb06.pdf).
The new step relative to that comparator is exact aggregation of low-score
keys and correction of high-score owners; it must not be described as a
new subset index. Newton interpolation and inclusion--exclusion are classical.
Partial aggregate materialization and compensating retrieval are also prior
art: [Drzadzewski--Tompa, Sections 5--6](https://link.springer.com/article/10.1007/s10115-015-0871-2)
was read through Algorithms 1--2 and the query steps in Sections 6.1--6.2
(PDF pages 23--25). It materializes large conjunction cells, rewrites tag
queries by inclusion--exclusion and retrieves members of smaller cells or
compensating differences. These are directly relevant components. Its
rewriting/member-retrieval costs must be reconciled with the threshold-score
and dynamic bounds here; no separation from that framework is established.
A matching specialization of that or other database/set-aggregation work
would classify this theorem as
a known consequence.

Reed--Solomon list decoding is a stronger baseline than a scan or a product
lookup for structured codeword queries; see
[Guruswami--Sudan](https://people.csail.mit.edu/madhu/papers/1998/gs.pdf).
List decoding alone is not an all-raw-bit aggregate guarantee, especially
when a query has many matching values. Conversely, our identity is not a
faster list decoder: it returns the specified sum and count, without
reporting all participating codewords. Full soft-decision/list-recovery
and range-aggregation comparisons remain necessary.

The candidate gives a complete constructive mechanism for efficient reads,
writes and independent retention in a more expressive regime than our own
earlier unique-neighbor theorem. It does not solve the full spectral-only
class, establish a learned encoder, or prove superiority over the strongest
published aggregation method. Preserve it for the precise priority audit;
do not start neural training or declare the original goal complete.

## A stronger coded comparator: soft interpolation and affine sums

The preceding Reed--Solomon witness needs a stronger comparator than either
a scan or our own full subset table. Classical soft interpolation, combined
with partial evaluation aggregates, gives the following deterministic bound.
This is a derived baseline, not an established new contribution. It applies
to a fixed polynomial codebook; the generic low-intersection admission class
in the preceding theorem is strictly broader. Neither result establishes a
learned-architecture advantage or completes the research goal.

### Model and complete cost statement

Fix integers k>=2 and 1<=h<k. Let p>=k be an odd prime,
n=p^k, and let the allowed keys be all degree-<k polynomials over F_p,
one-hot encoded by (a,f(a)) in d=p^2 raw bits. Keep any active subset of
these n possible keys, each with its own arbitrary signed fixed-grid
v-word value. A read accepts any raw d-bit vector z and returns the exact
weighted numerator and count, or their rational mean, at fixed radius p:

```text
weight(f,z) = max(0,p+1-dist(x_f,z)).
```

Empty reads have a distinguished empty result. All value coordinates,
occupancy, overwrites and deletions count. An inserted/replacement key is
supplied as raw bits: scan it, check one-hot rows, interpolate from k rows,
verify all p evaluations, then derive its coefficient address internally.
Reject malformed, noncodeword or duplicate keys before mutation. Admission
outside this fixed codebook is not supported. Coefficients are a lossless
description of a validated codeword, not a supplied latent encoding or a
compression of independent payloads.

In an exact finite-word model with adequate guard/address bits, the bounds
below count scalar word operations and field operations; field division,
larger payload arithmetic and bit-level costs must also be charged. With
fixed k,h, a precomputed field-inverse table and polylogarithmic precision,
the corresponding bit costs differ by polylogarithmic factors:

```text
state:                  O_k,h((n+C_h)*(v+1) + p^2) words,
complete raw read:      O_k,h(p^2 + p^(k-h-1)*(v+1)),
value overwrite/delete: O_k,h(p^h*(v+1)),
raw admitted key edit:  O_k,h(p^2 + p^h*(v+1)).
```

Here C_h=sum_{j<=h}binom(p,j)*p^j. For fixed k,h with 2h<=k, the
state simplifies to O_k,h(n*(v+1)), since C_h=O_h(p^(2h)). Without
that condition this fixed-k simplification is not asserted. Here n is
codebook capacity, even if few keys are active. Dense tables of that
capacity are charged. Initialize their zero state and metadata in
O_k,h((n+C_h)*(v+1)+p^2) steps. Loading n arbitrary raw keys costs
another O_k,h(n*(p^2+p^h*(v+1))) steps; there is no near-linear initial
loading claim. Updates are worst-case, with no random rebuild, horizon
counter or forgetting error. Optional external-slot ownership and duplicate
checks use two n-entry inverse arrays and are included.

For k=5,h=2, this is linear state and complete reads, raw admitted edits,
overwrites and deletes in O(n^(2/5)*(v+1)) scalar operations, up to the
stated bit/field costs. Full raw loading costs O(n^(7/5)*(v+1)). For
k=6,h=3, scalar reads instead cost O(n^(1/3)), while value writes cost
O(n^(1/2)*(v+1)); raw admission adds n^(1/3). These are different
read/write tradeoffs. Constants depend strongly on k and are not uniform
when k grows. In particular these bounds do not replace the preceding
growing-k subpolynomial theorem by treating a fixed-k constant as uniform.
An explicit growing-k comparison is given below; priority remains open.

### Store short evaluation sums, enumerate the missing fibers

For every consistent j-point evaluation constraint, 0<=j<=h, store its
active support count and exact current value sum. There are

```text
C_h = sum_{j=0}^h binom(p,j)*p^j = O_h(p^(2h))
```

possible cells. A codeword contributes to
A_h=sum_{j=0}^h binom(p,j)=O_h(p^h) cells. This proves the edit bounds;
recompute its p symbols in O_k(p) time if needed. Keep support counts
separately from signed sums: cancellation does not remove a live cell.
No size-k owner table is required.

Cells need no hash tables. For ordered distinct row IDs a_1<...<a_j,
the colexicographic rank sum_{i=1}^j binom(a_i,i) lies in
[0,binom(p,j)). Combine it with the base-p rank of their specified
evaluation values and an order offset. This indexes exactly C_h cells.
Precompute the O_h(p) binomial entries; all other metadata and arrays
fit in the state bound. A polynomial's coefficient tuple likewise has a
direct base-p address in [0,n). Deduplication uses fixed address marks and
a cleared touched list, rather than counters that grow with time.

Let T=weight(z), and let s_f count the selected pairs (a,f(a)) in z.
Rows may contain zero, one or several selected symbols. Exactly

```text
dist(x_f,z) = p+T-2*s_f,       g_T(s)=max(0,1-T+2*s).
```

If T>2p the result is empty. Otherwise first consider T<=6(k-1).
Use the preceding identity with b_j=Delta^j g_T(0), j<k, and
P(s)=sum_{j<k}b_j*binom(s,j). Enumerate the query's j-subsets with
nonzero coefficient. Inconsistent repeated evaluation rows contribute
zero. At j<=h fetch the cached cell. At h<j<k, solve its j independent
evaluation equations in k coefficients. Their solution is an affine fiber
of exactly p^(k-j) polynomials. Enumerate it, fetch active values by direct
address and sum them. Each missing cell therefore costs at most
O_k(p^(k-h-1)*(v+1)); no free aggregate oracle is assumed.

Every size-k subset of distinct rows specifies exactly one polynomial by
interpolation. Find these addresses, deduplicate them, and for each active
owner add (g_T(s_f)-P(s_f)) times its value and the matching count
correction. This gives exactly the preceding Newton sum, including arbitrary
large mixtures. There are at most sum_{j<=k}binom(6(k-1),j) subsets, a
constant only for fixed k. Raw input scanning costs p^2; owner distance
verification costs O_k(p) per candidate. This proves the low-T read bound.
When T>=2k-1 all b_j with j<k vanish, so this branch needs only the owner
interpolations, not the missing-fiber value scans.

### The other branch has at most two algebraic candidates

For T>6(k-1), build a nonzero Q(X,Y) of Y-degree at most two, vanishing
at every selected point, with (1,k-1)-weighted degree <s_0=ceil(T/2).
Here is an explicit quadratic-time construction, so the bound does not
hide a large linear-system solve or a general polynomial factorizer.

Start g_0=1,g_1=Y,g_2=Y^2. Order monomials by
(X-degree+(k-1)*Y-degree,Y-degree). At each selected (a,b), evaluate all
three discrepancies. If all vanish, continue. Otherwise choose the
nonzero-discrepancy pivot with least leading monomial. Using its old value,
subtract its appropriate field multiple from each other polynomial with
nonzero discrepancy, then multiply the pivot by X-a. Prior zeros persist
and the current point is now a zero of every polynomial.

Inductively the three leading Y-degrees remain 0,1,2: subtracting a
strictly smaller pivot cannot cancel another leading term. Multiplying the
pivot increases its leading X-degree by one. Thus if these leading
X-degrees are e_0,e_1,e_2, sum e_j<=T. If all three weighted degrees
were >=s_0, their sum would give

```text
sum e_j >= 3*s_0-3*(k-1) > T,
```

a contradiction. Choose a polynomial of least weighted degree as Q.
Its coefficient degrees are O_k(T); incremental evaluation and updates
take O_k(T^2) field operations and O_k(T) scratch words, even with
repeated a and different b. Selected identical pairs cannot repeat in
a raw bit vector.

For every near codeword, s_f>=s_0 and the univariate Q(X,f(X)) has
degree <s_0 and at least s_0 distinct evaluation roots. It is therefore
zero identically, so Y-f(X) divides Q. There can be at most two distinct
polynomial roots. This is the classical interpolation/score mechanism of
[Koetter--Vardy, Section III, Theorem 3](https://www.site.uottawa.ca/~zhcheng/RS_soft_decision.pdf),
specialized here to multiplicity-one raw supports. The primary PDF's score,
cost and factor argument were checked on pages 4--5; page 5 was visually
inspected. The explicit three-polynomial calculation above supplies our
particular threshold and cost rather than attributing the entire dynamic
aggregate theorem to that paper.

Root extraction is elementary over odd prime fields. Write Q=A*Y^2+B*Y+C.
If A=0 use exact division -C/B when B!=0; a nonzero constant-in-Y
polynomial has no roots. Otherwise form D=B^2-4*A*C. A polynomial
root implies D=(2*A*f+B)^2. If D=0 use its zero square root; otherwise
its degree must be even, and its leading coefficient must be a square.
Find a leading square root by scanning F_p, then recover the remaining
coefficients from highest to lowest by dividing by twice that leading
coefficient. Verify the entire square. At most the two exact quotients
(-B+H)/(2*A) and (-B-H)/(2*A) can be polynomial roots; retain only
zero remainders, degree <k and verified substitutions. This costs
O_k(T^2+p) field operations. Common X factors and repeated roots are
allowed. General algebraic root reconstruction is established prior art;
see [Roth--Ruckenstein, Section V, Proposition 5.2 and Section VI](https://www.site.uottawa.ca/~zhcheng/efficient.pdf).
No general factorization guarantee is needed for this quadratic case.

Finally check each candidate's actual raw distance, occupancy and current
value, and form its exact weighted numerator/count. Since T<=2p this
whole branch costs O_k(p^2+v+1). Combined with the low branch it proves
the complete read bound. Payload arithmetic is over signed integers, not
modulo p; the finite field is used only for key geometry.

### Precision, retention, and elementary lower bounds

For fixed k,h, low-branch subset counts and Newton coefficients depend
only on k. Each affine cell has at most n members. All integer aggregate
and correction magnitudes fit S+O_k(log(n+1)) bits per scalar if the
input scalar has S bits; the final count is at most n*(p+1). Address
width covers the full n-entry value/owner arrays and C_h cells, including
v. Preselect this precision; it never grows with overwrite count. The
output keeps a shared exact denominator rather than assuming a constant
cost floating division. These observations justify the finite-word model
and its stated polylogarithmic bit-cost qualification.

Distinct keys still have distance at least 2(p-k+1). For p>2k-2, every
raw perturbation of at most p-2k+1 bits around an active codeword returns
its current independent value at radius p. The whole pure-read family
contains an n-by-n identity. There is no global two-latent rank bound:
low-T queries may combine many independent values through jointly
constrained sums. For large enough p, high-T queries can also have two
independently changeable neighbors, for example the union of two disjoint
constant-codeword supports. The two-candidate bound is local to that
query branch, not a restriction on stored capacity.

An ordinary information bound requires n*S bits of state for n independent
S-bit scalar values: their pure reads distinguish every assignment. It
does not prove our write cost optimal. A separate exact-input bound matches
scalar read scaling up to word/precision factors. Define t(y)=1 except
t(p-1)=-(p-1), and store the legal scalar witness

```text
V_f = sum_{a in F_p} t(f(a)).
```

All n codewords are active. Evaluations at any two distinct rows are
independent and uniform for k>=2. The zero query returns zero. The
singleton raw query (a,y) selects exactly f(a)=y, all with weight two,
and returns t(y), which is never zero. Thus changing any one input bit
at the zero query changes the exact answer. Every deterministic exact
reader must inspect every packed input word there, requiring
Omega(p^2/W) word probes for W-bit words. The same conclusion holds for
zero-error randomized readers by fixing their coins on this finite input
family. The witness uses only O(log p) value bits. For k=5,h=2 and
k=6,h=3 scalar precision O(polylog n), our respective O(p^2) upper
bounds are therefore optimal up to word/precision factors for this
particular raw-input model. No v-factor optimality or dynamic write
lower bound follows. Exact raw admission has the analogous input bound:
flipping any bit of a valid one-hot key makes a row invalid.

These are elementary input and information bounds, not new lower-bound
techniques. They explain why cost should be compared to raw dimension
and retained value capacity, rather than a redundant dense matrix alone.

### The same growing-k regime as the generic candidate

The exact cell count also permits a stronger codebook comparator in the
preceding theorem's regime. Take k=floor(p/(4*ceil(log2 p))), h=k-1,
and sufficiently large p, so log n=Theta(p). In this paragraph expand
all dependence on k. Define

```text
F = sum_{j<k}binom(p,j),
B = sum_{j<=k}binom(6(k-1),j) <= 2^(6(k-1)).
```

Then C_h<=p^(k-1)*F=n*F/p. The usual binomial-sum estimate gives
F=exp(O(k*log(p/k+2)))=n^o(1). Low-branch reads have at most B subsets
and no missing affine fibers. Each subset lookup/address calculation
costs O(k); each size-k owner interpolation costs O(k^3), its raw
distance evaluation O(p*k), and its payload fetch O(v+1). Computing
Newton coefficients and the P(s) table costs O(k^2+p*k). Hence a
uniform conservative bound is

```text
read:  O(p^2+p*k+k^2 + B*(k^3+p*k+v+1)),
edit:  O(p^2+k^3+p*k + F*(k+v+1)),
state: O((n+n*F/p)*(v+1) + n*k + p*k + p^2) words.
```

The high branch still costs O(p^2+v+1): its coefficient degrees and
quadratic-root checks are O(p), and k<=p. Whole-state allocation plus
n raw insertions is n^(1+o(1))*(v+1). Guard widths follow the preceding
generic Newton bound with this F and B, giving O(S+polylog n) bits;
field and payload bit arithmetic remain charged. Thus reads and edits
are n^o(1)*(v+1), while state and preprocessing are
n^(1+o(1))*(v+1). More specifically the read subset factor is
2^O(k)=n^O(1/log(log n)), compared with the full-subset count over
O(p) selected points, whose elementary bound contains an additional
log(log(log n)) in that exponent. These are asymptotic bounds with
potentially very large finite constants, not measured speedups.

This comparison uses the same raw dimension, codebook, radius, full
independent values, spectral witness and noise retention as the generic
candidate's polynomial example. Classical decoding plus cached evaluation
sums already suffices there. It still does not admit arbitrary raw keys
satisfying only the low-intersection condition. A claim for the generic
class must therefore rest on that broader dynamic guarantee and a real
aggregation-priority gap, not solely on the polynomial example.

### The generic algebra is an intersection transform

A closer algebraic predecessor is
[Björklund--Husfeldt--Kaski--Koivisto, Theorem 1 and Sections 2.1--2.3](https://arxiv.org/abs/0809.2489).
Its full intersection transform returns the weighted sum of sets with
each possible query-intersection size. The ranked up/down zeta moments
are exactly our y_j=sum_X binom(|X intersect z|,j)*V_X. Their inverse
Pascal matrix recovers the intersection-size sums from these moments;
contracting them against an arbitrary score g gives the equivalent
finite-difference coefficients used here. The primary proof and Lemmas
3--5 were read on PDF pages 4--7. Its Theorem 1 charges complete
down-closures of the data and query families and gives an arithmetic
circuit; that statement alone is not our online compressed bound.

For low-intersection data, every occurring subset of size at least k has
one owner. Consequently its high-order moments equal
sum_{owners i:s_i>=k}binom(s_i,j)*V_i. Substituting this into the full
transform and grouping each owner's tail gives exactly
(g(s_i)-P(s_i))*V_i. The candidate is therefore a compressed application
of the published intersection-transform algebra plus unique-owner subset
indexing. The same identity holds for counts and signed vectors. The
additional exact controls exhaust three non-polynomial families and three
radii, comparing the full transform, direct intersection classes and this
compressed tail.

The potentially distinct statement is a compressed dynamic cost bound
under low-intersection admission, not a new transform or interpolation
identity. We have not identified a published theorem with exactly that
full bound, but that absence is not evidence that this elementary
specialization merits a paper. A substantive improvement, application or
matching obstruction is still needed to justify continuing this candidate
as a contribution rather than a baseline. Keep the partial-materialization
comparison alongside this much closer algebraic reference.

### Verification and decision

The new `affine_coded_baseline` report block checks exact full numerators
and denominators against independent raw Hamming scans. It exercises
polynomial division/square roots (including repeated roots, common factors,
nonsquares and no-root cases), incremental interpolation invariants,
affine fibers with repeated/conflicting rows, dense cell address ranges,
signed cancellations, invalid/duplicate admission, rekeys, overwrites,
deletion/reuse, noise retention and the every-bit sensitivity witness.
It also checks the full published intersection-transform algebra against
the generic compressed tail, separately from the coded comparator.
Complete small codebooks and sparsely occupied larger codebooks are
reported separately. Larger capacities are conceptual, with no claim that
their full arrays or all independent values were installed in the checker.
The Python tuple maps are correctness controls, not finite-word storage
measurements or runtime comparisons. The proof above accounts for the
dense arrays and full initialization explicitly.

This comparator substantially sharpens the benchmark for polynomial keys,
including the preceding growing-k example. The earlier full-subset
algorithm cannot claim a special advantage from that example alone.
It does not refute the generic candidate:
outside-codebook admitted keys need not have interpolation polynomials or
affine fibers. The remaining research question is whether the generic
low-intersection aggregation theorem gives a useful tradeoff beyond
published aggregate-materialization/set-query methods, or whether a
stronger generic construction or obstruction is needed. The composed
fixed-k bounds above also need a dedicated priority check before any
novelty claim. Preserve them as a baseline and proof artifact; publication,
learned usefulness and the original advantage goal remain unestablished.

## Exact exponential attention with signed query budgets

The subset mechanism extends to a substantially richer query function
than a score of the unweighted intersection size. The theorem below
computes full exponential attention for signed integer queries, with
every active value included. It is an exact structural result, not a
new exponential approximation or a learned-model result. Its multilinear
expansion is classical; a new useful advantage over the strongest
applicable methods and a publication contribution remain unproved.

### Finite model and uniform algorithm

There are at most n active raw binary keys in dimension d, each of weight
w, with pair intersections <k, 1<=k<=w. Each slot holds an arbitrary
independently overwritable signed fixed-grid v-coordinate payload. Write
its integer grid units as V_i. Use a known common dyadic step 2^(-g);
S charges both the unit-magnitude bits and the declared grid precision g.
The displayed formula is in grid units, with that fixed scale applied for
physical value units and charged when materializing output denominators.
A query is
all d raw integer coordinates q_a, with sum_a |q_a|<=B. B is fixed
before execution. Return the complete exact rational attention mean

```text
R(q) = sum_i 2^(sum_a q_a*x_ia)*V_i / sum_i 2^(sum_a q_a*x_ia),
```

or the canonical empty result if no keys are active. This is softmax with
logits ln(2)*<q,x_i>, expressed with exactly representable dyadic weights.
It is not unrestricted real-valued transformer attention. Queries have
no supplied neighbor address or encoded latent; scan every coordinate,
check integrality and the norm budget, and derive their nonzero support.
The theorem assumes this integer query representation; producing it from
semantic inputs or training an encoder is not free or proved here.

Keep the generic subset aggregates, support counts, size-k unique owners,
raw keys, sorted supports and per-key handles from the preceding theorem.
Admission is the complete low-intersection test, not membership in a fixed
polynomial codebook. With F=sum_{j<=k}binom(w,j), T=|supp(q)| and
Q_T=sum_{j<=k}binom(T,j), the following bounds include every payload word:

```text
L = k*ceil(log2(d+1)),
state:      O(n*(ceil(d/W)+w+v+1) + n*F*(L+v+1) + d+w*k) words,
read:       O(d + Q_T*(L+k+w*k+v+1)),
overwrite:  O(F*(v+1)),
raw edit:   O(d + F*(L+v+1)).
```

Deletion is included in the raw-edit bound; an ordinary value overwrite
uses cached handles. Initial loading charges n raw edits plus allocation
and zeroing of the full finite-capacity arenas, raw-key arrays, independent
values, metadata and scratch. Sorted supports cost n*w words, already
displayed. Field arithmetic and a learned codebook are unnecessary.
The finite-word model permits exact integer arithmetic on W-bit words;
choose W for the explicit guard/address bounds below and charge larger
bit arithmetic separately. The Python verifier is an exact map model,
not a measurement of this arena implementation or a runtime benchmark.

Put D=sum_{a:q_a<0}|q_a| and u_a=2^q_a-1. Scale every weight by the
common factor 2^D. For any subset A of the query support, its scaled
multilinear coefficient is the integer

```text
C_A = (-1)^(number of negative q_a in A)
      * 2^(D-sum_{a in A,q_a<0}|q_a|)
      * product_{a in A}(2^|q_a|-1).
```

For s_i=|supp(q) intersect supp(x_i)| define

```text
lambda_i = 2^(D+<q,x_i>),
P_i = sum_{A subset supp(q) intersect supp(x_i), |A|<k} C_A.
```

The product identity product_{a in supp(x_i)}(1+u_a) proves that
P_i=lambda_i whenever s_i<k. Therefore the full scaled sums are exactly

```text
den = sum_{A subset supp(q), |A|<k} C_A*c_A
      + sum_{i:s_i>=k}(lambda_i-P_i),
num = sum_{A subset supp(q), |A|<k} C_A*H_A
      + sum_{i:s_i>=k}(lambda_i-P_i)*V_i.
```

Enumerate subsets to fetch the aggregates and discover all high-order
owners, deduplicating by cleared address marks. There is no all-key scan
or omitted low-score mass. At most binom(T,k) owner occurrences arise.
Unlike the unweighted case, P_i generally depends on which coordinates
are shared, not only their number.

Compute P_i with elementary symmetric polynomials in O(w*k) operations.
Initialize E_0=2^D and E_j=0 for 0<j<k. For each shared coordinate a,
update E_j+=u_a*E_{j-1} in descending order of j. All E_j remain integers:
before introducing a negative coordinate of magnitude m, every E_j has
the common factor 2^(D-M), where M is the total magnitude of already
processed negative coordinates. Since D-M>=m, dividing E_{j-1} by 2^m
is exact. Multiply the resulting integer by -(2^m-1); positive coordinates
need no division. Then P_i=sum_{j<k}E_j. Compute lambda_i from its
nonnegative exponent D+<q,x_i>, which is at most B.

This proves the complete operation bounds. All active keys have strictly
positive true weights, including keys outside the owner list. Overwriting
any such value changes the exact output. The saved witness checks this
explicitly; the candidate list only locates corrections to cached sums.

### Fixed precision and an asymptotic class

Every |C_A|<=2^B and every lambda_i<=2^B. Partial elementary sums have
magnitude at most 2^B*F. With support counts <=n, a conservative scalar
guard is

```text
S+B+ceil(log2(n+1))+ceil(log2(F+1))+ceil(log2(Q_max+1))+O(1) bits,
Q_max = sum_{j<=k}binom(min(d,B),j).
```

This also bounds signed cancellation during the complete numerator/count
calculation. Counts govern cell lifetimes even if the value sum vanishes.
The final denominator is positive for a nonempty active state. Return
signed integer numerators and a shared integer denominator; canceling
the common 2^D does not require floating arithmetic. Address width covers
all arenas, n-slot arrays, query coordinates and v. No guard or metadata
counter grows with the number of overwrites or queries.

For w=Theta(log n), d=log^O(1)n, k=O(w/log w), B=O(w), and
polylogarithmic scalar precision, T<=B and the earlier binomial bounds
give F,Q_max=n^o(1). State and full raw preprocessing are
n^(1+o(1))*(v+1), and each complete read, admitted raw edit and value
overwrite costs n^o(1)*(v+1), including bit arithmetic up to the stated
polylogarithmic factors. This covers all future raw signed integer queries
within the budget, not only the stored keys or a supplied query list.
The integer step is essential: arbitrarily small real or fine-grid
coordinates can make a norm-bounded query dense. On a grid of step
2^(-b), the support bound becomes T<=B*2^b, and the claimed subset
cost must be recomputed; it need not remain subpolynomial. No free
transfer to dense real neural activations follows from quantization.

The finite factors matter. In the six saved parameter examples below,
F already exceeds n. The asymptotic notation does not certify useful
latency, small tables or a win over existing ML implementations.

### Generic admission is feasible without a hidden codebook

For a fresh independent uniform weight-w proposal and any fixed existing
weight-w key, a specified k-subset is contained in the proposal with
probability binom(w,k)/binom(d,k). Union over that existing key's
binom(w,k) subsets, then over at most n current keys, gives the conditional
per-proposal rejection bound

```text
rho <= min(1,n*binom(w,k)^2/binom(d,k))
     <= min(1,n*(e*w^2/(k*d))^k).
```

It is valid after any history if this proposal is freshly uniform and
independent of it. When the displayed untruncated bound is <1, resampling
until acceptance has conditional expected attempts at most 1/(1-rho)
and tail at most rho^t. Resampling is not a worst-case accepted-edit
bound. Every individual raw proposal is validated at the stated worst-case
cost; semantic data are not assumed to be uniform random proposals.
Correctness and retention are deterministic for every admitted history.
In an infinite run, this per-proposal statement is not an event that all
future proposals will be accepted.

For a batch of n independent proposals, the failure probability is at
most binom(n,2)*binom(w,k)^2/binom(d,k). For example let n=2^m,
w=8m, d=w^2 and k=ceil(6m/ceil(log2(m+1))). As m tends to infinity,
the rejection bound is n^(-5+o(1)) and the batch failure bound is
n^(-4+o(1)); simultaneously F,Q_max=n^o(1) at B=2w. This proves a
broad random family, not just an isolated codebook. Exact integer controls
at m=16,32,64,128,256,512 certify per-proposal rejection <=n^-2.
They also report F and Q_max bit lengths, rather than allocating their
very large conceptual tables. The hypergeometric/union-bound reasoning
is elementary and is not a new random-code capacity theorem.

### Retention, independent values, and the scope of a rank obstruction

Let q=b*x_t+e for a positive integer b and sum_a|e_a|<=eta, with this
query in budget. For any other active key,

```text
<q,x_t>-<q,x_i> >= b*(w-k+1)-eta = Delta.
```

Thus the total off-target/target weight ratio is at most
R=(n-1)*2^(-Delta). If |V_ij|<=M in actual value units, each output
coordinate differs from its current target value by at most
2*M*R/(1+R)<=2*M*R. If values lie on a known grid of step delta and
2*M*R<delta/2, ordinary rounding recovers every independent target
coordinate exactly. This survives all accepted key edits, overwrites,
deletions and slot reuse; there is no horizon or accumulated forgetting.
It does not claim fixed points or iterative Hopfield convergence.
The standard separation-to-softmax-concentration mechanism is already
part of [Ramsauer et al., Theorems 4--5](https://arxiv.org/abs/2008.02217).
The primary statements were read on PDF page 5. The cost theorem, signed
integer query model and grid accounting are the additional questions here,
not a discovery of exponential memory capacity or separation-based retrieval.

In the batch family above, b=1, eta=w/4, B=2w, M=1 and grid
delta=2^(-2m) satisfy the strict rounding inequality for sufficiently
large m. The six stated finite examples check that inequality exactly.
The verifier separately installs all 169 affine-polynomial keys over F_13,
weight 13, k=2, with independent integer values in [-3,3]. It exhausts
all signed integer perturbations of norm <=1 at four targets, including
states changed by earlier target overwrites, and checks pure reads after
each target's own independent overwrite. The conservative error bound is
63/128<1/2, sufficient for exact integer-grid rounding. A noncandidate
value overwrite still changes the unrounded rational attention output.

At queries q=x_i, let P be the normalized n-by-n attention matrix.
If gamma=2*(n-1)*2^(-(w-k+1))<1, then ||P-I||_infinity<=gamma,
where this norm is the induced maximum absolute row sum. For every
rank-deficient matrix A, choose Av=0 with ||v||_infinity=1. Then
||(P-A)v||_infinity>=1-gamma. Consequently a fixed feature method of
rank <n cannot approximate this operator uniformly on all independent
bounded values to error <1-gamma. In the fully installed 169-key example,
the exact row norm is 162/2129, giving lower bound 1967/2129.
This is an elementary rank obstruction, not a separation from dynamic
dictionaries, sparse feature states of dimension >=n, indexed correction
algorithms or value-dependent nonlinear methods. Independent grid-value
retention also gives the usual information bound requiring enough bits
to distinguish every value assignment; it does not prove our update cost
optimal. No two-latent value bottleneck is imposed by this architecture.
The operator-norm maximum is attained on sign vectors, so the obstruction
also covers finite value grids containing -1 and +1; it does not depend
on allowing an infinitely precise value witness.

### A simpler exact comparator: coordinate postings and a total

An ordinary inverted index supplies another important baseline. Maintain
the current total payload H, active count m, and a posting list of slots
for each coordinate. For query support S, traverse these lists, accumulating
c_i=sum_{a in S intersect x_i}q_a at touched slots. Untouched keys have
score zero. With D as above, exactly

```text
den = 2^D*m + sum_{touched i}(2^(D+c_i)-2^D),
num = 2^D*H + sum_{touched i}(2^(D+c_i)-2^D)*V_i.
```

This is the h=0 version of the correction idea. It returns the full
function, including untouched values, with no high-order subset tables.
Every signed integer query in the stated budget works; cancellation of
positive and negative coordinates can leave a touched slot with zero
correction, which is harmless. The exact guard is only
S+B+O(log(n+1)) scalar bits plus address width. The same direct dyadic
oracle checks both algorithms' complete numerators and denominators on
every valid test query.

If every coordinate has at most L active keys, at most T*L posting hits
and min(m,T*L) distinct values are visited. Complete read cost is
O(d+T*L+min(m,T*L)*(v+1)). Value overwrites cost O(v+1), updating H
and the independent slot. A linked-node arena with n*w nodes, d heads,
per-key deletion handles, free lists, bounded score counters and cleared
touch marks uses O(n*(w+v+1)+d) words. Sorted supports losslessly describe
the binary keys; supplying raw keys still costs d to validate. The subset
method may use the same lossless representation, so this is not a memory
comparison against a redundant dense key array. The verifier stores both
comparison states; it does not measure this standalone bound.

Generic low-intersection admission also works through postings. Before
mutation, check each proposed coordinate's load after excluding the old
slot. Count intersections at all slots reached through its w posting
lists, rejecting if any other slot reaches k. Keys outside those lists
have intersection zero. Ignore only the replacing slot; a duplicate is
therefore rejected by the same check. This completely validates both
promises in O(d+w*L+w+v) steps. Deletion costs O(w+v+1) using the
cached handles. Check before removal so rejected edits are atomic; accepted
edits update the lists, counts and H. There is no n-slot clearing pass or
unbounded epoch counter. Full initialization charges all input scans,
validation, allocations and pointer setup. A separate exact control tests
a load-only rejection that would pass the intersection condition.

The load cap is an additional promise. Low pair intersection alone does
not imply it. Conversely a random batch of n weight-w keys has coordinate
load mean mu=n*w/d; a union of standard binomial Chernoff bounds gives
Pr(max load>2*mu)<=d*exp(-mu/3), with integer rounding handled upward.
One may prove that elementary inequality directly from
E exp(t*load)=(1-w/d+(w/d)*exp(t))^n and Markov's inequality at
t=ln 2. At d=sqrt(n) times suitable logarithmic factors, w,B and v
polylogarithmic, this baseline has near-linear sparse-key/value state,
complete reads and raw admission in sqrt(n) times logarithmic factors,
and O(v+1) value overwrites. Constant k>=8 makes the simultaneous pair
intersection failure bound n^(-2+o(1)) in that wider dimension regime.
After admission both invariants and the retention bound are deterministic.
These statements assert neither an efficient semantic encoder nor that
the original polylogarithmic-dimension regime has the same load bound.

This is classical postings with an exact background correction, not a
new fundamental architecture. It shows why the subset factors must be
compared with load-sensitive methods and why a low-rank obstruction is
insufficient. Changing raw dimension changes input cost and the feasible
load, and must be stated when comparing storage, access and retention.

### Primary comparisons and decision

[Alman--Song, Theorems 3.7--3.8 and 4.6](https://arxiv.org/abs/2302.13214)
give polynomial-feature approximation and a conditional worst-case barrier
in their normalized real-input model. Their near-linear regime does not
directly cover the strongly scaled sparse example here. Conversely that
worst-case lower bound does not apply to every low-intersection promised
class, and its real-RAM approximation model differs from our exact dyadic
model. No contradiction or general softmax speedup follows.

[Aliakbarpour--Braverman--Yin--Zhang, Definition B.3 and Theorems C.3/I.5](https://arxiv.org/abs/2510.01643)
is a directly relevant sparse-exact plus polynomial comparator. The primary
PDF's Definition B.3 and Theorem I.5 were visually checked on pages 23
and 67. B.3 assumes sublinear total support of large input entries;
ordinary sparsity of every row does not imply it for n rows. I.5 removes
that assumption but its error includes epsilon_2*exp(3*B_entry^2), with
runtime proportional to epsilon_2^-4. These guarantees cannot be
advertised as entry-independent high-accuracy bounds by hiding that
factor. They also do not assert our dynamic exact-query guarantee.
Their methods or refinements might nevertheless cover this promised
class; a theorem mismatch is not proof that they cannot.
Definition I.1 also starts buckets from the minimum absolute input entry;
our sparse inputs have zeros. A zero-aware specialization is needed for
direct use, but that routine boundary issue is not an algorithmic separation.

The contemporary comparison is not complete. In particular,
[Gupta--Huang--Saha--Xu--Ye, ICLR 2026](https://proceedings.iclr.cc/paper_files/paper/2026/hash/8a01099096c85890b1d1aff3c6b4ea56-Abstract-Conference.html)
studies arbitrary-temperature attention with dimension-dependent bounds,
and [S3-Attention](https://arxiv.org/abs/2601.17702) uses learned sparse
features and inverted retrieval. Their primary abstracts were identified;
the full theorem/function/cost comparisons remain to be done. Neither is
treated as excluded by this iteration or as a matching exact dynamic
guarantee on the strength of its abstract.

Finally, the weighted multilinear subset expansion and cached conjunction
sums are established machinery. The new query function is richer than
the unweighted intersection transform, but this extension alone does not
establish a publication-worthy contribution. In particular classical
subset indexing and aggregate materialization remain applicable baselines.
Our own compressed index is already an obvious comparator for the new
application. The fixed-feature rank obstruction cannot substitute for
that strongest-method comparison.

The useful outcome is a precise exact-attention candidate with complete
precision, admission and retention accounting. Preserve it and its exact
certificates. Do not claim a new fundamental advantage, start training or
prepare a contribution paper on this extension alone. The next substantive
question is whether query/update cost can be improved in a useful finite
regime, beyond these established subset and sparse-correction mechanisms,
or whether this branch should be classified fully as a baseline.

## Capacity, admission, and the stronger all-subsets baseline

This iteration classifies the generic sparse-key exponential-attention
workload as a classical baseline. The pair-intersection admission promise
is unnecessary if all subsets of a sparse key are retained. Separate
query support and amplitude budgets permit exact high-temperature reads
without making every coordinate nonzero. This removes a real restriction
from the previous candidate, but does not establish a new architecture,
an advantage over the strongest applicable method, or a publishable result.

### Necessary packing bounds do not guarantee admission

For n distinct weight-w supports in [d], pairwise intersection less than k
implies n*C(w,k)<=C(d,k): each k-subset has at most one owner. Repeatedly
counting supports through a coordinate gives the stronger nested-floor
Johnson bound. Define J(d,w,0)=1 and

    J(d,w,k)=floor(d/w * J(d-1,w-1,k-1)).

Indeed a coordinate belongs to at most J(d-1,w-1,k-1) supports, and
counting all incidences gives the recurrence. It is only a necessary
bound; it supplies neither a family nor an efficient insertion procedure.
For the explicit materialized-subset recipe, a state of n keys has exactly
n*C(w,k) owner entries, while an arbitrary independent value overwrite
changes sum_{j<k} C(w,j) lower-order payload sums. These costs concern that
recipe, not every possible data structure. The certificate records the
smallest k not excluded by Johnson for sixteen (d,w) choices at capacity
2^20, and checks the nested-floor bound against the elementary packing
bound in 363 small integer cases. It does not certify existence or
retention for those sixteen choices. For example d=128,w=8 already needs
k>=5 by this necessary test, with 58,720,256 owner entries in that recipe.

A uniform family shows why free slots and low coordinate load do not
repair the admission promise. Let q be a prime, m>=2, v=q^m. Partition
d=(q-1)*v coordinates into q-1 groups, each an affine space F_q^m. Use
every affine line inside each group as a weight-q key. Any two keys
intersect in at most one coordinate; every within-group pair lies on
exactly one key. There are

    n=v*(v-1)/q keys,
    r=(v-1)/(q-1) keys through every coordinate.

Allow 2*n slots and a coordinate-load limit 2*r. Every fresh weight-q
insertion is still blocked: q points in q-1 groups include a within-group
pair, which already has an owner. Thus this is an insertion obstruction
despite spare resources. Its owned-pair fraction is
(v-1)/((q-1)*v-1), tending to 1/(q-1) as v grows. A small fraction of all
pairs can suffice to block every insertion.

This is an extremal classical packing construction. To see the lower
bound directly, form the graph of owned pairs in any maximal weight-w
pair packing. It has E=n*C(w,2) edges and no independent set of size w.
In a uniform random vertex ordering, the vertices earlier than all their
neighbors form an independent set, with expected size
sum_a 1/(deg(a)+1). Cauchy's inequality gives

    d^2/(2*E+d) <= sum_a 1/(deg(a)+1) <= w-1,
    n >= d*(d-w+1)/(w*(w-1)^2).

The affine construction attains equality with w=q. Balanced disjoint
cliques attain the graph bound, and exact pair designs inside them
partition their edges into keys. Novák's 1995 paper uses precisely this
Turán partition and exact-design mechanism for minimum maximal
quadruple packings: see Theorems 1.1, 2.1--2.3 and their proofs on printed
pages 368--369. We therefore treat the admission application as a transfer
of known extremal design facts, not a new packing theorem.
[Novák, *Packings of pairs with a minimum known number of quadruples*](https://www.dml.cz/handle/10338.dmlcz/126092).

Even replacing one slot can require extensive uniform resampling. Delete
one line in the construction. A weight-q proposal can then use j>=2
points from that deleted line, with the remaining q-j points in different
other groups. No other form is admissible. The exact number of accepted
proposals is

    A(q,v)=sum_{j=2}^q C(q,j)*C(q-2,q-j)*v^(q-j),
    p_accept=A(q,v)/C((q-1)*v,q).

This includes reinserting the removed key. In particular, rejecting that
same key cannot improve the acceptance probability. Since v>=q^2,
putting l=j-2 bounds A by
(q^2/2)*v^(q-2)*sum_l 1/(l!)^2 < 2*q^2*v^(q-2).
Also C((q-1)*v,q)>=((q-1)*(v-1))^q/q! and
q!<=q^q*2^(-floor(q/2)). The elementary bounds
(q/(q-1))^q<=6 and (v/(v-1))^q<=3 therefore give

    p_accept <= 64*(q/n)*2^(-floor(q/2)).

This is a conditional failure of uniform resampling in these states,
not a lower bound on intelligently chosen replacements or all admission
algorithms. It does not invalidate a theorem explicitly conditioned on
accepted updates. It also does not assert a jam family at polylogarithmic
dimension for the earlier candidate's parameters.

Retrieval can remain well conditioned in a jammed state. For inverse
temperature b, the symmetric kernel on these keys is
K(i,j)=2^(-b*(q-|x_i intersect x_j|)), with K(i,i)=1. Every key meets
a_hit=q*(r-1) others at one point, giving the exact off-diagonal row sum

    [a_hit*(2^b-1)+n-1]/2^(b*q).

Whenever this is less than one, strict diagonal dominance proves full
rank. The six certified finite families have this property at b=8.
The statement is not uniform in m at that fixed temperature. This
separates the admission obstruction from an automatic rank collapse.

### Exact attention for arbitrary distinct sparse keys

Here is the stronger classical comparator, including its dynamic costs.
Fix slot capacity N, raw dimension d, maximum key weight w<=d, a global
integer bias a>=0, integer query support bound T0<=d and amplitude bound
B. Every distinct binary key x_i with |x_i|<=w is admissible, including
the empty key. Every value V_i is an arbitrary independent v-dimensional
dyadic vector with declared grid and magnitude. Duplicate keys are
rejected; any other valid raw key can replace a slot. There is no
pair-intersection, codebook or coordinate-load constraint. The complete
query function is

    weight_i(q)=2^(<q,x_i>-a*|x_i|),
    R(q)=sum_i weight_i(q)*V_i / sum_i weight_i(q),
    q in Z^d, |supp(q)|<=T0, ||q||_1<=B.

All active values contribute. This is ordinary exponential attention
with a charged fixed length-bias coordinate on the keys. These Boolean
keys and bounded signed-integer queries are essential hypotheses; no
claim about unrestricted real transformer inputs follows.

Maintain counts c[t,A] and current value sums H[t,A] over keys of length
t containing A, for every subset A of every current key. Length strata
are essential: they permit the factor 2^(-a*t) without treating every
zero query coordinate as a nonzero feature. They also isolate every
independent value. If a key's support is X, its full cell (|X|,X) has
count one, so its sum equals exactly V_i. A separate copy of all N values
is unnecessary. Cached handles for every subset of a slot's key support
let an overwrite read this full cell, form the canonical value delta,
and change all its subset sums. Deletion subtracts the old value and
removes only zero-count cells. A zero signed sum does not mean that a
cell is empty. Admission validates the entire raw key and checks its
full cell for an exact duplicate before any mutation.

For S=supp(q), set u_j=2^(q_j)-1 and D=sum_{q_j<0}|q_j|. Distributing
the product gives, for a key X of size t,

    2^(<q,x>-a*t)
      =2^(-a*t)*product_{j in S intersect X}(1+u_j)
      =2^(-a*t)*sum_{A subset S, A subset X} product_{j in A} u_j.

Consequently querying all cells (t,A) with 0<=t<=w, A subset S, |A|<=t
computes the exact numerator and denominator. A common power-of-two
scale E=D+a*w makes every coefficient integral. Explicitly, the scaled
coefficient of c[t,A] or H[t,A] is

    C[t,A]=(-1)^(number of negative q_j in A)
            *2^(D+a*(w-t)-sum_{j in A,q_j<0}|q_j|)
            *product_{j in A}(2^|q_j|-1).

The exponent is nonnegative. Signed intermediate cancellation is exact;
the final denominator is positive whenever any slot is active. Empty
state returns no value. This is a direct all-subsets product expansion,
followed by the usual normalized kernel-attention aggregation. Both
mechanisms are established: see Polato--Aiolli Section 2 for the
all-subsets product and Katharopoulos et al. Sections 3.2--3.3,
Equations 3--5 and 9--12 for the attention sums. Negative feature
coefficients do not make the positive similarity weight negative.
[Polato--Aiolli, *Propositional Kernels*](https://pmc.ncbi.nlm.nih.gov/articles/PMC8391199/),
[Katharopoulos et al., *Transformers are RNNs*](https://arxiv.org/abs/2006.16236).

Let M bound the number of live (length,subset) cells and let

    M <= Mmax = min(N*2^w,
                   N+sum_{t=1}^w sum_{j=0}^{t-1} C(d,j)),
    Q(T)=sum_{t=0}^w sum_{j=0}^{min(t,T)} C(T,j)
         <=(w+1)*2^T.

The second M bound counts at most N full cells and every possible proper
subset separately within each length stratum. It need not be attained.
Use a balanced ordered dictionary, encoding a whole (t,A) tuple in
O(w*log(d+1)+log(w+1)) bits of a sufficiently wide word. It has O(M)
nodes and worst-case O(log(Mmax+1)) lookup. Per-slot cached cell handles
use at most N*2^w words; store sparse raw supports in O(N*w) words.
Tree deletion can transplant nodes rather than copy another live cell
into the deleted node, preserving all surviving cached handles. Charge
slot/free-node metadata and all counts and pointers as well as payloads.
With word width as specified below, conservative operation bounds are

    state: O(Mmax*(v+1)+N*2^w+N*w+d) words;
    complete raw read:
      O(d+Q(T)*(w+log(Mmax+1)+v+1)) operations;
    value overwrite: O(2^w*(v+1)) operations;
    raw admitted replacement:
      O(d+2^w*(w+log(Mmax+1)+v+1)) operations.

Deletion has the same conservative replacement bound without the raw
scan. Building a state costs its allocation plus the sum of all raw
admission costs; it is not free preprocessing. Scratch is
O(v+T0+w+1) guarded words with depth-first subset enumeration. The
word-operation model charges arithmetic on its explicitly wide words;
bit costs must additionally charge the chosen integer algorithms.
The Python certificate uses tuple maps and arbitrary integers checked
against declared guards; it is not a balanced-tree arena implementation
or a hardware benchmark.

For finite precision let the signed value grid units require at most s
bits, including sign, with the grid scale separately declared. State sums
fit in s+ceil(log2(N+1))+O(1) bits. Every scaled coefficient is at most
2^(B+a*w) in magnitude, and every partial numerator contains at most
Q(T0) terms bounded by N times this magnitude times the value-unit
bound. Thus

    s+B+a*w+ceil(log2(N+1))+ceil(log2(Q(T0)+1))+O(1)

bits suffice for scalar read intermediates. Choose the word width also
to cover tuple keys, coordinates, handles, slot indices and the total
arena address range. Output exact rational grid units with the declared
scale; output construction and any requested rounding are charged.
No fixed-width floating-point cancellation assumption is used.

Retention is possible for every distinct sparse key, rather than only
well-separated ones. Assume actual values are on the grid delta=2^-g
in [-1,1]. Query a target through q=2*a*x_t+e with integer ||e||_1<=eta.
For any other key x_i its score gap from the target is

    a*|x_t symmetric_difference x_i| + <e,x_t-x_i> >= a-eta.

Hence the total off-target weight relative to the target is at most
Roff=(N-1)*2^(-(a-eta)). Every coordinate's absolute output error is
at most 2*Roff/(1+Roff). Taking

    a=eta+g+ceil(log2(8*(N+1)))

makes that error strictly less than delta/4. Nearest-grid rounding
therefore returns the current arbitrary independent value exactly.
Bounds T0>=min(d,w+eta) and B>=2*a*w+eta cover every such integer
perturbation. Canonical deltas preserve the invariants after arbitrarily
many valid edits, with no finite history horizon or precision drift.

This is tolerance of small amplitude jitter around the scaled query
2*a*x_t. It is not a promise to correct binary-key bit flips or erasures.
Two admissible keys can be one bit apart; a shared corrupted binary
input cannot uniquely specify which arbitrary label was intended.
Increasing a cannot remove that ambiguity. There is also no theorem
here about learning a semantic sparse encoder.

The retention promise alone has a simpler dictionary baseline. Since
a>eta, thresholding each query coordinate at a recovers x_t exactly:
target-zero coordinates are below a and target-one coordinates are
above a. A sparse-key dictionary can then return its current value in
O(d+log(N+1)+v) guarded-word operations, with O(N*(w+v+1)) state.
It does not compute the complete exponential-attention function on all
admissible queries. Thus exact rounded retention tests by themselves
cannot establish an access or robustness advantage over dictionaries;
the complete query function needs a separate reason to be useful.

For example w=O(log(N)/log(log(N))), T0=O(w), d=polylog(N),
eta=O(w), g=O(log(N)) and s=polylog(N) give a=O(log(N)),
B=O(w*log(N)). Then state and preprocessing are N^(1+o(1))*(v+1)
words/operations and complete reads and updates cost N^o(1)*(v+1),
including polylogarithmic bit-operation factors. This statement admits
every distinct raw key of weight at most w. Neither a specially selected
codebook nor the earlier random-admission bound is required. It is a
sparse classical feature baseline, not a new separation from it. In
particular k=w in the older fixed-weight formula already permits all
distinct weight-w keys; that endpoint is not a new identity either.
Installing N distinct values additionally requires sum_{j<=w} C(d,j)>=N.
For example d=ceil(log2(N))^3 and
w=ceil(log2(N)/log2(log2(N))) satisfy this for sufficiently large N;
an arbitrary smaller declared universe need not have N distinct keys.

The finite costs need not all be enormous. At N=2^20,d=128,w=4,T0=8,
the exact arithmetic resource bound is 16 subset cells per largest-key
update, 303 cell lookups per read, and at most 358,020 proper-subset
payload cells plus N independent full cells, or 1,406,596 payload cells.
Tuple nodes, handles, sparse key storage, raw input scanning and scalar
bit widths are additional. These are combinatorial counts, not bytes,
measured memory, measured latency or an equal-cost learned-recall win.
No million-slot state was allocated. The six large parameter examples in
the preceding section exposed expensive constants for those examples;
they did not prove that all sparse parameter choices were impractical.

### Contemporary primary-method audit and exact controls

The previously pending contemporary sources have now been read beyond
their abstracts. The official ICLR 2026 version of
[Gupta--Huang--Saha--Xu--Ye, *Subquadratic Algorithms and Hardness for Attention with Any Temperature*](https://proceedings.iclr.cc/paper_files/paper/2026/file/8a01099096c85890b1d1aff3c6b4ea56-Paper-Conference.pdf)
gives a static approximate batch theorem at constant dimension d:
Theorem 1.1 has time tilde-O(n^(2-1/d)*polylog(B_entry/epsilon)).
Theorem 1.2 gives the corresponding constant input-rank r statement
with an nd term. The range-searching and polynomial-feature construction
in Section 3 and Appendix A has dimension-dependent feature costs;
the constants cannot be treated as uniform for growing d. Its conditional
hardness theorem concerns worst-case inputs, not every sparse promised
class. These results neither assert the dynamic exact latest-value
guarantee above nor establish that our specialization beats their best
applicable specialization. The audit covers these theorem statements and
construction, not every gradient or hardness extension in the paper.

[S3-Attention, arXiv:2601.17702v2](https://arxiv.org/abs/2601.17702v2),
Sections 3.3--3.6, Algorithm 1 and Appendix D.6--D.7, learns Top-k sparse
autoencoder features, maintains CPU postings, retrieves text windows,
fuses retrieval signals and reruns the model on a compressed context.
Its constant GPU-memory claim is chunk-local; its CPU index grows with
context length, the indexed layers and feature sparsity. Retained source
text, posting hits, score initialization over the context, convolution,
window selection and the final model pass must all be charged. This is
not an exact all-independent-value sum with latest-value overwrite
guarantees. Conversely our Boolean-key theorem does not prove superior
learned semantic retrieval or end-to-end language-model performance.
The different query specifications are a comparison boundary, not an
advantage certificate.

Reproduce the controls with

    python3 -B scripts/run_timing_restart_checks.py
    ./build.sh
    python3 tests/test_smoke.py

The new `capacity_frontier_baselines` report block checks 10,076 complete
reads against independent direct biased dot-product sums, including
scaled numerator and denominator. It installs every raw key of weight
at most w in four finite universes (d,w)=(5,2),(6,3),(7,3),(8,4), including
the empty key and all length strata. It checks 336 norm-one integer-noise
queries at targets before and after their own value overwrites,
duplicate/out-of-range/too-heavy admission rejection, twelve malformed
queries, zero-sum live cells, deletion and reuse. Six prime affine-line
families certify pair ownership, balanced load, the attained minimum
maximal bound and exact kernel row sums. The four smaller families
exhaust all 25,654 raw insertion proposals and the same proposals after
one deletion, independently checking the replacement formula. The larger
two use complete pair coverage plus the uniform proof, without claiming
exhaustive proposal enumeration. A blocked support (0,1,9) in the
24-key, d=18,w=3 jam is admitted by the full-subsets baseline as a 25th
key, and both an old and the new value are retrieved after rounding.

The certificates support the elementary proofs; they do not establish
novelty. Preserve the earlier certificates and classify this workload's
factorization, generic admission and retention result as a baseline.
The active research goal remains unmet: the next candidate needs a
useful primitive or a nonroutine lower bound beyond this stronger
baseline, with explicit noise and query semantics. Merely replacing a
codebook, changing the temperature or exploiting an uncharged resource
does not qualify. No further neural training or contribution manuscript
is justified by this iteration alone.
