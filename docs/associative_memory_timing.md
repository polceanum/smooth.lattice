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
  the primary abstract identifies a broader existing generic robustness
  reduction and a blinking-adversary refresh model. Only that abstract was
  inspected in this pass. Its exact theorem/model comparison remains pending;
  neither equivalence nor a gap is established by the abstract.

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
