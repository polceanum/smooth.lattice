# Mathematical research on exponent stream interleaving

## Current working assumption and deferred version check

On 2026-10-17 the user instructed us to assume, for ongoing work, that the
relevant results are unchanged between the following versions of
Elizalde--Galvan, *Triangular partitions: enumeration, structure, and generation*:

- Published article: *Advances in Applied Mathematics* 163 (2025), 102807,
  [DOI 10.1016/j.aam.2024.102807](https://doi.org/10.1016/j.aam.2024.102807).
  The public page is readable after the human check, but the full text is paywalled.
- Inspected free preprint: [arXiv:2312.16353v1, 26 December 2023](https://arxiv.org/abs/2312.16353v1).
  This has not been verified as identical to the published article.

**Deferred check for the user:** compare Section 4's Mobius results,
Lemma 7.1's subpartition recurrence, Section 8.2's pyramidal discussion, and
new references or acknowledgements. Record whether any change affects our
interval formula, slice obstruction, counting recurrence, or novelty assessment.
The version assumption allows work to continue; it does not establish novelty
or make the draft submission-ready.

The user also instructed: **never send any emails or messages.** No author
query or other outbound contact is permitted. The old correspondence draft
is retained solely for the user's inspection and is retired from the workflow.

## Research overview

The main objective is a general structural or quantitative contribution to
fixed-prime rank/unrank, with a theory of occurrence transformations retained
where it supplies a useful mechanism or an independently relevant theorem.
The three-group classification is a worked example, not the research endpoint.
Propositions 14--16 below add general reconstruction and transient laws, an
obstruction to recovering three-prime order from finite pairwise patterns,
and a generic exact-unranking baseline derived from established theory.
Proposition 17 gives an integer representation of an entire ordered prefix
with a cutoff score linear in its length for fixed dimension. It is now
included as a supporting lemma in the quantitative manuscript. The paper
claims no new general weight-compression method or faster unranking algorithm.
The preserved structural note concerns Propositions 18--20: a formula for the Mobius function of
pyramidal-partition intervals in every dimension, obtained by a precise
reduction to Edelman's region-poset theorem. It addresses a specific question
in the 2023 Elizalde--Galvan preprint. Proof and finite checks are recorded;
the published-version and subsequent-literature audit is still incomplete.
Proposition 19 supplies an infinite family refuting the coordinate-slice
recognition conjecture in Galvan's 2023 thesis. Every two-dimensional
coordinate section is triangular, but a balanced pair of included and
excluded points obstructs a single global cutting plane. This is a direct
local-to-global obstruction, rather than a claim about prime distribution.

The quantitative goal is recorded in
`results/research/quantitative_goal.json`. The final sections below record a
multiscale proof for corner cuts in every fixed dimension. Candidate
22 gives a uniform finite sieve and cumulative bound. Candidate 23 adds a
same-dimension injectivity argument at a prescribed size and gives the matching
pointwise lower bound. Closed constants, exact lattice-patch witnesses and
the internal proof audit are recorded. The consolidated assessment at the end
finds a substantive improvement over the inspected published bounds and
supports the main lower-count novelty claim under the documented primary
comparison. No referee assessment or journal acceptance is claimed. The
structural note and its elementary consequences do not satisfy this goal.

**2026-10-18 source update:** the user supplied the full Vasiliev--Pavlov
(2009) paper. All eleven pages have now been read and visually checked.
Its strictly admissible finite orderings are exactly the ordered-prefix
class in Candidate 24. Its algorithms and tables are prior art; it states
neither our proposed asymptotic counting orders nor a prime-height bound.
The detailed comparison and exact table checks are recorded at the end of
this note. The specific source-access blocker is resolved. The final
internal contribution and submission assessment is recorded below.

The current manuscript also includes the complete ordered-prefix count,
linear integer-score certificate and sharp worst-case prime-height proof.
Its title is *Corner cuts and monomial prefixes: sharp counts and prime
realizations*. The main pointwise lower count remains its central result;
the additional theorems give a quantitative consequence for the established
finite-ordering class and an explicit arithmetic size guarantee.

**Final internal assessment, 2026-10-18:** the main matching lower count
supports a defensible submission in discrete geometry or enumerative
combinatorics under the documented primary-literature comparison and the
user's final-version assumption. The detailed acceptance audit at the end
separates this contribution judgment from independent refereeing or journal
acceptance. The source-access issue is resolved; it is not a remaining gate.

The first bounded study is complete: for the reconstructed rounded-occurrence
rule, Propositions 10--13 below give a partition reduction, a prefix-period
lifting lemma, a complete classification of periodic frequency partitions for
base b>=3 and word length n<=b^2, and the resulting universal word periods.
For 2<=n<=b+2 every input eventually has period 2; for b+3<=n<=b^2 every
input eventually has period 6. This is a proved statement about the candidate
rule, not a recovery of the user's original infinite construction.

The literature comparison identifies established background and rules out direct
equality with two cited partition operations. Novelty of the restricted theorem
remains unresolved. The decision is to **hold publication and specialized cycle
catalogues pending that comparison**, while investigating the general ordering
problem rather than extending one rounding rule case by case.
The smooth-number algorithms remain a separate computational tool. Prime
applications require additional arithmetic information and uniform estimates.

The propositions below have elementary proofs. They establish a starting point;
novelty requires a literature comparison beyond these derivations. In
particular, an exact reformulation of the next-prime problem is
not itself progress on locating the next prime.

## Existing mathematical connections

For two primes p and q, ordering p^a q^b is equivalent to ordering a+b theta,
where theta=log(q)/log(p). If the coordinates are shifted from a,b>=0 to
a+1,b+1>=1, every weight receives the same additive constant. Thus the first
exponent stream plus 1 is precisely the signature sequence of theta in
Kimberling's terminology. This identifies the two-generator object with an
established family, rather than merely an analogous visual pattern.

Kimberling studies signature sequences, trimming operations, and arrays of
positions in [Fractal sequences and interspersions, 1997](https://combinatorialpress.com/ars-articles/volume-045-ars-articles/fractal-sequences-and-interspersions/).
The specific log(3)/log(2) example is recorded directly by its contributors in
[OEIS A191475](https://oeis.org/A191475).

The usual valuation streams v_p(n), indexed by all positive integers, are
related but differently ordered objects. A recent exposition explicitly
discusses their self-similarity and a sieve built from them:
[Tillman, How Prime Factors Form Fractals, 2025 preprint](https://arxiv.org/abs/2502.15743).
That work does not establish a new bound for the rank maps considered here.

## Prior-art boundary (checked 2026-10-16)

The following comparisons concern the stated definitions, not a claim to have
exhausted the literature. Page numbers refer to the PDFs linked here.

| Source and location | Established connection | Boundary of the comparison |
| --- | --- | --- |
| [Kimberling, 1997](https://combinatorialpress.com/ars-articles/volume-045-ars-articles/fractal-sequences-and-interspersions/); [author's explicit signature definition](https://faculty.evansville.edu/ck6/integer/fractals.html) | Ordering a+b log(q)/log(p) gives precisely the two-prime exponent signature after adding 1. | Exact identification of an existing object; no novelty claim for the two-prime streams. The publisher PDF could not be retrieved for a page-by-page audit in this pass. |
| [Kimberling--Brown, 2004](https://cs.uwaterloo.ca/journals/JIS/VOL7/Kimberling/kimber67.pdf), pp. 2--4, Theorem 1.1 and Corollary 1.2; pp. 17--18, Corollary 4.6 | Their rank array orders the same weighted pairs. Reciprocal weights transpose it. The formula a(i,j)=a(i,0)+a(0,j)+ij-1 gives known two-dimensional rank structure. | Their partial-complement operator t acts on increasing position sequences. It is not the occurrence operator T used here; the paper's convergence conjecture must not be transferred to T or F_b. |
| [Levine, Fractal Sequences and Restricted Nim](https://arxiv.org/pdf/math/0409408), pp. 5--10, Proposition 3.2, Theorems 3.3 and 3.6 | Characterizes Kimberling fractal sequences through Restricted Nim, interspersions and subadditive triangles. | Deleting first occurrences is part of this definition. Deleting zeroes and subtracting 1 is a different condition; the single-coordinate law alone does not establish membership in Levine's class. |
| [Kimberling, 2022](https://cs.uwaterloo.ca/journals/JIS/VOL25/Kimberling/kimber16.pdf), Section 4 and Example 7, pp. 3--4 | Transposition of a regular position array defines a dual sequence; dualizing twice returns the original. | On an infinite lattice word with every symbol occurring infinitely often, this dual is exactly T, after changing between zero- and one-based labels. The smooth-stream occurrence partner is therefore an instance of established duality. |
| [Amini--Jahangiri, arXiv:1612.09481](https://arxiv.org/pdf/1612.09481), Section 3, p. 4; Theorem 12, p. 5 | Uses inclusive occurrence counts to construct a companion sequence and claims a characterization of doubly fractal sequences by signatures. | Recorded as a preprint; neither its claimed characterization nor its construction is needed in our proofs. It is additional evidence against presenting ordinary occurrence duality as new. |
| [Hopkins--Jones, 2006](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1r80/pdf), Section 2, pp. 2--4 | Studies the finite dynamical graph of Bulgarian solitaire D_1(lambda)=sort(ell,lambda_i-1), including known cycle descriptions. | Useful methods and terminology; its map preserves residual rows, unlike our grouped transpose. |
| [Hopkins, 2009](https://math.colgate.edu/~integers/j6proc/j6proc.pdf), Section 2 and Lemma 1, pp. 2--3; Theorems 3 and 5, pp. 7--10 | D_k(lambda)=sort(lambda'_1,...,lambda'_k,lambda_i-k), omitting nonpositive parts. It interpolates between Bulgarian solitaire and conjugation. | We combine entire blocks of columns. The witness below proves that our map is not any D_k on the same partition space. No conjugacy or reduction to all variants has been excluded. |
| [Eriksson--Jonsson--Sjöstrand, arXiv:1703.07099](https://arxiv.org/pdf/1703.07099), Section 1, pp. 2--3, Definition 1 and Theorem 1 | A sigma-solitaire move is sort(sum sigma(lambda_i),lambda_i-sigma(lambda_i)). The well-behaved case imposes monotonicity on sigma and h-sigma(h). | The same witness excludes direct equality with a single sigma-solitaire move. Their stability/limit-shape results do not automatically classify our cycles. |
| [Etessami--Stewart--Yannakakis, 2014](https://arxiv.org/pdf/1304.5429), Theorem 3.1, Propositions 3.5--3.7, pp. 7--10 | Explicit lower bounds for nonzero linear forms in integer logarithms and polynomial-time logarithm approximation give polynomial-time comparison of products of powers when the number of factors is fixed. | This already supplies a generic comparison/precision baseline. Generalizing to arbitrary primes alone would not establish a new theorem. |
| [Barvinok--Pommersheim, 1999](https://library.slmath.org/books/Book38/files/barvinok.pdf), Algorithm 5.2, pp. 112--113 | Counts lattice points of rational polytopes in polynomial time for fixed dimension. | Together with certified rationalization of logarithmic weights, this gives the generic unranking baseline derived in Proposition 16. It is not an algorithm introduced by this project. |
| [Frank--Tardos, 1987](https://doi.org/10.1007/BF02579200), theorem as restated in [Etscheid--Kratsch--Mnich--Roglin, Theorem 1, p. 4](https://arxiv.org/pdf/1507.03439) | Polynomial-time integer weight compression preserves signs on bounded integer vectors. | Small integer weights and short encodings are established tools. Proposition 17 concerns a prefix-dependent linear cutoff, but a distinct bound alone does not establish a substantial advance. The restatement was read; the original paper has not been fully audited. |
| [Onn--Sturmfels, Cutting Corners, 1999](https://www3.risc.jku.at/research/theorema/Groebner-Bases-Bibliography/gbbib_files/publication_536.pdf), Section 3, Lemma 3.1, pp. 35--36; [author-uploaded text](https://www.researchgate.net/publication/2640136_CUTTING_CORNERS_Shmuel_Onn) | A smooth exponent prefix, regarded as a set, is a corner cut. Their criterion separates maximal included points from minimal excluded points by linear inequalities. | The separation LP in Proposition 17 uses this established geometry and adds order constraints. Recognition by an LP is not new. The RISC PDF fetch failed; the author-uploaded text of the criterion and surrounding argument was inspected. |
| [Wagner, On the Number of Corner Cuts, 2002](https://doi.org/10.1016/S0196-8858(02)00014-3); [primary preprint](https://citeseerx.ist.psu.edu/document?doi=4f8e4d8de830cb1d5a398beee165d10531b58dc2&repid=rep1&type=pdf), Observation 4.4 and Theorem 4.5 | Uses lattice geometry and a product bound to encode corner cuts and bound their number in fixed dimension. | These are close precedents for a compact encoding. The indexed preprint text of the introduction, Section 2 and upper-bound argument was inspected; direct PDF retrieval failed. A full comparison of integer normals, ordered prefixes and subsequent work remains open. |
| [Elizalde--Galvan, arXiv:2312.16353v1](https://arxiv.org/pdf/2312.16353), Section 8.2, p. 37 | Reports experimental evidence for Mobius values in {-1,0,1} for three-dimensional pyramidal partitions. | Proposition 18 addresses that specific observation in every dimension. The author's [publication list](https://math.dartmouth.edu/~sergi/research.html) identifies the journal version as Adv. Appl. Math. 163 (2025), 102807, 45 pp.; that final text has not yet been fully inspected. No claim that the question is still open in 2026 follows just from reading the preprint. |
| [Edelman, 1984](https://doi.org/10.1090/S0002-9947-1984-0737888-6), Lemmas 1.2--1.3 and Theorem 1.11, pp. 619--622; [author-uploaded text](https://www.researchgate.net/publication/243079827_A_partial_order_on_the_regions_of_dissected_by_hyperplanes) | Computes the Mobius function of a hyperplane-arrangement region poset by a facial criterion and codimension. | This is the main established theorem used in Proposition 18. The candidate contribution is the finite interval representation and its explicit cutting-layer interpretation, not the region-poset theorem. General interval homotopy is not claimed. |
| [Galvan, Properties of triangular partitions and their generalizations, 2023 thesis](https://upcommons.upc.edu/bitstreams/639d30d2-4361-48cb-be64-a3bc821699a1/download), Conjectures 4.1 and 4.5, printed pp. 51 and 60 | Proposes recognition by triangular coordinate-plane slices and a three-dimensional Mobius rule using joins of two or three atoms. | The full 85-page PDF was retrieved; the relevant definitions, conjectures, conclusion and references were inspected, including visual checks of printed pp. 51, 59--60. Proposition 19 contradicts Conjecture 4.1 as stated; the quadrilateral example contradicts Conjecture 4.5 as stated. This establishes the comparison with the thesis, not priority over all later work. |
| [Galvan, Reports@SCM 9 (2024), 31--40](https://publicacions.iec.cat/repository/pdf/00000506/00000065.pdf), Section 7 and references | States the convex-hull separation characterization in all dimensions, unbounded numbers of addable/removable cells, and an enumeration congruence; cites Pilaud (2023) for a different higher-dimensional failure of a planar characterization. | Section 7 contains neither the coordinate-slice criterion nor a Mobius formula. Convex separation is established background, not a new criterion introduced in Proposition 19. This short summary does not replace inspection of the final 2025 journal article. |

Here is a concrete comparison that does not rely on names. For base 3,

    H_3((9)) = (6,2,1).

On the single row (9), every Hopkins D_k gives (9-k,1,...,1), with at
most one part greater than 1. Every sigma-solitaire move gives at most two
positive parts. Thus neither family contains H_3 as the same one-step map
under the identity identification of partitions. This is a limited distinction,
not a proof that H_b or Proposition 12 has never appeared elsewhere.

Unresolved comparisons include alternative encodings, compositions or factors
of published operations, and a complete reference/citation-chain audit. The
1995 numeration-basis paper and 1999 card-sorting paper identified in the
author's bibliography are not claimed as fully reviewed. There is currently
**no verified literature gap** on which to base a publication claim.

## Definitions

Let P={p_1,...,p_k} be distinct primes, and let

    s_1=1 < s_2 < ...

be their smooth integers. Write e(n) for the full exponent vector of s_n and
E_i(n) for its i-th coordinate. In this note use the value-based count

    C_P(x) = #{s in S_P : s <= x},

with C_P(x)=0 for x<1 and S_empty={1}. This differs from the log-height notation
in theory.md; C_P(exp(T)) is the log-height count there.

Define the insertion map

    R_i(n) = C_P(p_i s_n).

It is the rank occupied by the n-th term after multiplying it by p_i.

## Exact identities for the insertion maps

**Proposition 1.** Each R_i is strictly increasing, and

    e(R_i(n)) = e(n) + unit_i,
    R_i(R_j(n)) = R_j(R_i(n)).

The image of R_i consists exactly of the positions with E_i>0. Consequently,
deleting those positions with E_i=0 and subtracting unit_i from every remaining
full vector reproduces the entire exponent-vector sequence.

**Proof.** Multiplication by p_i is strictly increasing. Unique factorization
gives the statement about exponents and the image. The two compositions both
give the rank of p_i p_j s_n. Division by p_i identifies the ordered image with
the original sequence. This proves a statement about full vectors, stronger
than the corresponding statement about a single exponent stream.

**Proposition 2.** The displacement of a copy is a lower-dimensional count:

    R_i(n) - n = C_{P without p_i}(p_i s_n).

Moreover, if Z_i(N) counts zeroes among the first N terms of E_i, then

    Z_i(N) = C_{P without p_i}(s_N).

**Proof.** Partition the smooth integers at most p_i s_n into multiples of p_i
and those with i-th exponent zero. There are exactly n of the former, by
division by p_i, and the latter form the indicated face. Apply the same face
description at s_N for the second identity.

Subtracting Proposition 2 at consecutive ranks gives an exact local question:

    R_i(n+1) - R_i(n)
      = 1 + #{face values in (p_i s_n, p_i s_{n+1}]}.

Thus insertion gaps are governed by lower-dimensional values inside scaled
intervals between consecutive full-dimensional values. This identity specifies
what a useful interleaving theorem must control; it does not evaluate those
interval counts by itself.

## What self similarity alone cannot determine

There is considerable freedom even before combining several streams. Choose
any set of zero positions containing position 1, with infinitely many nonzero
positions left. Let Z(n) count the chosen zero positions at most n. Then the
recursion

    E(n) = 0                         at a chosen zero position,
    E(n) = 1 + E(n-Z(n))             otherwise

defines a unique stream: the argument n-Z(n) is smaller than n whenever it is
used. At the m-th nonzero position the value is 1+E(m), so deleting zeroes and
subtracting 1 reproduces E exactly. The zero positions of an individual
self-similar stream can therefore be prescribed arbitrarily under these
conditions. The extra arithmetic constraints come from combining the streams
and their ordering.

**Proposition 3.** The vector trimming rule and commuting insertion maps do not
characterize prime logarithms, or even orders defined by one real weight per
coordinate.

**Proof.** For any positive rationally independent weights w_1,...,w_k, order
nonnegative exponent vectors by their scalar product with w. There are finitely
many vectors below any fixed weight. Translation by unit_i preserves order and
maps all vectors onto those with positive i-th coordinate. The proof of
Proposition 1 therefore applies without prime numbers.

There is an even simpler arithmetic control: replacing every generator p_i by
p_i^m for a common positive integer m preserves the complete exponent ordering,
because it raises every sequence value to the same power. For m>=2 all the
generators are composite, although every exponent pattern remains identical.

For a stronger counterexample, order pairs (a,b) first by a+b and then by a,
both increasingly. Every initial segment is finite, and coordinate translation
preserves the order, so trimming and commutation still hold. If positive real
weights u,v reproduced all these comparisons strictly, (0,1)<(1,0) would imply
v<u. But (d,0)<(0,d+1) for every positive integer d would imply
u<(1+1/d)v for every d, hence u<=v, a contradiction.

This last example is a degree lexicographic monomial order. Self-similar
trimming is an order-preserving translation phenomenon. A prime application
must also use the numerical weights or another arithmetic constraint.

## Why two primes are simpler than three

Extend Z_i to Z_i(0)=0 and define Q_i(n)=n-Z_i(n), including Q_i(0)=0.
This counts positions at most n whose i-th exponent is positive. Division by
p_i identifies those values, in order, with the first Q_i(n) smooth values.

**Proposition 3a.** For distinct coordinates i and j,

    Z_i(n) - Z_i(n-Z_j(n))
      = C_{P without p_i and p_j}(s_n).

The right side is symmetric in i and j. Also,

    Q_i(Q_j(n)) = Q_j(Q_i(n)).

**Proof.** Among the first n values, the number with i-th exponent zero and
j-th exponent positive is Z_i(Q_j(n)): divide those values by p_j. Subtract
this number from Z_i(n) to count the values with both exponents zero. These
are precisely the smooth values for the remaining primes at most s_n.
Similarly, Q_i(Q_j(n)) counts those with both exponents positive, by division
by p_i p_j, which proves commutation.

For two primes the remaining face contains only 1, so the identity reduces to

    Z_i(n) - Z_i(n-Z_j(n)) = 1.

For three primes p_i,p_j,r, it instead becomes

    Z_i(n) - Z_i(n-Z_j(n)) = #{a>=0 : r^a <= s_n}.

Equivalently, the right side is 1+floor(log(s_n)/log(r)); the exploratory
script evaluates counts with exact integers, not floating logarithms.

This makes the dimensional change explicit. The individual trimming law is
unchanged. The coupling of two zero masks now carries a nonconstant term from
the third coordinate. With four or more primes the same term is an entire
lower-dimensional smooth count. A comparison of this identity with published
signature-sequence and interspersion results is needed before treating it as
a new result.

## An elementary quantitative baseline

Fix P and k>=2. Put

    w_i = log(p_i), W = product_i w_i, sigma = sum_i w_i,
    A_N = (k! W N)^(1/k).

**Proposition 4.** As N tends to infinity with P fixed,

    R_i(N)-N = alpha_i N^(1-1/k) + O_P(N^(1-2/k)),
    Z_i(N)  = alpha_i N^(1-1/k) + O_P(N^(1-2/k)),

where

    alpha_i = k w_i / (k! W)^(1/k).

**Proof.** Let F(T)=C_P(exp(T)) and V(T)=T^k/(k! W). Unit boxes anchored at
the lattice points counted by F(T) cover the positive simplex of height T and
are contained in the simplex of height T+sigma. Thus, for T>=0,

    V(T) <= F(T) <= V(T+sigma).

At T_N=log(s_N), F(T_N)=N, so

    max(0,A_N-sigma) <= T_N <= A_N.

Apply the analogous bounds in dimension k-1 to the two face counts in
Proposition 2. Replacing T_N by A_N+O_P(1), or T_N+w_i by the same quantity,
changes their leading polynomial term by O_P(A_N^(k-2)). Their leading
coefficient is 1/((k-1)! product_{j!=i} w_j), giving the result.

In particular, the relative number of zeroes tends to zero for fixed P, and

    Z_i(N) / Z_j(N) -> log(p_i) / log(p_j).

The streams therefore retain relative weight information. For the ordinary
integers, by contrast, the proportion with v_p(n)=0 tends to 1-1/p. Passing from
a fixed smooth sequence to all integers is not a uniform consequence of the
fixed-dimensional asymptotic. Even the elementary bound above has a height
uncertainty sigma, which grows with the generator set.

## The precise connection to primes

**Proposition 5.** If P consists of the first k primes, its first missing
positive integer is p_{k+1}, and

    p_{k+1} = min{n >= 1 : s_n > n}.

**Proof.** The least missing integer cannot be composite: both its proper
factors would be smaller and present, and their product would be present. It
is therefore the first prime outside P. The sequence agrees with the integers
before this gap and exceeds its rank at the gap.

**Proposition 6.** For y=p_k and integer y<=x<=2y,

    pi(x) = k + x - C_P(x).

**Proof.** If a composite n<=2y had a prime factor q>y, write n=q m with m>=2.
Then n>2y, a contradiction. Thus all composites in the range are P-smooth,
and the missing integers at most x are exactly the primes in (y,x].

For 0<=H<=y, the deficit (y+H)-C_P(y+H) consequently counts primes in
(y,y+H]. An independent positive lower bound for this deficit would bound
the next-prime gap. Merely substituting the deficit for the prime count does
not supply such a bound. The challenge is uniformity as k grows, rather than
an estimate whose hidden constants depend on a fixed P.

## Counting occurrences and returning after two passes

The occurrence transformation replaces the first occurrence of each value by
0, its second by 1, and so on, with counts taken across the whole sequence:

    T(a)_n = #{m<n : a_m=a_n}.

For example, T(0000)=0123 and T(0123)=0000. This is a different operation from
deleting zeroes and subtracting one, although the smooth exponent streams
have both properties.

**Proposition 7.** For every finite or infinite sequence a,

    T(T(T(a))) = T(a).

Consequently no sequence has a minimal return period greater than 2 under
this particular transformation. Every input enters a cycle after at most
one pass; an arbitrary input need not itself be in the cycle. For example,

    011 -> 001 -> 010 -> 001 -> ... .

**Proof.** Call a word a lattice word when, in every prefix, its occurrence
counts satisfy count(0)>=count(1)>=count(2)>=... . The output b=T(a) always
has this property: each occurrence of r+1 in b comes from an original value
whose r-th occurrence index has already appeared. Thus there are at least
as many r's as r+1's in each prefix.

Now let b be any lattice word. At a position with b_n=r, let c be the number
of earlier r's, so T(b)_n=c. Just before this position, each value j<r has
appeared at least c+1 times, since appending r must preserve the lattice
inequalities. Each j>=r has appeared at most c times. Every value j that
has appeared at least c+1 times contributes exactly one earlier c to T(b).
There are exactly r such values. Hence T(T(b))_n=r=b_n. Applying this result
to b=T(a) proves the identity.

Equivalently, a nonnegative integer word returns after two passes exactly
when it is a lattice word. These words can be represented by growing Young
diagrams: at step n place n in the next cell of row a_n. The transformed
value is the column index of that cell. A second pass exchanges row and
column again. For background on the lattice-word and standard-tableau
correspondence, see
[Chen, Yang and Zhang, 2016](https://arxiv.org/abs/1601.05863).
The direct proof above establishes the occurrence identity without needing
an external theorem or assuming that it is new.

An infinite sequence cannot have period 1: its first entry would have to be
0; its second entry would have to be 0 or 1, and the transformation exchanges
those two possibilities. Its return period, when it returns at all, is
therefore exactly 2. A return after 4 or 6 passes is possible, but a return
after 2 then necessarily occurs as well. Longer minimal cycles would require
a change to the operation, such as an additional reordering or relabelling.

**Proposition 8.** Every individual smooth exponent stream E_i returns after
two passes, for any finite prime set P. Its transformed partner is

    T(E_i)(n) = C_{P without p_i}(s_n / p_i^E_i(n)) - 1.

**Proof.** Put a=E_i(n) and u=s_n/p_i^a. Earlier occurrences of exponent a
correspond exactly to face values v<u, because multiplication by p_i^a
preserves their order. Their count is the face rank of u minus 1. Fixing
u instead, its appearances have exponents 0,1,2,... in that order, so the
occurrence index of this face-rank label is a. Thus T(T(E_i))=E_i.

For two primes the remaining face consists of the powers of the other prime,
so the operation simply exchanges the two exponent streams. For three or
more primes it exchanges an exponent with the rank of a whole remaining
exponent vector. The return law survives; the intermediate sequence no
longer gives just one other prime exponent.

## Longer cycles when occurrence counts are restricted

Restricting the input alphabet to a subset of the natural numbers does not
change Proposition 7: T depends only on equality of entries. Encoding the
output count c by the c-th member f(c) of an increasing allowed set also
cannot produce a longer cycle. Such f is injective, so T(f(a))=T(a), and
the modified transformation F=f composed with T still satisfies F^3=F.
For example, writing counts as 0,2,4,6,... is only an injective relabelling.

A different way to enforce allowed values is to round each ordinary occurrence
count to an allowed number. This can merge distinct counts and permits longer
cycles. The leading reconstruction uses powers of a base, as follows.

For base 3, allow 0,1,3,9,27,... and round every ordinary occurrence count
down to the largest allowed value. There is an exact six-cycle:

    010111 -> 001113 -> 010110 -> 001111 -> 010113 -> 001110 -> 010111.

The user reports that this example looks familiar. The operation remains a
candidate reconstruction, rather than a confirmed recovery of the original
infinite sequence.

**Proposition 9.** Let b>=3 be any integer, and let F_b apply T and then round
each count down to 0 or a power of b. Write 1^r below for a block of r ones,
not exponentiation, and put

    A_x = (0,1,0, 1^(b-1), x),
    B_x = (0,0,   1^b,     x).

Then these words of length b+3 form the minimal six-cycle

    A_1 -> B_b -> A_0 -> B_1 -> A_b -> B_0 -> A_1.

**Proof.** In these words all occurrence counts lie between 0 and b.
Rounding sends 0 to 0, each count 1,...,b-1 to 1, and b to b.
Counting in A_1 gives the initial outputs 0,0,1, followed by b-1 ones
and a final b; this is B_b. Counting in B_b gives 0,1, then 0 and b-1
ones, and finally 0; this is A_0. The same count in A_0, B_1, A_b,
and B_0 gives B_1, A_b, B_0, and A_1 respectively. The A words differ
from the B words at their second entry, and the three final entries
0,1,b are distinct. Hence all six words are distinct and no earlier
return occurs.

The family includes prime bases 3,5,7,..., but the argument also works for
composite bases. This is a fact about a restricted counting operation; it
does not establish a property specific to prime numbers. These words also
do not satisfy the earlier delete-zeroes-and-subtract-one law: for example,
the positive entries of A_1 would all become zero after subtracting one.
Recovering an infinite sequence with an appropriate fractal rule remains
a separate problem.

Other nearby reconstructions give longer cycles as well. One uses the allowed set {0,2,4,6,...} and

    F(a)_n = 2 ceil(T(a)_n / 2).

In words, count previous occurrences as usual, then round that count up to
the nearest even integer. It has the exact finite cycle

    020222 -> 002224 -> 020220 -> 002222 -> 020224 -> 002220 -> 020222.

These six words are distinct, so the first return is after six passes.
For example, the ordinary counts of 020222 are 001123, which round to
002224. The zero-based convention is retained: a first occurrence maps to 0.

Another natural restriction is to round counts down to squares, using
G(a)_n=floor(sqrt(T(a)_n))^2. It gives the finite six-cycle

    0011111 -> 0101114 -> 0011110 -> 0101111 -> 0011114 -> 0101110 -> 0011111.

These are additional candidate reconstructions of the remembered variation.
The exact workbench checks every transition and the minimal period of both
examples, and checks the base-power templates for bases 3 through 31.

## Partition reduction of the reconstructed operation

The user confirmed T (count all strictly earlier equal entries). The additional
rule q_b below was introduced during reconstruction and remains unconfirmed:

    q_b(0)=0,
    q_b(c)=b^j when b^j <= c < b^(j+1), for c>=1,
    F_b(a)=q_b(T(a)), entry by entry.

All results in this section concern finite words and eventual cycles. They do
not assert that the initial word itself returns, or define an infinite fractal.

For a finite word a, let lambda(a) be its positive symbol frequencies, sorted
in nonincreasing order. Let lambda'_j be the height of column j of its Ferrers
diagram, using j>=1. For any fixed count map q, define a labelled histogram
and its sorted partition by

    L_q(lambda)_s = sum_{r>=0 : q(r)=s} lambda'_(r+1),
    H_q(lambda) = sort(the positive entries of L_q(lambda)).

Only finitely many columns and output labels occur at each fixed word length.

**Proposition 10 (factor map).** For every finite word a and count map q,

    frequency_histogram(q(T(a))) = L_q(lambda(a)),
    lambda(q(T(a))) = H_q(lambda(a)).

**Proof.** A symbol with multiplicity m contributes exactly the occurrence
indices 0,...,m-1, independently of its positions. There are lambda'_(r+1)
symbols contributing index r. Combine those indices with equal q(r), retaining
their labels for L_q and sorting their positive totals for H_q.

For q=q_b this transposes the diagram and combines columns in the blocks

    {1}, {2,...,b}, {b+1,...,b^2}, {b^2+1,...,b^3}, ... .

**Proposition 11 (all-prefix period lifting).** Let a have length n>=1.
For each j=1,...,n let p_j be the minimal eventual period of the partition
orbit starting from lambda(a_1,...,a_j) under H_q. Then the minimal eventual
period of a under F=q composed with T is exactly

    lcm(p_1,...,p_n).

**Proof.** All these orbits are eventually periodic: after one pass each word
uses only labels in the finite set q({0,...,n-1}), and partitions of each mass
are finite in number. F commutes with taking a prefix. Each prefix partition
is therefore a factor of the word orbit, so p_j divides its eventual period.

Conversely choose a time after every prefix partition has entered its cycle,
and put L=lcm(p_1,...,p_n). By Proposition 10, the labelled histogram of each
prefix at the next time is determined by that prefix's current partition.
All these histograms therefore repeat after L steps. The histogram of prefix
j minus that of prefix j-1 is one unit at exactly the label of the j-th entry.
Thus the entire word repeats after L steps, proving the reverse divisibility.
This argument concerns eventual periods; it does not equate transient lengths.

Keeping only the partition of the whole word is insufficient. For example,
the six-cycle at base 3 has full-word frequency partition period 3. Its
length-2 prefix supplies a factor 2, giving word period lcm(2,3)=6.

## Complete classification before the second power threshold

Write sort for sorting positive parts, and omit zero parts. Throughout this
section b>=3 is an integer and 1<=n<=b^2. The maximum occurrence index is
n-1<b^2, so only the output labels 0,1,b can appear. If lambda has mass n
and ell parts, put

    R(lambda) = sum_i max(lambda_i-b,0).

The factor formula becomes

    H_b(lambda) = sort(ell, n-ell-R(lambda), R(lambda)).             (A)

In particular every partition enters the set of partitions with at most three
parts in one step. This justifies exhaustive checking of that reduced set;
it does not require enumeration of every unrestricted partition of n.

**Proposition 12 (all periodic partitions).** The following lists every
periodic partition and its transition under H_b. Families outside n<=b^2
are simply absent. Cycles are directed, and all displayed tuples are sorted.

| Mass n | Periodic partitions and transitions |
| --- | --- |
| 1 | (1) is fixed. |
| 2 | (2) <-> (1,1). |
| 3<=n<=b+2 | sort(n-2,2) is the only periodic partition, and is fixed. |
| b+3 | sort(b+1,2) -> sort(b,2,1) -> sort(b,3) -> sort(b+1,2). |
| b+4 | sort(b+1,3) -> sort(b+1,2,1) -> sort(b,3,1) -> sort(b+1,3). |
| b+5<=n<=2b+3 | Put s=n-3 and h=s-b. Exactly S_t=sort(3,t,s-t), for 1<=t<=h-1, are periodic, with S_t -> S_(h-t). |
| 2b+4<=n<=3b+2 | Put s=n-3 and d=s-2b. Exactly S_t=sort(3,t,s-t), for d<=t<=b, are periodic, with S_t -> S_(b+d-t). |
| 3b+3<=n<=b^2 | sort(3,2b,n-2b-3) is the only periodic partition, and is fixed. |

The two reflection families have period 1 when t equals its partner and
period 2 otherwise; count each pair only once. The two displayed 3-cycles
have three distinct partitions. Thus the only possible partition periods in
this range are 1,2,3, and period 3 occurs only at n=b+3 or b+4.

**Proof.** Formula (A) gives the small masses and the two displayed 3-cycles
directly. We show the list is exhaustive.

First, a one-part partition (n) can only be preceded by (1,...,1). A periodic
state has at most three parts, so n<=3; (3) leads to the fixed state (2,1).
This leaves only the one-part states already listed for n=1,2.

For n>=3 a periodic two-part partition must therefore be an image of a state
with two or three parts. Since (A) has two positive entries only when R=0,
it must be U=sort(n-2,2) or V=sort(n-3,3), with V used only for n>=4.
For n<=b+2, U is fixed and V leads to U. For n=b+3 or b+4, direct substitution
gives the displayed cycles; at n=b+4 the remaining candidate U maps to
sort(b,2,2) and then to V, so U is not periodic.

For n>b+4 neither U nor V is periodic. To see this without a search, write

    H_b(U) = sort(2,b,n-b-2),
    H_b(V) = sort(2,b+1,n-b-3).

If n<=2b+2, the next image after H_b(U) is V. If n>=2b+3, that image has
three parts including 3; at n=2b+3 it is S_1 in the first reflection family,
and for n>2b+3 every state with three parts including 3 stays three-parted,
since s=n-3>2b makes R positive. Similarly, for n<=2b+3, the next image
after H_b(V) is Y=sort(3,n-4,1)=S_1. It lies in the first reflection family
(h=n-b-3>=2). For n>2b+3 that image again stays three-parted. Thus these
paths cannot return to U or V.

It remains to classify cycles consisting entirely of three-part states.
Each such state is the image of another three-part state, so it contains 3.
Write it uniquely as S_t=sort(3,t,s-t), removing one copy of 3 and choosing
1<=t<=floor(s/2). Because 3<=b, its transition is

    R = max(t-b,0) + max(s-t-b,0),
    H_b(S_t) = sort(3,s-R,R).                                    (B)

If s<=b+1 then R=0, so it leaves the three-part set. If b+2<=s<=2b, put
h=s-b. Exactly 1<=t<h give R>0; on these states (B) gives t -> h-t.
The other parameters leave the three-part set. This proves the first family.

If 2b<s<3b, put d=s-2b. On d<=t<=b, (B) gives t -> b+d-t. Parameters
t>b map to d; parameters t<d map to a value greater than b and then to d.
Thus every remaining state enters this reflection family, proving the second.

Finally, if s>=3b, every t<b maps to a parameter at least b. All parameters
t>=b map to the constant min(s-2b,2b), which is itself at least b. The
resulting state is the single fixed partition sort(3,2b,s-2b). This exhausts
the possibilities and completes the proof.

**Proposition 13 (universal eventual word periods).** For every starting word
a of length 1<=n<=b^2, independently of its entries, F_b has minimal eventual
period

    1,  if n=1;
    2,  if 2<=n<=b+2;
    6,  if b+3<=n<=b^2.

**Proof.** Proposition 12 bounds every prefix partition period by 1,2,3.
The mass-2 prefix always has partition period 2. Below b+3 all other prefix
partitions have period 1, and Proposition 11 gives 2. At and above b+3,
the prefix of mass b+3 has period 3 for every initial partition, since the
classification gives only one cycle at that mass. Proposition 11 now gives
lcm(2,3)=6. A length-1 word maps to the fixed word (0).

This proves a family rather than just the examples of Proposition 9. For
base 3, a constant seed of length 6 reaches a six-cycle after three passes;
it does not itself return. The actual periodic seed 010111 returns after
exactly six passes. Outside the classified range, the constant base-3 word
of length 64 has transient 10 and minimal eventual period 12. Its full cycle
and all prefix partition periods are saved as an exact finite witness, showing
why the restriction n<=b^2 cannot simply be dropped.

## General reconstruction, rather than a special cycle rule

Propositions 10 and 11 already apply to **any fixed count map q and every
finite word length**. They do not require monotonicity, powers, prime bases,
three frequency groups, or the bound n<=b^2. The following sharpens that general
statement. It remains a structural derivation, with no novelty claim attached.

**Proposition 14 (one-step reconstruction and transients).** Let

    Pi(a) = (lambda(a[:1]),...,lambda(a[:n]))

be the entire chain of prefix frequency partitions. The next word F_q(a)
is determined by Pi(a): its j-th label is the unique unit increment between
L_q(lambda(a[:j])) and L_q(lambda(a[:j-1])), with the empty histogram at j=0.
If mu is the minimal word transient and nu is the maximum of the minimal
transients of its prefix partition orbits, then

    nu <= mu <= nu+1.

**Proof.** Proposition 10 gives both labelled histograms, whose difference
counts just the j-th output entry. This reconstructs the whole next word.
Projection from a word orbit to its prefix partitions shows nu<=mu. After
nu steps all prefix partitions are on their cycles. Their tuple repeats with
period L=lcm(p_1,...,p_n). Reconstruction therefore puts the word on an
L-periodic orbit by time nu+1. Proposition 11 makes L its minimal period,
giving the upper bound. Both bounds can occur: for q(c)=c, word 0000 has
nu=mu=0, whereas word 011 has nu=0 and mu=1.

This is a complete finite reduction of the occurrence dynamics to prefix
partition dynamics, including reconstruction of the periodic words, not just
a bound on their periods. It is not a closed classification of every H_q
orbit, nor a claim that computing those orbits is efficient in all parameters.

The numerical labels of q are also inessential. If q'=f composed with q for
an injective relabelling f, then for every t>=1,

    F_(q')^t(a) = f(F_q^t(a)).

Indeed T(f(w))=T(w), so induction gives the identity. Eventual directed cycles
are therefore carried bijectively to one another by f. Initial transient
lengths may differ. Only which occurrence indices q merges determines this
eventual dynamics; using primes as labels supplies no additional arithmetic
information to this operation.

## A general obstruction to extending pairwise patterns

**Proposition 15 (finite pairwise data are insufficient).** For every finite
M>=1 there exist two ordered triples of distinct primes, with the same first
and second primes, such that all three corresponding two-prime sequences
have identical exponent-vector prefixes of length M, but their three-prime
exponent vectors differ at rank 5.

This concerns exponent patterns. The actual prime values are available to an
arithmetic algorithm and distinguish the two instances; they are not assumed
to be part of the pattern-only data.

**Proof.** We may assume M>=2, since that also implies the case M=1.
Start with normalized positive real weights

    (1, alpha, beta),  alpha=sqrt(2), beta=1+alpha.

For a pair of these weights, no nonzero integer two-coordinate linear form
vanishes: alpha, beta and beta/alpha are irrational. Among vectors in the
finite box {0,...,M-1}^2 there are consequently finitely many strict pairwise
comparisons and a positive minimum margin. A sufficiently small perturbation
of beta both upwards and downwards preserves all those comparisons.

In the three-coordinate order, however, beta crosses 1+alpha. With a small
enough perturbation, the first four vectors in both cases are

    (0,0,0), (1,0,0), (0,1,0), (2,0,0),

and the next two are (0,0,1) and (1,1,0) in opposite orders. The strict
inequalities 1<alpha<2 ensure that 2*alpha and 3 are larger than beta and
1+alpha, so no other vector intervenes.

These open conditions can be realized by primes, not merely real weights.
Fix such a perturbation epsilon>0. For arbitrarily large K, Bertrand's
postulate supplies primes in each of

    [exp(K),2 exp(K)], [exp(K alpha),2 exp(K alpha)],
    [exp(K(beta-epsilon)),2 exp(K(beta-epsilon))],
    [exp(K(beta+epsilon)),2 exp(K(beta+epsilon))].

Use the first two primes in both triples, and one of the latter two as the
third prime. Their logarithms divided by K approximate the target weights
with error at most log(2)/K. Taking K sufficiently large preserves every
strict comparison above and makes the four intervals disjoint. Thus the
triples use distinct primes and straddle the product of their shared first
two primes.

Finally, the first M vectors of a positive-weight two-coordinate sequence
lie inside {0,...,M-1}^2: an exponent >=M would already have at least M
strictly smaller pure powers. The preserved box comparisons therefore give
the same prefixes of length M. This proves the statement.

For a concrete exact witness, take

    P_minus = (101,683,68963), P_plus = (101,683,68993),
    68963 < 101*683 = 68983 < 68993.

Every corresponding pair has the same first **10,000 exponent vectors**.
The two full triples have respectively (0,0,1) and (1,1,0) at rank 5.
The script verifies primality by trial division, compares the pair prefixes
using exact integer ordering, and checks short prefixes with the independent
pointer generator. The evidence includes hashes of the shared pair words.

This explains a limitation of a natural extension strategy: no fixed amount
of pairwise pattern data works uniformly over prime triples. It does not rule
out adaptive methods, algorithms using the numerical primes, or recovery from
the entire infinite pairwise streams. The elementary separating-wall argument
is not presented as a publication claim.

## A generic exact-unranking baseline from existing theory

The relevant broad problem remains the repository's original one: given P
and rank N, return the exponent vector, without constructing the enormous
integer itself. There is already a theoretical generic route in fixed
dimension. The following is our deduction combining established results;
neither cited paper is being attributed this exact smooth-unranking statement.

**Proposition 16 (standard-results corollary).** For each fixed k>=1, exact
unranking over k input distinct primes is computable in deterministic time
polynomial in their binary encoding length and log(N+1), returning the
exponent vector in binary. Constants and polynomial degrees may depend on k.
This is an existence/bit-complexity baseline, not an implemented kernel or a
claim of practical speed.

**Derivation.** N=1 returns the zero vector. For N>=2 let H=N-1. Every
coordinate of the N-th vector is at most H, because a_i>=N would put the
N distinct pure powers p_i^0,...,p_i^(N-1) strictly below it. Thus selection
may be restricted to the bounded exponent box {0,...,H}^k.

Unique factorization makes every comparison of distinct vectors strict.
An explicit linear-form-in-logarithms bound gives a simultaneous positive
separation for all their differences, whose coefficients have absolute value
at most H. For fixed k the number of bits needed to resolve this separation
is polynomial in the prime encoding length and log(N+1). Theorem 3.1 and
Propositions 3.5--3.7 of
[Etessami--Stewart--Yannakakis](https://arxiv.org/pdf/1304.5429) provide the
bound and the required polynomial-time rational log approximation.

Consequently choose a common dyadic scale Q and positive integer weights A_i
so that A_i/Q approximate log(p_i) to less than delta/(4kH), where delta is
a certified lower bound for that simultaneous separation. Their encoding
lengths are polynomial for fixed k. For every distinct a,b in the box,

    abs((a-b) dot (log(P)-A/Q)) < delta/4.

Hence the integer scores A dot a are distinct and give exactly the same
order as the products of prime powers throughout the box. This is one
simultaneous rationalization, not separately answering every comparison.

For an integer threshold t, count points of the rational polytope

    0 <= a_i <= H,  A dot a <= t.

[Barvinok's fixed-dimension counting theorem](https://doi.org/10.1287/moor.19.4.769)
does this in polynomial time in its encoding length. Binary search over
0<=t<=H sum A_i finds the N-th score in polynomially many calls. At that
score there is one integer point. Each coordinate is recovered by binary
searching its upper bound within the score-equality polytope, again using
counts. The number and encoding size of all queries remain polynomial for
fixed k. This completes the deduction.

Neither the separation bound nor fixed-dimension counting promises a cheap
algorithm as k grows. Irrational weights cannot simply be handed to a
rational-polytope algorithm without the separation step. These are precisely
the qualifications that prevent an apparent generic solution from hiding
its difficult part.

## Ordered prefix compression: first candidate under the active goal

**Proposition 17 (proved bound, novelty unverified).** Let w be a positive real
k-vector whose scalar products with distinct points of Z_{>=0}^k are distinct.
Let a_1=0,...,a_N be its first N lattice points in increasing scalar-product
order. There is a positive integer k-vector A such that these are exactly
the first N points under A, in exactly this order, with no ties through a_N.
Writing S=A dot a_N, one can arrange

    N-1 <= S <= C_k N,       C_k = k (k!)^2 2^(k-1).

For N=1 use A=(1,...,1) and S=0. The assertion concerns the order of the
exponents, not equality of integer scores and smooth-number values. Ties
after the prefix are permitted. In particular it applies to w_i=log(p_i)
for any distinct primes, independently of their sizes.

**Proof.** Suppose N>1, put T=w dot a_N and
D={a_1,...,a_N}. A coordinate i is active when w_i<=T. Let r be the number
of active coordinates, so 1<=r<=k. Inactive coordinates vanish throughout D;
for the moment restrict to the active coordinates. Put

    B_i = floor(T/w_i) >= 1.

Then D is a downset in the box 0<=a_i<=B_i. The simplex
{x>=0 : w dot x<=T} is covered by the N unit cubes based at its integer
points: for each x in the simplex, floor(x) belongs to D. Comparing volumes,

    product_i B_i <= T^r / product_i w_i <= r! N.                 (F)

Let E be the set of coordinatewise minimal points outside D. It is finite;
every point outside D dominates an element of E. Every b in E has
b_i<=B_i+1, since otherwise b-e_i is already outside D. Consider the LP

    x_i >= 1,
    (a_{j+1}-a_j) dot x >= 1           (1<=j<N),
    (b-a_N) dot x >= 1                (b in E).                 (G)

The finitely many left sides are strictly positive at w, so a sufficiently
large multiple of w satisfies (G). The feasible polyhedron has a vertex:
minimize sum_i x_i, using x_i>=1 to make each bounded sublevel set compact,
and choose a vertex of the compact optimal face. At such a vertex choose r
linearly independent tight constraints, giving an integer matrix M with
M x=1. Every entry in column i has absolute value at most 2 B_i.

By Cramer's rule A_i=|det M| x_i are positive integers. The determinant
expansion for the matrix with column i replaced by ones gives

    A_i <= r! 2^(r-1) product_{j!=i} B_j.

All constraints remain satisfied after this scaling, since |det M| is a
positive integer. Consequently the included points have strictly increasing
integer scores, and every excluded point has score greater than S. Moreover,

    S <= sum_i B_i A_i
      <= r r! 2^(r-1) product_i B_i
      <= r (r!)^2 2^(r-1) N <= C_k N.

Restore each inactive coordinate with weight S+1. Any point using one of
them has score greater than S. There are N distinct nonnegative integer
scores starting at zero, which gives S>=N-1. This completes the proof.

The primitive integer vector obtained by dividing the active weights by
their common divisor also satisfies the strict order and the same bounds.
The proof only needs strict ordering through the included points and the
cutoff; global rational independence is a sufficient hypothesis.

**What this supplies.** For fixed k, the entire ordered prefix has a
description consisting of N and k positive integers with O_k(log(N+1))
bits, independent of prime bit lengths. The cutoff grows linearly in N;
the lower bound N-1 shows that its order of growth cannot be improved.
This does not assert an optimal dimension constant. Nor does it certify
agreement with an input prime set merely from the short description:
that agreement needs to be established separately.

**Known ingredients and current limits.** The separation LP is an ordered
refinement of Onn--Sturmfels' corner-cut criterion. Determinant bounds for
rational vertices and general integer weight compression are established
tools. Wagner's product-bound encodings of corner cuts are particularly
close prior art. We have not established that the linear cutoff bound for
the full ordered prefix is absent from this literature or follows only
with a new argument. Treat Proposition 17 as a proved internal lemma until
that comparison is resolved, rather than advertising a new compression theorem.

The implementation takes an explicitly enumerated prefix, builds its minimal
excluded border, and enumerates rational LP vertices exactly. It is a small
verification oracle, not a random-access algorithm. A search cap raises an
error and never counts as evidence of infeasibility. Thus the existence result
does not yet improve the project's core rank/unrank task or pass the impact gate.

## Pyramidal partitions: a concrete question from recent literature

Here a pyramidal partition is a finite set of exponent vectors cut off by a
positive hyperplane. It is a different object from the frequency partitions
used for the rounded-occurrence dynamics. The smooth-number connection is
direct: forgetting the order within a prefix leaves precisely such a cut,
with weights log(p_i).

Let K_k be the poset, under inclusion, of finite sets

    D(w,t) = {a in Z_{>=0}^k : w dot a < t},       w_i > 0,

whose cutting hyperplane contains no lattice point. Include the empty set.
Using a weak cutoff instead gives the same collection of finite sets, since
the cutoff can be moved slightly. Translating every coordinate by 1 identifies
this convention with the positive-coordinate convention of the cited papers.
For A<=B, the interval [A,B] is finite. Define mu(A,A)=1 and

    mu(A,B) = - sum_{A<=D<B} mu(A,D)              (A<B).

This is the incidence Mobius function of a containment poset, not the
arithmetic Mobius function on positive integers.

Section 8.2 of [Elizalde--Galvan's 2023 preprint](https://arxiv.org/pdf/2312.16353)
reports evidence that mu only takes values -1, 0, 1 in dimension 3. The
following argument supplies a formula in all dimensions, conditional only
on the correctly stated published theorem of Edelman. The remaining research
question is priority and publishability of this application, not a missing
uniform proof.

**Proposition 18 (application of Edelman's theorem).** Let A<B in K_k and
U=B\A. Then mu(A,B) is nonzero if and only if there are w_i>0 and t with

    w dot a < t       for every a in A,
    w dot u = t       for every u in U,
    w dot b > t       for every b outside B.                    (H)

When (H) holds,

    mu(A,B) = (-1)^(1 + affine_dimension(U)).                    (I)

Consequently mu takes values only in {-1,0,1}, in every dimension.
Condition (H) says that U is one entire lattice layer on a positive cutting
hyperplane, with A below it and the rest of the complement above it.

**Finite interval representation.** Write E=min(Z_{>=0}^k\B), the finite
minimal excluded border of B. In R^{k+1}, with coordinates (w,t), take the
central hyperplane arrangement consisting of

    w_i = 0                  (1<=i<=k),
    h_a(w,t)=w dot a-t = 0    (a in B union E).

This arrangement is essential: the coordinate normals and any one h_a
span all k+1 coordinates. There are regions R_A and R_B corresponding to
the strict cuts A and B. They have positive w_i, negative h_a for included
points and positive h_a for excluded points. Relative to R_A, the separating
hyperplanes of R_B are exactly {h_u=0 : u in U}.

Order regions by inclusion of their separating sets, as in Edelman's
definition. A region in [R_A,R_B] keeps all coordinate signs positive, all
h_a for a in A negative, and all h_b for b in E positive; only signs for U
can change. Its positive weights therefore define a finite cut D containing
A. Every point outside B dominates a point of E, so no such point belongs
to D. Thus A<=D<=B. Conversely every such cut D defines exactly one of
these arrangement regions. Both directions preserve inclusion, proving

    [A,B] in K_k is isomorphic to [R_A,R_B] in the region poset.  (J)

This is an interval of a finite arrangement, not an assertion that the whole
infinite poset is a lattice or that arbitrary restrictions preserve Mobius values.

**Application of the facial criterion.** Edelman's Theorem 1.11 says that,
relative to a base region R_A, the value at R_B is nonzero precisely when
its separating set is the set of all arrangement hyperplanes containing a
nonempty face F of the closed base region. The value is then (-1)^codim(F).

Here the separating set is exactly the hyperplanes for U. A point in the
relative interior of such an F has equality for U and the strict R_A signs
for every other hyperplane, giving positive w, the inequalities for A and E,
and hence (H). Conversely a point satisfying (H) has precisely that zero
set and belongs to the relative interior of a face of the base region. To
see it is in the closure, add epsilon times any point of R_A and let epsilon
decrease to zero; the unchanged strict signs persist and all U signs agree
with R_A. Thus the facial criterion is exactly (H).

The normals of the hyperplanes for U are (u,-1). Their rank is
1+affine_dimension(U), so this is codim(F). The face has positive dimension,
since positive w gives a nonzero ray in this central arrangement. Equation
(I), the value range, and the uniform proof now follow from (J).

**Computing a value without enumerating its interval.** It suffices to impose
the inequalities for maximal elements of A and minimal excluded elements of
B. Fix u_0 in U, eliminate t=w dot u_0, and solve

    w_i > 0,
    (u_0-a) dot w > 0             (a maximal in A),
    (b-u_0) dot w > 0             (b minimal outside B),
    (u-u_0) dot w = 0             (u in U).

Because there are finitely many homogeneous strict inequalities, feasibility
is equivalent to replacing every >0 by >=1 after a common positive scaling.
It is a rational LP together with an exact rank calculation. Given explicit
point lists for A and B, the border has at most k|B| candidates (or one when
B is empty), so the problem has polynomial encoding size and can be solved
by a polynomial-time LP algorithm. This statement does not require fixed k.
It concerns explicit point-list input, not binary-encoded gigantic prefixes.
The supplied small-case verifier uses exact Fourier--Motzkin elimination;
that implementation is not claimed to have polynomial running time.

**Evidence and scope.** The verifier exhaustively enumerates every downset
inside five specified simplex tops in dimensions 1--4, recognizes all corner
cuts by strict rational separation, and computes mu by its defining recurrence
for every comparable pair. All 27,100 intervals agree with (H)--(I). The
simplex tops have (dimension, degree) equal to (1,6), (2,4), (3,2), (3,3),
(4,2); this is not an enumeration of all partitions of a given mass. Exact
integer cutting witnesses, masks and interval values are saved in
`results/research/pyramidal_poset.json`. Separate vertex enumeration also
checks recognition on 233 nonempty downsets in the smaller fixtures. Three
targeted controls distinguish a complete line, a complete plane and a partial
line whose missing lattice point makes the value zero. No floating-point
tolerance is used. The standalone source `paper/pyramidal_mobius.tex` presents
the interval representation and value formula as a short research draft with
its incomplete priority audit clearly stated.

A separate exact witness connects the structure back to prime prefixes.
Put A={2x+3y+5z<11} and B={2x+3y+5z<=11} in Z_{>=0}^3. Their difference
consists of (0,2,1), (1,3,0), (3,0,1), (4,1,0), a parallelogram. Of its 16
subsets, precisely the two opposite pairs fail strict separation, so [A,B]
has 14 elements and four atoms. Direct Mobius recurrence gives -1 at the
top. Every three atoms already have a proper upper bound, so the top is
not the join of two or three atoms. This illustrates why affine dimension
of an entire layer, rather than a bound on how many atoms generate its top,
is the useful invariant.
The full thesis comparison now confirms that this example contradicts
Conjecture 4.5 as stated on printed p. 60. That conjecture concerns arbitrary
intervals, not only intervals starting at the empty partition. No claim is
made that we first observed the poset's failure to be a lattice: that failure
was already reported in the later Elizalde--Galvan preprint.

For primes (3,5,13), A is the first 20 exponent vectors and B is the first
24. The consecutive smooth values at ranks 20--25 are

    243, 325, 351, 375, 405, 507.

The verifier checks exact heap enumeration and, independently, the products
on maximal elements of A and the minimal excluded border of B: their
maximum and minimum are respectively 243 and 507. The four layer products
lie strictly between them. The saved border certificate excludes the entire
infinite complement by coordinatewise monotonicity. This witness is not
included in the 27,100-interval simplex total above.

The interval representation gives a concrete connection between the project's
exponent geometry and an explicit question in an established literature.
The value formula is an application of an older theorem; this attribution
must remain prominent in any paper. A short research note may be justified
if the precise application has not already been recorded and the geometric
criterion is useful to researchers studying these intervals. Neither an
unsuccessful search nor the computations alone establishes that case.

## Coordinate sections versus one global cutting plane

**Proposition 19 (uniform slice obstruction).** For each integer m>=1 there
is a three-dimensional downset D_m with 4m+5 points such that all its
coordinate-plane sections are corner cuts, while D_m is not separated from
its complement by any affine hyperplane. Embedding it in higher dimensions
gives the same failure of recognition by coordinate two-dimensional sections
for every k>=3.

**Proof.** Define D_m by its z-heights, with columns x=0,1 and rows y=0,1,2:

    H_m = ((2m+2, m),
           (m+1,  1),
           (1,    0)).

All other heights are zero; a height h includes z=0,...,h-1. Decreasing
heights make a downset and their sum is 4m+5. The following table supplies
positive integer weights c and twice a threshold 2t for every nonempty
coordinate section. Remaining coordinates retain their original order.
Each section is exactly {v>=0 integer : c.v<t}.

| Fixed coordinate | Weights c | 2t |
| --- | --- | --- |
| x=0 | (2m+1,2) | 8m+5 |
| x=1 | (2m-1,2) | 4m-1 |
| y=0 | (2m+3,2) | 8m+5 |
| y=1 | (m,1) | 2m+1 |
| y=2 | (1,1) | 1 |
| z=0 | (3,2) | 11 |
| 1<=z<m | (1,1) | 3 |
| z=m | (2,1) | 3 |
| m<z<=2m+1 | (1,1) | 1 |

For example, at x=0 the inequality (2m+1)y+2z<(8m+5)/2 allows
2m+2,m+1,1 values of z at y=0,1,2 and none for y>=3. The other rows
follow by the same integer inequalities. Odd twice-thresholds exclude
lattice points on the cutting lines. All other sections are empty and
use c=(1,1), t=-1/2.

The included points a=(0,0,2m+1), d=(1,1,0) and excluded points
b=(1,0,m), c=(0,1,m+1) satisfy a+d=b+c. A strict affine inequality
including the first pair and excluding the second would give contradictory
inequalities for the same sum. Their segments meet at
(1/2,1/2,m+1/2), so this is the usual convex-hull obstruction to separation.
No new convex separation theorem is claimed.

For k>3 use D_m times {0}^{k-3}. The same four points obstruct separation.
Every coordinate two-dimensional section is empty, one of the original
sections, a finite initial axis segment, or a singleton. Each is a corner
cut; an axis segment uses a sufficiently large weight on the other axis.
This statement does not require all higher-dimensional proper sections to
be corner cuts. QED.

For m=1 the height matrix is ((4,1),(2,1),(1,0)) and the set has nine points.
No global minimality claim is made. This family refutes Conjecture 4.1 on
printed p. 51 of Galvan's thesis, with the positive-coordinate convention
translated by (1,1,1). It identifies an obstruction to combining separately
valid two-coordinate descriptions, an issue directly relevant to the original
attempt to extend two-prime interleaving. It gives no new prime estimate.

The exact artifact checks m=1,...,32, including 1,280 nonempty sections.
It verifies every table certificate against the included points and minimal
excluded border, independently solves each section's strict rational LP,
independently finds the global cutting LP infeasible, and checks the four
balanced points. Full witnesses are saved in `results/research/pyramidal_poset.json`.
The proof, rather than those finite checks, covers every parameter. Priority
relative to later literature and the value of a combined short paper remain
open; a counterexample to a thesis conjecture alone does not guarantee impact.

## A recurrence over removable exposed faces

**Proposition 20 (consequence of Proposition 18 and ordinary Mobius inversion).**
Let B be a nonempty corner cut, P=conv(B), and I(B) the number of corner cuts
contained in B, including the empty cut. Call a nonempty face F of P admissible
if some positive w exposes it as the maximum-score face on B and

    t=max_{a in B} w.a < min_{e in minimal excluded border(B)} w.e.

Put U_F=B intersect F. Then

    I(B)=1+sum_{admissible F} (-1)^dim(F) I(B minus U_F).       (J)

Each B minus U_F is a corner cut. These are precisely the lower endpoints
with nonzero proper mu(A,B). There are at most sum_{j=1}^{min(k,n)} binom(n,j)
terms, where n=|B|. For k=3 there are at most 6n; when conv(B) has dimension
three and v vertices, at most 6v-10 terms contribute.

**Proof.** Proposition 18 identifies a nonzero difference U with the entire
maximum-score face F=conv(U) of B under a positive weight vector, strictly
below the excluded border. Conversely such a face gives exactly the layer
criterion, and a threshold just below its score cuts off B minus U_F. This
is a bijection. Standard Mobius inversion gives

    sum_{A<=B} mu(A,B) I(A)
      = sum_{D<=B} sum_{D<=A<=B} mu(A,B) = 1.

Isolating A=B and substituting Proposition 18 proves (J). Each admissible
face has an affine basis of at most k points of B and equals P intersect its
affine span, proving the binomial bound. For a full-dimensional polytope
in dimension three, Euler's relation gives v+e+f=2e+2<=6v-10. When |B|>1,
P itself is not admissible: a positive functional takes distinct values on
the origin and a nonzero point of the downset. A polygon has 2v proper
nonempty faces and a segment has two. A singleton has one admissible face.
This proves the uniform 6n bound. QED.

All contributing terms can be found without enumerating the subcuts: enumerate
affinely independent subsets S of at most k coordinatewise maximal points
of B, form U=B intersect aff(S), deduplicate, and test the strict layer LP
for (B minus U,B). Every admissible face has such a basis, because its layer
is coordinatewise maximal. Affine rank and rational LP give polynomial time
in the explicit input for fixed k. This is not a result for a prefix given
only by its binary rank, nor a bound on all states reached by recursive use
of (J), nor a claimed faster overall enumeration algorithm.

This consequence connects the interval theorem to the existing subpartition
counting problem. It generalizes the two-dimensional recurrence in
[Elizalde--Galvan, preprint Lemma 7.1](https://arxiv.org/html/2312.16353v1#S7),
using removable exposed faces. Mobius inversion and polyhedral face counting
are standard; this is a consequence of the candidate application, not a
separate new inversion theorem. Its prior-art gate is the same as Proposition 18.

For B={x+y+z<=2}, there are three admissible vertex layers, three edge layers
and one two-dimensional layer. Their deletion leaves respectively 32,18,9
contained corner cuts, yielding I(B)=1+3*32-3*18+9=52. The edge and face
layers include all three or six lattice points, not only their vertices.
The verifier constructs terms from affine bases, compares them with every
nonzero proper entry of the complete Mobius columns, and checks (J) against
direct contained-cut counts for every recognized top in the saved fixtures.

## Active goal and research iteration policy

The user requested an active goal on 2026-10-16. Its acceptance condition is a
substantive proved contribution, a documented comparison with the closest
primary literature, reproducible evidence, and an honest submission-ready
paper explaining its mathematical use. A draft alone does not complete the
goal. Publication and significance cannot be promised in advance.

The primary target remains general fixed-prime rank/unrank and the geometry of
its exponent order. The occurrence-transform theory is a secondary route only
when it answers a natural mathematical question beyond a single rounded rule.

The first iteration investigates Proposition 17. Its proof mechanism is the
interaction of the prefix's downset volume and an integer vertex of its order
cone. The associated exact checks test ordering, ties, exclusion of the entire
complement via its minimal border, active and inactive coordinates, and both
determinant and cutoff bounds. They also include an infeasible ordered downset.

The next decisions are explicit:

1. Compare the complete ordered-prefix bound with corner-cut encodings,
   integer representations of monomial orders, and weight compression. If it
   is a routine known consequence, record that and do not build a paper around it.
2. If a gap survives, seek a constructive consequence that addresses an existing
   problem: for example, computing or verifying a compressed prefix description
   without enumerating all N terms. Any claimed saving must be compared with
   the generic baseline in Proposition 16, including dimension and bit complexity.
3. Continue toward a paper only if the result and its consequence together
   warrant one. Otherwise retain the lemma as background and redirect the goal.

No larger search framework or additional cycle catalogue is needed for this
iteration. Prime-distribution claims remain outside the proved results.

The second iteration makes a narrower research decision without narrowing the
goal: retain Proposition 17 as background, and prioritize Proposition 18's
answer to the specific pyramidal-poset observation. The finite interval
isomorphism and the single-layer formula have a full argument and exact
checks. Before treating a draft as submission-ready, inspect the journal
version and subsequent work, and audit Edelman's theorem
against the application. If the connection is already established, cite it
and redirect; do not present a corollary as a new arrangement theorem. No
messages have been sent to authors or journals.

The third iteration retrieved the complete thesis, checked the relevant
definitions and conjectures in context, and proved Proposition 19's infinite
family. The draft now treats the local separation question together with the
interval formula. Next assess whether either answer is already recorded in
the final journal version or later research, and whether the combined results
warrant a short note. Preserve the explicit known-theorem attribution and
the unresolved priority gate; two disproved old conjecture statements do not
alone establish a new publishable contribution in 2026.

The fourth iteration audits the arrangement argument against Edelman's
full Theorem 1.11, including its unbounded-base case, and derives Proposition
20's recurrence and bounds. It supplies an explicit counting use for the
candidate formula while retaining attribution to standard Mobius inversion.
The publisher API also proved inaccessible; it is not evidence that the
result remains open. A later-version comparison and an external assessment
of significance are still needed before submission.

On 2026-10-17 the user superseded the access blocker for ongoing work with
the provisional version assumption recorded at the start of this note.
The final-version comparison is now a deferred pre-submission check, not a
reason to halt research. The no-messages instruction rules out contacting
authors or other researchers as a way to complete that check.

## Research decision after the bounded study

The feasibility gate is passed: the restricted problem has a uniform proof,
an explicit classification and a consequence for every finite input in its
range. The elementary partition reduction and prefix lemma explain the
mechanism. Primality plays no role in this theorem.

The novelty gate is unresolved. The comparison above rules out equality with
specific cited moves, but it does not rule out an older equivalent formulation,
a conjugate system or a direct consequence of a more general theorem. No
claim of an open problem, new prime estimate or publication-ready result follows.

**Decision: pursue the general ordering problem; hold specialized cycle
catalogues and publication claims.** The primary target is a structural or
certification theorem for rank/unrank across arbitrary finite prime sets, with
explicit dependence on dimension, rank and arithmetic separation. A useful
advance must simplify or strengthen what the existing logarithm-comparison
and lattice-counting framework already supplies. Proving generic existence
again, or adding more cycle examples, does not pass this test.

The occurrence-transform strand now uses the general reconstruction, period
and transient statements as its organizing results. The three-group theorem
is a worked example. Develop this strand further only if it gives a stronger
general classification of a natural class, contributes to the ordering problem,
or survives specialist scrutiny as independently useful mathematics. More
generality in the definition alone is insufficient evidence of significance.

Proposition 15 supplies a specific constraint on the search: a uniform solution
must capture multi-coordinate arithmetic comparisons beyond any fixed stock
of pairwise patterns. Proposition 16 supplies a positive baseline and separates
fixed-dimension tractability from the growing-dimension problem. The immediate
research question is what the translation/self-similarity structure contributes
that is absent from these general methods, for example a provably smaller
certificate or a recurrence with controlled cost and precision. A claim must
name that improvement and compare against the relevant published theorem.

Novelty checks of grouped conjugation and the new structural derivations remain
necessary before a paper claim. No external messages have been sent. If the
only surviving results are elementary reformulations or known special cases,
park the paper effort while retaining the useful software and exact evidence.

The original insertion-discrepancy question remains separate:

    D_i(N)=R_i(N)-N-alpha_i N^(1-1/k).

A future prime application must improve the elementary estimate in Proposition
4 or supply another arithmetic constraint, with explicit dependence on the
number and sizes of generators. Rephrasing prime counts through smooth-number
deficits is insufficient. No work on kernels or broad discrepancy searches is
required to finish the bounded combinatorial study recorded here. The next
research stage has the broader objective above; it does not presume a prime
distribution consequence or authorize replacing the core kernels with an
unimplemented theoretical algorithm.

## Exact exploratory checks

Run from the repository root:

```bash
python3 scripts/run_interleaving_research.py --ranks 20000
```

The dependency-free script uses exact integer heap generation, checked against
an independent exact pointer generator for a short prefix. It checks full-vector
trimming, insertion commutation, the displacement identity, invariance under
common powers, the coupled zero-mask identity, arbitrary scalar zero masks,
and a degree lexicographic control. It also checks the occurrence return law
and its face-rank interpretation for every coordinate, and checks T^3=T on
all equality patterns of lengths 1 through 10. It records the finite
six-cycles for the restricted-count candidate rules above. A separate ordinary sieve
checks the first-gap and prime-count identities for the first 1 through 16
primes. The JSON output is saved at
`results/research/exponent_interleaving.json`, with the command and source state.

The bounded partition study additionally checks all 1,385,262 partitions with
at most three parts, over every mass n=1,...,b^2 and every base b=3,...,16.
It compares the geometric-column step with formula (A), traverses every basin,
and checks the observed directed cycles against the closed classification,
including their minimal periods. The saved JSON contains the full cycle list
at every mass, state counts, basin counts and maximum reduced transient.

The independent word checks cover all equality patterns of lengths 1--7 for
each base, geometric rounding thresholds and neighbours, and constant-word
witnesses at the length boundaries. There are 16,241 base-power word orbits,
with direct counting checked against labelled Ferrers-cell totals and the
partition factor, then checked against the all-prefix LCM. Another 729 words
check reconstruction, the lifting lemma, the one-step transient bound and
injective relabelling for identity, constant and parity count maps. The
base-power word cases also check reconstruction and the transient bound. Saved
word witnesses separate transient length, minimal period, whether the initial
seed is periodic, and the actual cycle states. No cycle is inferred merely
from a repeat observed before an arbitrary iteration limit.

The exact pairwise-pattern obstruction additionally checks primality of all
four fixture primes, all 10,000 exponent-vector positions of each corresponding
pair, independent pointer prefixes, and the first differing three-prime rank.
Proposition 16 is a deduction from the cited algorithms and bounds, with no
Barvinok kernel implemented or benchmarked in this repository.

The prefix-compression checks cover 230 prefixes of lengths 1--20, 32, 64 and
100 from ten generator sets, including dimensions 1--4, close primes, inactive
coordinates, composite prime powers, and both triples in Proposition 15.
Each saved witness includes its integer weights, cutoff, axis bounds, minimal
excluded border and exact LP basis. Independent additive-heap enumeration
checks the whole ordered prefix and the next score, while direct border
inequalities certify exclusion of the entire complement. An incompatible
ordered downset is rejected as infeasible. These checks support Proposition
17's proof; they do not establish novelty or provide a faster random-access
method.

These are finite consistency checks of the proved propositions. They do not
certify the fast floating-point kernels or establish a new discrepancy bound.

## Quantitative growth investigation, 2026-10-17

**Status: complete candidate argument recorded for audit; priority and
publication significance unresolved.** This section is not a claim that the
goal has been achieved. It does not assert the conjectured exact-size
three-dimensional lower bound for every n.

Write C_d(n) for the number of n-point corner cuts in the nonnegative
d-dimensional integer lattice, and A_d(n)=sum_{j=1}^n C_d(j).
The [published baseline](https://doi.org/10.1016/S0196-8858(02)00014-3) gives
C_d(n) between Omega_d(n^{d-1} log n) and O_d((n log n)^{d-1}).
Wagner's [thesis](https://doi.org/10.3929/ethz-a-004708408), Chapter 3,
Lemma 3.5, also gives C_{d+1}(n)>=A_d(n). The Chapter 3 proof and remarks
were read in the [full-text mirror](https://paperzz.com/doc/7535624/k-sets-and-applications),
with allowance for OCR errors in subscripts and inequality signs. His
[survey](https://alpergur.xyz/wagner_April13_survey.pdf), Section 4.5,
explicitly retains the gap and points to weaknesses of lifting. The 2023
Elizalde--Galvan preprint cites the growth work but supplies no improved
higher-dimensional growth theorem. This does not establish the current
strongest bound or novelty in 2026. The independently obtained bounds
credited to Remond in Wagner's abstract remain an unresolved comparison.

### Candidate 21: a multiscale cumulative lower bound

The proposed conclusions are

    A_3(n) = Theta(n^3 (log n)^2),
    C_d(n) = Omega_d(n^{d-1} (log n)^2) for every fixed d>=4.

The second statement would improve the lower logarithmic exponent in the
inspected baseline, if not already known. The first would also show that
the logarithmic exponent in its three-dimensional upper bound cannot be
uniformly reduced. A pointwise lower bound C_3(n)>=c n^2(log n)^2 remains
unproved by this argument.

**Step 1: many primitive triangles in three separated boxes.** Fix
delta=1/100. Let Q=LMH tend to infinity with
Q^{1/4}<=L,M,H<=Q^{1/2}, all integers. Put s_X=floor(delta X) and define

    B1 = [L,L+s_L] x [0,s_M] x [0,s_H],
    B2 = [0,s_L] x [M,M+s_M] x [0,s_H],
    B3 = [0,s_L] x [0,s_M] x [H,H+s_H],

intersected with the integer lattice. Each box has at least delta^3 Q
points. For a uniformly chosen triple p in B1, q in B2, r in B3, let
v=q-p, u=r-p. We claim that a fixed positive fraction satisfy

    gcd((v cross u)_1,(v cross u)_2,(v cross u)_3)=1.       (21.1)

Here is an elementary uniform sieve argument. Fix p_x. The independent
integers v_x,u_x lie in the same interval of length ell_L=s_L+1, are
nonzero, and have absolute value at most 2L. The probability they share a
prime factor is at most

    sum_{primes l<=2L} (1/l+1/ell_L)^2
      <= 3/4 + O_delta(log L/L).

The constant 3/4 follows from
sum_{j=2}^infinity 1/j^2 < 1/4+integral_2^infinity x^{-2}dx = 3/4.
The error bound uses the harmonic sum and at most 2L summands.
Thus gcd(v_x,u_x)=1 has probability at least 1/8 once Q is sufficiently
large, uniformly in p_x.

Condition on any such x coordinates and fix p_y,p_z. Set

    D_y=v_x u_y-u_x v_y,  D_z=v_x u_z-u_x v_z,
    P=floor(Q^{1/8}).

For every prime l, at least one of v_x,u_x is invertible modulo l. Counting
one of the two remaining y coordinates, and independently one of the two
z coordinates, gives

    Pr(l divides D_y and D_z)
       <= (1/l+1/ell_M)(1/l+1/ell_H).

The union bound over l<=P is at most 3/4+o(1), uniformly: the errors are
O_delta(log P*(1/M+1/H)+P/(MH)). For larger primes first fix all y
coordinates. The box positions imply D_y is positive and at most 3LM.
It has at most log(3LM)/log P<=16 distinct prime divisors larger than P
for sufficiently large Q. For each, the conditional probability of also
dividing D_z is at most 1/P+1/ell_H. Their total contribution tends
uniformly to zero. Therefore, conditionally on the coprime x pair, the
probability that gcd(D_y,D_z)=1 is at least 1/8 for sufficiently large Q.
This implies (21.1). At least 1/64 of all triples are primitive, giving
at least delta^9 Q^3/64 primitive triples. No unproved distributional
assertion about primes is used; primes only index divisibility events.

**Step 2: injectivity of the primitive triples.** A primitive triple
determines the plane a.x=h with a=v cross u and h=a.p. To see that v,u
form a basis for the lattice parallel to this plane, choose an integer z
with a.z=1. The matrix [v,u,z] has determinant 1. Every integer vector
d with a.d=0 is consequently s v+t u with integers s,t.

After scaling coordinates by L,M,H, the y,z submatrix of [v,u] differs
from the identity in maximum row-sum norm by at most 2delta. Its inverse
has norm at most (1-2delta)^{-1}. If two lattice points of the plane
belong to the same B_i, their scaled y,z coordinates differ by at most
delta each. Their integer coefficient difference (s,t) thus has norm at
most delta/(1-2delta)<1, so they coincide. Each B_i contains exactly one
point of this plane, and two different primitive ordered triples cannot
give the same plane. This proves injectivity, rather than assuming that
most triples define different planes.

**Step 3: from planes to distinct corner cuts.** The normalized row
matrix of p,q,r differs from the identity in maximum row-sum norm by at
most 3delta. If w is the solution of p.w=q.w=r.w=1, then

    |L w_1-1|, |M w_2-1|, |H w_3-1|
        <= 3delta/(1-3delta) < 1/10.

In particular all coefficients are positive, h is positive, and each
primitive triple gives a different interior vertex of the arrangement
of planes b.w=1 in the open box

    Omega = (9/(10L),11/(10L)) x
            (9/(10M),11/(10M)) x
            (9/(10H),11/(10H)).

Take the finite arrangement from all nonnegative integer b with
b_i<=ceil(10 L_i/9), where (L_1,L_2,L_3)=(L,M,H). Every other lattice
point has b.w>1 throughout Omega. Each chamber therefore corresponds
to a different finite corner cut {b:b.w<1}.

For completeness, in any dimension the number of chambers of a finite
hyperplane arrangement restricted to an open convex domain is at least
the number of its interior vertices. Add distinct hyperplanes in a fixed
order. The increase in chambers is the number of chambers of the
previous arrangement restricted to the newly added hyperplane. Assign
each interior vertex to the last hyperplane containing it. It is a
vertex of that restricted arrangement: its incident normals span the
ambient space, so two-dimensional sections in dimension three, and
the corresponding lower-dimensional sections generally, have full
rank. Induction on dimension bounds the assigned vertices by the new
chambers. The one-dimensional case counts points and intervals. This
argument allows coincident restricted hyperplanes and all other affine
dependencies; duplicate full hyperplanes are ignored.

Consequently Omega yields at least delta^9 Q^3/64 distinct cuts. Every
one is contained in the box b_i<10 L_i/9 and has fewer than 4Q points
for sufficiently large minimum side length. Its maximal coordinate on
axis i lies between (4/5)L_i and (6/5)L_i for large L_i.

**Step 4: sum over disjoint aspect ranges.** Choose r maximal with
Q=4^r<=n/4, so Q>n/16. Take

    L=4^i, M=4^j, H=4^k,
    i+j+k=r, i,j,k>=ceil(r/4).

There are binom(r-3ceil(r/4)+2,2)=Theta(r^2) choices. All sides satisfy
the conditions of Step 1. Every cut produced has size at most n. If
two scale triples differ, some side lengths differ by at least a factor
of four. The corresponding axis-maximum intervals are disjoint, so
the two collections of cuts are disjoint. Summing Step 3 proves

    A_3(n) >= c n^3 (log n)^2

for all sufficiently large n, with an absolute c>0. The matching upper
bound follows by summing the inspected published upper bound for C_3.
Wagner's lifting lemma then gives C_4(n)>=A_3(n), and iterating that
lemma gives the claimed lower bound in every fixed dimension d>=4.

### What this argument would and would not establish

If the proof survives audit, the three-dimensional upper exponent is
sharp on average. In particular, no uniform bound
C_3(n)=O(n^2(log n)^{2-epsilon}) with fixed epsilon>0 is possible.
The exact-size three-dimensional lower bound at every n remains a
separate problem. The four-dimensional consequence is uniform in n.

The proof uses standard divisibility sieving, primitive lattice bases,
arrangement counting, and the published lifting lemma. A comparison with
grid threshold-function counts, primitive-frame counts and later corner-cut
work is needed to decide whether their combination is already established
or immediate from a stronger known theorem. No novelty inference is made
from unsuccessful searches. No prime-distribution consequence or faster
rank/unrank algorithm is claimed.

The exact checker `scripts/run_corner_cut_growth_checks.py` checks 243,512
ordered triples in three finite fixtures. Their 115,286 primitive triangles
give pairwise distinct planes within each fixture and satisfy the rational
parameter-box bounds. Nonprimitive triples do collide. An independent
strict-feasibility enumeration of all 256 sign patterns of an eight-plane
arrangement finds 48 chambers and 18 interior vertices in its parameter
box. These finite checks can expose mistakes in the local construction;
they do not prove the sieve's uniform asymptotics or establish priority.

Reproduce with:

```bash
python3 -B scripts/run_corner_cut_growth_checks.py
```

Exact witnesses, commands and source hashes are saved in
`results/research/corner_cut_growth_checks.json`. Candidate 22 below extends
this construction and records the subsequent internal proof audit. A broader
comparison with primary literature remains necessary before a novelty claim.

### Candidate 22: uniform growth in every fixed dimension

Candidate 21 extends to every fixed d>=2. The proposed conclusions are

    A_d(n) = Theta_d(n^d (log n)^{d-1}),
    C_d(n) = Omega_d(n^{d-1} (log n)^{d-2}) for d>=4.

The exact-size consequence follows from Wagner's lifting lemma applied to
A_{d-1}, not from assuming a pointwise lower bound in dimension three. It
leaves one logarithmic factor between the inspected exact-size upper and
lower bounds in every dimension d>=3. The argument below gives all
quantifiers and a finite large-prime bound; novelty is a separate question.

#### Supporting lemma: primitive simplices in polynomially comparable boxes

Fix d>=2, put delta=1/(100d), and let L_1,...,L_d be positive integers with

    m=min_i L_i,  max_i L_i <= m^2.

For each i let B_i be the integer points of the box whose jth coordinate
interval is [L_j,L_j+floor(delta L_j)] if j=i, and
[0,floor(delta L_j)] otherwise. At least a fraction 8^{-(d-1)} of the
ordered choices p_i in B_i have primitive difference matrix

    V = [p_2-p_1 ... p_d-p_1],

meaning that the gcd of its (d-1)-rowed maximal minors is 1, provided
m>=2^{2S_d}. Here S_d is any integer satisfying

    S_d >= d^2+4(d-1),
    E_d(S_d) <= 1/8,

where, with a=100d and b=4d-3,

    E_d(s) = 2a(s+1)/2^{2s} + a^2/2^{3s}
             + b(1/2^s + a/2^{2s}).                    (22.1)

Such S_d exists because E_d(s) tends to zero and is decreasing for s>=1.
All the bounds in its definition are rational, so the threshold is explicit
and can be checked without numerical logarithms.

There is also a closed choice, valid for every d>=2:

    S_d=4d^2,  m>=2^{8d^2}.                            (22.1a)

Indeed, a<=2^{4d}, b<=2^{2d}, and 4d^2+1<=2^{3d}. Each inequality
holds at d=2 and is preserved on increasing d by one. Therefore

    E_d(4d^2) <= 2^{1+7d-8d^2} + 2^{8d-12d^2}
                + 2^{2d-4d^2} + 2^{6d-8d^2}.

The four exponents decrease for integers d>=2. At d=2 their values are
-17,-32,-12,-20, so the displayed sum is less than 2^{-10}<1/8.
Also 4d^2>=d^2+4(d-1). The threshold is deliberately coarse; no finite
search over dimensions is needed for this choice.

**Proof.** Use induction on the number r of boxes, 1<=r<=d, keeping the
same delta throughout. For r=1 the difference matrix has no columns and
its single maximal minor is the empty determinant 1. At each subsequent
step condition on any primitive projected configuration of the first
r-1 points in the first r-1 coordinates. Write U for its
(r-1)-by-(r-2) difference matrix and a_0 for its cofactor normal; a_0 is
primitive. Coordinates not in this projection remain independent and
uniform on their original intervals.

Let u be the projected difference from the first point to the new point
p_r, and D=det[U u]. This determinant is nonzero for every allowed choice.
Indeed, the normalized row matrix of the first r-1 projected points is
within (r-1)delta of the identity in maximum row-sum norm. Its cutting
normal is positive and within (r-1)delta/(1-(r-1)delta) of the vector of
ones. The normalized new point has every coordinate at most delta,
and lies strictly on the origin side of their affine span. Also,

    |D| <= (r-1)! 2^{r-1} (max_i L_i)^{r-1}
         <= 2^{d^2} m^{2(d-1)}.                       (22.2)

For any prime l, U has rank r-2 modulo l. The congruence D=0 modulo l
is a nontrivial linear equation in the independent projected coordinates
of p_r, since a_0 is primitive. If ell is the shortest coordinate-interval
length, its probability is at most 1/l+1/ell. Conditional on D=0, u lies
in the column span of U modulo l. There is a unique vector lambda with
u=U lambda. If z is the last row of the old columns and u_r is the last
entry of the new column, the full matrix fails to have rank r-1 modulo l
exactly when

    u_r = z lambda (mod l).                            (22.3)

Conditional on everything except the last coordinate of p_r, this selects
one residue class. Thus its conditional probability is at most 1/l+1/ell_r.
Every interval length is larger than delta m.

Put s=floor((log_2 m)/2) and P=2^s, so
2^{2s}<=m<2^{2s+2}. The union bound for primes l<=P is at most

    sum_{l prime, l<=P} (1/l+a/m)^2
      <= 3/4 + 2a(s+1)/2^{2s} + a^2/2^{3s}.          (22.4)

This uses sum_{j>=2} j^{-2}<3/4 and H_P<=1+log P<=1+s.
For l>P, first fix all projected coordinates. Every bad prime divides
the fixed nonzero integer D. From (22.2), the number of its distinct
prime divisors larger than P is at most

    log_2 |D| / s
      <= 4(d-1)+(d^2+4(d-1))/s <= b.

For each such divisor, (22.3) still bounds the conditional probability
by 1/P+a/m. The total large-prime contribution is at most
b(2^{-s}+a 2^{-2s}). This is a finite bound before taking any limit;
it does not interchange an infinite union with natural density.

Equations (22.1)--(22.4) bound the conditional probability of a
nonprimitive extension by 3/4+E_d(s)<=7/8 whenever s>=S_d.
The previous projected sides still satisfy max L_i<=(min L_i)^2,
so the induction is uniform for all projected configurations. Each of
the r-1 extensions succeeds with conditional probability at least 1/8.
More explicitly, let E_r be the event that the first r points projected
onto the first r coordinates form a primitive frame. Conditioning on
E_1,...,E_{r-1} imposes conditions only on the previous projected coordinates.
After fixing those coordinates, the rth old coordinates and the new point
still have their product distributions. Thus
Pr(E_r | E_1,...,E_{r-1})>=1/8, and
Pr(E_1 intersect ... intersect E_d)>=8^{-(d-1)}. This event is contained
in E_d. No independence of the events E_r themselves is assumed.
This proves the stated fraction. QED.

The primitivity of the projected configuration is essential: without it,
the top block can lose more than one rank modulo a prime, and varying the
last coordinate need not repair it. The checker includes an explicit
four-dimensional negative control with this failure.

#### Geometric injection and the cumulative count

For one primitive choice let a be the cofactor normal of V. There is an
integer z with a.z=1, so [V z] has determinant 1. Consequently the
columns of V form a basis of the full parallel integer lattice
{x in Z^d:a.x=0}; this is the usual primitive-lattice criterion.
The submatrix of V in coordinates 2,...,d, after coordinate scaling,
differs from the identity by at most (d-1)delta in maximum row-sum norm.
Two points of the same B_i on the affine span would have an integer
coefficient difference of norm at most

    delta/(1-(d-1)delta) < 1.

They therefore coincide. The supporting hyperplane determines the unique
point in each B_i and hence determines the entire selected simplex.

The normalized row matrix of the d points differs from the identity by
at most d delta. The unique solution of p_i.w=1 is positive and satisfies

    |L_i w_i-1| <= d delta/(1-d delta) < 1/10.

Thus every selected primitive simplex gives a different interior vertex
of the arrangement b.w=1 in

    Omega_L = product_i (9/(10L_i),11/(10L_i)).

For completeness, any finite affine-hyperplane arrangement in an open
convex domain has at least as many chambers as interior vertices. Delete
duplicate hyperplanes and add them in a fixed order. Adding H creates
exactly as many new chambers as there are chambers of the restricted
arrangement in H intersected with the domain: each old chamber is convex
and is split if and only if it meets H. Assign each vertex to its last
incident hyperplane H. The normals of the incident hyperplanes span the
ambient space; the nonzero normal of H extends to a basis among them.
Consequently the assigned vertex is an interior vertex of the restriction
to H. Induction on dimension bounds the number assigned to H by the
number of new chambers. In dimension one, k distinct interior points
give k+1 chambers. Summing proves the assertion, including all degeneracies.

Only points with b_i<=ceil(10L_i/9) can change membership in Omega_L;
the arrangement is finite. Every chamber gives a different corner cut
of size at most (3/2)^d Q for sufficiently large minimum side, where
Q=product_i L_i. Each B_i has at least delta^d Q points, so the lemma
and injection give at least

    8^{-(d-1)} delta^{d^2} Q^d                         (22.5)

distinct cuts in this parameter box. The constants depend only on d.
In the cardinality bound, each coordinate has fewer than 10L_i/9+1
possible values. This is at most (3/2)L_i for L_i>=3, an explicit
condition already implied by (22.1a).

Choose R maximal with Q=4^R<=n/2^d. For each tuple of integers

    r_1+...+r_d=R,  r_i>=ceil(R/(d+1)),

set L_i=4^{r_i}. These sides lie between Q^{1/(d+1)} and
Q^{2/(d+1)}, so their maximum is at most their minimum squared and
the lemma applies for all sufficiently large n. There are

    binom(R-d ceil(R/(d+1))+d-1,d-1)
        = Theta_d((log n)^{d-1})

such tuples. All cuts have at most n points. Their maximal coordinate
on axis i lies between (4/5)L_i and (6/5)L_i for large L_i, so no cut
is counted in two different scale tuples. Since Q>n/(4*2^d), summing
(22.5) gives A_d(n)=Omega_d(n^d(log n)^{d-1}). Summing Wagner's
published upper bound supplies the matching upper bound.

#### A fully explicit lower bound

Set

    R_0(d)=4d^2(d+1),  n_0(d)=2^d 4^{R_0(d)},
    kappa_d = 8^{-(d-1)} (100d)^{-d^2}
              / ((d-1)! [4(d+1)]^{d-1} 2^{d(d+2)}).

Then the construction proves the concrete inequality

    A_d(n) >= kappa_d n^d (log_4 n)^{d-1}
    for every integer d>=2 and n>=n_0(d).              (22.6)

Here is the constant calculation. With R maximal as above, R>=R_0(d),
so every side is at least 4^{R/(d+1)}>=2^{8d^2}. Write
h=R-d ceil(R/(d+1)). Since R>=2d(d+1),

    h>=R/(d+1)-d>=R/[2(d+1)].

There are binom(h+d-1,d-1)>=h^{d-1}/(d-1)! scale tuples.
Moreover R>=log_4 n-d/2-1>=(log_4 n)/2 at this threshold, and
Q>n/2^{d+2}. Substitute these three inequalities into (22.5) to
obtain (22.6). No decimal approximations enter this calculation.
The constants make the quantifiers checkable; their numerical size is
not a claim of an effective small-input enumeration method.

#### Consequences and the boundary of the claim

Wagner's lifting lemma gives the displayed exact-size lower bound in
dimension d>=4. For d=3 it recovers only the previously known
Omega(n^2 log n). The exact-size three-dimensional gap is still open
within this investigation. In every fixed dimension, however, an upper
bound O(n^{d-1}(log n)^{d-1-epsilon}) for a fixed epsilon>0 would
contradict the cumulative lower bound.

More precisely, the cumulative result and the pointwise published upper
bound imply that there are c_d,eta_d>0 such that, for every large N,
at least eta_d N integers j<=N satisfy

    C_d(j) >= c_d j^{d-1}(log j)^{d-1}.

To see this, write A_d(N)>=a_d N^d(log N)^{d-1} and
C_d(j)<=b_d j^{d-1}(log j)^{d-1} for large arguments, absorbing the
finitely many smaller ones into an error. If c_d<a_d/4, the terms
below this proposed threshold contribute at most
c_d N^d(log N)^{d-1}. Every remaining term is at most
b_d N^{d-1}(log N)^{d-1}, forcing a positive fraction of them.
This remains a statement about many sizes, not every size.

The new exact checks cover dimensions 2--4 exhaustively in stated
finite boxes and dimensions 5--6 by deterministic samples. Independent
modular Gaussian elimination checks the determinant criterion and the
one-residue-class extension condition. The values of S_d and (22.1)
for d=2,...,12 are checked with rational arithmetic. These checks support
the argument; the proof above, not sampling, supplies its uniformity.

There is a direct interpretation in an established application. For a
cut D let alpha_D=sum_{x in D} x. For fixed size n, different cuts have
different alpha_D, and each alpha_D is a vertex of their convex hull.
To prove both assertions, choose positive weights strictly cutting off D.
Every removed point of D has smaller weight than every replacement
point outside D. Hence w.alpha_D<w.alpha_E for every other n-point
set E, in particular every other corner cut. Thus alpha_D is uniquely
exposed and cannot equal alpha_E. Conversely, every vertex of this finite
convex hull is one of its generating points. Its vertex count is exactly
C_d(n).

This is the standard corner-cut-polytope dictionary, not a new theorem.
[Berstein et al., 2010](https://doi.org/10.1007/s10463-010-0291-8),
[free journal PDF](https://www.ism.ac.jp/editsec/aism/pdf/062_4_0673.pdf),
Definition 7 and Theorem 2, identify these polytopes and the algebraic
fans of generic experimental designs. The candidate bound would therefore
give Omega_d(n^{d-1}(log n)^{d-2}) vertices and models for fixed d>=4.
This quantifies an existing enumeration problem. It neither gives a
faster enumeration algorithm nor changes the meaning of generic design
used in that source.

#### Prior-art audit added for the general proof

Primitive matrix density and the lattice criterion are established
background. [Maze--Rosenthal--Wagner, 2011](https://doi.org/10.1016/j.laa.2010.11.015)
give the density product; [Micheli--Weger, 2019](https://doi.org/10.1137/18M1177093),
[free text](https://arxiv.org/pdf/1803.08785), explain a defect in the earlier
proof and prove the density result using a controlled sieve. Their
Propositions 3--4 give the standard minor and modular-rank criteria;
Corollary 17 gives the rational-integer density. We claim neither a new
matrix criterion nor a new unrestricted density formula. A density in
centered expanding cubes alone does not supply the uniform estimate in
our moving anisotropic boxes; the finite argument above supplies that
estimate directly. Whether a stronger published counting theorem already
implies the complete multiscale result remains to be checked.

The [Zamaraeva--Zunic threshold-function paper](https://arxiv.org/html/2007.03986v3)
uses primitive segments and injectivity in dimension two. Its stated
counting results concern finite planar grids and intersections of two
halfplanes, rather than corner cuts of bounded cardinality in arbitrary
dimension. This establishes a difference in definitions, not novelty of
our counting conclusion. The broader grid-threshold and lattice-incidence
literature remains part of the priority audit.

The independent Remond result is acknowledged as a communication in
Wagner's thesis, immediately after (3.2); that passage does not cite a
separate multidimensional proof. Its acknowledged bounds are the same
baseline. An unavailable independent proof must not be silently treated
as inspected, but it also need not prevent analysis of the stated bounds.

[Balko--Cibulka--Valtr](https://arxiv.org/pdf/1703.04767), Theorems 3--5,
study lattice coverage by subspaces, with bounds in terms of successive
minima. These cover-number estimates do not themselves count cuts by
included cardinality or assert the moving positive-intercept estimate
proved here. Their incidence applications are relevant background, not
an identified prior proof of (22.6).

#### Internal adversarial proof audit, 2026-10-17

This is an internal audit by the research agent, not an independent human
referee report. The following checks pass for the written argument:

| Step | Audit question | Resolution |
| --- | --- | --- |
| Dimension induction | Does conditioning destroy the uniform coordinates needed for the next step? | The events through r-1 depend only on old projected coordinates. Fixing those leaves the new point and old rth coordinates uniform; the joint-event calculation is explicit above. |
| Controlling determinant | Could the finite large-prime bound use D=0? | Near-identity projected rows give a positive cutting normal; the new projected point is strictly below their affine span. Thus every allowed D is nonzero. |
| Modular extension | Can failure impose more than one residue on the new last coordinate? | The old projected frame has full column rank modulo every prime. Its unique lambda gives exactly (22.3). The nonprimitive negative control demonstrates why the hypothesis matters. |
| Sieve tails | Was a limiting density exchanged with an infinite union? | No. Small primes use a finite union. All remaining bad primes divide one fixed nonzero D, with at most b large factors, before the last coordinate is averaged. |
| Uniformity in d | Is the threshold specified only by computation? | The closed choice S_d=4d^2 follows from four elementary exponential bounds for all d>=2. |
| Duplicate planes | Do distinct primitive simplices ever produce the same plane? | A primitive frame is the full parallel integer-lattice basis. The scaled inverse bound forces each plane to meet each B_i in only its selected point. |
| Chambers | Do dependencies or domain boundaries invalidate vertex counting? | All selected vertices are interior. The last-incident-hyperplane assignment and dimension induction work for arbitrary degeneracies in an open convex domain. |
| Duplicate cuts | Could disjoint parameter boxes still represent the same cut? | Different scale tuples differ by a factor of at least four in one L_i, while each cut's maximum axis coordinate is in (4L_i/5,6L_i/5). Those observable ranges are disjoint. |
| Cardinality and constants | Are all constructed cuts within the stated n and threshold? | The coordinate bound gives size at most (3/2)^d Q<n; (22.6) records every scale, floor and threshold loss explicitly. |
| Exact-size inference | Does the argument prove the missing dimension-three pointwise lower bound? | No. Lifting yields the new displayed exact-size bound only for d>=4. The positive-proportion statement is the precise averaged consequence. |

The internal proof pass is complete, with no unresolved logical step
identified in this argument. That is not a guarantee against an undiscovered
error. The main research objective is still incomplete: the counting theorem
needs a defensible priority comparison and a manuscript whose contribution
survives that comparison. The available sources support comparison with the
inspected baseline, not a claim that the result is the strongest known in 2026.

### Quantitative manuscript and additional primary comparisons

The separate working manuscript is `paper/corner_cut_growth/main.tex`,
*Cumulative growth of corner cuts in fixed dimension*. It contains the full
uniform proof, explicit constants, the arrangement lemma with degeneracies,
and a self-contained proof of Wagner's lifting step. The structural note in
`paper/pyramidal_mobius.tex` remains separate. The quantitative manuscript is
a research draft; its existence does not establish novelty or submission
readiness, and the priority-status paragraph states the remaining checks.

Two further primary comparisons clarify the boundary of the contribution.
[Babson--Onn--Thomas](https://arxiv.org/pdf/math/0207135), Proposition 3.2
and Theorem 3.5, give an upper vertex bound for the Hilbert zonotope and its
refinement of state polyhedra. Section 5 identifies the generic-point case
with corner cuts. This is the known universal-refinement and algorithmic
framework; an upper count of a refinement supplies no lower count of its
coarser corner-cut fan. We do not claim their framework as new.

[Chen--Feng--Liu--Wu, 2024](https://doi.org/10.1007/s11424-024-1313-6),
[inspected preprint v2](https://arxiv.org/pdf/2105.05383v2), Theorem 1 and
Sections 1.3 and 2, give finite lower probability bounds when independent
rows from one common interval are appended to a prescribed primitive matrix.
The interval length bounds the prescribed entries. Our free coordinates
are in translated, unequal intervals whose shortest length can be much
smaller than the largest entry, so those stated hypotheses do not directly
give our window estimate. Prime-wise rank analysis and determinant bounds
are known techniques; the multiscale bounded-cardinality count remains the
contribution to assess. No new general primitive-matrix density theorem is
claimed.

The [Ghilezan--Pantovic--Zunic encoding paper](https://doi.org/10.1109/TNN.2007.891678)
was checked at the authored-abstract level. It studies moment codes for
partitions of a finite cube by parallel hyperplanes, with bit rate depending
on the cube side, dimension and number of planes. A fixed-cube encoding
statement by itself neither controls included cardinality nor sums disjoint
aspect windows. The full-text comparison is still not complete; it is
explicitly identified as such in the draft rather than represented as a
finished priority check.

### Additional attribution and the exact prime-base connection

[Müller, *Corner cuts and their polytopes*, 2003](https://ftp.gwdg.de/pub/misc/EMIS/journals/BAG/vol.44/no.2/b44h2mue.pdf),
printed pp. 323--333, was inspected in full. Lemma 2.6 explicitly proves
that the moment vectors of corner cuts are exposed vertices, so the
polytope dictionary above must also be attributed to that source. Section 1,
p. 324, directs enumeration to the planar result and Wagner's bounds.
Theorems 3.2--3.3 give a planar vertex order; Theorems 4.9, 5.2 and 5.6
concern pointedness and adjacency. None of those statements is a
bounded-cardinality asymptotic count. This resolves that specific comparison;
it does not settle all priority questions.

**Supporting proposition (prime realization, no separate novelty claim).**
Every nonempty finite corner cut D in dimension d is the exponent set of a
smooth-number prefix for some ordered tuple of distinct primes. The primes
can all exceed any prescribed B. In particular, there are infinitely many
realizing tuples.

**Proof.** Write D={x:w.x<1} with strict separation. Set

    g=1-max_{x in D} w.x>0,  H=max_{x in D} sum_i x_i.

Temporarily relabel coordinates so w_1<=...<=w_d. Choose T with
exp(T w_1)>max(B,2) and d H log(2)/T<g. Bertrand's postulate supplies
distinct primes in the disjoint intervals

    2^{i-1} exp(T w_i) < p_i < 2^i exp(T w_i).

For real endpoints apply the integer statement to their floors. Equal
weights cause no problem because adjacent intervals remain disjoint.
With v_i=log(p_i)/T, we have 0<v_i-w_i<d log(2)/T. Thus every inside
point has v.x<1, while every outside point has v.x>=w.x>1. This latter
inequality controls the entire infinite complement, not just a finite box.
The smooth products below exp(T) therefore have exponent set exactly D.
Unique factorization makes distinct exponent vectors give distinct values.
Undo the coordinate relabeling. Since B was arbitrary, infinitely many
tuples exist. For the classical prime theorem, see
[Erdős, 1932](https://users.renyi.hu/~p_erdos/1932-01.pdf), pp. 197--198.

Let S_{d,B}(n) count distinct length-n exponent-prefix sets obtained over
ordered tuples of distinct primes greater than B. Count shapes once, not
their realizing tuples. The proposition and the elementary converse give

    S_{d,B}(n)=C_d(n),
    sum_{j=1}^n S_{d,B}(j)=Theta_d(n^d (log n)^{d-1}).

Consequently Candidate 22 counts precisely the finite shapes available to
the original smooth-number project when the prime bases vary. For one fixed
prime tuple there is only one prefix per length. This does not give a new
prime-distribution theorem or improve random-access complexity. The
realization step is elementary background supporting the quantitative
contribution, not a substitute for it.

Six exact fixtures in `scripts/run_corner_cut_growth_checks.py` include
equal weights, permuted labels, unused coordinates and the one-point cut.
For integer weights a and score bound s, set b=2^{dH}+1 and choose primes
in (2^{i-1}b^{a_i},2^i b^{a_i}) after sorting weights. Every inside product
is at most 2^{dH}b^s<b^{s+1}, and every outside product is at least
b^{s+1}; actual interval choices make the separator strict. The checker
verifies primality by trial division, exact integer boundary products,
and an independently generated heap prefix. These fixtures falsify finite
mistakes; the proposition's proof gives the general statement.

### Candidate 23: prescribed size through primitive normal directions

**Status: new complete candidate proof, with exact falsification checks;
priority and submission review are incomplete.** This extends Candidate 22
in the same dimension, without Wagner lifting. If the argument survives
further review, it gives, for every fixed d>=2,

    C_d(n)=Theta_d(n^{d-1}(log n)^{d-1}).                 (23.1)

In particular, the proposed dimension-three conclusion is the pointwise
bound C_3(n)=Theta(n^2(log n)^2), not an average or a subsequence statement.
It matches Wagner's inspected upper bound. No assertion that this is new
over all subsequent literature is made yet.

#### 1. Count directions, rather than heights

Use Candidate 22's axis boxes and primitive frames. Orient the primitive
cofactor normal a positively, and put H_a=a.p_1. If P is the matrix of the
point rows, then H_a=det(P)>0. After column scaling, P differs from the
identity by a matrix of maximum row-sum norm at most d delta=1/100.
Its eigenvalues therefore lie within 1/100 of one; its determinant is
positive, and the following coarse bounds hold:

    2^{-d} Q <= H_a <= 2^d Q,
    9/10 < a_i L_i/H_a < 11/10.                         (23.2)

The normal bounds were already proved in Candidate 22. A primitive frame
V=[p_2-p_1 ... p_d-p_1] is a basis of the entire integer kernel of a.
Different selected simplices give different planes (a,H_a). At most
2^d Q positive integer heights accompany one fixed primitive normal.
Thus the rho_d Q^d planes of Candidate 22 supply at least

    2^{-d} rho_d Q^{d-1} distinct primitive normals,
    rho_d=8^{-(d-1)}(100d)^{-d^2}.                      (23.3)

#### 2. Define an n-point cut for each integer normal

Order the lattice points first by a.x, then lexicographically, and let
D_n(a) be the first n points. This is a corner cut. To verify that assertion
without treating tied points as a strict cut, let h be the least integer
with F_a(h):=#{x>=0:a.x<=h}>=n. All points at score at most h have each
coordinate at most n-1: the n points on a minimum-weight axis already
give h<=(n-1)min_i a_i. The lexicographic tie order on this finite box is
represented by an integer vector u with digits in base n, reversing the
coordinate significance as needed. Taking w=M a+u, with M greater than
max_{a.x<=h} u.x, preserves the a-levels up to h and strictly orders the
boundary level. Positivity prevents any higher a-level from entering.
A threshold between successive w-values gives D_n(a). For n=1 the set
{0} is immediate.

In particular, all points with a.x<h belong to D_n(a), none with a.x>h
do, and a is a weak separator:

    a.x <= a.y for every x in D_n(a), y outside D_n(a).  (23.4)

#### 3. Two integer patches force injectivity at a fixed size

Set K=4^{d+1}. Suppose a and b arise in the same aspect box, and
D_n(a)=D_n(b). We prove a=b whenever the boundary height h for a obeys

    h/H_a >= 16d^2(K+1).                               (23.5)

Since a is primitive, choose an integer z with a.z=1. Reduce it modulo
the kernel basis V so that

    z=p_1/H_a+V c,  |c_j|<=1/2.                        (23.6)

This reduction is possible because z-p_1/H_a lies in the real kernel.
Each row of V has absolute sum less than 2(d-1)L_i, hence
|z_i|<=(d+1)L_i<=2dL_i.

Let y_i=(h-1)/(d a_i), so a.y=h-1. Reduce the integer point (h-1)z
modulo V about y to obtain an integer x_0 with

    a.x_0=h-1,  x_0-y=V beta,  |beta_j|<=1/2.

Consider all integer points in the two patches

    x_0+V k,  x_0+2z+V k,   k_j in {-K,...,K}.          (23.7)

They are nonnegative. Indeed x_0 differs from y in coordinate i by
less than dL_i, V k has absolute coordinate at most 2KdL_i, and
2z has absolute coordinate at most 4dL_i. Meanwhile (23.2) and (23.5)
give

    y_i/L_i > (h-1)/(2dH_a) > 2Kd+5d.

The first patch has a-score h-1 and is entirely inside D_n(a); the
second has score h+1 and is entirely outside, regardless of how the
score-h tie was resolved. Put S=sum_j |b.v_j|, where v_j are the
columns of V. Weak separation by b and the extrema on the two patches
give

    K S <= b.z <= b.p_1/H_a + S/2.                     (23.8)

For b in the same aspect family, (23.2) yields H_b/H_a<=4^d and

    b.p_1/H_a < 2*4^d,

because sum_i (p_1)_i/L_i<=1+d delta. Hence

    (K-1/2) S < 2*4^d < K-1/2.

But S is a nonnegative integer. Thus S=0, b annihilates every column
of V, and b is parallel to a. Both normals are primitive and positive,
so a=b. This is the new fixed-size injectivity step. It never assumes
that C_d(n) is monotone or converts a cumulative count by averaging.

#### 4. Enough height and observable aspect ranges

Let

    D_d=[64d^2(K+1)]^d,
    G=(d! n/Q)^{1/d},

and take n/Q>=D_d, m>=2^{8d^2}. Covering lattice points by their
half-open unit cubes gives the elementary volume bounds, for t>=0,

    t^d/(d! product_i a_i) <= F_a(t)
        <= (t+sum_i a_i)^d/(d! product_i a_i).          (23.9)

Use F_a(h-1)<n<=F_a(h), together with (23.2), to obtain

    (9/10)G-(11/10)d/m < h/H_a
        < (11/10)G+1/H_a.

In particular (4/5)G<h/H_a<(6/5)G. The chosen D_d makes the lower
bound larger than 16d^2(K+1), proving (23.5).
If M_i is D_n(a)'s maximum axis coordinate, then
h/a_i-2<M_i<=h/a_i, so

    G L_i/2 < M_i < 2G L_i.                            (23.10)

G is common to every aspect tuple having the same Q. Consequently
scale tuples differing by a factor of at least four in an L_i produce
disjoint observable ranges in (23.10), even after padding to size n.

#### 5. Uniform pointwise count

Put R_0=4d^2(d+1) and n_star(d)=D_d 4^{R_0}. For n>=n_star(d),
take R maximal with Q=4^R<=n/D_d. As in Candidate 22, take all tuples
sum_i r_i=R with r_i>=ceil(R/(d+1)), and L_i=4^{r_i}. They satisfy
the primitive-window hypotheses and number at least

    R^{d-1}/([2(d+1)]^{d-1}(d-1)!).

Also Q>n/(4D_d) and R>=(log_4 n)/2. For the latter closed inequality,
d^2<=4^{d-1} gives D_d<=4^{2d^2+4d}, while
R_0>=2d^2+4d+2. Combining (23.3), fixed-size injectivity, and
disjointness (23.10) proves

    C_d(n) >= lambda_d n^{d-1}(log_4 n)^{d-1},
    lambda_d = rho_d /
      [2^d (4D_d)^{d-1} [4(d+1)]^{d-1}(d-1)!].         (23.11)

Wagner's upper bound then gives (23.1). All constants depend only on
d; they are deliberately coarse, not practical enumeration thresholds.

The exact checker constructs Bezout increments, reduces both affine cosets
with rational Gaussian elimination, and checks nonnegative patch corners
and reversed b-orders for 2,280 ordered pairs of distinct normals in
dimensions 2--6. Witnesses contain the actual two lattice points and both
normal scores. These are checks of the fragile injection step, not a
computational proof of the uniform count. The source proof also needs a
fresh independent logical review before the manuscript can be regarded as
ready; the literature gate remains open.


#### Additional primary citation checks

[Kalka--Kuhnt, arXiv:2005.04051v2](https://arxiv.org/pdf/2005.04051),
Section 2.2, printed p. 6, explicitly repeats Wagner's two orders for the
number of corner cuts. The paper studies numerical statistical fans;
that background passage supplies no improved corner-cut lower count.
[Maruri-Aguilar--Wynn, arXiv:1207.2968](https://arxiv.org/pdf/1207.2968),
Section 6.2, printed pp. 21--22, describes the known state-polytope and
corner-cut-model interpretation and refers cardinality to the planar
paper and Wagner. These are inspected passages, not assertions of a
complete audit of all results in those papers. A public citation-index
listing was used only to locate sources, not as proof of novelty.


#### Fresh internal audit of Candidate 23

The prescribed-size argument was re-derived from the definitions after it
was written, without assuming that the cumulative result implies it.
No unresolved logical gap was identified. This is an internal review by
the research agent, not an independent human referee report.

| Step | Adversarial check | Resolution |
| --- | --- | --- |
| Primitive windows | Are unobserved coordinates still uniform after earlier rank events? | All preceding events depend only on old projected coordinates; new coordinates retain their independent distributions. The finite small-prime and large-prime bounds are uniform in each fixed configuration. |
| Counting planes | Can a primitive plane meet a selected box twice? | Its frame is the full integer kernel basis, and the scaled inverse bound forces the coordinate difference to have an integer coefficient vector of norm less than one. |
| Counting directions | Does using H_a hide a multiplicity? | Divide the plane count by at most 2^d Q integer heights per primitive normal, then fix one witnessing simplex. D_n(a) depends on the normal alone; its auxiliary height need not be unique. The draft now states this explicitly. |
| Partial ties | Is a lexicographic prefix always a strict corner cut? | The boundary score is at most (n-1) min_i a_i. Every relevant exponent is a base-n digit vector, so w=n^d a+(n^{d-1},...,1) strictly implements the prefix, including a partial boundary level. |
| Coset reduction | Can a unit-level integer increment or centered lattice point fail to exist? | Primitivity gives a.z=1. Subtracting integer combinations of V centers any real coefficient vector without changing integrality or its level. The construction works for every integer boundary height. |
| Patch signs | Can the centered construction produce negative coordinates? | The explicit condition h/H_a>=16d^2(K+1) exceeds the sum of all centered, patch and increment displacements in every coordinate. |
| Separation | Does a partial tie affect patch inclusion? | Patches have scores h-1 and h+1. All points of the first lie inside the prefix and all points of the second lie outside, for every tie choice at h. |
| Integer rigidity | Does the second separator require generic or real weights? | Both counted normals are integer and primitive. Equal prefixes imply K S<=b.z<=b.p_1/H_a+S/2, where S=sum |b.v_j| is integer and b.p_1/H_a<2*4^d<K-1/2. Hence S=0 and the normals agree. |
| Prescribed n | Could a boundary-height estimate only hold for selected ranks? | F_a(h-1)<n<=F_a(h) holds for every rank. Unit-cube volume bounds and n/Q>=T_d give the required height uniformly. No monotonicity of C_d(n) is assumed. |
| Aspect duplication | Can padding to size n merge different aspect tuples? | All tuples have the same Q and G=(d!n/Q)^{1/d}; the observable axis ranges (G L_i/2,2 G L_i) are disjoint when one L_i changes by a factor of four. |
| Quantifiers | Do floor operations or the large constants restrict n to powers? | R=floor(log_4(n/T_d)) is valid for every integer n>=n_star. The composition and logarithm bounds were checked algebraically with the displayed closed thresholds. |

The original patch fixtures used h equal to an integer multiple of H_a,
which samples only one boundary residue. They now use offsets
0, 1, floor(H_a/2), H_a-1 and H_a+17. All 11,400 ordered-normal-pair
and height combinations and 7,290 pairs of patch corners pass in dimensions
2--6. Actual integer points, residues and reversed scores are retained.
The 200 direct tie-resolution checks include 150 partial boundary levels
and 50 complete levels. These computations target finite failure modes;
the uniform proof is the argument above, and priority remains unestablished.

#### Fixed-cardinality comparison and the average-completion consequence

The finite-grid comparison can be made exact, rather than resting on a
difference in terminology. Put G_n={0,...,n}^d. The n-point corner cuts
are exactly the n-sets of G_n containing the origin. Here an n-set means
a subset strictly separable from its finite-grid complement. The forward
implication follows because an n-point order ideal has coordinates at most
n-1. Conversely, write a strict separator as w.x<t. Inclusion of the
origin implies t>0. If w_i<=0, all n+1 axis points 0,e_i,...,n e_i would
be included. Thus every w_i>0. Each n e_i must be excluded, since its
inclusion would again include n+1 axis points. Therefore n w_i>t, and
every nonnegative lattice point outside G_n also has score greater than t.
This proves the converse. Wagner's thesis Section 3.2 uses a smaller
hyperbolic-cross ground set; finite reduction is established background.

[Haukkanen--Merikoski, 2013](https://doi.org/10.1016/j.dam.2012.08.012),
[inspected v2](https://arxiv.org/pdf/1110.3566v2), Section 1, equations
(1)--(2), and Section 2, equation (9), count all threshold functions on
an m-by-r rectangle. Their exact primitive-segment formula sums over
all included cardinalities; their asymptotic is
t(m,r)=6(mr)^2/pi^2+O(m r^2) for 2<=m<=r. In particular the total on
G_n in dimension two has order n^4, whereas the selected origin class
C_2(n) has order n log n. A lower count for the total does not bound
this particular cardinality class from below. This resolves the direct
comparison with that theorem, without claiming to exclude every theorem
in threshold-function research. Sections 1--2 and 5, including the theorem
statements and their parameter definitions, were inspected; the full
Section 4 error-term proof was not re-audited.

The related Farey-diagram lead needs an explicit correction in our
literature record. Khoshnoudirad's
[2015 article](https://doi.org/10.13069/jacodesmath.79416), Definitions
2.8--2.9, counts integer height patterns on a fixed rectangular domain,
not n-point orthant cuts. More significantly, its proposed improved
upper bound must not be treated as established: the
[2016 published erratum](https://doi.org/10.13069/jacodesmath.00924),
[journal PDF](https://dergipark.org.tr/en/download/article-file/267848),
Sections 4.4--5, printed p. 121, leaves an additional term uncontrolled
and states that this method does not improve the known Farey-vertex
upper bound. We inspected the definitions and stated bounds of the
2015 text, and the abstract, definitions and corrected conclusion of the
2016 text. Neither corrected statement supplies a prescribed-size
corner-cut lower bound. We do not use the withdrawn upper estimate.
The downloaded full text remains local and excluded from the public
repository.

There is also a direct consequence for the question in
[Wagner's thesis, Section 3.5](https://doi.org/10.3929/ethz-a-004708408),
under "Weakness of the lower bound". For every nonempty (d-1)-dimensional
corner cut T of size at most n, let F_n(T) count its n-point completions
in dimension d along the coordinate slice x_d=0. Every n-point cut has
one such nonempty slice, giving the established identity

    sum_{1<=|T|<=n} F_n(T) = C_d(n).

The number of possible slices is A_{d-1}(n), so their mean number of
completions is C_d(n)/A_{d-1}(n). Candidate 23 therefore gives

    mean F_n(T) = Theta_d(log n)  for every fixed d>=2.

For d>=3 use the matching count in dimensions d and d-1 and sum the
latter; for d=2 use A_1(n)=n. This is precisely the mean logarithmic
growth suggested in Wagner's discussion. The restriction, fibre identity
and lifting viewpoint are not new. This consequence explains what the
new matching count would settle; it is not a substitute theorem or an
independent novelty claim. It does not give logarithmically many
completions for each T: if |T|=n then F_n(T)=1.

The manuscript now includes this corollary and the exact finite-grid
comparison. Its main theorem and proof constants are unchanged. These
source-level comparisons are completed for the stated results; the
overall priority assessment and final submission review remain open.

The closest recognition comparison inspected in this pass is
[Coeurjolly--Sivignon--Dupont--Feschet--Chassery, 2005](https://doi.org/10.1016/j.dam.2005.02.022),
[author-hosted published text](https://perso.liris.cnrs.fr/florent.dupont/publications/2005-DAM.pdf).
Equations (4)--(5), Proposition 6 and Theorems 8 and 11 describe feasible
parameters of a digital plane using its two leaning polygons. Theorem 12
bounds the faces of one such parameter polyhedron under additional
hypotheses. The stated results and their proofs in Sections 3--4 were
inspected. Boundary-point recognition is established background; the
remaining contribution to assess is the uniform prescribed-cardinality
injection and its multiscale count, not that general geometric principle.

### Consolidated contribution assessment

**Decision: prepare the matching-count argument for mathematical review as
a focused discrete-geometry paper.** The case for doing so is the proved
quantitative improvement over a precisely identified published problem,
rather than the visual origin of the project or the absence of search hits.

| Requirement | Evidence and conclusion |
| --- | --- |
| Same mathematical object | The definition in Wagner's published abstract is exactly the finite orthant cut used here. Counting sets, coordinate labels and prescribed included cardinality agree. |
| Quantitative advance | Wagner's lower and upper orders are n^(d-1) log n and n^(d-1)(log n)^(d-1). Candidate 23 proves the latter order from below for every sufficiently large n in each fixed d. The gain over that lower bound is (log n)^(d-2); in d=3 it closes one logarithmic factor. |
| Complete uniform mechanism | The finite sieve controls all primes by a fixed nonzero determinant. The two-patch argument controls collisions at the prescribed cardinality, including partial ties. Observable extents separate aspect ranges. Every threshold and constant is displayed for all d>=2. |
| Closest later counting sources | Muller's full 2003 structural paper sends counting to Wagner; Kalka--Kuhnt Section 2.2 explicitly repeats his two orders; Elizalde--Galvan v1 Section 8.2 refers to his growth work. The last comparison uses the user-authorized assumption for the 2025 journal version. These are source statements, not a claim that citation indexes are exhaustive. |
| Known ingredients | Restriction, lifting, the fibre identity, primitive-matrix methods, boundary-point recognition and the corner-cut/model dictionary are attributed to their earlier sources. The claimed advance is the uniform matching count obtained from the construction. |
| Mathematical consequence | The theorem determines the logarithmic average completion multiplicity suggested by Wagner and the order of the known vertex/model count. The prime-realization proposition gives its exact interpretation as a count of smooth-prefix shapes over varying prime bases. |
| Size of the claim | This is a substantive specialist counting result if new. It gives no leading constant, growing-dimension estimate, prime-distribution theorem or fixed-prime random-access improvement. Its proof constants are deliberately coarse. |
| Remaining uncertainty | This comparison supports a favorable research assessment, but it is not an independent priority adjudication or referee report. Overall novelty and submission readiness remain unverified, so the active goal is not complete. |

The Kalka--Kuhnt source was retrieved directly for this assessment. Its
arXiv-v2 PDF has 24 pages and explicitly prints both Wagner orders on
p. 6. The title page is dated 27 November 2024; this internal date is
recorded separately from the arXiv version identifier, rather than used
to infer a journal publication date. Its Wagner bibliography entry has a
volume typo; our reference uses the primary Wagner record, volume 29.

The anonymous manuscript at `paper/corner_cut_growth/main.tex` now states
the actual advance directly, cites Wagner's upper theorem and Remond's
independent attribution, separates the inspected Elizalde--Galvan preprint
from the journal citation, and discloses the use of OpenAI Codex in proof
development, literature comparison, checks and drafting. It does not invent
an author, affiliation, independent review or acceptance. Earlier draft
status statements in this chronological notebook are retained as history.
The final contribution and submission assessment remains an open task;
the changes here complete the consolidated internal assessment, not that
external verification.

### Broader extremal problem: a bounded feasibility audit

The follow-up suggested after the matching-count assessment is the maximum
number of reduced Groebner bases of an arbitrary n-point configuration.
This is an existing problem, not a new definition introduced here.
[Onn--Sturmfels, 1999](https://doi.org/10.1006/aama.1999.0645),
[author-hosted text](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Onn/Onn4.pdf),
printed pp. 43--47, define F(n,d) by maximizing over n-point configurations
in d-dimensional affine space and over the underlying fields. Corollary
5.4(a) gives F(n,d)=O_d(n^{2d(d-1)/(d+1)}). Generic configurations give
C_d(n) bases. Thus Candidate 23 would give the starting range

    Omega(n^2 (log n)^2) <= F(n,3) <= O(n^3).

The lower estimate depends on the candidate proof and the known generic
dictionary. The upper estimate is a verified published bound, not a claim
that every subsequent improvement has been excluded. The printed source
pages, including the full proof of Proposition 5.6, were visually inspected.

The generic case is not an exact extremizer. Onn--Sturmfels Proposition
5.6 proves C_2(7)=8<F(7,2)=10. Any proposed transfer from the generic
corner-cut fan to every point ideal must respect this counterexample.
In particular the generic fan is not a universal common refinement.

#### Exact reproduction of the failed transfer

The existing growth checker now reproduces that known example over F_29.
Reduce the two equations for their final point modulo 29: alpha=17 and
beta=10 satisfy both equations. The seven distinct points are

    (0,0), (1,1), (2,2), (3,4), (5,7), (11,13), (17,10).

All fifteen partitions of seven are checked. Their evaluation determinants
vanish exactly for the three staircase moments (7,4), (6,6), and (4,7),
as in the published real example. The twelve remaining evaluation bases
have ten positive-weight minimizing moments; the unrestricted fifteen
staircases give eight. The calculation checks all 46 intervals between
positive pairwise moment crossings. At each sampled rational weight it
independently compares determinant-based minimization, direct monomial
rank selection, and heap generation of the generic corner cut.

The weights (30,19) and (42,29) give the same strict seven-point corner cut

    {(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(2,0)},

but give initial staircases with moments (3,9) and (5,6), respectively,
for this point ideal. Consequently a single generic corner-cut cone can
split. The full sets, determinant residues and exact weights are saved in
`results/research/corner_cut_growth_checks.json`, under
`nongeneric_point_ideal_control`. Reproduce with

```bash
python3 -B scripts/run_corner_cut_growth_checks.py
```

This is a reproduction and a negative control for a proposed proof route.
It is not a new counterexample or an asymptotic improvement. Enumerating
all fifteen partitions suffices here; no broader search framework is added.

#### Later primary comparisons and their limits

* [Babson--Onn--Thomas, 2003](https://doi.org/10.1016/S0196-8858(02)00509-2),
  [inspected preprint](https://arxiv.org/pdf/math/0207135), Corollary 2.5,
  extends the same exponent to every ideal of colength n, including
  nonradical ideals. Theorem 2.4 identifies initial staircases with state
  vertices; Theorem 3.5 gives the Hilbert-zonotope refinement; Theorem 4.2
  gives the polynomial-time algorithm. The latter is distinct from
  Theorem 3.5. The statements and the state-polyhedron, bound and greedy
  algorithm arguments were read. This supplies an existing general
  framework, not an improvement obtained from Candidate 23.
* [Zhang--Stigler, 2020](https://doi.org/10.1007/978-3-030-42687-3_9),
  [inspected arXiv v1](https://arxiv.org/pdf/1907.01080v1), Section 5,
  Theorem 4 and equations (8)--(13), specializes the moment-volume bound
  to n points in F_p^d. Written with our parameter names, its asymptotic
  expression is

      O_d((p^2 floor(n/p)+(n mod p)^2)^{d(d-1)/(d+1)}),

  using the complementary point count when n exceeds half the grid.
  For p>n the displayed expression is exactly the old asymptotic order
  n^{2d(d-1)/(d+1)}. It therefore does not improve the unrestricted
  extremal bound. Pages 15--17 of v1 were read and visually checked.
  We retain big-O constants rather than treating the preprint's displayed
  unit-constant formula as an exact universal numerical estimate.
  The [NSF-hosted later text](https://par.nsf.gov/servlets/purl/10177092)
  has Theorem 8, uses a containing hypercube instead of the v1 simplex,
  and states the bound asymptotically. Its indexed primary theorem
  passage was inspected, but full retrieval failed. The two versions
  are not represented as having received the same reading scope.
* [Dimitrova--He--Robbiano--Stigler, 2020](https://doi.org/10.1142/S0219498820500875),
  [inspected v2](https://arxiv.org/pdf/1704.03680v2), Theorem 4.8 and
  Corollary 4.11, show that complementary point subsets inside a grid
  have the same Groebner fan. The hypotheses, proof of Theorem 4.8,
  and the radical-grid specialization were read. Replacing a small
  configuration by its grid complement can increase the point count,
  but it preserves the number of bases. It cannot amplify this number
  to improve an asymptotic lower construction.

The lattice-ideal direction also has established constructions.
[Sabariego--Santos](https://arxiv.org/html/0705.3631), Theorems 1.2--1.3,
give arbitrarily many coherent minimum-distance diagrams for certain
triple-loop networks. Their displayed family has size
N=k(k-1), k=q^2+q+2, with q not divisible by three, and 3(q+2)
such diagrams: this particular lower
construction has order N^{1/4}. These are graded initial ideals of a
lattice ideal, a restricted family. The introduction and stated
theorems were inspected, not the full later proof. Neither the lower
estimate nor its being smaller than our generic count bounds the total
number of arbitrary-order bases of those ideals. It is not a new
construction available to claim here.

#### Decision from this pass

The broader extremal problem is a credible high-upside research target.
This pass does not establish its strongest current bounds or a new proof
mechanism. Direct generic-fan transfer fails exactly; grid-complement
amplification is ruled out by a known theorem. The finite-field estimate
does not itself close the unrestricted gap.

Do not expand the manuscript around a conjectural generalization or turn
the counterexample into a substitute result. A successful next step must
either bound the additional initial staircases caused by evaluation
dependencies, with a uniform improvement on the general upper bound, or
construct an infinite family whose growth exceeds the generic order.
In either case, the relevant later-bound audit must be completed before
claiming novelty. Candidate 23 and its original contribution assessment
remain intact; the active goal remains incomplete.

### Candidate 24: full ordered prefixes have the same growth order

This is a complete internal upper-bound argument and a consequence of
Candidate 23. It concerns an established family of finite monomial
orderings. Its priority comparison, especially with Vasiliev--Pavlov
(2009), remains incomplete. It is not yet being added as a new claim in
the manuscript or used to declare the active goal complete.

Let O_d(n) count distinct ordered tuples (x_1,...,x_n) consisting of the
first n points of Z_{>=0}^d under x -> w.x, as the positive real vector w
varies with rationally independent coordinates. Count the entire ordered
exponent tuple, with labeled coordinates. This differs from C_d(n),
which forgets the internal order. The assumption makes every scalar
product distinct; a finite prefix strictly realized by positive weights
can also be perturbed to such a vector, since its comparisons and minimal
excluded border impose only finitely many strict inequalities.

**Proposition 24.1.** For every fixed d>=2,

    C_d(n) <= O_d(n) = O_d(n^(d-1) (log(n+1))^(d-1)).

Consequently Candidate 23 gives

    O_d(n) = Theta_d(n^(d-1) (log n)^(d-1)).

The upper argument below also reproves the known corner-cut upper order.
Wagner's upper theorem retains its original attribution. No improvement
over that upper order is claimed.

**Proof of the upper bound.** The case n=1 has one tuple. For n>1, fix
one prefix and witnessing w, and put t=w.x_n>0. Its active coordinates
are J={i:w_i<=t}; write r=|J|>=1. Every point in the prefix has zeroes
outside J. Restrict to these r coordinates and define

    b_i=t/w_i >= 1,
    L_i=2^floor(log_2 b_i),
    Q=product_i L_i=2^R.

There are exactly n lattice points of w.x<=t in this coordinate face.
The simplex {x>=0:w.x<=t} is covered by the unit cubes based at those
points: floor(x) has no greater score. It follows that

    product_i b_i / r! <= n.

On the other hand all the n lattice points lie in their coordinate
bounding box, giving

    n <= product_i(floor(b_i)+1) <= 2^r product_i b_i.

Since L_i<=b_i<2L_i, these inequalities imply the shell bounds

    n < 4^r Q,                 Q <= r! n.                 (24.1)

For fixed r, the integer R therefore takes at most
2r+ceil(log_2(r!))+2 values, all at most log_2 n+log_2(r!).
For each R there are binomial(R+r-1,r-1) nonnegative integer tuples
(log_2 L_1,...,log_2 L_r). Thus the number of possible aspect tuples
is O_r((log(n+1))^(r-1)), rather than the count of all tuples with
product merely bounded above. The lower bound in (24.1) is essential
to this logarithmic exponent.

For a fixed aspect tuple, every point of the prefix lies in

    G_L = product_i {0,...,2L_i-1}.

Every comparison between two points of this box changes only on a
central hyperplane w.v=0 for a nonzero integer difference vector v with
|v_i|<=2L_i-1. Up to the sign of v there are at most

    (product_i(4L_i-1)-1)/2 <= 4^r Q <= 4^r r! n          (24.2)

such hyperplanes. Repeated hyperplanes only reduce this count. An
arrangement of m central hyperplanes in R^r has at most

    2 sum_{j=0}^{r-1} binomial(m-1,j) = O_r(m^(r-1))

chambers. This follows by adding one hyperplane at a time: the number
of new chambers is at most the chamber count of its restricted central
arrangement in dimension r-1. The formula follows by induction from
one hyperplane and dimension one, and is an upper bound even with
degeneracies.

Within each chamber the total order on G_L is fixed. Whenever a weight
with this chart witnesses a global n-point prefix, that prefix is the
first n points of G_L in this total order. Hence one chamber contributes
at most one ordered prefix for the fixed chart. A chamber need not
contain a weight realizing this chart, and the same prefix may occur
in several charts; both issues only overcount in an upper bound.
Rational independence ensures the witnessing weights lie off every
comparison hyperplane.

Equations (24.1)--(24.2) therefore give
O_r(n^(r-1)(log(n+1))^(r-1)) orders for each active coordinate set J.
There are binomial(d,r) choices of J. Summing over 1<=r<=d proves
the asserted upper bound. This active-face step handles unused
coordinates; applying a full-dimensional volume estimate without it
would be invalid when some b_i<1.

Every corner cut admits a generic positive separator after a sufficiently
small perturbation, so forgetting order is a surjection from these
ordered prefixes to the n-point cuts. This proves C_d(n)<=O_d(n).
Candidate 23 supplies the matching lower order. No numerical check is
used in either uniform argument.

An ordered prefix determines every initial set of sizes 1,...,n, so this
also controls the common refinement of the generic corner-cut fans over
those sizes. It is not a universal refinement for arbitrary point ideals:
their greedy basis construction may reject early monomials because their
evaluation vectors are dependent. The counterexample in the preceding
feasibility audit still applies to that broader problem.

**Consequence 24.2 (average ordering multiplicity).** For an n-point
corner cut D, let m(D) be the number of different orders in O_d(n)
whose underlying set is D. Then

    1 <= (sum_{|D|=n} m(D))/C_d(n) = O_d(n)/C_d(n) = O_d(1).

The first inequality holds for every n and the last bound for all
sufficiently large n in each fixed dimension. This is a uniform mean
over the finite set of cuts, not a bound on m(D) for every cut and not
an assertion that O_d(n)=C_d(n). In dimension two at n=3 there are four
ordered prefixes and three underlying shapes: the shape {0,e_1,e_2}
has two possible orders. The finite checker retains such witnesses.

**Consequence 24.3 (prime bases).** Every ordered prefix counted by
O_d(n) can be realized over distinct ordered prime bases above any
prescribed B. To strengthen the preceding prime-realization proof,
normalize its witnessing positive weights so that the separating
threshold is one. In addition to the positive cutoff gap g, let

    epsilon=min_{j<n} w.(x_{j+1}-x_j)>0

when n>1. Choose the Bertrand-interval scale T large enough that every
included product's normalized logarithmic perturbation is less than
min(g,epsilon/2). All perturbations are nonnegative coordinatewise,
so excluded products remain above the threshold and consecutive
included products remain in the prescribed order. For n=1 only the
cutoff condition is needed. The disjoint-interval construction already
ensures the primes are distinct and exceed B. Thus the number of full
smooth exponent-prefix words over varying labeled prime bases greater
than B equals O_d(n), not just the unordered count C_d(n). This is an
elementary use of the same arithmetic input, not a prime-distribution
theorem or a fixed-prime random-access improvement.

#### Exact falsification checks

`ordered_prefix_chart_checks` in the existing growth checker verifies
1,080 strict prefixes in ambient dimensions two through six. Every
effective dimension one through six occurs, including prefixes confined
to a single axis. All chart containments, integer forms of the volume
inequalities, dyadic product bounds and comparison-direction estimates
pass. A further 132 cases test product-shell widths in effective
dimensions one through twelve. The original uniform lower-bound checks
remain separate.

Four additional prime-realization fixtures check the entire prescribed
order by exact integer products and an independent heap, including
opposite orders on the same shape and an unused coordinate. These are
reported separately from the original six unordered realization checks.

For dimension two the critical comparison slopes are coprime a/b>0.
The first change in the prefix occurs at n=(a+1)(b+1)/2: at weights
(b,a), direct counting gives exactly (a+1)(b+1)/2-1 lattice points
strictly below the tied points (a,0) and (0,b). Accordingly

    O_2(n)=1+#{(a,b)>=1:gcd(a,b)=1, (a+1)(b+1)<=2n}.

This planar family and the threshold rule were already recorded by
Franklin T. Adams-Watters in [OEIS A163823](https://oeis.org/A163823)
(2009). They are not new claims here. The new finite fixture sweeps
every such slope for lengths 1--60, checks the full heap-generated
orders and all 261 first-crossing ranks at length 60, and agrees with
the contributor's first twelve values. It records all sixty computed
counts, ending in 262. These are checks of finite failure modes, not
evidence for the uniform theorem's priority.

Reproduce both fixtures with

```bash
python3 -B scripts/run_corner_cut_growth_checks.py
```

#### Newly identified close prior art and outstanding comparison

[Vasiliev--Pavlov, 2009](https://doi.org/10.1134/S0361768809020030),
*Enumeration of finite monomial orderings and combinatorics of universal
Groebner bases*, explicitly studies finite intervals of monomial orders,
their enumeration, and several admissibility classes. Its primary
publisher abstract and reference list were inspected; its full definitions
and quantitative statements have not been compared. An authored summary
of Pavlov's 2012 thesis also reports an established correspondence with
signature prefixes. The finite-ordering family and that correspondence
must not be presented as introduced by this project. The thesis summary
was inspected through a text mirror; it is a locator for a full primary
comparison, not a substitute for one.

[Vietri, 2002](https://doi.org/10.1023/A:1021236722756),
*Combinatorial aspects of total weight orders over monomials of fixed
degree*, studies orders on a fixed homogeneous layer and also discusses
prime-induced orders. Its authored abstract, available through the
Sapienza institutional record, was inspected. A layer of fixed total
degree is a different domain from the first n monomials across all
degrees, but full results must be checked before drawing a stronger
priority conclusion. In particular the broad connection between primes
and monomial weight orders is established research.

The next source comparison should decide whether Proposition 24.1's
quantitative upper bound is already stated or immediate there. If it is,
credit it and use the matching lower count only for the consequence;
do not convert a known upper argument into an additional novelty claim.
The complete dyadic proof is saved now so that the mathematical claim
can be assessed independently of the source-access question. The main
manuscript and Candidate 23's novelty gates have not been relaxed.

#### Precise Vietri comparison and prime-realization attribution

The complete primary text of [Vietri (2002)](https://doi.org/10.1023/A:1021236722756)
is publicly readable in the [RISC bibliography copy](https://www3.risc.jku.at/research/theorema/Groebner-Bases-Bibliography/gbbib_files/publication_455.pdf).
The inspected scope is the definitions and discussion in Sections 1--2,
the enumeration argument in Section 3, and Proposition 6.2 with its proof
and Corollary 6.3 in Section 6 (printed pp. 259--260). This is not an
assertion that every proof in Sections 4--5 and 7 has been audited.

**The prime-realization bridge is established prior art.** Proposition
6.2 states that any finite family of strict product inequalities with
nonnegative real exponents, satisfiable by increasing real bases greater
than one, is satisfiable by increasing prime bases. Its proof uses
Bertrand's postulate after a sufficiently large power scaling. That
scaling can be increased arbitrarily, so all primes can exceed a
prescribed bound B. The extra lower-bound requirement here is therefore
also covered by that proof, rather than being a new prime theorem.

Here is the exact specialization, including the infinite-complement
issue. Let D={a_1,...,a_n} be a corner cut with a prescribed strict
weight order, and let E be the coordinatewise minimal elements of its
complement. The set E is finite: it is contained in {0,...,n}^d,
since a minimal excluded point has every immediate predecessor in D.
Choose a generic positive witnessing weight w, sort its distinct
coordinates temporarily, and put X_i=exp(w_i)>1. Impose

    X^(a_j) < X^(a_(j+1))       for 1 <= j < n,
    X^(a_n) < X^b              for every b in E.

Apply Vietri's proposition and undo the coordinate permutation. Every
excluded point dominates some b in E; since the resulting primes are
all greater than one, all excluded products exceed the final included
product. Thus both the unordered realization proposition in the
manuscript and Consequence 24.3 for entire ordered prefixes follow from
Vietri's finite-inequality theorem. For n=1 the first family is empty.
The manuscript now credits Proposition 6.2 explicitly and retains its
short direct proof. Existing exact fixtures remain checks of our
implementation, not novelty evidence for this bridge.

**Exact relation to the fixed-degree domain.** Vietri's set D_h^k consists
of all (k+1)-tuples of nonnegative integers whose sum is h, with weights
in increasing coordinate order. His beta-linearizations are the
restrictions of weight orders to that whole homogeneous layer. There is
an explicit homogenization link to our problem, rather than merely a
difference in terminology. Map

    a in Z_{>=0}^d, |a|_1 <= h  ->  (h-|a|_1,a) in D_h^d.

Choose r_0>0 and r_i=r_0+w_i, sorting the positive w_i as needed.
The layer score becomes

    r_0 h + sum_i w_i a_i.

Consequently layer orders contain all weight orders on the entire
finite total-degree ball. Every point in an n-term prefix has total
degree at most n-1: a coordinate-increasing path to it contains
|a|_1+1 earlier-or-equal points. Taking h=n-1 therefore allows the
full layer order to determine the ordered n-prefix. This is a
comparison via a much larger domain, not a direct identification of
h with n or of the two counting functions.

Section 3 enumerates the layer problem in three variables. Its count
has quadratic order in h; transferring that order by h=n-1 supplies
only O(n^2) for planar ordered prefixes, rather than the sharper
O(n log n) of Proposition 24.1. For three-dimensional prefixes the
homogenization uses four variables; Vietri explicitly states in the
introduction (printed p. 240) that his approach does not enumerate the
case k>=3. These inspected statements do not supply Candidate 23's
prescribed-cardinality lower bound. This is a bounded comparison of
specified results, not proof of originality over all finite-ordering
literature. The Vasiliev--Pavlov full quantitative comparison remains
pending; its publisher abstract is insufficient to close that issue.

**Impact assessment.** Candidate 23 remains the substantive candidate:
a uniform matching cardinality count in every fixed dimension, with
consequences for known corner-cut polytopes, generic design models and
smooth-prefix sets. Candidate 24 is a natural extension to complete
finite words and would give bounded mean internal-order multiplicity,
subject to its own priority check. The prime-realization proposition is
now specifically attributed background. The broader arbitrary-point-ideal
problem still has the previously documented obstruction and no new
asymptotic mechanism here. Finish the proof and priority assessment of
the focused counting contribution before expanding that harder target.
No prime-distribution or fixed-prime unranking improvement is claimed.

#### Global monomial orders: exact finite-prefix scope

A further primary comparison is now complete for
[Vasiliev (2003), *Monomial Orderings, Young Diagrams and Groebner Bases*](https://www3.risc.jku.at/research/theorema/Groebner-Bases-Bibliography/gbbib_files/publication_2491.pdf),
printed pp. 397--401. All four numbered theorems, their proofs and the
conclusion were read in the published text. The paper establishes the
Dickson-set and universal-basis framework and a linear-algebra membership
criterion; it does not state a quantitative count by prefix cardinality.
Its final remark says that the weakly admissible extension is outside
that paper's scope. This completes the comparison with that earlier
paper, not with the later Vasiliev--Pavlov (2009) enumeration paper.

Pavlov's authored thesis introduction and conclusion, reproduced in the
[public text mirror](https://www.dissercat.com/content/polinomialnye-modeli-osnovannye-na-diagrammakh-yunga),
explicitly report the correspondence between strictly admissible finite
orders and signature prefixes, and cite OEIS A163823 and Kimberling.
Those passages have been inspected. They do not contain the detailed
quantitative statements needed for the 2009 comparison. The mirror's
abstract-viewer PDF was also checked and returned Not Found. No inference
of novelty is drawn from this source-access failure. The author's
[public thesis announcement](https://shmat-razum.blogspot.com/2013/02/blog-post.html)
provides another locator for the full thesis; that copy has not been read.

The exact relationship with the global class can nevertheless be proved
from established primary work. Robbiano's *On the Theory of Graded
Structures*, [MPI/SFB preprint 85-17](https://archive.mpim-bonn.mpg.de/2939/1/preprint_1985_17.pdf),
Theorem 2.5 and its proof, printed pp. 23--25, represent an ordered
rational exponent group by finitely many lexicographically compared
linear forms. The accompanying term-order remark applies this to
monomials. The equivalent classification is Theorem 5 of
[*Term Orderings on the Polynomial Ring*](https://doi.org/10.1007/3-540-15984-3_321).
The latter theorem was also read in a reproduction of the author's text.
This representation and the resulting finite-weight equivalence are
known background, not new claims here.

**Finite-prefix lemma (known-method deduction).** Let x_1,...,x_n be the
first n exponent vectors of any global monomial order on d variables:
a total additive order with zero least. There is a positive weight
vector, rationally independent after a sufficiently small perturbation,
whose entire ordered n-prefix is exactly this tuple. Conversely every
positive rationally independent weight vector defines a global monomial
order. Thus Candidate 24's O_d(n) counts exactly the finite n-prefixes
of all global monomial orders, with labeled variables.

**Proof.** The underlying set D is a finite coordinate downset. Its
minimal excluded border E is finite, and every excluded point dominates
an element of E. Extend the monomial order to the additive rational
group by cancellation. By Robbiano there are real linear forms
u_1,...,u_s such that the first nonzero entry of

    (u_1.v,...,u_s.v)

is positive whenever v is positive in this ordered group. Apply this
to the finite family

    v=e_i,
    v=x_(j+1)-x_j  (1<=j<n),
    v=b-x_n       (b in E).

For all sufficiently small epsilon>0, every polynomial
sum_l epsilon^(l-1) u_l.v is positive: its first nonzero coefficient
is positive and the family is finite. Consequently
w=sum_l epsilon^(l-1)u_l has positive coordinates, strictly orders the
included tuple and puts all E above x_n. Positivity handles the entire
infinite complement. The strict constraints define an open neighborhood,
which contains vectors with rationally independent coordinates. Conversely
positive weights are additive and have only finitely many nonnegative
lattice points below any fixed score, so rational independence gives a
well-ordered total monomial order. This proves both directions.

This is a finite-prefix statement. A single real weight vector need not
represent a given global monomial order at every length simultaneously;
lexicographic and graded orders are important examples. It also does not
identify all weakly admissible finite tableaux with this class, or extend
the bound to arbitrary point ideals whose greedy basis search can reject
early monomials. The existing seven-point negative control still applies.

The existing Proposition 17 therefore applies to all these prefixes.
There is no need for a second weight-compression construction: the saved
LP/Cramer proof already gives positive integers A_i and S=A.x_n with

    n-1 <= S <= c_d n,     c_d=d(d!)^2 2^(d-1).

For an active coordinate e_i belongs to D, so A_i<=S. Inactive
coordinates were assigned S+1. In particular

    max_i A_i <= c_d n+1.                                      (24.3)

For n=1 use all A_i=1. Both the cutoff bound and its existing exact
oracle are retained; the generic weight-compression problem is not being
presented as new.

**Quantitative prime realization (supporting corollary).** For every such
ordered n-prefix and B>0, it is enough to use distinct prime bases with

    B < p_i < 2^[q(c_d n+1)+d],
    q=d n+ceil(log_2(max(B,2))).                                (24.4)

In particular, for fixed d and fixed B their binary lengths are O_d(n^2).
This is a deduction combining Proposition 17 with the established
finite-realization method; novelty of this quantitative formulation has
not been established. It is not a new prime-distribution estimate.

**Proof.** Sort the integer weights A_i temporarily. Bertrand gives a
prime in each disjoint interval

    2^(q A_i+i-1) < p_i < 2^(q A_i+i),        1<=i<=d.

The lower endpoints exceed B. Every included exponent vector x has
total degree at most n-1, since a coordinate-increasing path to x lies
in D. Its product is at least 2^(q A.x), and is strictly less than
2^[q A.x+d(n-1)], except that the zero vector has product exactly one.
Since q>d(n-1), an integer gap of at least one between successive
included scores preserves their strict product order. Every excluded
vector has integer score at least S+1, and hence product at least
2^[q(S+1)], above all included products. The zero vector is included,
so there is no excluded zero exception. This proves the whole-prefix
claim, and (24.3) gives (24.4). Undo the coordinate permutation.

This supplies a finite prime-size budget; it gives neither an efficient
search for a realizing tuple nor an algorithm for recovering the prefix
of a supplied prime tuple without enumeration. Its proof uses classical
arithmetic and previously saved geometry. Its role is an explicit scope
and size consequence, not a substitute for the sharp counting contribution
or for the outstanding priority comparison.

**New falsification checks.** `check_monomial_prefix_realization` in the
existing interleaving script checks 90 prefixes of lengths 1, 2, 3, 5,
8 and 13 in dimensions two through four. Inputs come from integer matrix
orders (lexicographic, reversed variable lexicographic, graded lexicographic,
graded reverse lexicographic and elimination), rather than logarithms of
prime tuples. The original exact LP oracle constructs scalar certificates;
an independent scalar-score heap verifies the entire ordered prefix and
strict cutoff. The singleton and all effective dimensions one through
four occur. Exact matrices, exponent tuples, integer weights and next
scores are recorded. These checks test the finite/global and inactive-face
failure modes. They are not evidence of novelty and are not a replacement
for either uniform argument.

Reproduce with

```bash
python3 -B scripts/run_interleaving_research.py --ranks 20000
```

The main manuscript remains focused on Candidate 23; Candidate 24's
quantitative priority assessment is still open. The finite monomial-order
family is established mathematics. The next decisive source is still the
full 2009 Vasiliev--Pavlov paper, rather than its abstract, bibliography or
an inference from inaccessible copies.

### Candidate 25: sharp worst-case prime size for finite-prefix realization

**Status: complete new internal quantitative argument; priority unverified.**
This sharpens the quadratic bit bound (24.4). It uses additional arithmetic
information, namely a published short-interval prime estimate, and gives a
concrete prime-size consequence. It is independent of Candidate 23's lower
count. Neither the use of prime bases nor the finite-weight bridge is new.
The pending finite-ordering comparison remains necessary before assessing
this formulation as a publication contribution.

Fix d>=2. For an ordered n-prefix W of a global monomial order, let
h_B(W) be the least value of log(max_i p_i) over ordered tuples of distinct
primes p_i>B whose smooth exponent prefix is exactly W. Logarithms in
this section are natural. Let

    H_d(n,B)=max_W h_B(W).

The finite-prefix lemma shows that this is the same class of words counted
by O_d(n). For n>=2 and B>0, the proposed sharp order is

    H_d(n,B)=Theta_d(n log(max(B,2))).                           (25.1)

The same worst-case order holds if only the underlying n-point corner-cut
set is prescribed. Thus every admissible finite pattern is realizable
with primes having O_d(n) bits when B is fixed, and a linear number of
bits is necessary for some patterns. The constants may depend on d;
there is no claim of optimal constants, fast construction or a statement
about primes in a fixed input tuple.

#### Arithmetic input and exact source scope

Dusart, [*Estimates of Some Functions over Primes without R.H.*,
arXiv:1002.0442v1](https://arxiv.org/pdf/1002.0442v1), Proposition 6.8,
printed p. 8, states that for x>396738 there is a prime in

    x < p <= x(1+1/[25 (log x)^2]).                              (25.2)

The statement and its complete displayed proof were read directly in
Section 6.1.3. This is a passage-level comparison, not a full audit of
all twenty pages or the earlier analytic estimates used in that proof.
The construction below uses x>=exp(30), so none of the small-threshold
computational transitions is needed here. We use this established
unconditional theorem, not an assumed prime-gap conjecture or a new
prime-distribution estimate.

#### Proof of the upper bound

For n>=2 take the integer certificate supplied by Proposition 17 and the
finite-prefix lemma. Write its weights as A_i, its last included score
as S, and put c_d=d(d!)^2 2^(d-1). Then

    n-1 <= S <= c_d n,       1 <= A_i <= S+1 <= c_d n+1,
    D = {x in Z_{>=0}^d : A.x <= S},

and all included scores are distinct in the prescribed order. Temporarily
sort the coordinates so m=A_1<=A_2<=...<=A_d.

**Axis-collision observation.** For every two distinct coordinates i,j,

    lcm(A_i,A_j)>S.                                           (25.3)

Otherwise the two distinct pure-axis points of that common score are both
in D, contradicting the strict included order. In particular, for j>=2,

    m A_j >= lcm(m,A_j)>S>=n-1,    A_j>=m,
    A_j^2>=m A_j>n-1.                                         (25.4)

There is at most one coordinate weight at or below sqrt(n-1). This is
why it is sufficient to anchor a single prime exactly.

Put

    a=max(30,d,log(max(B,2))+1),
    delta=1/[25 a^2(n-1)].

Bertrand supplies a prime p_1 in (exp(a m),2 exp(a m)). Define

    T=log(p_1)/m,       a<T<a+log(2).

For j=2,...,d put

    X_j=exp(T A_j)(1+delta)^(j-2).

Then X_j>=exp(30). By (25.4),

    1/[25 (log X_j)^2] < 1/[25 a^2(n-1)] = delta.

Apply (25.2) separately to each X_j, obtaining

    X_j < p_j < X_j(1+delta).                                  (25.5)

These intervals are disjoint, including when some A_j agree:
the next lower endpoint is at least the previous upper endpoint.
Also X_2>=p_1, so p_1,...,p_d are distinct increasing primes above B.
Set eta_i=log(p_i)-T A_i. The anchor gives eta_1=0, and (25.5) gives

    0<eta_j<(j-1)log(1+delta)<(d-1)delta   for j>=2,

where the last strict bound may be replaced by <= at j=d without
changing the following estimates. More simply, eta_j<=(d-1)delta
for all j. Every included x has |x|_1<=n-1, so

    0 <= eta.x < (d-1)/(25 a^2) < a < T.                        (25.6)

For any successive included x,y, A.y-A.x>=1; hence

    log(p^y)-log(p^x) >= T-eta.x > 0.

Every excluded y has A.y>=S+1 and eta.y>=0, whereas every included
x has log(p^x)<TS+T. Thus all excluded products exceed every included
product, controlling the infinite complement as well as internal order.
Undoing the coordinate permutation realizes W.

Finally, (25.5) and (25.6) imply the explicit bound

    log(max_i p_i) < (a+log 2)(c_d n+1)+1.                      (25.7)

Here the residual (d-1)delta is less than one because a>=d and n>=2.
Since a=O_d(log(max(B,2))), this proves the uniform upper order in
(25.1). For n=1 any d distinct primes above B suffice; repeated Bertrand
intervals give log(max p_i)<=log(max(B,2))+O_d(1).

#### Matching lower bound

Consider the admissible prefix

    W=(0,e_1,2e_1,...,(n-1)e_1).

Any realizing tuple must satisfy p_j>p_1^(n-1) for every j>=2, because
e_j is excluded. Therefore

    h_B(W) >= (n-1)log(max(B,2)).                               (25.8)

For n>=2 this is at least n log(max(B,2))/2. The same axis example
works when only the underlying set is prescribed. Equations
(25.7)--(25.8) prove the claimed sharp worst-case order. The singleton
has the corresponding order log(max(B,2)).

#### What changes, and what remains open

Vietri's Proposition 6.2 already supplies qualitative prime realization
for arbitrary finite inequality systems. The argument here bounds the
worst-case height uniformly in prefix length and proves that its order
is necessary. The prefix-dependent linear integer-score bound and the
axis-collision observation allow one prime to be fixed exactly while
Dusart controls the others. No assertion that Dusart's prime estimate,
Robbiano's representation, finite prime realization or generic
weight-compression methods are new is made.

Combining (25.7) with Candidates 23--24 would show that allowing primes
of only O_d(n) bits already realizes the entire family of n-term words
and shapes with the asserted sharp counting order. This is a finite
universality statement under an explicit size budget. It does not
improve fixed-prime unranking complexity or guarantee a fast enumeration
of the realizing tuples. The formulation's priority and its significance
as an addition to the main paper remain to be established.

The exact compression fixtures now check (25.3) for every pair of weights
in the existing 230 cases, and the 90 matrix-order fixtures independently
check (25.3), (25.4), scalar prefix agreement and the rational error budget
in (25.6). This is a check of failure modes, not verification of (25.2)
by experiments and not the proof of the uniform theorem. No large-prime
search framework has been added. The main manuscript has not been
expanded while the finite-ordering priority comparison remains pending.

#### Candidate 25: fresh internal proof audit and published input

A fresh internal review checked the following obligations in the full
prime-height argument. This is an audit by the research agent, not an
independent human referee report. No unresolved logical step was found.

| Obligation | Verification |
| --- | --- |
| Scope of the integer certificate | Robbiano's finite comparison argument handles all global monomial-order prefixes; Proposition 17 applies after choosing a generic positive separator. It does not give one scalar representation for the entire infinite order. |
| Active-face volume bound | The closed score-T simplex is covered by exactly n based unit cubes in its active face. Its volume bounds the intercept product. Inactive axes are excluded from this volume calculation. |
| LP existence and Cramer bound | The finite strict border and consecutive comparisons become feasible unit-margin inequalities after scaling. Minimizing the positive coordinate sum gives a compact optimal face and a vertex. Column bounds 2B_i and the product estimate give S<=c_d n. |
| Bound for every restored weight | An active axis point e_i lies in D, so A_i<=S. An inactive weight is S+1. Thus max A_i<=c_d n+1, including unused coordinates. |
| Axis-collision rigidity | A common axis score <=S would put two distinct lattice points in D with equal integer scores. Therefore lcm(A_i,A_j)>S. Sorting yields A_j^2>n-1 for every j>=2. Ties after S are permitted and irrelevant. |
| Exact anchor | Choose p_1 in (exp(a m),2exp(a m)) and define T=log(p_1)/m. Then a<T<a+log2 and the anchor error is exactly zero, rather than an uncontrolled error on a frequently repeated small-weight coordinate. |
| Interval range and distinctness | Every X_j>=exp30, and its logarithm exceeds a sqrt(n-1). The short-interval upper width is strictly smaller than delta. Ordered endpoints are disjoint even when larger integer weights agree; every selected prime is above B. |
| Internal order and infinite complement | All coordinate errors are nonnegative and their total on an included point is <T. Each included score gap is at least one. Every excluded point has integer score >=S+1, so the argument covers the infinite complement, not just the finite border. |
| Height bound, arbitrary B and singleton | The residual logarithmic error is <1; a=O_d(log(max(B,2))) gives the uniform upper order. The singleton is handled separately with Bertrand intervals. |
| Matching lower family | Excluding e_j from the pure-axis n-prefix forces p_j>p_1^(n-1). This applies to both words and underlying sets and gives the same order uniformly in B. |

The arithmetic input also has a directly inspected published form.
Dusart, [*Explicit estimates of some functions over primes*](https://doi.org/10.1007/s11139-016-9839-4),
*The Ramanujan Journal* 45 (2018), 227--251, Proposition 5.4, printed
p. 242, states that for every real x>=89693 there is a prime

    x < p <= x(1+1/(log x)^3).

The statement was read in the indexed text of the
[published PDF reproduction](https://piyanit.nl/wp-content/uploads/2020/10/art_10.1007_s11139-016-9839-4.pdf),
including the surrounding Section 5.2 and Corollary 5.5. Direct retrieval
of the complete PDF timed out, so this is a statement-level comparison,
not a claim to have read the complete published proof or all analytic
prerequisites. Publisher metadata and the author's publication list agree
on the article. For every x>=exp30,

    1/(log x)^3 < 1/[25 (log x)^2].

Consequently the published proposition supplies the exact weaker input
used in Candidate 25 throughout its construction. The explicit constants
and the proof of (25.1) remain unchanged. The original 2010 statement and
its displayed proof stay recorded as an additional inspected source;
publication dates and numbering are not conflated.

The [author's publication list](https://www.unilim.fr/pages_perso/pierre.dusart/Publications.html)
also links a [correction revised 26 February 2018](https://www.unilim.fr/pages_perso/pierre.dusart/Recherche/correctif_RJ.pdf).
Both pages of that correction were read. It corrects two typographical
errors in Theorem 3.5 and gives PARI/GP code for reproducing Table 1.
It does not state a replacement for Proposition 5.4 or Corollary 5.5.
No claim that all analytic dependencies have been independently rederived
is made; the published interval theorem is an established input.

#### Required primary-text comparison before expanding the paper

The detailed Vasiliev--Pavlov (2009) comparison has remained open across
three consecutive saved research turns, visible in public commits
`1047a97`, `b89d9ef` and `d6cd1cf`. The exact source is

> N. N. Vasiliev and D. A. Pavlov, *Enumeration of finite monomial
> orderings and combinatorics of universal Groebner bases*,
> Programming and Computer Software 35(2) (2009), 79--89,
> [doi:10.1134/S0361768809020030](https://doi.org/10.1134/S0361768809020030).

The [RISC bibliography record](https://www3.risc.jku.at/research/theorema/Groebner-Bases-Bibliography/details.php?details_id=3375)
was inspected directly: it lists this paper and its publisher link but
provides no attached full text. The public publisher abstract and the
thesis introduction/conclusion establish that the finite-ordering family
and signature-prefix correspondence are prior art. They do not supply
the definitions, theorems and quantitative arguments required to decide
whether our bounds are already stated or immediate consequences there.
The earlier Vasiliev paper and the precise Vietri statements have been
compared separately; they cannot substitute for the missing 2009 text.

A readable copy of this paper, or Pavlov's full thesis treatment of these
finite orderings, is now the required input for the next priority gate.
The comparison must check the order classes, whether the domain is an
initial interval or an arbitrary finite set, the counting parameter,
any asymptotic estimates and the complexity analysis of its enumeration
algorithms. It must also check any quantitative prime-realization bounds.
If those results cover a candidate, credit them and reassess the actual
remaining contribution without weakening the original acceptance criteria.

The main matching-count manuscript, the ordered-prefix proof and the
prime-height proof are preserved. The new candidates are not being added
to the manuscript while this concrete source comparison remains open.
Source access failure is not evidence of novelty, and further isolated
examples would not answer the missing comparison.

### Full Vasiliev--Pavlov comparison: source gate resolved, 2026-10-18

The user supplied a readable copy of the published eleven-page paper,
*Enumeration of finite monomial orderings and combinatorics of universal
Groebner bases*, Programming and Computer Software 35 (2009), 79--89,
[doi:10.1134/S0361768809020030](https://doi.org/10.1134/S0361768809020030).
Every page, all nine propositions, eight tables, algorithm figures,
implementation comments and references have been read; the page renders
were also inspected. The supplied PDF has SHA-256
`9b7379754cca81198d84ce21ced054de59fdce485070e0cf786f25a77873da20`.
The source PDF and extracted text are local reading materials and are not
being added to the public repository.

#### Definitions and the exact common class

The paper uses *degree n* for the number of monomials in a finite chain,
not for their total polynomial degree. Definition 4, p. 80, imposes
divisor closure together with the relevant finite ordering axiom.
Definition 5, on the same page, calls a finite ordering *strictly
admissible* when it is the beginning of a total admissible monomial
ordering. This is exactly Candidate 24's family, with coordinates labeled.
The finite-weight consequence of Robbiano already recorded above proves
the equivalence with strict positive-weight prefixes. This equivalence
and the existence of enumeration algorithms are established background.

The distinction from the paper's wider finite *admissible* class matters.
Remark 1 gives the chain

    1 < x < x^2 < y < y^2.

It satisfies the finite local condition but cannot be an initial interval:
x<y forces xy<y^2, yet xy is missing. Thus Tables 5--6 must not be used
as counts of our ordered prefixes. The relevant counts are Tables 7--8
on p. 88. Forgetting the order in a strictly admissible prefix gives a
corner cut; every corner cut has at least one such strict internal order.
The set count C_d(n) and the word count O_d(n) remain different parameters.

| Paper location | Actual result or computation | Boundary for this project |
| --- | --- | --- |
| Definitions 1--5 and Proposition 1, p. 80 | Four ordering classes, their inclusions and the distinction between local finite admissibility and a global initial interval | Credit the object and terminology; do not claim to introduce the ordered-prefix family. |
| Sections 3--4, Propositions 2--3, pp. 80--82 | Recursive divisor-border enumeration and convex-hull tests for larger classes | These classes include many words outside the strict prefix class. Their table sizes do not give its asymptotic count. |
| Section 5, Propositions 4--5 and Figure 7, pp. 82--83 | Weight cones, finite linear comparisons and a reduced constraint set for finite admissible orderings | Finite inequalities, weight feasibility and cone-based enumeration are prior art. No bit-size or uniform prefix-cardinality estimate is stated. |
| Section 6, Propositions 6--7 and Figure 8, pp. 83--84, 86 | Multiplication/front constraints and a recursive enumeration procedure for strictly admissible finite orderings | This concerns our exact word class. It supplies an enumeration method, not the sharp asymptotic count or a prescribed-size lower-bound construction. |
| Section 7 and Tables 7--8, pp. 84--85, 88 | Numerical strict-prefix counts in dimensions two and three; discussion of partition enumeration | The values agree with the exact checks below. No all-n asymptotic theorem is stated for strict prefixes or corner-cut shapes. |
| Sections 8--9, Propositions 8--9, pp. 85--89 | Finite descriptions of Groebner-fan cells, pseudobases and extended universal bases | These known algebraic constructions do not state a quantitative corner-cut or strict-prefix cardinality bound. No new fan algorithm is claimed here. |
| Section 10, p. 89 | Common Lisp implementation for arbitrary dimension, exact integer arithmetic and normalized integer normals | Exact arithmetic is prior implementation practice. No uniform runtime, integer-weight-height or realizing-prime-size theorem is supplied. |

There is a reversed inequality in the opening prose of Section 6: its
literal statement about excluded monomials conflicts with Definition 5.
The comparison uses the explicit initial-interval definition and the
published strict-prefix tables. It does not rely on that sentence or
claim to validate every line of the published pseudocode.

#### Independent finite identification checks

The existing research checker now includes
`vasiliev_pavlov_table_checks()`. It compares **all 24 entries of Table 7**,
including n=50,100,200,...,700, with the exact planar crossing formula

    O_2(n) = 1 + #{(a,b)>=1 : gcd(a,b)=1, (a+1)(b+1)<=2n}.

All agree. This reuses the established signature-prefix crossing rule;
it is not a new planar counting claim.

For dimension three, the checker exhausts every positive-weight strict
prefix through length 10 by extending the minimal excluded border.
The exact unit-margin feasibility oracle includes both consecutive
comparisons and the final excluded border. Each accepted certificate is
checked by an independent additive heap for the full order and strict
cutoff. An infeasible LP is rejected only after exhausting its vertex
bases; a search cap raises an error rather than silently rejecting it.
Every genuine strict prefix is reached because every next point is in
the minimal excluded border of its predecessor. Counts for n=3,...,10 are

    9, 21, 39, 69, 123, 177, 243, 339,

agreeing with all eight corresponding Table 8 entries. Across lengths
2--10, 1,023 accepted certificates are heap-checked and 2,856 infeasible
extensions are rejected. The paper's Remark 1 chain is separately
rejected as a strict prefix. These checks identify the parameter and
family; they do not prove the asymptotic results, reproduce the larger
three-dimensional table entries, or establish novelty.

Reproduce with

    python3 -B scripts/run_corner_cut_growth_checks.py

The report's `vasiliev_pavlov_table_checks` field records the exact
published/computed counts, canonical word-list hashes, representative
integer certificates, method and limitations. No new search framework
or mandatory dependency has been added.

#### Consequences for the proposed contribution

* **Candidate 23 is retained.** The full paper does not state a matching
  asymptotic corner-cut count or an argument yielding the multiscale
  prescribed-cardinality lower bound. The comparison does not preempt
  the proposed improvement over Wagner's inspected bounds.
* **Candidate 24 is retained with explicit attribution.** It would give
  the asymptotic order of the exact strict-prefix family enumerated in
  this paper. The class, finite scalar representation, feasible cones and
  recursive enumeration are prior art. The candidate contribution is the
  sharp uniform counting order and its bounded mean order multiplicity.
* **Candidate 25 is retained.** The paper contains no prime-realization
  theorem, prescribed prime-height bound or asymptotic integer-certificate
  size estimate. The known qualitative prime-realization theorem remains
  credited to Vietri. Our candidate is the sharp worst-case prime-bit
  budget, using established prime-interval estimates.

The exact missing-primary-text comparison is now complete and its
source-access blocker is resolved. This is favorable evidence for
continuing the quantitative paper, not a certificate of originality
against all literature or an independent referee report. Do not mark the
publication goal achieved solely because this paper does not preempt the
claims. The next task is the final contribution and manuscript assessment,
using the completed source comparison and unchanged acceptance criteria.

### Integrated quantitative manuscript, 2026-10-18

The source-access comparison is resolved, so Candidates 24--25 and their
supporting lemmas are now integrated into `paper/corner_cut_growth/main.tex`.
The resulting paper is titled *Corner cuts and monomial prefixes: sharp
counts and prime realizations*. This is a manuscript development step,
not a claim of journal acceptance or an independent referee assessment.

The argument now has the following complete dependency structure:

| Manuscript label | Statement and role | Attribution and scope |
| --- | --- | --- |
| `thm:main` | Sharp pointwise set count C_d(n)=Theta_d(n^(d-1)(log n)^(d-1)), including the prescribed-size primitive-normal injection | Wagner's upper bound and the known planar order are credited. The asserted advance is the matching lower bound in d>=3. |
| `lem:finiteweight` | Ordered positive-weight prefixes equal finite initial intervals of global monomial orders | Known consequence of Robbiano. The proof uses finite unit, consecutive-order and minimal-excluded-border comparisons. The finite prefix is essential. |
| `thm:orders` | The full ordered-prefix count O_d(n) has the same sharp growth order | The family and its enumeration are credited to Vasiliev--Pavlov. A bounded-width dyadic product shell and central arrangements give the uniform upper bound; the main set theorem supplies the lower bound. |
| `lem:linear-score` | Positive integer weights with n-1<=S<=c_d n and max A_i<=S+1, c_d=d(d!)^2 2^(d-1) | Supporting LP/Cramer/product-bound lemma, with no separate claim of a new general compression method. The proof handles unused coordinates explicitly. |
| `thm:height` | Worst-case minimum log(max p_i) equals Theta_d(n log(max(B,2))) for n>=2, B>0, for ordered words and underlying sets | Qualitative prime realization is Vietri's. The upper bound uses Dusart's published Proposition 5.4, the integer certificate and an exact anchor; a pure-axis word gives the matching lower bound. |
| `cor:budget` | Imposing the explicit linear-bit prime budget leaves both complete pattern families unchanged | A consequence of the two counts and the prime-height theorem, not a prime-distribution estimate or efficient construction algorithm. |

The earlier qualitative prime-realization proof is replaced in the paper
by its attribution and the stronger quantitative argument. The separate
structural note, original occurrence rule and smooth-number CLI remain
unchanged. No additional search framework or mandatory dependency is
introduced.

#### Internal audit during integration

The newly typeset proofs were checked against the saved arguments, with
the following boundary cases made explicit:

* A finite global-order prefix is controlled by finitely many differences,
  including all unit vectors and the minimal excluded border. Positivity
  then controls the infinite complement; no scalar representation of the
  whole infinite order is asserted.
* The ordered count restricts to active coordinate faces before applying
  volume bounds. Its product lies in a bounded-width dyadic shell, which
  supplies exactly the logarithmic exponent d-1. Chart intersections and
  repeated hyperplanes only overcount.
* In the certificate LP, the excluded border is taken in the active face.
  Compact sublevel sets give an optimal vertex. Cramer scaling preserves
  every unit-margin inequality. Active weights are at most S because
  their unit vectors are included; inactive weights are restored as S+1.
* Axis collisions at lcm(A_i,A_j)<=S are impossible because included
  scores are strict. Thus all except possibly the smallest weight exceed
  sqrt(n-1), even when some coordinates are unused.
* The smallest-weight prime is anchored exactly. All other logarithmic
  errors are positive and their total on any included point is below one
  scaled integer-score gap. Distinctness of primes holds for equal larger
  weights, and positivity of the errors controls every excluded point.
* Dusart's published interval is used only at x>=exp(30), where its
  1/(log x)^3 width is strictly smaller than 1/[25(log x)^2]. The source
  statement and publisher metadata were rechecked; its full analytic
  proof is an established input, not independently rederived here.
* The theorem is uniform in B>0 for fixed d; n=1 is handled separately.
  The same pure-axis family proves the word and set lower bounds.

The primitive-window, affine-patch, height-to-rank and cross-aspect steps
of the main set-count argument were also checked for consistency with the
new consequences. The two new theorems do not alter that proof. No
unresolved mathematical step was identified in this internal pass.

The cited Robbiano source is the 1985 EUROCAL chapter, with the inspected
ordered-group formulation in MPI/SFB 85-17, Theorem 2.5, pp. 23--25.
Its 2005 online date is not used as the original publication date.
The Dusart citation is the 2018 published paper, not the differently
numbered 2010 preprint. Both sources retain their earlier recorded
reading-scope qualifications.

The current finite evidence includes 230 exact integer-score fixtures,
90 matrix-order fixtures and the 32 published-table comparisons with
1,023 accepted strict-prefix certificates. The chart, lattice-patch and
primitive-sieve checks remain supporting falsification evidence. None of
these numerical checks proves a uniform asymptotic theorem or priority.

The original goal is not reduced to manuscript preparation. Its final
contribution and submission assessment must evaluate the complete set
of statements against the strongest inspected primary work and the
original acceptance criteria. The full-source comparison and the newly
integrated proofs make that assessment possible; typesetting alone does
not establish its outcome.

### Wagner's stronger flag bound and the ordered extension

The final comparison re-read the author's thesis Sections 3.4--3.5,
especially Lemma 3.11, Corollary 3.12, Theorem 3.15 and the subsequent
flag argument. The product-bounded difference-matrix encoding controls
the sets and face flags of the single size-n corner-cut polyhedron.
We explicitly credit this stronger scope. The main lower-count result
does not improve that upper order.

Ordered words specify initial sets at all lengths through n. They cannot
be identified with maximal face flags at their final-cut vertex. One
exact planar example at n=50 has row lengths

    (10,9,8,7,6,4,3,2,1)

as x increases from 0 to 8. The five weights

    (60,49), (84,71), (112,97), (144,127), (18,17)

give five distinct strict words with precisely this set. Its moment
vertex is (130,155), whereas a vertex of a full-dimensional planar
polyhedron has two incident maximal face flags. A direct finite-box
sort checks all 51 leading scores for each weight, independently of
the heap used in the crossing sweep. The 50th/51st scores are

    (485,490), (691,704), (933,948), (1211,1228), (157,158).

Thus the naive identification preserving the final vertex is false.
This does not exclude a different deduction of the same asymptotic
upper bound from established methods. The manuscript keeps its direct
dyadic proof and does not treat the flag distinction as an independent
publishable result. The exact shape, weights, scores and word hashes
are saved under `planar_ordered_prefix_checks.flag_vertex_control`.

### Final internal contribution and acceptance assessment

The research objective is evaluated on the main prescribed-cardinality
theorem, with the ordered and arithmetic results as additional quantified
consequences. The argument is uniform for each fixed dimension and does
not rely on the isolated checks as proof.

| Original requirement | Evidence and assessment |
| --- | --- |
| A substantive quantitative advance arising from smooth-exponent geometry | `thm:main` gives the sharp order C_d(n)=Theta_d(n^(d-1)(log n)^(d-1)) for every fixed d>=2. In d=3 it gains one logarithmic factor over Wagner's lower bound; in general it gains (log n)^(d-2). The upper bound retains its original attribution. |
| A complete uniform proof, including duplication control | The primitive-window lemma, division by bounded integer height, prescribed-size two-patch injection and disjoint observable aspect ranges are all proved in the manuscript. Closed constants and thresholds specify the quantifiers. Partial boundary ties and inactive dimensions are handled. |
| Comparison with the strongest relevant primary statements | The author thesis proof and stronger flag scope, full Muller paper, later Kalka--Kuhnt bound statement, Elizalde--Galvan v1, full Vasiliev--Pavlov paper, Vietri's exact realization theorem, primitive-matrix extensions and digital-plane/threshold comparisons are recorded with their reading scopes. The retained same-object bounds in later primary sources are positive comparison evidence; an empty search result is not used to establish novelty. |
| A defensible novelty claim for the main result | The main lower bound improves the same cardinality problem in the inspected lineage and supplies the missing prescribed-size mechanism. No compared theorem gives that lower order. This supports presenting it as a new result in a submission, with the comparison made explicit. It is a research judgment based on the specified primary audit, not a claim to know every unpublished or unseen work. |
| A useful consequence for an existing problem | The main theorem determines the mean coordinate-slice completion multiplicity raised in Wagner's discussion. The ordered theorem quantifies the established strict finite-ordering class, and the prime-height theorem gives a sharp size budget that realizes every such pattern. Known polytope/model correspondences are credited, not claimed as new applications. |
| Arithmetic information for any prime application | The prime-height proof uses the published unconditional Dusart interval, exact anchoring and the collision property of the integer certificate. Its matching axis lower bound makes the worst-case order optimal. No new prime-gap estimate or fixed-prime speedup is claimed. |
| Reproducible supporting evidence | The three existing research scripts save exact witnesses, commands and source hashes. They check determinants and modular ranks, both affine patches, strict tie handling, exact compression LPs, independent scalar orders, published finite counts and the flag-control counterexample. The saved reports state their finite scope. |
| An honest complete manuscript | The TeX source contains every new proof and all supporting lemmas, exact prior-art attribution, defined counting parameters, limitations and an explicit AI-assistance disclosure. Its PDF is compiled and rendered for layout review; the publication record contains its exact hashes. |
| Preserve the original toolkit and separate structural findings | The CLI, C++ kernels, build and smoke sources, and separate Mobius manuscript are unchanged. Existing side results do not substitute for the quantitative objective. |
| Observe the final-version assumption and communication constraint | Both Elizalde--Galvan versions remain linked, with the final-version check explicitly deferred under the user's assumption. No email, author query, message to another task or journal submission has been sent. |
| Save only reviewed public artifacts on main | The publication review checks an explicit file list, source hashes, anonymous PDF metadata and private-file exclusions before committing and verifying the main-branch push. Source-paper PDFs, correspondence and local reading records are excluded. |

**Contribution judgment:** the complete prescribed-size count provides
a substantive answer to an existing discrete-geometric counting problem.
The general ordered result and the optimal prime-size guarantee make a
coherent extension of the original exponent-sequence investigation.
This is sufficient mathematical content for a focused submission; there
is no need to advertise a new general density theorem, a prime-generation
method or an unranking improvement.

**Publication judgment:** the manuscript can be put before referees with
the stated contribution and explicit literature boundary. Journal
acceptance, independent verification and broader recognition have not
occurred and are not promised. The anonymity and assistance disclosure
are intentional. Actual external submission is not authorized by the
standing prohibition on outbound messages.

The final validation state and the requirement-by-requirement audit are
recorded in `results/research/quantitative_goal.json` and
`results/research/corner_cut_growth_audit.json`. Their completion decision
must be made only after inspecting the current reports, compiled artifact,
source hashes and reviewed public Git state.

## Higher-impact research plan: universal algebraic complexity

### Ambition and precise question

The next objective is a general theorem explaining how algebraic relations
change the complexity of monomial orderings. The smooth-prefix theorem is
the generic geometric starting point. The larger class consists of arbitrary
finite-dimensional polynomial quotients, including ideals of distinct points
and ideals with multiplicities. This is a change in the mathematical research
scope requested by the user; the toolkit's fixed-prime rank/unrank task remains
unchanged.

For a field K and I contained in K[x_1,...,x_d], put

    length(I) = dim_K(K[x_1,...,x_d]/I),
    g(I) = #{in_prec(I) : prec is a global monomial order},
    G_d(n) = max g(I), over all fields K and ideals of length n,
    F_d(n) = max g(I(P)), over all fields and n distinct K-rational points.

Keep G and F separate: F_d(n) <= G_d(n). This counts different monomial
initial ideals, equivalently initial staircases in the state-polyhedron
framework. It does not count the polynomials in a single basis or the time
required to compute a basis. No generic maximality assumption is made.

The long-term target is the sharp asymptotic growth of G_d(n) in each fixed
dimension, together with a mechanism explaining when dependencies increase
the generic count. Begin with d=3; use d=2 as a proof diagnostic. A theorem
for every n-point configuration is also a substantial target, but its radical
scope must be explicit rather than silently promoted to all colength-n ideals.

The current manuscript and the inspected primary bounds give the starting
comparison

    C_d(n) <= F_d(n) <= G_d(n) = O_d(n^(2d(d-1)/(d+1))),
    C_d(n) = Theta_d(n^(d-1)(log n)^(d-1)).

Thus the currently checked three-variable range is

    Omega(n^2(log n)^2) <= F_3(n) <= G_3(n) <= O(n^3).

These are verified starting statements, not a certification that later work
has not improved the unrestricted bounds. The first milestone must establish
the strongest relevant primary baseline. The candidate equality G_d(n) of
the same order as C_d(n) is a question to test, not an adopted conjecture.

### Established framework and first comparison gate

* [Onn--Sturmfels, Cutting corners (1999)](https://doi.org/10.1006/aama.1999.0645),
  Corollary 5.4 and Proposition 5.6, define the point-configuration problem
  and already show that generic configurations are not exact maximizers:
  C_2(7)=8 whereas F_2(7)=10. The previous exact reproduction is a mandatory
  negative control for any transfer argument.
* [Babson--Onn--Thomas (2003), inspected preprint](https://arxiv.org/pdf/math/0207135),
  Corollary 2.5, extends the displayed upper exponent to all colength-n
  ideals. Theorem 3.5 supplies a common Hilbert-zonotope refinement;
  Theorem 4.2 supplies a fixed-dimension polynomial-time algorithm.
  The projected-matroid and greedy-basis framework is already present.
  Rephrasing it, or providing another universal list without improving a
  proved bound, is not the proposed contribution.
* [Zhang--Stigler, inspected v1](https://arxiv.org/pdf/1907.01080v1),
  Section 5, supplies a finite-field specialization. The sixteenth pass
  additionally inspects the later NSF deposit's Theorem 8 and its proof
  through indexed pages 15--16. Its bound recovers the old exponent when
  p>n. Fixed field size is not the maximum over all fields; a complete
  publisher-version comparison is not claimed.
* [Dimitrova--He--Robbiano--Stigler, inspected v2](https://arxiv.org/pdf/1704.03680v2),
  Theorem 4.8 and Corollary 4.11, give the grid-complement invariance.
  Taking a complement inside a grid cannot amplify the basis count.

The new comparison must trace later uses of the extremal bound and its
proof ingredients, including lattice-polytope vertex bounds, state polytopes,
universal monomial-order lists, and algorithms for enumerating Groebner fans.
For each potentially stronger result record its exact parameter, field,
radicality, dimension and input model. Inspect the actual theorem and relevant
proof. Search snippets and failure to find a paper do not establish novelty.
The current comparison is partial; no new open-problem claim is certified here.

### Proof programme

1. **Locate the loss in the general bound.** Reconstruct the established
   state-polyhedron upper proof and determine which lattice-polytope estimates
   ignore constraints imposed by multiplication in a finite algebra. Record a
   specific inequality that a new lemma would have to improve. Do not assume
   the generic corner-cut fan refines the fan of every ideal; the known
   seven-point example rules this out.

2. **Use a common exact model for dependencies.** Represent a length-n
   quotient by commuting multiplication matrices M_1,...,M_d and the cyclic
   vector representing 1. Residues of monomials are M_1^a_1...M_d^a_d v.
   Greedy linear independence under a monomial order recovers the standard
   monomial basis. This is established finite-algebra machinery, not a new
   lemma. Any proposed lower-bound construction must actually define an
   algebra: check commutation, cyclicity, dimension, and the claimed radical
   or point-configuration properties. Arbitrary vector matroids cannot be
   silently substituted for realizable polynomial quotients.

3. **Test two proof routes before expanding computation.** For an upper
   bound, try to charge extra state vertices or subdivisions to minimal
   monomial-residue dependencies and use their multiplication constraints to
   improve the general count. Any structural parameter must be explicitly
   defined independently of g(I). For a lower bound, seek a formula for an
   infinite family whose dependencies create provably many distinct positive
   normal cones, with exact witnesses and control of duplicate staircases.
   Tensor, union, deformation or complement operations need a proved
   fan/count transformation before being treated as amplification mechanisms.

4. **Demand an asymptotic theorem for a broad class.** The leading target
   is a sharp order for G_3(n), followed by general fixed d if the mechanism
   supports it. An intermediate success must substantially advance the
   strongest checked bound for G_3 or all-point F_3. A power improvement
   O(n^(3-epsilon)), epsilon>0, over the currently checked universal upper
   bound is an example, conditional on the literature gate. An explicit
   family growing as Omega(n^(2+epsilon)) would similarly separate it from
   the current generic order, but must also improve the strongest relevant
   published lower bound. A theorem giving sharp growth with the same
   exponent and the correct logarithms is another valid major milestone.
   Small examples, an extra finite-field case, or a renamed known result do
   not meet this gate.

5. **Establish a consequence beyond recounting the theorem.** Seek an
   improved bound or a proved limitation for universal order selection or
   Groebner-fan enumeration. A smaller number of outputs alone does not give
   a faster algorithm: such a claim requires an enumeration procedure, an
   input/output model and a complete complexity analysis. If an algorithmic
   consequence does not follow, identify a separate precise geometric or
   algebraic consequence. Do not advertise an implication for prime
   distribution or fixed-prime random access without an arithmetic theorem.

6. **Prepare a separate honest manuscript only after these gates.** Supply
   complete uniform proofs, exact prior-art comparisons, independently
   implemented finite checks and certificates, and explicit unresolved
   points. Preserve the current corner-cut manuscript as its own result.
   Internal assessment is distinct from independent refereeing or acceptance.
   No correspondence, messages or submission is authorized.

### Bounded first milestone and decision rule

The first milestone consists of three recorded passes: a theorem-level
baseline audit, an attempted upper-bound lemma, and an attempted infinite
lower-bound construction. Each proof attempt must state its hypotheses and
the precise obstacle or falsifying witness if it fails. The existing
seven-point split and grid-complement theorem must be respected.

After those passes, record whether there is a viable uniform mechanism,
whether the target is already settled, or which specific obstruction remains.
Do not expand searches in n, characteristic or dimension without a conjecture
and a discriminating test. A failed attempt is not evidence that the full
problem is intractable. An inconclusive feasibility report does not complete
the theorem goal; it triggers an explicit assessment and a proposed redirect.

The arithmetic route remains a possible later redirect, not simultaneous
open-ended work. It would first need a precise fixed-prime theorem exceeding
the existing log-separation plus fixed-dimension lattice-counting baseline,
with a concrete rank/unrank consequence. Merely changing prime bases, adding
cycle examples or demonstrating fixed-dimension tractability is insufficient.

The plan is saved in `results/research/structural_complexity_goal.json`.
Its status records research progress, not an assertion of mathematical impact.
This planning addition changes no proof or executable input of the completed
manuscript; earlier computational reports retain their original source hashes.

### First feasibility pass: two rigorous route obstructions

The previous goal turn completed the requested goal registration and saved the
plan. That is progress on setup, not a new mathematical advance. This pass
adds two uniform negative controls and an independent exact quotient checker.
Neither control meets the new theorem goal, and neither is presented as a new
publishable result. They rule out specific proof mechanisms.

#### Upper-bound attempt: shrinking the full state-polytope volume fails

The proof of Babson--Onn--Thomas Corollary 2.5 encloses all staircase moments
in [0,n^2]^d and applies Andrews's lattice-polytope vertex bound. It is natural
to ask whether multiplication constraints force a smaller full volume. They
do not: already a generic point configuration attains order n^(2d).

Here is a uniform exact witness. Fix d>=2, r>=1, and n=binom(r+d,d).
In a single generic n-point ideal all corner cuts are initial staircases,
by the established generic-point correspondence. In particular, the d pure
axis chains and the degree simplex

    Delta_r = {a in N_0^d : a_1+...+a_d <= r}

are initial staircases. The axis moments are A e_i with A=n(n-1)/2.
By counting weak compositions, the degree-simplex moment is B(1,...,1),
where B=nr/(d+1)=binom(r+d,d+1). The simplex formed by these d+1 moment
points has volume

    A^(d-1) (A-dB) / d!.

Indeed its edge matrix is A Id - B 11^T, whose determinant is
A^(d-1)(A-dB). For fixed d, r=Theta_d(n^(1/d)), hence dB/A tends to zero
and the displayed volume is Theta_d(n^(2d)). All of these points belong to
the convex hull P of initial-staircase moments for the same ideal. The
existing cube gives the matching volume upper order. The pure axis chains
and Delta_r are separated sets even though Delta_r has internal degree
ties; no internal ordering assertion is needed.
There is also a closed lower estimate: Delta_r contains 0 and the dr nonzero
pure-axis points, so n-1>=dr and dB/A<=2/(d+1). Its moment-simplex volume
is therefore at least ((d-1)/(d+1)) A^d/d! for every r>=1.

Consequently a universal o(n^(2d)) bound for the full volume is false,
even just along these infinitely many exact cardinalities and even for
radical point ideals. Applying the same vertex-versus-volume inequality
to a supposedly smaller whole polytope cannot improve the exponent. A
successful argument needs finer control of the positive-normal boundary,
its local caps, or the algebraic restrictions on its vertices. This does
not rule out an improved vertex count or a different use of volume.

The exact checker enumerates the simplex points independently of the moment
formula, then computes the determinant by rational elimination in dimensions
2--6 at r=1--6. The finite checks support the displayed uniform derivation;
they do not prove the asymptotic statement by extrapolation.

#### Lower-bound attempt: disjoint socle relations cannot amplify enough

Let D be a finite monomial order ideal of size N in three variables, and let
S be its coordinatewise maximal points, of size m. Choose disjoint subsets
S_1,...,S_q of S, each of size at least two, and nonzero coefficients. Set

    J_D = (x^a : a outside D),
    f_j = sum_{a in S_j} c_{j,a} x^a,
    I = J_D + (f_1,...,f_q).

Every nonconstant monomial multiple of each f_j lies in J_D, since its
support is in maximal points. Thus the additional relation space is exactly
the q-dimensional span of the f_j, and length(I)=n=N-q. Disjoint supports
make the relations independent over any field. For a generic positive
weight w, the initial staircase is exactly

    D minus {argmax_{a in S_j} w.a : 1<=j<=q}.

To see this directly, a polynomial in I has, modulo J_D, a scalar linear
combination of the f_j. Its possible leading monomials in D are the maxima
of the disjoint supports. They cannot cancel between different relations;
all their nonconstant multiples are already in J_D. This proves the formula
for the entire ideal, not merely for a proposed generating set.

For one relation this gives at most m initial staircases. For all disjoint
relations together it gives O(m^2) in three variables. Each maximum is a
vertex of P_j=conv(S_j). A polytope with s vertices in dimension at most
three has at most 3s edges, with the lower-dimensional cases included.
The normal fans of all P_j therefore have a common refinement by at most
3m central edge hyperplanes. On the positive weight simplex these are at
most 3m lines, making at most 1+(3m)(3m+1)/2 regions. Every possible
initial staircase is realized by a full-dimensional region. Other global
monomial orders introduce no extra staircase: the finite strict comparisons
for their reduced bases admit a positive scalar weight, as in the established
finite-order representation.

For completeness, m=O_d(N^((d-1)/d)) has a simple chain-packing proof in
every fixed dimension. Choose h=ceil(N^(1/d)). There are at most h^(d-1)
maximal points in [0,h-1]^d, because a vertical line contains at most one
point of an antichain. For each coordinate i, maximal points with a_i>=h
give disjoint h-point downward chains parallel to that coordinate. Their
number is at most N/h. Taking the union over i yields

    m <= h^(d-1) + dN/h = O_d(N^((d-1)/d)).

In the present construction q<=m/2<=N/2, so N<=2n. Therefore

    g(I) = O(n^(4/3))  (three variables, disjoint socle relations).

This is much smaller than the generic n^2(log n)^2 order. Arbitrarily many
disjoint relations of this type cannot provide the desired lower family.
The bound does not cover arbitrary overlapping relation spaces on S, or
relations whose multiples survive in the quotient. Those remain possible
mechanisms. The socle and convex-hull bounds are elementary tools, not new
general compression or matroid theorems.

#### Independent finite verification

Run

```bash
python3 -B scripts/run_structural_complexity_checks.py
```


The script checks 30 volume fixtures and 18 socle-quotient fixtures. For each
quotient it constructs multiplication matrices, checks every commutator on
a basis, and verifies cyclic monomial evaluation on D and its border. Exact
polygon clipping enumerates the positive-weight chambers of all within-support
comparison planes. A separate line-intersection recurrence checks the region
count, including concurrency. Each resulting initial staircase is independently
recovered by greedy linear algebra over Q, F_2, F_3 and F_29. The complete
finite run has 398 chambers and 1,592 greedy checks. It saves defining order
ideals, relation supports and integer weight witnesses, with source hashes,
in `results/research/structural_complexity_checks.json`.

#### Prior-art update and next discriminating question

The unrestricted bound remains the strongest inspected statement, not a
claim to have excluded every subsequent improvement. Additional sources
already constrain possible algorithmic claims:

* [Maruri-Aguilar, Universal Groebner Bases for Designs of Experiments
  (2005)](https://rendiconti.dmi.units.it/volumi/37/05.pdf), Theorem 5.1
  and Section 6, give the planar Hilbert-zonotope normal fan directly from
  primitive vectors in [0,n-1]^2. The construction removes the logarithmic
  overhead of the previous planar universal-order construction. This is
  established prior art; generating that same list is not a new consequence.
  The theorem statement and indexed proof/discussion passages were inspected,
  not the entire paper or all of its proof details.
* [Deza--Manoussakis--Onn, inspected 2015 preprint](https://optimization-online.org/wp-content/uploads/2015/12/5266.pdf),
  Theorem 3.3 and its adjacent propositions, bound projections of arbitrary
  matroid polytopes in terms of a fixed coordinate range p, using Euler
  zonotopes. Their general matrices need not be multiplication-compatible
  monomial residue data. The theorem statement, parameter definitions and
  displayed construction were inspected, not the full paper.
* [Hausbrandt--Bachtler--Ruzika--Schaefer, inspected 2023 preprint](https://arxiv.org/pdf/2310.05147),
  Section 2, recalls the one-parameter matroid breakpoint bound
  O(M r^(1/3)) for ground-set size M and rank r, with references to Dey and
  Eppstein. The article's principal problem is interdiction, a different
  problem. This excerpt is a guide to the original sources, not a new
  two-parameter or quotient-algebra bound.

The next proof question is whether overlapping, multiplication-compatible
dependencies admit a useful bound on the additional positive-normal boundary
pieces, or an explicit realizable family forcing faster growth. The two
controls above show why shrinking the full volume and adding independent
socle gadgets do not answer that question. The strongest-current-bound audit
and the unrestricted asymptotic theorem remain unfinished. The goal stays
active, and these auxiliary controls do not change the completed manuscript.

### Second feasibility pass: interacting relations and a conditional count

The preceding turn made progress by proving and testing two specific route
obstructions. This pass strengthens one obstruction to all socle relation
spaces and proves a conditional counting lemma that permits non-generic
dependencies. These are supporting results for the unrestricted objective;
neither the objective nor its acceptance criteria is replaced by this class.
Priority of the conditional result is not established.

#### All socle relation spaces have a small moment budget

Keep a finite order ideal D in d variables, of size N, and its maximal-point
antichain S, of size m. For each coordinate i, the downward coordinate chains

    {a-t e_i : 0<=t<=a_i},  a in S,

are disjoint subsets of D: two maximal points cannot have the same other
coordinates. Consequently

    sum_{a in S}(a_i+1) <= N,
    sum_{a in S} a_i <= N-m.                         (SC1)

Let R be any r-dimensional linear subspace of the span of the socle
monomials x^a, a in S, over any field. No disjoint-support assumption is made.
Then I=J_D+R is an ideal, because every nonconstant monomial multiple of R
lies in J_D. Its length is n=N-r. Every standard monomial basis consists
of the fixed set D minus S and a basis of the vector configuration S/R.
The latter statement follows either by row reduction of R in the chosen
monomial order or by greedy independence of the residue columns.

Translate the state moments by the fixed sum of D minus S. The remaining
moments are sums of subsets of S. By (SC1), each coordinate lies between
0 and N-m<=n. Thus their convex hull is an integer polytope in [0,n]^d.
The established projected-matroid/state-polyhedron correspondence identifies
each initial staircase with a distinct positive-normal vertex of this
polytope. Andrews's lattice-polytope vertex bound, used in the original
general proof, now gives

    g(I) = O_d(n^(d(d-1)/(d+1))).                    (SC2)

If the affine dimension is smaller, apply the bound in that dimension; its
exponent is no larger. In three variables (SC2) is O(n^(3/2)). The earlier
O(n^(4/3)) bound remains stronger for disjoint one-relation blocks. In every
fixed dimension, the exponent in (SC2) is strictly less than d-1. Therefore
even arbitrary overlapping socle relations cannot beat the existing generic
corner-cut order. A viable larger lower construction must use relations
whose nonconstant multiples also survive.

The antichain budget and application of the known vertex theorem are
elementary supporting tools; no separate novelty claim is made for (SC2).
The checker constructs six genuine quotient algebras with overlapping
relation spaces. It checks multiplication, cyclicity and the moment budget,
then compares descending relation-row reduction with ascending residue-column
greediness on every exact weight chamber. Characteristics are kept separate:
for one 32-dimensional quotient the counts are 18 over Q and F_29, 10 over
F_2, and 14 over F_3. Reducing a rational configuration modulo a prime does
not generally preserve its entire fan.

#### Conditional count from rectangular spanning

Let A=K[x_1,...,x_d]/I have dimension n. Suppose an integer kappa>=1 has
the following property:

    Every monomial box B(L)={a:0<=a_i<L_i} of volume
    L_1...L_d >= kappa*n has residues spanning A.     (RS)

This condition is stated in terms of quotient ranks, independently of the
number of initial ideals. It allows dependencies among particular n-element
staircases. Then

    g(I) <= W_d(M),  M=d^d*kappa*n,                 (RS1)

where W_d(M) denotes the number of geometric ordered prefixes. In particular,
using the completed ordered-prefix theorem,

    g(I) = O_d((kappa*n)^(d-1)
                 (log(kappa*n+1))^(d-1)).          (RS2)

Proof. Choose a positive rationally independent weight w and a threshold t
strictly between its M-th and (M+1)-st lattice scores. Put

    L_i=floor(t/(d*w_i))+1.

The box B(L) lies strictly inside the t-cut. Moreover

    M <= product_i(floor(t/w_i)+1)
      <= d^d product_i(floor(t/(d*w_i))+1)
       = d^d |B(L)|.

Thus |B(L)|>=kappa*n, and (RS) shows that the first M monomials already
span A. Greedy independence on this ordered word recovers its n standard
monomials. For fixed I, identical words give identical standard bases, so
the number of bases is bounded by the number of words. Every global
monomial initial ideal is represented by such a generic positive weight,
by the established finite-order representation. This proves (RS1), and
the previously proved word count proves (RS2).

There is a precise consequence for universal order selection in this class.
For fixed d and kappa, choose one order for each length-M geometric word.
This list, independent of I, recovers every initial ideal of every quotient
satisfying (RS). Its cardinality is O_d,kappa(n^(d-1)(log n)^(d-1)). The
class contains generic point ideals, since a box of volume at least n
contains an n-element staircase and every such staircase is basic in the
generic case. A single generic ideal already has C_d(n) distinct initial
ideals. Hence any universal list for this class has the matching lower
order. The minimum list size is therefore

    Theta_d,kappa(n^(d-1)(log n)^(d-1)).            (RS3)

This is a cardinality result, not a new enumeration-runtime claim. Known
universal-order constructions and the prior-art boundary still need to be
compared before treating it as a paper contribution.

The condition is strictly broader than the full generic staircase model:
the published seven-point configuration over F_29, with ten initial
staircases rather than the generic eight, satisfies (RS) with kappa=1.
The checker verifies all 49 boxes with sides at most seven. Each coordinate
axis gives a seven-element basis, so any larger box also spans; this finite
check proves (RS) for that exact configuration. The resulting positive
weight tests complete rank by indices 8, 7 and 7, respectively, within the
guaranteed length 28. The first case explicitly contains a rejected monomial.

#### Why rank delay alone cannot settle the unrestricted problem

The parameter kappa cannot be bounded uniformly for all point configurations.
For n distinct diagonal points (t_j,...,t_j), the rank of B(L) is

    min(n, L_1+...+L_d-d+1),

by ordinary Vandermonde independence. Deficient boxes have sum_i L_i <=
n+d-2; their largest volume is attained by balanced side lengths and has
order n^d. The smallest integer kappa in (RS) consequently has order
n^(d-1), even though this ideal has exactly d initial staircases for n>=2.
All linear relations x_i-x_j eliminate every variable except the least one
in the order, leaving its powers 0 through n-1. Thus large rank delay does
not imply a large fan.

Even choosing a different representing weight for each staircase cannot
make rank completion uniformly linear in n. Set

    D_L={0} union {j e_i:1<=j<=L,1<=i<=d},  n=1+dL.

Over Q, the same points form a Cartesian configuration. Its ideal has a
Groebner basis for every monomial order given by x_i*x_j (i!=j) and the
univariate polynomials product_{j=0}^L(x_i-j). Their leading monomials
are x_i*x_j and x_i^(L+1) for every global order, so D_L is the unique
standard staircase. Indeed, the displayed leading monomials leave exactly
n standard monomials, while the polynomials vanish on n distinct points;
the matching dimensions prove the Groebner and point-ideal assertions.
Translating each coordinate by L+1 preserves all
initial ideals and makes every point coordinate nonzero. Translation
preserves leading monomials because its additional terms are divisors of
the original terms.

Every geometric word containing D_L reaches score at least L*max_i w_i.
All points with coordinates at most floor((L-1)/d) have strictly smaller
score. Every representing word therefore has length at least

    (floor((L-1)/d)+1)^d = Omega_d(n^d).

This is a radical point-ideal family with g(I)=1 and no monomial evaluating
to zero identically on its points. It rules out extending (RS2) simply by
asserting a short representative word for every initial staircase. The
proof needs a way to separate a large forced staircase from genuine
variation; it cannot rely solely on the time at which rank becomes n.

#### Additional primary comparison and remaining scope

[Gallardo--Schmidt, arXiv:2205.15238v2](https://arxiv.org/html/2205.15238v2),
Section 6.2, Definitions 6.3--6.4, Theorem 6.5 and its proof, and Lemmas
6.6--6.7, have now been inspected. They use the established generic
corner-cut correspondence to compute selected flat limits and GIT indices.
Those statements do not give an unrestricted extremal count in our sense;
the GIT linearization walls are a different parameter problem. The rest of
the paper has not received a full proof audit. The primary SMU dissertation
page is readable, but direct access to its download returned HTTP 403; its
counting chapter is not claimed as read.

The new exact tests cover 517 rectangular boxes, nine geometric-prefix
checks, six overlapping-socle quotients and nine translated Cartesian
rank-delay fixtures. The report retains the earlier controls, with separate
field-dependent witnesses. The unrestricted asymptotic advance, full
priority assessment and target manuscript are still missing. The goal
remains active. The next construction must have surviving multiples of
relations; the next upper argument must control variation despite a large
fixed staircase.

### Third feasibility pass: binomial quotients and the boundary of the method

**Status:** a uniform bound for a substantial restricted class follows from
existing tools. This is not the unrestricted theorem, not a sharpness claim,
and not an established novelty claim. The larger goal remains active.

#### A count for every finite binomial quotient

Let I be any binomial ideal of colength n in K[x_1,...,x_d], where d>=2
is fixed and K is any field. Monomial generators are allowed. Then

    g(I) = O_d(n^(d-1)(log(n+1))^((d-1)^2)).       (BQ1)

In three variables this is O(n^2 log^4(n+1)), compared with the inspected
unrestricted bound O(n^3). The assertion covers nonradical quotients and
relations whose nonconstant multiples survive. The generic corner-cut
lower bound does not automatically apply to this subclass, so (BQ1) is
not asserted to be optimal.

Here is the complete reduction. By
[Eisenbud--Sturmfels, *Binomial ideals*, Proposition 1.1](https://arxiv.org/pdf/alg-geom/9401001),
fixing any reference monomial order gives a standard monomial basis
b_1,...,b_n in which every monomial residue is zero or a nonzero scalar
multiple of exactly one b_j. This statement uses neither radicality nor
an assumption about the characteristic.

Every member a of an n-element staircase has all its divisors in that
staircase. Consequently all possible standard monomials belong to

    V_n^d = {a in N^d : product_i(a_i+1)<=n}.

Write N=|V_n^d|. The elementary recurrence

    |V_n^d| = sum_{k=1}^n |V_floor(n/k)^(d-1)|

and the harmonic-sum bound give N<=n(1+log n)^(d-1). This is the established
hyperbolic-cross ground set used by Babson--Onn--Thomas, not a new ground-set
reduction.

Partition the nonzero residues on V_n^d into the n classes S_j that map to
the lines K*b_j. Every class is nonempty, since the reference standard
basis lies in V_n^d. Zero residues are loops and can be discarded. A set
of n residue columns is a basis precisely when it chooses one monomial
from each S_j. Therefore its moment polytope is exactly

    B(I) = sum_{j=1}^n conv(S_j).                 (BQ2)

For a weight generic on all within-class comparisons, the unique
minimum-weight residue basis chooses the minimum-weight member of each
class. The usual finite-comparison representation of a global monomial
order by a positive weight applies here. Thus each monomial initial ideal
gives a positive exposed vertex of (BQ2). Distinct coherent standard
bases have distinct moment sums: if two had the same sum, the weight
selecting one uniquely would make the other equally good. In particular,

    g(I) <= f_0(B(I)).

The total number of input vertices in (BQ2) is at most N. The needed
published geometric result is
[Weibel, *Maximal f-vectors of Minkowski sums of large numbers of polytopes*,
Corollary 2 and Section 5](https://arxiv.org/pdf/1002.0155): in fixed dimension
d, a Minkowski sum with M input vertices in total has O_d(M^(d-1))
vertices. The displayed proof treats full-dimensional summands in general
orientations; its worst-case conclusion suffices here.

For completeness the dimension issue can be removed directly. Ignore
point summands, which only translate the result. If fewer than d summands
remain, the product bound already gives O_d(M^(d-1)). Otherwise add a
small full-dimensional simplex to each summand. This refines the normal
fan and multiplies the input-vertex budget by at most d+1. An arbitrarily
small generic perturbation preserves the finitely many strict supporting
witnesses for existing vertices. Apply the published bound to the
perturbed full-dimensional sum. This proves the required O estimate even
with lower-dimensional or coincident summands. Combining it with the
bound on N proves (BQ1).

This proof is short because its ingredients are already substantial
published theorems. It is recorded as their consequence. An unsuccessful
search for this exact combined formula would not make it a new paper.

#### Actual algebra models and an exact obstruction to automatic extension

The checker now uses the kernels of

    K[x,y,z] -> K[t]/(t^n-epsilon),
    x |-> t,  y |-> c_y*t^a,  z |-> c_z*t^b,

where epsilon is zero or one and c_y,c_z are nonzero. These kernels are
generated by x^n-epsilon, y-c_y*x^a and z-c_z*x^b, so their colength is
exactly n. All selected fixtures have a+1<n, making x*(y-c_y*x^a) a
relation between nonzero monomial residues. They are not confined to the
socle. The epsilon=0 models are nonradical; epsilon=1 models remain valid
when the characteristic divides n, when they too are nonradical.

For lengths 4,5,6,7,9 and two models per length, the script independently
constructs commuting multiplication matrices with cyclic vector 1. It
compares their recursive monomial evaluation with the closed residue-class
formula. All positive-weight chambers of the within-class comparison
arrangement are covered by exact rational polygon clipping, with the
number of chambers independently checked by line intersections. In every
chamber, Gaussian independence over Q, F_2, F_3 and F_29 is compared with
classwise minimum selection. Exact standard staircases, moment sums and
integer weight witnesses are saved. The tests support the reduction; the
uniform proof above is what establishes the bound.

The difficulty for all ideals is already visible for four ordinary
distinct points, not just exotic nonradical schemes. Over F_101 take

    P = {(t,t^2,t^3): t=1,2,3,4}.

The five monomials 1,x,y,z,y^2 all lie in V_4^3. Their evaluation columns
satisfy

    -24*1 + 50*x - 35*y + 10*z - y^2 = 0.       (BQ3)

Every four of these columns are independent: the first four are a
Vandermonde basis, and each coefficient in the displayed circuit is
nonzero modulo 101. Thus there are five distinct nonzero residue lines
of rank four. Choosing one from each line is not a basis. The checker
verifies all five determinants exactly modulo 101.

This refutes an automatic extension of the parallel-class proof to all
point ideals. It does not refute every possible Minkowski decomposition.
At the level of all monomials, the existence of a residue basis in which
each monomial has zero or one nonzero coordinate is equivalent to
binomiality: the forward implication supplies monomial or two-term
relations reducing every monomial to its basis representative, which
span the kernel; the reverse implication is Proposition 1.1 above. A
general theorem therefore needs to handle genuine multi-term circuits,
not just rename the parallel classes. No equivalence is asserted from
testing only the finite ground set V_n^d.

#### Primary comparisons checked in this pass

- [Boffi--Logar, arXiv:1608.08004](https://arxiv.org/pdf/1608.08004),
  introduction and Section 5, Propositions 5.1--5.3 and adjacent proofs,
  give an explicit planar full-rank lattice-ideal description. General
  border bases form a larger collection than initial staircases. That
  planar classification must not be repackaged as a new application of
  (BQ1). The rest of the paper has not received a full proof audit here.
- [Morita, arXiv:0712.2892](https://arxiv.org/pdf/0712.2892), introduction,
  Theorem 4.2 and its proof, and Corollary 4.3 and adjacent proof,
  establish the known link between the normalization of a finite
  abelian-group orbit Hilbert scheme and its Groebner fan. The inspected
  statements are a geometric correspondence, not the asserted new
  unrestricted counting theorem. No new G-Hilbert-scheme result is claimed.
- [Nemesch--Ruzika--Thielen--Wittmann, 2025 survey](https://doi.org/10.1007/s10898-025-01512-6),
  parameter definitions and Section 4.4.1, was checked to distinguish
  single-parameter bounds from the required two independent normalized
  weights in three variables. The inspected multi-parameter graphic
  bound does not improve the present unrestricted lattice bound. This
  survey is a guide to original results, not a proof of a new count.
- [Eppstein, *Geometric lower bounds for parametric matroid optimization*](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf),
  Theorem 7, Corollary 2 and conclusions, concern one-parameter lower
  constructions. A quotient-algebra realization with a controlled
  colength cost is not supplied by that result. Multiplying two
  breakpoint counts does not by itself construct a realizable
  three-variable ideal.
- [Hassin--Tamir, 1989](https://www.tau.ac.il/~hassin/matroids_89.pdf),
  primary metadata and abstract only, uses a pair of objective data for
  convex optimization in a plane. Its terminology must not be mistaken
  for a theorem on two independent weight parameters in our problem.
  No algorithmic result is imported from an unread proof.

#### Decision and next proof target

The binomial reduction provides a positive control with surviving
multiples, and isolates multi-term circuits as the missing phenomenon.
It does not change the acceptance criterion. There is still no improved
uniform bound for all point ideals or all colength-n ideals, no stronger
infinite lower family, and no completed target manuscript.

The next upper-bound attempt should explicitly account for overlapping
multi-term circuits under commuting multiplication. The next lower-bound
attempt must realize such circuits at a colength cost low enough to
exceed n^2 log^2 n in three variables. An abstract matroid arrangement,
a binomial subclass count, or a large integer weight certificate alone
will not meet that test. Keep the goal active and retain (BQ1) as a
comparison tool; do not launch a separate binomial search framework or
paper merely because this bound is easy to state.

#### A realizable next model: truncating degree classes

There is a concrete way to move beyond the binomial case without losing
control of realizability. For distinct positive numbers t_1,...,t_n and
positive integer exponents A_1,...,A_d, take

    P_A = {(t_j^A_1,...,t_j^A_d):1<=j<=n}.

Equal integer degrees A dot a give equal evaluation columns. Columns with
r distinct degrees have rank min(n,r).
To prove this, a nonzero polynomial with r nonzero terms has at most r-1
distinct positive zeros: divide by its smallest monomial and induct on
r using Rolle's theorem on the derivative. A dependence among at most n
distinct-degree columns would have n positive zeros and contradict this
bound. More than n distinct columns contain n independent ones. Thus the
residue matroid is the rank-n truncation of the degree-class partition
matroid. Its dependencies include circuits between n+1 distinct classes,
so this is genuinely outside independent rank-one class selection.

For a positive weight, select the cheapest representative in each degree
class and then the n cheapest distinct classes. This is an exact greedy
description for every choice of positive sample parameters, not an
assumption about generic numerical coordinates. Taking A_1=1 guarantees
distinct points and a cyclic x_1 multiplication operator.

The bounded control A=(1,6,9), n=7, t_j=j has 37 initial staircases over Q,
where the generic three-variable count is 36. The script enumerates all
86 seven-box staircases and tests strict positive feasibility by exact
rational clipping. For a candidate D with seven distinct degrees, the
inequalities say that every unselected member of a represented class is
more expensive than its selected representative, and every unrepresented
class is more expensive than all selected representatives. These conditions
are necessary and sufficient for the truncated-class greedy rule. Every
feasible weight is independently checked by rational Gaussian elimination
on the actual point evaluation matrix. Exact points, staircases and
integer weight witnesses are retained in the report.

This demonstrates extra variation in a simple characteristic-zero point
family; it is not evidence of a stronger asymptotic exponent. The known
Onn--Sturmfels nongeneric phenomenon already precludes treating a finite
increase as a new principle. Before increasing search size, the next task
is a symbolic count or a uniform obstruction for this truncated-class
family. It offers a realizable test of multi-term dependencies, while the
unrestricted acceptance criterion remains unchanged.

### Fourth feasibility pass: fixed curves cannot supply the required separation

**Status:** the fixed-curve gate now has a proof, including a matching
infinite family. Its growth is below the existing unrestricted generic
lower bound, so it does not satisfy the goal. The results below use the
established lattice-polytope estimate; novelty is not asserted.

#### An anisotropic version of the established moment argument

For any n-element staircase D in a box with integer sides q_1,...,q_d,
write mu(D)=sum_{a in D} a. Columns parallel to coordinate i have lengths
ell<=q_i, so

    mu_i(D) = sum_columns ell*(ell-1)/2
            <= n*(q_i-1)/2.                         (MB1)

Take any subcollection of the initial staircases of a fixed colength-n
ideal. Each of their moment sums remains a vertex of their convex hull:
a positive weight selecting that staircase uniquely also selects its
moment uniquely among any smaller collection. This is the established
state-polytope correspondence, with no generic-point hypothesis. If the
subcollection is contained in the same box, (MB1) and the lattice-polytope
vertex estimate used in
[Babson--Onn--Thomas, proof of Corollary 2.5](https://arxiv.org/pdf/math/0207135)
give a bound

    O_d((product_i [n*(q_i-1)/2+1])^((d-1)/(d+1))).  (MB2)

The added ones make the dimension issue harmless. Add the integral unit
cube to the moment polytope. Its normal fan refines the old one, so it has
at least as many vertices; it is full-dimensional and has volume bounded
by the product displayed in (MB2). Apply the usual volume bound there.
Consequently the weaker, simpler expression

    O_d((n^d * product_i q_i)^((d-1)/(d+1)))

always suffices. A finite explicit cover by boxes allows summing this
bound over the boxes. This is an anisotropic reuse of a known proof, not
a new unrestricted vertex theorem.

There is also a useful layer estimate in three variables. If a
subcollection has a_j<D throughout, then its moment vertices occupy at
most floor(n*(D-1)/2)+1 integer layers perpendicular to coordinate j.
In each layer they remain vertices of a two-dimensional lattice
polytope in a square of side at most n^2. The two-dimensional version of
the same vertex theorem bounds that layer by O(n^(4/3)). Thus this
subcollection has

    O(n^(7/3)*D)                                  (MB3)

initial staircases. This can be stronger than using only the full
three-dimensional volume. In particular, if I contains a degree-D
polynomial whose coefficients of x^D,y^D,z^D are all nonzero, then for
a positive weight with a unique largest coordinate its leading monomial
is that coordinate's D-th power. Three groups therefore cover all
initial ideals, proving g(I)=O(n^(7/3)*D). The condition D=O(n^(1/3))
would give O(n^(8/3)) for this class of ideals.

The condition is not automatic, even for radical point ideals. For the
Cartesian axis configuration D_L already used above, a polynomial of
degree at most L vanishing on the points restricts to a univariate
polynomial with L+1 roots on each axis. These restrictions are zero, so
all its pure top-degree coefficients vanish. Translation preserves the
top homogeneous part, giving the same obstruction with all point
coordinates nonzero. Here n=1+3L and g(I)=1. A general proof must also
handle mixed leading monomials and forced staircases; it cannot assert
the small-D condition for every point set.

#### Every finite quotient on a fixed monomial curve

Let A=(A_1,...,A_d) be a vector of positive integers, with d>=2. Suppose
I has colength n and contains the primitive pure-power relations

    x_j^(A_i/gcd(A_i,A_j)) - x_i^(A_j/gcd(A_i,A_j))

for every i!=j. In particular this holds for every finite quotient of
the monomial-curve ring K[t^A_1,...,t^A_d], in every characteristic.
Nonzero scalar twists of the relations give the same leading-power
argument. The ideal I itself need not be binomial or radical.

For a positive weight generic on the finitely many ratio ties, let i be
the unique minimizer of w_i/A_i. The displayed relation has leading
monomial x_j^(A_i/gcd(A_i,A_j)) for every j!=i. Hence the standard
staircase has coordinate bounds

    q_i=n,
    q_j=min(n,A_i/gcd(A_i,A_j)) for j!=i.

Every monomial initial ideal has a weight avoiding these ratio ties:
its finite set of strict defining comparisons remains true on an open
set. There are only d choices of i. Apply (MB2) to each group of moment
vertices to obtain the explicit uniform estimate

    g(I) <= C_d * n^(d-1)
            * sum_i (product_{j!=i}
                     min(n,A_i/gcd(A_i,A_j)))^((d-1)/(d+1)).       (MC1)

All constants in (MC1) depend only on d. In particular,

    g(I)=O_A(n^(d-1))                               (MC2)

when the curve is fixed. No Vandermonde or positivity assumption is used
for this upper bound. It covers all its finite quotients, not just the
positive point configurations tested in the preceding pass.

The dependence on A must not be suppressed when A varies with n. With
every capped ratio as large as n, (MC1) returns the old unrestricted
exponent 2d(d-1)/(d+1). Thus it is a quantitative restriction on a proposed
construction, not an unrestricted improvement in disguise.

The checker verifies 97 existing binomial and positive-point witnesses
against these sector boxes and the stronger moment bound (MB1). It first
perturbs each integer weight to remove ratio ties, preserving every
previously strict comparison on V_n. For the point witnesses the perturbed
basis is checked independently by rational Gaussian elimination.

#### The exponent is attained on the affine rational normal curve

For each fixed d>=2, let I_n be the ideal of any n distinct positive real
samples on

    (t,t^2,...,t^d).

Then

    g(I_n)=Theta_d(n^(d-1)).                       (MC3)

Here is an explicit lower construction, valid for every sufficiently
large n. Write n=(d+1)m+e, 0<=e<=d, take m>=16d, and put
h=floor(m/(8d)). Independently choose integers r_i in [m,m+h] for
1<=i<d, and set r_0=n-sum_{i=1}^{d-1}r_i. Thus r_0-r_i lies in
[7m/8,m+d]. With x_0=1, define the staircase

    D_r = union_{i=0}^{d-1} {x_i*x_d^k:0<=k<r_i}.

It has n members and is closed under division, since r_0>=r_i.
Choose integer weights

    w_d=d,    w_i=d*(r_0-r_i)+i for 1<=i<d.

For 1<=i,j<d, write i+j=q*d+s, 0<=s<d. The toric relations

    x_i*x_j - x_s*x_d^q

strictly reduce the product on the left. Indeed its weight is at least
(7/4)*d*m+2, while the right side has weight at most
d*m+d^2+2d-1 (or just d when s=0). The chosen m ensures strictness.
These reductions leave exactly the d rays x_s*x_d^k. Their images are
t^(s+d*k), one for each nonnegative integer degree, so they are
independent in K[t] and give the unique normal representatives there.
This also proves the claimed reduction system, without an assumption
about a generic toric Groebner basis.

The threshold T=d*r_0-1/2 selects exactly D_r from these rays: the last
selected monomial on each ray has weight at most d*r_0-1, and the first
unselected monomial has weight at least d*r_0. The selected degrees are
distinct. The positive-root argument in the previous pass makes their
evaluation columns independent on the n positive samples, so D_r is
exactly the minimum-weight residue basis. Monomials outside the rays
have a strictly cheaper representative of the same degree and cannot
alter this conclusion.

Different allocations give different staircases, providing
(h+1)^(d-1)=Omega_d(n^(d-1)) witnesses. Apply (MC2) with A=(1,...,d)
for the matching upper bound. This proves (MC3) for all sufficiently
large n, including every remainder e, and for every choice of distinct
positive samples. For a rational, completely explicit family take t_j=j.

The checker verifies 726 allocation certificates in dimensions 2--6,
including nonzero remainders. Each record gives the allocation ranges,
exact first and last certificates, a digest of the full deterministically
generated list, quadratic-reduction margins and cut thresholds. Four
additional evaluation-rank checks are performed modulo the proved prime
1000003. A nonzero integer evaluation determinant modulo that prime
certifies a nonzero determinant over Q. This is a check of those rational
examples, not a uniform assertion about arbitrary finite fields.

This is a sharp infinite restricted-family control. It is below
C_d(n)=Theta_d(n^(d-1)(log n)^(d-1)), and supplies no new lower bound
for F_d(n) or G_d(n). In particular, the earlier finite excess 37 versus
36 cannot justify a claim of asymptotic amplification on a fixed curve.

#### Prior-art distinction: a class region can contain different bases

An accessible publisher copy of
[Katoh--Tokuyama, *K-levels of concave surfaces*](https://link.springer.com/content/pdf/10.1007/s00454-001-0086-z.pdf)
has now been inspected at Theorem 12, Corollary 13, Proposition 15 and
their proofs, the set-equivalence definition, and Section 4. Proposition
15 bounds connected halfspace-deformation classes preserving which sets
are intersected. Corollary 13 bounds a single concave-surface level.
Neither statement, by itself, counts all representatives selected inside
the intersected degree classes. The proof of their general Lovasz lemma
has not been fully audited here; its statement and role were inspected.

The existing seven-point curve control makes the distinction exact.
Weights (131,22,27) and (181,26,45) select the same degree set

    {0,6,9,12,15,18,21}

but different initial staircases. Along their straight weight segment,
every pairwise-order change on V_7 occurs at one of

    0, 1/26, 7/38, 13/50, 1/2, 29/34, 1.

Exact checks at all seven endpoints and the six open intervals show that
the selected degree classes remain the same and the gap to every
unselected class remains positive. Choosing the halfspace threshold in
that gap gives a continuous deformation. Thus even one connected class
region contains two different bases. The 37 bases use 34 distinct degree
sets, but the number 34 alone is not asserted to count connected regions.
An additional refinement argument is required to convert a colored-level
bound to an initial-ideal bound.

The 2025 primary work of
[Bae--Oliver--Papadopoulou, full version arXiv:2504.06960](https://arxiv.org/pdf/2504.06960),
Section 6.1, Lemma 33, Corollaries 34--35 and their adjacent arguments,
was inspected alongside the conference version. It distinguishes level
complexity and refined color diagrams, and still cites the older
single-level bound. Its precise convex-distance or position hypotheses
must be checked before importing its sharper Voronoi counts. No
unrestricted polynomial-quotient count is inferred from these results.

[Conca--De Negri--Rossi (2007)](https://msp.org/ant/2007/1-3/ant-v1-n3-p01-p.pdf),
abstract and introduction only, studies Cohen--Macaulay initial ideals
of the homogeneous projective rational normal curve and their relation
to contracted ideals. Its parameter and ideal differ from the varying
colength-n finite affine quotients above. That distinction does not
establish novelty for (MC3); the full relevant theorem comparison is
still pending, and the curve result is not being promoted to a paper.

#### Decision for the full goal

The previous turn was progress: it produced the binomial reduction,
actual cyclic quotient checks and a realizable nonbinomial control,
saved in commit 57318ab. This turn resolves the promised fixed-curve
asymptotic gate and identifies an exact gap in a possible literature
transfer. It does not prove the required unrestricted advance.

Do not extend the fixed-A enumeration. A stronger lower construction
must vary A with n while respecting (MC1), or use a different realizable
multi-term dependency mechanism. For the upper route, (MB2) identifies
what a relation-based box cover would need to save; arbitrary mixed
leading monomials do not supply the pure coordinate bounds used here.
Keep the original goal active, with the unrestricted theorem, further
consequence and target manuscript still outstanding.

### Fifth feasibility pass: how many variables can propagate a relation?

**Status:** a proposed lower mechanism with one surviving coordinate is
now sharply bounded. This is a construction-level result, not the
unrestricted theorem. It uses the established state-polytope estimate;
novelty is not asserted. The preceding turn was progress, with the
fixed-curve gate and exact certificates saved in commit af20074.

#### A gluing shortcut fails even for a simple mixed relation

The relation xy=0 suggests treating the (x,z) and (y,z) components
separately. Their initial ideals alone need not determine the answer.
For any field K and integer m>=1, take

    I=(xy,x+y,z^m)=(x+y,x^2,z^m).

The two elimination ideals are (x^2,z^m) in K[x,z] and (y^2,z^m)
in K[y,z]. Each has a unique initial ideal. But I has two initial
staircases: {1,x} times {1,z,...,z^(m-1)} when w_y>w_x, and the same
set with y replacing x when w_x>w_y. The displayed generators reduce
triangularly and leave exactly 2m standard monomials, proving the claim.
Weights (1,2,1) and (2,1,1) are exact witnesses. The connecting relation
x+y contains information absent from both component fans. This refutes
the proposed common-refinement shortcut; it is not a novelty claim.

#### A uniform relative-socle bound

Let D be any finite monomial order ideal of size N in d variables, and
let J_D be its monomial complement ideal. Fix a set U of s coordinates,
1<=s<=d, and define the relative socle support

    S_U={a in D : a+e_i is outside D for every i in U}.

Let R be a linear subspace supported on S_U, stable under multiplication
in K[x]/J_D by the remaining variables. Multiplication by every variable
in U is already zero on this support, so J_D+R is an ideal. The stability
hypothesis is essential: arbitrary relation rows cannot be substituted
for a polynomial quotient. Write n=N-dim R>=1.

Then

    g(J_D+R)=O_d(n^((2d-s)*(d-1)/(d+1))).          (RSoc1)

To prove this, put q=|S_U|. For each i in U, the downward i-chains from
the distinct members of S_U are disjoint. Two members on one such chain
would contradict maximality in coordinate i. Therefore

    sum_{a in S_U}(a_i+1)<=N,
    sum_{a in S_U}a_i<=N-q<=n.                    (RSoc2)

Every initial staircase contains the same fixed set D minus S_U.
The other standard monomials are selected from S_U by the residue
matroid of span(S_U)/R. Subtract the moment of the fixed set. In the s
specified coordinates the remaining moment is at most n by (RSoc2).
In every other coordinate it is at most the total moment of an
n-element staircase, hence at most n(n-1)/2. The variable moment
polytope thus lies in

    [0,n]^s times [0,n^2]^(d-s).

Adding the integral unit cube handles lower-dimensional cases. Its
volume is O_d(n^(2d-s)); the lattice-polytope vertex bound in the
[Babson--Onn--Thomas proof of Corollary 2.5](https://arxiv.org/pdf/math/0207135)
gives (RSoc1). The constant depends only on d, not on N, R, the field
or the chosen monomial order ideal. This recovers the full-socle result
at s=d and permits genuinely nonzero multiples in the other coordinates.

For d=3 the consequences are more informative when the layer argument
is included:

| Coordinates killing the relation support | Uniform upper bound |
| --- | --- |
| Three | O(n^(3/2)) |
| Two | O(n^2) |
| One | O(n^(7/3)) |

The last row improves the volume-only O(n^(5/2)): the bounded moment
coordinate has at most n+1 integer values. Each slice is a planar
lattice polytope in a square of side n^2, hence has O(n^(4/3)) vertices.
The other two rows use (RSoc1). These statements concern this specified
relation-support construction. They do not classify arbitrary ideals
containing a mixed polynomial.

#### The one-surviving-coordinate exponent is attained

The quadratic upper bound in the second row is sharp along an explicit
infinite family. Let L>=2, let F_L(z)=product_{j=1}^L(z-j) over Q, and set

    J=(x^4,x^2*y,x*y^2,y^4,F_L(z)),
    I_L=J+(x*y-z^L*y^3, x^3-z^(2L)*y^3).         (RSoc3)

The two-variable monomial reference diagram is

    E={1,x,x^2,x^3,y,y^2,y^3,xy},

with maximal monomials y^3,xy,x^3. Their multiples by x and y vanish
in the reference ring, but their z multiples need not vanish. Modulo J,
this support is a free rank-three module over B=Q[z]/F_L. Its quotient
in (RSoc3) is B, under

    y^3 |-> 1,   xy |-> z^L,   x^3 |-> z^(2L).

The kernel is a B-submodule of dimension 2L, so I_L has length 8L-2L=6L.
The other five monomials of E, each times z^k for 0<=k<L, remain fixed
independent standard columns. The reference J has the same diagram
E times [0,L) for every monomial order: its generators are monomials
in x,y and one monic univariate polynomial in z. Normal-form reduction
preserves the relative support. Thus the proof of (RSoc1) applies
verbatim to this product reference, giving g(I_L)=O(L^2).

For the lower bound, evaluate B at z=1,...,L. The 3L support monomials
have columns t^k, t^(L+k), t^(2L+k), 0<=k<L. Every L of them are
independent by the positive-root/Vandermonde argument already proved.
Therefore the varying part of the standard basis is exactly the L
cheapest support monomials, from three arithmetic progressions.

Here are integer weights for quadratically many different choices.
Write L=4m+e, 0<=e<=3, take m>=8, and put h=floor(m/8). Independently
choose a,c in [m,m+h], and set b=L-a-c, u=b-a, v=b-c. Both u and v
are positive. In the progression order y^3,xy,x^3 choose counts (a,b,c)
and weights

    w_x=3u+6v+5,  w_y=6u+3v+4,  w_z=9.

Relative to the middle offset w_x+w_y, the outer offsets are 9u+3
and 9v+6. The threshold

    T=w_x+w_y+9b-1/2

selects exactly a,b,c monomials from the three progressions. Each last
selected score is at most w_x+w_y+9b-3, while each first unselected
score is at least w_x+w_y+9b. Different allocations give distinct
standard staircases. There are (h+1)^2=Omega(L^2) choices, proving

    g(I_L)=Theta(L^2)=Theta(length(I_L)^2).       (RSoc4)

The reference now contains F_L(z) rather than a monomial z^L, which
allows the positive evaluation proof. This extension is justified by
its order-independent reference diagram and support-preserving normal
forms; it is not silently assumed for arbitrary reference ideals.

#### Exact checks and field scope

The checker constructs full commuting multiplication matrices and a
cyclic vector for L=2,3,4,5. Companion-polynomial remainders implement
(RSoc3) independently of its positive evaluation coordinates. In every
positive-weight comparison chamber, descending relation-row reduction
is checked against ascending Gaussian selection over Q, F_2, F_3 and
F_29. Rational selection is also checked against the separate positive
evaluation matrix and the L-cheapest-column rule. Full matrix greedy
selection is checked once for every realized state in every field.

The exact numbers of initial staircases for these fixtures are:

| L | Q | F_2 | F_3 | F_29 |
| --- | --- | --- | --- | --- |
| 2 | 5 | 3 | 3 | 5 |
| 3 | 8 | 4 | 4 | 8 |
| 4 | 12 | 3 | 8 | 10 |
| 5 | 14 | 5 | 6 | 14 |

The relation module has dimension 2L and the quotient length is 6L in
all these characteristics. The characteristic-zero lower proof uses
positive distinct samples and is not promoted to a uniform finite-field
lower assertion. The checks cover 896 row/greedy comparisons, 112 full
matrix states, and 116 explicit lower-allocation certificates at larger
L. Saved witnesses include the actual reference diagrams, F_L
coefficients, relation matrices and positive integer weights.

#### Decision and remaining mathematical gap

Adding z levels to these socle-supported constructions can produce a
quadratic number of states, but cannot exceed the existing generic
n^2 log^2 n lower order. This rules out that proposed amplification
mechanism with a proof rather than a larger search.

Within this support-based approach, the remaining candidate has only
one annihilating coordinate and a relation module stable under two
surviving coordinate actions. Its proved O(n^(7/3)) ceiling leaves room
for a power improvement over the current generic lower bound; no family
achieving that improvement has been constructed. Any proposed module
must satisfy both multiplication actions, not just specify unrelated
matroid rows. Do not expand the two-annihilator enumeration.

The unrestricted upper problem also remains open. The gluing example
above shows why component fans alone do not settle mixed-leading cases.
Further primary-source access attempts for the NSF-hosted later
finite-field chapter did not retrieve its full text; its improved
finite-field theorem is still not represented as proof-audited. This
access issue does not halt independent work and is not evidence that
the unrestricted bound is new. Keep the full goal active; the
qualifying unrestricted advance, further consequence and target
manuscript are still missing.

### Sixth feasibility pass: bounded connectors also limit sparse plane deformations

The next test uses genuine first-order deformations, rather than choosing
unrelated relation rows. Both surviving multiplication actions can be
nonzero. Nevertheless, the entire hook family below has only O(n^2)
initial ideals. This is another restricted obstruction, not the
unrestricted advance required by the goal.

#### A counting lemma that keeps the connecting information

Let a finite matroid have ground set

    E = E_1 disjoint-union ... disjoint-union E_p disjoint-union C,
    |C|=k,

and suppose deleting C leaves the direct sum of the restrictions to
the E_i. Element costs are linear functions of the weight parameters.
Suppose the joint arrangement of comparisons *within* the E_i has R
open regions. Then the number of distinct minimum bases realized by
generic weights is at most

    R product_{j=1}^k(p+j) = R (p+k)!/p!.       (BC1)

This uses the full internal element orders, rather than just the
component minimum-basis fans. The gluing counterexample in the fifth
pass therefore does not contradict it.

Proof. Fix an internal-order region. The minimum basis after deleting
C is fixed, by the greedy algorithm in each summand. Insert the k
elements of C in a fixed sequence. Given the previous minimum basis B
and a new element c, the usual matroid insertion rule either appends c
if the rank increases, or removes a maximum-cost member of the unique
circuit in B union {c}. Removing c itself means retaining B.

For completeness, this rule follows from the greedy prefix-rank
criterion. At any cost threshold, the members of a minimum basis below
the threshold span all old elements below that threshold. If c is in
the span of B, it enters this span exactly when all old members of its
fundamental circuit have entered. Thus replacing the most expensive
old circuit member by c gives precisely the required rank increase
between their two costs when c is cheaper; otherwise B stays optimal.
If c is outside the old span it must be appended.

At insertion j, condition on a particular previous basis. The most
expensive circuit member from each E_i is fixed by that block's
internal order. Together with the previous j-1 inserted elements and
c, there are at most p+j possible elements to remove. The rank-increase
case has a single outcome. Induct on j and sum over the R regions.
Generic perturbations avoid ties; every uniquely exposed basis is
represented in an open region. This is an elementary consequence of
the existing greedy/circuit rule, with no novelty claim for that rule.

#### The full hook deformation class and its exact length

Take S=K[y,z], arbitrary field K, and

    M=(y^r,yz,z^s),  r,s>=2,
    D={1,y,...,y^(r-1),z,...,z^(s-1)},
    B=S/M,  m=|D|=r+s-1.

An S-linear map phi:M -> B is determined by the three tails
f_0=phi(y^r), f_1=phi(yz), f_2=phi(z^s). Its compatibility conditions
are exactly

    z f_0 = y^(r-1) f_1,    z^(s-1) f_1 = y f_2  in B.       (HD1)

Indeed, the two adjacent least-common-multiple syzygies generate the
syzygy module of these three plane monomial generators: the remaining
pair syzygy is their telescoped combination. Comparing the independent
y and z chains in B shows that (HD1) is equivalent to

    f_0 in span{y,...,y^(r-1),z^(s-1)},
    f_1 in (y,z)B,
    f_2 in span{y^(r-1),z,...,z^(s-1)}.          (HD2)

These are 2m freely chosen scalar parameters, over every field. Define
the actual three-variable ideal

    I_phi=(x^2, y^r-x f_0, yz-x f_1, z^s-x f_2).             (HD3)

It has basis D disjoint-union xD and length n=2m. Here is a direct
proof, so flatness is not inferred from an unchecked finite deformation.
Multiplying the three relations by x gives xM contained in I_phi.
Then M^2 is contained in I_phi too: for u,v in M, replace u by
x phi(u), and its product with v is zero modulo xM. Consequently D
and xD span the quotient.

For independence, suppose a+xb is in I_phi with a,b in span(D).
Reducing modulo x gives a in M, hence a=0. Write a representation using
coefficients A_i+xB_i modulo x^2 and denote the three generators of M
by g_i. The constant coefficient gives sum A_i g_i=0. The x
coefficient, reduced modulo M, then gives

    b = -sum A_i f_i = -phi(sum A_i g_i) = 0 in B.

Thus b=0 as well. The quotient is free of rank m over K[x]/(x^2).
In this note, 'first-order deformation' means this algebra over the
dual numbers; it does not mean substituting an ordinary nonnilpotent
parameter into an arbitrary perturbation.

If a constant 1 is added to f_1 while the other tails are held fixed,
both compatibilities fail. The resulting relations force
x y^(r-1)=x z^(s-1)=0, so the claimed free basis would be false.
The checker includes this negative control.

#### The three connecting columns

The monomial reference ideal

    J=(x^2,xM,M^2)

is contained in I_phi. Let D_0 be the standard diagram of M^2;
|D_0|=3m+1. Its reference domain is D_0 disjoint-union xD.
All relations between its residue columns are supported on

    T=(D_0 minus D) disjoint-union xD.

The plain D columns are independent fixed summands. The variable
residues in T all lie in xB, with relation rows a-x phi(a) for
a in D_0 minus D. These rows are closed under multiplication by both
y and z and killed by x in the reference quotient.

The column x*1 is a coloop, since none of the tails has a constant
term. Remove it, and reserve the three connecting columns

    C={y^r,yz,z^s}.

The other columns split into

    E_Y = {xy^i:1<=i<r}
          union {y^i:r+1<=i<2r}
          union {y^i z:2<=i<=r},

    E_Z = {xz^j:1<=j<s}
          union {z^j:s+1<=j<2s}
          union {y z^j:2<=j<=s}.

By (HD2), E_Y has residues only in span{xy,...,xy^(r-1)},
and E_Z only in span{xz,...,xz^(s-1)}. Their identity columns span
these disjoint spaces, in every characteristic, regardless of which
tail coefficients vanish. Deleting C therefore leaves a direct sum.
There are p=2 blocks and k=3 connectors in (BC1).

Each block consists of three exponent progressions of length r-1 or
s-1, with respective offsets w_x, 0 and w_z (or w_y for E_Z).
Within E_Y, only comparisons between progressions can change. Each
of its three pairs has at most 2r-3 distinct comparison hyperplanes;
the same count with s holds for E_Z. In particular, their joint
internal arrangement has at most H=6(r+s) hyperplanes. Normalizing
w_x+w_y+w_z=1 gives an open triangle, so it has at most

    R <= 1+H(H+1)/2

regions. Every global monomial initial ideal is represented by
positive weights on the finite reference comparisons, as in the
earlier state-polytope reduction. Applying (BC1), and restoring the
fixed columns, proves uniformly

    g(I_phi) <= 60[1+H(H+1)/2] = O(n^2),  n=2(r+s-1).        (HD4)

No genericity of coefficients or characteristic-zero assumption is
needed for this upper bound. Although balanced hooks have reference
bounding boxes of volume of order n^2, escaping the earlier box-volume
obstruction, their bounded number of connecting columns still prevents
an asymptotic improvement over the generic n^2 log^2 n lower order.

#### Reproduction and relation to existing deformation theory

The exact checker uses the explicit integer tails

    f_0=sum_{i=1}^{r-1}(i+1)y^i+(r+3)z^(s-1),
    f_1=sum_{i=1}^{r-1}(i+2)^2 y^i
        +sum_{j=1}^{s-1}(j+3)^3 z^j,
    f_2=(s+5)y^(r-1)+sum_{j=1}^{s-1}(j+5)^2 z^j.

It checks all positive comparison chambers for (r,s)=(2,2),(2,3),(3,3)
over Q, F_2, F_3 and F_29. Descending relation-row reduction is
compared with ascending Gaussian residue selection. Full commuting
multiplication matrices and cyclicity on the reference domain and its
border independently check that the rows come from the claimed
quotient. Full-algebra selection is checked for every realized state;
the checker also verifies the two block subspaces, the internal
arrangement count and the bound of 60 states per internal-order region.
These finite checks support (HD4); the proof above supplies its uniform
scope. The (3,3) fixture has both surviving actions nonzero over Q.

First-order and formal deformations of plane monomial ideals are
established mathematics. [Oszer, 2025](https://doi.org/10.1093/qmath/haaf032)
provides explicit Hilbert--Burch families: the inspection covered
Definition 3.1, Theorems 5.4--5.5 and their proofs, and Examples 6.1--6.2.
The formal parametrization does not justify arbitrary finite
perturbations; the examples explain this limitation. We use only the
dual-number family proved directly above. No new deformation theory or
published priority for (HD4) is asserted; the full paper was not audited.

The comparison with [Eppstein, 1998](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf)
was extended to Section 4, Theorems 4--5 and 7, Corollaries 1--2 and
their construction arguments. Those large one-parameter examples use
transversal matroids built from convex chains of line arrangements.
They do not themselves realize a commutative quotient algebra with
the required exponent weights and controlled colength. A realization
and its size cost remain necessary before importing that lower bound.

#### Decision and next gate

Stop enlarging hook parameters or changing their coefficients in an
attempt to beat the generic lower order: (HD4) rules out the entire
family. Also reject the initially considered dense three-plateau
module: its reference box has volume proportional to its quotient
length, so the previously proved moment-box estimate already gives
O(n^2), without needing another enumeration.

Within the first-order approach, the next meaningful test is a plane
monomial ideal whose number of corners grows with its colength. First
determine how its syzygies limit the possible connecting residue
columns, and whether a deletion/direct-sum argument still limits the
count. Construct a symbolic compatible family before increasing any
enumeration. This is a feasibility direction, not an asserted lower
bound. The unrestricted upper-bound and primary-priority work remain
necessary alternatives. Keep the full goal active: no qualifying
unrestricted asymptotic advance, further consequence or new target
manuscript has been obtained.

### Seventh feasibility pass: interface rank, rather than corner count

The previous turn made mathematical progress: commit e611ebf proved the
uniform hook obstruction and saved independently checked quotient models.
This pass tests the proposed escape through a growing number of corners.
A rank-sensitive version of the deletion argument rules out a sparse
family with arbitrarily many corners as well. The unrestricted target
remains unchanged and unproved.

#### A linear-matroid interface lemma

Partition the columns of a finite vector configuration into E_1,...,E_p.
Write V_i=span(E_i), and set

    W_i = V_i intersect sum_{j!=i} V_j,
    q = sum_i dim(W_i).

Suppose their joint internal element-order arrangement under linear
costs has R open regions. Count the unique minimum bases realized at
weights with distinct element costs, as holds for generic weights on
the monomial columns here. Their number is at most

    2^q R.                                                   (IR1)

The intersections are with the sum of *all* other spaces. Pairwise
intersection dimensions alone would not suffice: three distinct lines
in a plane give a counterexample to that replacement.

Proof. Fix an internal-order region. Let B_i be its fixed greedy basis
of E_i. A column outside B_i is dependent on cheaper columns in its own
block and therefore can never be selected by global greedy selection.
This follows because, after every prefix of a greedy scan, the selected
columns span all columns scanned so far, including the rejected ones.

Scan B_i in its fixed order and project its prefixes to V_i/W_i.
Mark a column private when it increases the projected rank. Every marked
column is forced into the global greedy basis: dependence on previously
scanned global columns would, modulo the span of the other blocks,
contradict that projected rank increase. There are exactly
dim(V_i)-dim(W_i) such private columns. The remaining dim(W_i) local
basis columns are the only optional columns from block i. Thus all
variation occurs among at most q fixed columns in this internal region,
and there are at most 2^q possible subsets. This proves (IR1).

This proof uses only linear algebra and the existing matroid greedy
criterion. It does not establish priority for the lemma. Unlike (BC1),
it permits an unbounded number of columns in the coupling, provided
the interaction spaces have bounded dimension.

There is also a constructive consequence. In each internal region,
compute the fixed local bases, private columns and optional columns
using rank calculations. Refine the region only by pair comparisons
between the optional columns. There are at most binom(q,2) additional
hyperplanes; in three variables the normalized region is a convex
polygon, so at most

    1+binom(q,2)(binom(q,2)+1)/2

open cells are needed. On every such cell the global greedy result is
fixed: all private columns are forced and only the optional-column
order can change the basis after contracting that fixed set. Choosing
one positive rational weight per cell gives a per-instance order list
of size O((q+1)^4 R) recovering all initial ideals. This is a statement
about the stated representation and exact rank/cell operations, not a
new general bit-complexity bound. No coefficient-independent order list
for the entire class is inferred.

#### A sparse staircase with arbitrarily many corners

For k>=2 and r,s>=2k+1, take

    M=(y^r, y^k z, y^(k-1) z^2, ..., y z^k, z^s) in S=K[y,z].

Its diagram D consists of the two long pure-power arms and all mixed
monomials of total degree at most k. It has

    m=|D|=r+s-1+k(k-1)/2,

and k+2 minimal generators. Write B=S/M and define

    A_Y=span{y^i:k<=i<r},  A_Z=span{z^j:k<=j<s},
    W=Soc(B)=span{y^(r-1),z^(s-1),y^i z^(k-i):1<=i<k}.

Label the generators g_0=y^r, g_i=y^(k+1-i)z^i for 1<=i<=k,
and g_(k+1)=z^s. Every S-linear map phi:M -> B has exactly the tails

    f_0,f_1 in A_Y+W,
    f_i in W for 2<=i<=k-1,
    f_k,f_(k+1) in A_Z+W.                                    (SC1)

These tail choices are independent. Their total dimension is
2r+2s+(k-2)(k+1)=2m, including the overlaps of A_Y or A_Z with W
when taking the displayed spans.

Here is a direct classification proof. The adjacent syzygies are

    z f_0=y^(r-k)f_1,
    z f_i=y f_(i+1) for 1<=i<k,
    z^(s-k)f_k=y f_(k+1).                                   (SC2)

As before, they generate all syzygies. Since r-k>=k+1, the right
side of the first equation has only pure y terms, whereas the left
side has a positive z exponent. Both must vanish separately. Thus
f_0 lies in Ann_B(z)=A_Y+W, and f_1 has no pure y term of degree
less than k, including the constant. The last equation similarly
gives f_(k+1) in A_Z+W and excludes pure z terms of degree less
than k from f_k.

For a homogeneous component of degree h<k, the middle equations are
equations in the ordinary polynomial ring: B contains every monomial
of degree at most k. They imply

    (f_i)_h = y^(k-i) z^(i-1) H,
    deg(H)=h-(k-1).

This is zero for h<k-1. For h=k-1 the first endpoint restriction
forces the scalar H to be zero. In degree at least k, the only
nonsocle monomials in B are pure arm powers. The middle equations
force all nonsocle z powers out of f_1,...,f_(k-1) and all nonsocle
y powers out of f_2,...,f_k. This proves (SC1). Conversely, every
tail in (SC1) makes both sides of each equation in (SC2) zero, so
all those choices are valid over every field.

For every such phi, the ideal

    I_phi=(x^2, g_0-x f_0, ..., g_(k+1)-x f_(k+1))             (SC3)

has length n=2m and free basis D over K[x]/(x^2), by the elementary
flatness proof in the sixth pass. The reference J=(x^2,xM,M^2) is
contained in I_phi. Again its plain D columns are independent fixed
summands, and all variable relations have residues in xB.

#### Three blocks with bounded interface dimension

Remove the variable coloops x times (D minus (A_Y union A_Z union W)),
where the spans stand for their monomial bases, and remove zero columns.
Partition the remaining columns into the following blocks; zero columns
can also remain in the displayed blocks without affecting the proof:

    E_Y = {xy^i:k<=i<r-1}
          union {y^i:r+1<=i<2r}
          union {y^i z:k+1<=i<r+k},

    E_Z = {xz^j:k<=j<s-1}
          union {z^j:s+1<=j<2s}
          union {y z^j:k+1<=j<s+k},

    E_C = {g_0,...,g_(k+1)} union x Soc(B).

Every unassigned variable column is zero. Indeed, every positive
multiple of an intermediate f_i is zero, and the only possibly
nonzero positive multiples of f_0,f_1 are their pure y multiples;
the other endpoint tails behave symmetrically. This also proves that
the displayed arm columns have residues in xA_Y and xA_Z respectively.

The two arm residue spaces are disjoint. The core residue space has
y-arm components contained in

    span{xy^(r-1), x(f_0)_Y, x(f_1)_Y},

of dimension at most three, and z-arm components in the analogous
three-dimensional space from f_k,f_(k+1). All its other components
lie in the mixed socle. Hence the interface dimensions in (IR1) are
at most 3,3,6 for E_Y,E_Z,E_C respectively. In particular q<=12,
independently of k, r, s, coefficients and characteristic.

It remains to bound internal orders, not just component basis fans.
Each arm block consists of three parallel exponent progressions, so
at most 6r or 6s distinct internal comparison hyperplanes suffice.
In the core, the k mixed generators lie on the line y+z=k+1 and
the k-1 mixed socle identity columns lie on the parallel line
x=1,y+z=k. Comparisons within these progressions contribute one
hyperplane, and between them at most 2k-2. The four remaining core
columns are y^r,z^s,xy^(r-1),xz^(s-1); comparing them with the
progressions and with each other adds at most 8k+2. Thus the joint
arrangement has at most

    H=10(r+s+k+1)=O(m)

hyperplanes, and at most R<=1+H(H+1)/2 positive normalized regions.
Applying (IR1) proves uniformly

    g(I_phi) <= 2^12 [1+H(H+1)/2] = O(n^2).                  (SC4)

The constructive consequence supplies O(n^2) positive rational orders
for each specified quotient, with a constant independent of the number
of corners. These are restricted consequences of the proved interface
lemma, not the requested unrestricted advance or an established
novelty claim.

For example, r=s=k^3 gives length of order k^3 and an unavoidable
staircase bounding box of volume of order k^6, since the fixed plain
D already contains both long arms and x is a variable coloop. Thus
the uniform conclusion is not just the dense-box obstruction: the
number of corners grows while the box-to-length ratio is unbounded.

#### Exact controls and next decision

The checker uses k=2,3,4 with r=s=2k+1, giving m=10,16,23 and
4,5,6 minimal generators, and the asymmetric control (k,r,s)=(3,11,13)
with m=26. It verifies the full adjacent-syzygy matrix
has tangent dimension 2m over Q,F_2,F_3,F_29; the displayed free
tail coordinates satisfy all equations. Explicit integer tails have
coefficient (i+2)(a_y+1)^2+(a_z+2)^3 at each allowed monomial a
of f_i. Independent relation rows, commuting multiplication matrices,
cyclicity including the reference border, and full-algebra greedy
selection check the quotient model. Selected positive integer weights
also check every private/optional column assertion in (IR1), and recover
the same global basis after contracting the forced set and scanning
only the optional columns. A separate three-line configuration verifies
the failure of the pairwise-intersection replacement over all four fields.

These are bounded exact weight controls, not exhaustive fan counts.
The report saves the actual tails, relation matrices, blocks, fields,
weights and private/optional witnesses. The symbolic proof supplies
the uniform assertion; no asymptotic count is fitted to the samples.

The later Stigler--Zhang chapter remains accessible as a primary NSF
record and an indexed Theorem 8, but its full PDF request again failed.
Its theorem has not been promoted to a full proof audit. The publisher's
[chapter record](https://doi.org/10.1007/978-3-030-42687-3_9) confirms
the finite-field scope. This does not settle the strongest-current-bound
audit, and failed retrieval is not evidence of novelty.

Increasing the corner count alone is therefore not enough: (SC4)
covers an unbounded-corner sparse class. A lower construction through
this route must have unbounded interaction dimension or substantially
more internal comparison hyperplanes, in addition to two surviving
coordinate actions. The next feasibility test should target that
specific missing property, rather than enumerating larger instances
of (SC3). The unrestricted upper problem and its priority audit remain
active alternatives. The full goal remains incomplete.

#### A concrete next test: thick hooks, with a substitution control

There is an explicit compatible family whose natural arm interface
can grow. For a,b>=1, R>=2a+1, S>=2b+1, take

    M=(y^R,y^a z^b,z^S),
    m=Rb+Sa-ab,  B=K[y,z]/M.

The two adjacent syzygies separate into disjoint arm supports, giving
the full tail classification

    f_0 in Ann_B(z^b),
    f_1 in (y^a,z^b)B,
    f_2 in Ann_B(y^a).                                     (TH1)

Indeed z^b f_0=y^(R-a)f_1 has disjoint z>=b and y>=R-a>=a
supports, so both vanish; the other endpoint behaves symmetrically.
The two annihilators constraining f_1 are both the complement of the
bottom a-by-b rectangle. The three dimensions are Rb, m-ab, Sa,
whose sum is 2m. Every choice in (TH1) is valid and gives the flat
length-2m quotient (x^2,y^R-xf_0,y^a z^b-xf_1,z^S-xf_2).

For example f_0=y^a+z^(S-b), f_1=0, f_2=z^b is compatible.
In a natural arm/core partition placing the top z-arm identity columns
in the core, those columns span an a*b dimensional subspace also
spanned by pure multiples of the z^S relation: for 0<=i<a,0<=j<b,

    phi(y^i z^(2S-2b+j))=y^i z^(S-b+j).

These monomials belong to the reference diagram of M^2. Thus this
partition's interaction dimension is unbounded when a*b grows. This
does not prove a large fan or preclude a different useful partition.

In fact a necessary negative control rules out a tempting inflation.
If an ideal I in K[x,Y,Z] is merely extended by Y=y^a,Z=z^b,
its quotient length is multiplied by a*b but its number of initial
ideals is unchanged. The larger ring is free with basis y^i z^j,
0<=i<a,0<=j<b. Its monomial residue relations split into these
classes. Each class has the same greedy problem at weights
(w_x,a*w_y,b*w_z), with a fixed additive cost offset. Consequently
the initial ideal is the corresponding free extension of the old
initial ideal, and positive weight scaling is bijective. This proves
the equality of fan counts without a novelty claim.

The displayed interface witness is itself such an inflation when
a divides R and b divides S, so it must not be presented as a lower
amplification. The meaningful next test is the full (TH1) family with
tails mixing the residue classes, together with a symbolic count or
a new uniform obstruction. Growing interface rank alone is insufficient.

### Eighth feasibility pass: a thickness threshold and a convolution constraint

The preceding turn made progress in commit 2f92871: the interface lemma
proved a uniform bound for an unbounded-corner sparse class and supplied
the compatible thick-hook candidate. This pass treats the full (TH1)
tail class, including coefficients that mix all residue classes. It
proves a parameter-sensitive upper bound and identifies an exact
commutativity obstruction to importing a support-matching lower model.
Neither conclusion is the requested unrestricted advance.

#### Use the interface lemma polynomially when the interface grows

The constructive part of (IR1) already implies, in three variables,

    g <= A(q) L,
    A(q)=1+binom(q,2)(binom(q,2)+1)/2=O((q+1)^4),             (IR2)

where L is the number of joint internal-order regions. Indeed, in each
region only q optional columns remain after contracting the private
columns. Their at most binom(q,2) comparison lines cut the normalized
convex weight region into at most A(q) cells. This also bounds the
number of bases, irrespective of the larger bound 2^q. No new
parametric-matroid theorem is imported for this elementary observation.

For M=(y^R,y^a z^b,z^S) as in (TH1), split its diagram into the
bottom rectangle D_* and the disjoint outer arms

    U={y^i z^j:a<=i<R,0<=j<b},
    V={y^i z^j:0<=i<a,b<=j<S}.

Let U_top be the a-by-b rectangle at the right end of U, and V_top
the a-by-b rectangle at the top of V. The full tail classification is

    f_0 in span(U union V_top),
    f_1 in span(U union V),
    f_2 in span(U_top union V).

The xD_* columns are fixed coloops in the variable residue matroid.
After removing them, define three blocks:

    E_Y = x(U minus U_top)
          union {y^(R+i)z^j:a<=i<R,0<=j<b}
          union {y^(a+i)z^(b+j):a<=i<R,0<=j<b},

    E_Z = x(V minus V_top)
          union {y^(a+i)z^(b+j):0<=i<a,b<=j<S}
          union {y^i z^(S+j):0<=i<a,b<=j<S},

    E_C = x(U_top union V_top)
          union {g_l y^i z^j:l=0,1,2,0<=i<a,0<=j<b}.

Here g_0=y^R,g_1=y^a z^b,g_2=z^S. These disjoint column sets lie
in the reference diagram of J=(x^2,xM,M^2); E_C has 5ab columns.
Any remaining variable column is zero. For example, multiplying the
f_1 tail by y^i with i>=a kills its V component, whereas multiplying
by z^j with j>=b kills its U component. Both together kill the tail.
For f_0 and f_2 the endpoint annihilators give the corresponding
one-arm statement. This proves the claimed partition without assuming
independent coefficient rows.

The arm residue spaces lie in x span(U) and x span(V) and are
disjoint. The core's U projection is contained in

    x span(U_top)
    + span{x y^i z^j (f_0)_U, x y^i z^j (f_1)_U:0<=i<a,0<=j<b},

of dimension at most 3ab. Its V projection has the same bound, using
f_1 and f_2. Therefore the three interface dimensions are at most
3ab,3ab,6ab, and

    q<=12ab.                                                (TH2)

The internal comparison arrangement remains linear in m even when
the interface grows. In E_Y the x differences have three possible
values, y differences have at most 4R+1 values and z differences at
most 4b-1. Thus at most 60Rb distinct hyperplanes suffice. Similarly
E_Z uses at most 60Sa. The core is the union of five a-by-b exponent
rectangles. Each ordered pair of rectangles supplies at most
(2a-1)(2b-1) difference vectors, with fixed x difference, hence at
most 100ab altogether. Consequently

    H<=60Rb+60Sa+100ab=60m+160ab<=120m,
    L<=1+H(H+1)/2,

because m=Rb+Sa-ab>=3ab. Applying (IR2), with n=2m, proves

    g(I_phi)=O(n^2 (ab)^4).                                  (TH3)

The independent one-annihilator moment-layer estimate also applies:
all variable relations in this reference are killed by x and closed
under y,z. Combining the two bounds gives, over every field and for
every compatible set of tails,

    g(I_phi)=O(min{n^(7/3), n^2 (ab)^4}).                     (TH4)

The constants are absolute for this three-variable class. In
particular, fixed thickness cannot beat the existing generic order,
even with arbitrary mixed-residue coefficients. More generally,
ab=O(sqrt(log n)) cannot give an asymptotic order larger than
n^2(log n)^2. A proposed lower bound Omega(n^(2+epsilon)) would require
ab=Omega(n^(epsilon/4)); testing the ceiling n^(7/3) therefore requires
at least ab=Omega(n^(1/12)). These are necessary conditions, not
existence assertions.

#### A perfect support matching can still give a forced dependence

One possible lower strategy would replace the convolution columns by
independent entries with the same zero pattern, then import a hard
transversal-matroid example. Commutativity obstructs that step even
with generic polynomial coefficients.

Take a=2,b=3,R=5,S=7, and let

    F=1+2y+3z,  G=5+7y+11z,
    f_0=y^2 F,  f_1=y^2 G,  f_2=0.

These satisfy (TH1), giving an actual length-46 quotient. The six
plain monomials

    g_0, y g_0, z g_0, g_1, y g_1, z g_1

have xB residues x y^2 times F,yF,zF,G,yG,zG. All their products
fit inside the arm, so truncation removes no coefficient. With rows
1,y,z,y^2,yz,z^2 after removing the common x y^2 factor, the matrix is

    [1 0 0  5  0  0]
    [2 1 0  7  5  0]
    [3 0 1 11  0  5]
    [0 2 0  0  7  0]
    [0 3 2  0 11  7]
    [0 0 3  0  0 11].

Its support has a matching of size six, but its actual rank over Q
is five. The exact relation is

    5F+7yF+11zF-G-2yG-3zG=GF-FG=0.                         (TH5)

The saved control also gives a nonzero 5-by-5 minor of determinant 48
in its explicitly recorded row/column order. Thus the rank is exactly
five. More generally GF-FG=0 holds for symbolic coefficients too;
the displayed nonzero minor shows that rank five is generic for this
matrix pattern. The support-matching model would instead have generic
rank six. This is a classical Koszul relation, used as a realizability
control, not a new syzygy theorem or a lower-bound result.

The relevant established matrix literature must keep its hypotheses.
[Hutchinson--Smarandache--Trumpf, 2006](https://arxiv.org/abs/cs/0607089),
Definition 2.8 and the statements of Theorems 3.2--3.3, concern
superregular lower triangular Toeplitz matrices and a specified
block construction. Those statements were inspected; their external
proof dependencies were not audited. They do not justify replacing
the present two-polynomial bivariate convolution by independent
nonzero entries. No result from that paper is used in (TH2)--(TH5).

#### Exact controls and the remaining lower-bound gate

Four mixed-tail fixtures use (a,b,R,S)=(1,2,3,5),(2,2,5,5),
(2,3,5,7),(2,2,9,11). The longer asymmetric arms in the last fixture
test the proof outside the minimum-arm examples. Every allowed tail
coefficient is an explicitly saved positive integer, generated by a
fixed-seed standard-library procedure in sorted monomial order. No
genericity is inferred from those choices. Their supports include all
ab power-substitution residue classes over Q.

The checker verifies the full tangent dimension 2m, adjacent syzygies,
commuting multiplication matrices, cyclicity on the reference domain
and its border, the three block spaces, the interface bounds and the
linear internal-hyperplane budget. At selected strict positive integer
weights it compares descending relation rows with Gaussian residue
selection and then with the full algebra matrices over Q,F_2,F_3,F_29.
It does not enumerate the full fan or extrapolate a count from samples.

A separate 2-by-2 power-substitution control lifts every checked basis
of a length-10 hook quotient to the corresponding length-40 quotient
at rescaled weights, in all four characteristics. The exact uniform
free-extension proof from the preceding pass remains the justification
for equality of full fan counts. The integer Koszul matrix, matching,
null vector and nonzero minor are saved separately.

The next lower construction must therefore grow its thickness at the
required rate, avoid pure power inflation, and respect polynomial
syzygies in the convolution blocks. A support-matching or arbitrary
matroid count alone is insufficient. The concrete next proof target
is a symbolic rank/count description of these two-polynomial blocks,
before any larger enumeration. No family with the required new growth
has been obtained. Keep the unrestricted upper route, strongest-primary
baseline comparison and full goal active; this pass does not establish
publication novelty or complete the requested theorem.

#### The untruncated two-polynomial block can be solved exactly

The symbolic rank/count subproblem has an exact answer before boundary
reduction. This gives a sharper next gate than treating all convolution
blocks as unexplained: the remaining difficulty is specifically the
quotient boundary or a more complicated syzygy module.

Let F,G be nonzero polynomials in S=K[t_1,...,t_d], any field and
dimension. Remove their common factor and call the resulting coprime
polynomials F',G'. Multiplication by the common factor is injective
in S, so this removal preserves every linear dependence of their
multiple columns. Let A=P intersect N^d and B=Q intersect N^d,
where P,Q are bounded convex polytopes. Consider the finite vector
configuration consisting of the two tagged sets

    {t^u F:u in A},  {t^v G:v in B}.

There is no reduction modulo another ideal. Define the finite set

    K_0={h in N^d:h+Newt(G') subset P,
                     h+Newt(F') subset Q},  h_0=|K_0|.       (UC1)

Then its rank is exactly

    |A|+|B|-h_0.                                           (UC2)

Proof. Every relation UF'+VG'=0 is uniquely (U,V)=(HG',-HF'),
by coprimality in the polynomial UFD. The Newton-product identity
Newt(HG')=Newt(H)+Newt(G') says that HG' has support inside P
exactly when every exponent of H satisfies the first containment in
(UC1); similarly for the second domain and F'. Thus H ranges over precisely the
monomial vector space spanned by K_0. This proves the kernel dimension
and the rank. Convex, lattice-complete multiplier domains are essential
to this argument; arbitrary sets with holes are not silently included.

The full coherent-basis description is equally explicit. Assign the
first copy cost vectors p+(0,u) and the second q+(0,v), with disjoint
tags and fixed offsets p,q. Form the tagged Newton polytope

    C=conv({p+(0,g):g in supp(G')}
           union {q+(0,f):f in supp(F')}).

At a generic weight, let c be its unique maximizing tagged term in C.
The leading positions of the whole relation space are exactly

    {c+(0,h):h in K_0}.                                    (UC3)

Indeed, the winner between the two components of (HG',-HF') is the
same tagged term c for every H; its other factor is the leading term
of H. Taking all monomial H in K_0 gives precisely (UC3). Descending
row reduction and the complementary greedy-basis criterion then say
that the coherent basis consists of all input columns except this
translated set. There are at most |supp(F')|+|supp(G')| possible bases.
If h_0>0, distinct exposed tagged terms give different bases; if
h_0=0 there is the single all-column basis.

There is an exact polytope consequence. Write T for the sum of all
input cost vectors and S_0=sum_{h in K_0}(0,h). The projected matroid
basis polytope is

    T-S_0-h_0 C.                                           (UC4)

For every generic weight the minimizing basis has the point
T-S_0-h_0 c by (UC3). These support functions determine the polytope;
the formula also covers h_0=0. Thus the Newton polytope determines
all coherent bases and their normal fan without enumerating minors.
For negative cost weights this is a finite vector-configuration
statement, not an assertion about global monomial orders on an
infinite polynomial ring.

This is derived from the classical two-generator syzygy, the
Newton-product identity and the existing complementary greedy rule.
It is not represented as an established novel theorem. More
importantly, multiplication by a nonzero polynomial need not be
injective after quotient reduction, and truncated products need not
have the Newton sum above. Consequently (UC1)--(UC4) do not settle
the thick-hook boundary problem or the unrestricted goal.

The exact checker takes F=1+y+z,G=1+2y+3z, coprime in all four
checked fields, with three pairs of rectangular multiplier domains.
It checks (UC2)--(UC4) by direct Gaussian selection at positive and
negative finite cost weights. One asymmetric pair has ranks 15 over
Q and 14 over F_2: the support change in G changes the erosion domain,
and the formula correctly accounts for this difference. The next
symbolic target is the truncated version of these blocks, where
boundary annihilation supplies additional relations. The larger
theorem goal and its original acceptance criteria remain active.

### Ninth feasibility pass: boundary amplification is real, but independent arms remain quadratic

The preceding turn made progress in commit 6f44807: it proved the
thick-hook thickness bound and solved the untruncated two-polynomial
rank and basis-polytope problem. This pass supplies an exact boundary
counterfamily and an infinite quadratic hook family. Both have complete
count proofs. They isolate a real mechanism, but their growth does not
improve the unrestricted generic lower order or complete the goal.

#### A factorial Toeplitz determinant

For L>=1 put

    c_L=(L-1)!,
    A_L(z)=c_L sum_{j=0}^{L-1} z^j/j! in Z[z].

Work over characteristic zero, or a field of characteristic p>L-1.
For a,b>=0 with a+b<=L, take the identity columns indexed by
0,...,a-1 and the convolution columns A_L,z A_L,...,z^(b-1)A_L
in K[z]/(z^L). Modulo the identity columns, restrict to rows
a,...,a+b-1. Their determinant is

    det [c_L/(a+i-j)!]_{0<=i,j<b}
       = c_L^b product_{j=0}^{b-1}j! / product_{i=0}^{b-1}(a+i)!,       (FB1)

with a negative-factorial entry interpreted as zero. To prove this,
multiply row i by (a+i)!/c_L. Column j becomes the falling-factorial
polynomial (a+i)_j, of degree j with leading coefficient one. Its
determinant is the Vandermonde determinant at a,...,a+b-1, namely
product_{j=0}^{b-1}j!. All factors in (FB1) are nonzero in the stated
characteristics. This also proves independence of every such pair of
prefixes, not just of the final full basis.

The matrices are diagonal rescalings of the classical lower Pascal
matrix. The determinant argument is included explicitly; neither a
new total-positivity theorem nor a generic-minor assumption is claimed.

#### An actual quotient with L+1 initial ideals

Consider

    I_L=(x^2,z^L-x A_L(z)) in K[x,z].                         (FB2)

Its length is 2L over every field, including the smaller characteristics
excluded from the count assertion. The reference monomial ideal
(x^2,x z^L,z^(2L)) is contained in I_L. The monomials 1,z,...,z^(L-1)
are independent fixed columns. The variable columns are

    x z^j -> x z^j,
    z^(L+j) -> x z^j A_L(z) mod x z^L,  0<=j<L.

Their residue rank is L. The usual direct first-order basis argument,
or the explicit commuting cyclic matrices in the checker, gives the
claimed length and validates these normal forms.

For positive weights, each column family is an arithmetic progression.
The first L columns in their merged order consist of an identity
prefix and a convolution prefix, hence are independent by (FB1).
Thus every coherent basis is exactly such a prefix allocation. If b
convolution columns are selected, the resulting standard staircase is

    {1,z,...,z^(L+b-1)} union {x,xz,...,x z^(L-b-1)},

and its initial ideal is

    (x^2,x z^(L-b),z^(L+b)),  0<=b<=L.                      (FB3)

The redundant generator at either endpoint can be omitted. Every
allocation occurs: choose (w_x,w_z)=(4b+1,2). The last selected
score is strictly below the first unselected score. There are no
additional initial ideals on a comparison boundary, since its finite
term comparisons can be perturbed to generic positive weights.
Consequently, uniformly in L,

    g(I_L)=L+1,  length(I_L)=2L.                              (FB4)

The positive coherent state points are explicitly

    (mu_x,mu_z)=(L-b,L(L-1)+b^2),  0<=b<=L.

Their consecutive fan boundaries are w_x=(2b+1)w_z. Boundary
reduction has created a lower convex lattice chain, rather than just
another translate of the untruncated Newton polytope.

Indeed, retain the same two tagged multiplier sets but do not reduce
products modulo z^L. Their rank is 2L-1, and (UC1)--(UC4) give only
two positively exposed bases, separated by w_x=w_z. For L>=3, weights
(5,2) and (9,2) lie in the same open untruncated Newton cell, yet select
different truncated bases. The number of boundary bases is unbounded
while the tagged Newton polygon has at most three vertices. This
disproves a fan-preserving extension of (UC3)--(UC4) to this quotient;
it is not an asymptotic improvement for G_2 or G_3.

#### Two independent boundary arms make the quadratic hook bound sharp

For t>=2, define the three-variable ideal

    H_t=(x^2,yz,
         y^(t+1)-x y A_t(y),
         z^(t+1)-x z A_t(z)).                               (FB5)

It is a flat first-order deformation of the hook
(y^(t+1),yz,z^(t+1)), with length 4t+2. Its variable residue
configuration splits into a coloop x and two disjoint t-dimensional
arm spaces. Each arm is the one-variable convolution problem above,
with a common monomial shift. Any unassigned reference column is zero.

The y arm can choose b convolution columns, 0<=b<=t, and the z arm
can independently choose c, 0<=c<=t. All choices are simultaneously
realized by positive integer weights

    w_x=(4b+1)(4c+1),
    w_y=2(4c+1),  w_z=2(4b+1).

The two ratios are exactly those in the one-arm proof. Conversely,
independence of the residue spaces means global greedy selection is
the union of the two local selections and the fixed coloop. Therefore

    g(H_t)=(t+1)^2=Theta(length(H_t)^2),                     (FB6)

over characteristic zero or p>t-1. Its coherent moment grid is

    (2t+1-b-c, t(t+1)+b^2, t(t+1)+c^2).

This proves sharpness of the earlier O(n^2) hook bound along lengths
4t+2. The interface dimension is zero after removing fixed summands
and loops. Thus boundary growth and a quadratic count can already
occur without the genuine interaction required for a stronger lower
construction. This family remains below n^2(log n)^2.

#### Exact checks and characteristic scope

The checker covers L=1,...,6 for (FB2), and t=2,3,4,5 for (FB5).
It embeds (FB2) in three variables by adding y to the ideal, preserving
the initial-ideal count.
It verifies every determinant (FB1) with a+b<=L, constructs actual
commuting cyclic quotient matrices, compares relation rows with
Gaussian residue selection, and checks the full quotient once per
realized state. The one-arm internal comparison intervals are exactly
the 2L positive intervals cut by w_x=k w_z, 1<=k<2L. For two arms,
the product of the corresponding 2t intervals covers the whole
positive fan because their residue spaces are independent. Cross-arm
comparisons do not affect the basis. This is an exhaustive fan
verification for these fixtures, not a numerical sample extrapolation.

The exact counts over Q,F_2,F_3,F_29 respectively are:

| L | one-arm counts |
| --- | --- |
| 1 | 2, 2, 2, 2 |
| 2 | 3, 3, 3, 3 |
| 3 | 4, 2, 4, 4 |
| 4 | 5, 3, 2, 5 |
| 5 | 6, 2, 3, 6 |
| 6 | 7, 3, 4, 7 |

The two-arm counts are the squares of the separately computed one-arm
counts at L=t, including the smaller characteristics. The uniform
formulas (FB4) and (FB6) are only asserted under their stated field
hypotheses. The integral defining ideals retain their claimed lengths
in all characteristics; smaller-prime fan changes are saved rather
than hidden. The checker also saves the same-Newton-cell witnesses,
integer coefficients, prefix determinants and exact basis witnesses.

#### Existing module theory and the next mathematical gate

[Naldi--Neiger, 2020](https://arxiv.org/abs/2002.06404) already studies
syzygies in finite-dimensional polynomial quotients, including the
Hermite--Pade form of this boundary problem. The inspection covered
the setup, Schreyer orders, Theorem 4.3 and its proof, Lemma 5.2 and
its proof, and Algorithm 4 with the Proposition 5.5 complexity argument.
It computes a basis for a prescribed order; that operation count does
not count all orders. No new Pade engine or algorithm claim is made.

Independent one-variable boundary arms are now understood sharply
enough to stop enlarging their search. A next lower construction must
have genuinely coupled two-variable boundary relations, growing
thickness as required by (TH4), and more growth than the quadratic
surface above. Its geometric-prefix ranks need a symbolic proof;
arbitrary support matching and one-order algorithm timings do not
supply one. The unrestricted upper route remains an alternative.

The working unrestricted upper bound is still the primary-proof-checked
Babson--Onn--Thomas exponent. Its status as the strongest current bound
has not been certified. The remaining literature-priority audit, a
qualifying unrestricted asymptotic advance and the target manuscript
remain unfinished. Keep the original goal active.

### Tenth feasibility pass: collisions preserve a fan, and sparse prefixes have Hermite defects

The immediately preceding goal turn restated the plan without changing an
artifact or establishing a new mathematical fact. It is classified as no
progress. The last substantive commit, 08625d5, proved boundary amplification
and the sharp quadratic count for independent arms. This pass gives two
structural constraints on further constructions, with uniform proofs and
independent exact controls. Neither is the qualifying unrestricted asymptotic
advance, and neither is being promoted to a separate paper claim.

#### Positive collisions preserve the whole monomial residue matroid

Fix A=(A_1,...,A_d) in positive integers with A_1=1. Over R, let

    q(t)=product_j (t-s_j)^(m_j),  s_j>0 distinct,
    sum_j m_j=n,
    I_(A,q)=(x_i-x_1^(A_i) for i>=2, q(x_1)).              (PC1)

The quotient is R[t]/(q), of dimension n, including when roots have
multiplicity. A monomial x^a maps to t^(A.a). Its normalized Hermite
evaluation column has entries

    binomial(A.a,r) s_j^(A.a-r),  0<=r<m_j,                 (PC2)

with entries zero when r>A.a. For every finite monomial set E,

    rank(E mod I_(A,q))=min(n, #{A.a:a in E}).              (PC3)

Here is a proof that includes multiplicities. A nonzero real polynomial
with r nonzero monomial terms has at most r-1 positive roots counted with
multiplicity. Divide by its lowest power of t, which changes no positive
zero. The derivative has r-1 terms. Rolle's theorem, including the
multiplicity lost at each root and the intervening roots, proves the
assertion by induction on r. Consequently a nonzero combination of at
most n distinct powers cannot vanish to total order n at the positive
nodes. Any at most n distinct degree columns are independent. Equal
degrees give equal columns, proving (PC3). Rank is unchanged by passing
between polynomial remainders and the Hermite coordinates (PC2).

In particular, all ideals (PC1) with fixed A,n have the same monomial
residue matroid, hence the same coherent standard basis for every generic
positive weight and the same monomial initial ideal for every global
monomial order. Their labeled positive Groebner fans coincide as cone
decompositions. This does not identify the coefficients in their reduced
Groebner bases or their nonmonomial weight initial ideals.

There is an explicit flat collision. Set

    q_epsilon(t)=product_{j=1}^n (t-1-j*epsilon).

The corresponding algebra is free of rank n over R[epsilon], since
q_epsilon is monic and x_1 maps to t. For epsilon>0 its fiber consists
of n distinct positive curve samples; at epsilon=0 it is the single
curvilinear point R[t]/((t-1)^n). All fibers for epsilon>=0 have exactly
the same positive fan by (PC3). Thus increasing multiplicity by this
collision cannot amplify the number of initial ideals. This is a
specific fan-preserving family, not a transfer theorem from arbitrary
point configurations to arbitrary local algebras.

The interpolation ingredient has established prior art. In
[Gonzalez-Serrano--Maximenko, arXiv:2312.15680](https://arxiv.org/abs/2312.15680),
Definition 4.1 and Theorem 1.1 give the generalized confluent Vandermonde
matrix and its Schur-polynomial determinant quotient. The inspected
proof includes Proposition 5.1, Proposition 5.4 and the first proof of
Theorem 1.1 in Section 5. The positive-node nonsingularity needed here
is a consequence of that classical determinant structure. No new
confluent Vandermonde or Chebyshev-system theorem is claimed. The
[published article](https://doi.org/10.1080/03081087.2025.2464639) explicitly
acknowledges an equivalent earlier formula; the preprint comparison is
sufficient to rule out claiming the determinant identity anew.

The checker takes A=(1,6,9), n=7, and four positive divisors: seven
simple nodes, one node of multiplicity seven, and two mixed multiplicity
profiles. It builds commuting cyclic quotient matrices independently
from the monic polynomial q, compares direct jets with polynomial
remainders, checks all 86 seven-element staircases, and checks every
one of the existing 37 coherent-basis witnesses in both representations.
The complete count stays 37, including at the single local point; the
generic count at this length is 36. The all-divisor equality follows
from (PC3), not from checking four cases.

#### Which abstract algebras admit the generic corner-cut fan?

Let K have characteristic zero, let B be an n-dimensional commutative
unital K-algebra, and fix d>=1. The following existence statements are
equivalent:

1. B is monogenic, meaning B=K[t]/(q) for a monic polynomial q of degree n.
2. There is a generating d-tuple in B for which every n-element monomial
   staircase is a basis.
3. There is a generating d-tuple whose kernel has the labeled generic
   corner-cut fan: under each global monomial order its standard basis
   is the first n monomials of that order.

In addition, when B is monogenic, a nonempty Zariski open subset of B^d
has properties 2 and 3. This is a generic-embedding statement, not a
claim that every embedding has this fan or that this fan is maximal.

For n>=2 set A_i=n^(i-1). On the finite hyperbolic cross V_n^d, every
coordinate is at most n-1, so a -> A.a is injective. For a scalar c in K
take the coordinate images

    x_i |-> (t+c)^(A_i) mod q.                             (ME1)

The first coordinate is t+c and generates B. To prove simultaneous
independence, choose any distinct degrees lambda_1<...<lambda_n from
this finite ground set. Write the remainder of (t+c)^lambda_j modulo
q in the basis 1,t,...,t^(n-1). Its coefficient in row r is

    binomial(lambda_j,r)c^(lambda_j-r)
       + terms of strictly smaller c-degree,  0<=r<n.      (ME2)

Indeed, the univariate binomial terms of degree below n do not reduce.
Terms t^k with k>=n have c-degree lambda_j-k, strictly below
lambda_j-r; their remainder coefficients are independent of c. A
negative apparent degree in (ME2) denotes a zero leading term.

The determinant therefore has exact c-degree

    E=sum_j lambda_j-n(n-1)/2,

and leading coefficient

    product_{i<j}(lambda_j-lambda_i)
       / product_{r=0}^{n-1} r! != 0.                      (ME3)

This is the falling-factorial Vandermonde formula. Finitely many
nonzero determinant polynomials exclude finitely many values of c.
Since K is infinite, one choice of c makes every n-subset of V_n^d
independent. Every n-element staircase lies in this cross. Also the
first n monomials of any global order form a staircase; their
independence forces them to be the greedy standard basis. This proves
1 implies 2 and 3.

For either converse, the axis staircase
{1,x_1,...,x_1^(n-1)} must be a basis. In statement 3 it occurs for
an elimination order in which every power of x_1 precedes each other
variable. Hence x_1 alone
generates B, proving monogenicity. The case n=1 is immediate. Finally,
each finite-ground determinant is a polynomial in the coordinates of
a tuple in B^d. Their simultaneous nonvanishing defines the asserted
nonempty open subset. The axis determinant also ensures that all its
tuples generate B.

This recovers a broad generic fan within every fixed monogenic algebra,
including nonsplit algebras and algebras with repeated roots. It follows
that even the abstract algebra K[u]/(u^n), with intrinsic embedding
dimension one, has both a one-initial-ideal embedding

    x_1 |-> u, x_i |-> 0 for i>=2,

and an embedding with C_d(n) initial ideals. An explicit instance of
the latter, supported at the origin, is

    x_i |-> (1+u)^(n^(i-1))-1 mod u^n.                    (ME4)

Before subtracting one, (ME3) is exact at c=1 when q=t^n. Subtracting
one changes the columns of any staircase by a unit triangular divisor
transformation and preserves its determinant. Equivalently use the
previously proved translation invariance of all monomial initial
ideals. This preserves staircase bases; it does not preserve the
whole finite monomial matroid. In (ME4) four monomials with positive
degree can all lie in the codimension-one maximal ideal and therefore
be dependent when n=4. The checker includes this distinction.

The argument uses established Vandermonde and open-minor machinery.
[Braun--Pokutta, Example 2.6](https://arxiv.org/abs/0912.1502) already
records independence of all monomial sets of the right size for
algebraically generic reduced points.
[Kreuzer--Robbiano, Theorem 2.7 and proof](https://arxiv.org/abs/1004.1044)
express common border-basis charts by nonvanishing change-of-basis
determinants; Proposition 4.1 and Theorem 4.2 describe the classical
shape representation and open parameter neighborhoods at radical
points. Those statements and their adjacent proofs were inspected.
They are not asserted to provide this particular fixed-algebra
formulation verbatim. The elementary consequence above is saved as a
structural constraint, without an established novelty claim.

Exact controls check every maximal minor on four small hyperbolic
crosses, including q=(t-1)(t^2+1) and q=(t^2+1)^2. In both nonsplit
examples c=1 has a saved zero minor whereas c=2 succeeds: a single
arbitrary shift must not be treated as generic. They also check every
staircase determinant of (ME4) for d=2,3 and 2<=n<=7, and d=4,
2<=n<=5, plus exact forward differences for (ME3). A characteristic-two
failure at n=4 is deliberate. The uniform embedding statement here is
over characteristic zero; no extension to small characteristic is made.

#### A sparse coupled boundary block has a genuine geometric-prefix defect

The one-variable prefix determinant from the ninth pass cannot simply
be extended by assuming that every sparse two-variable geometric prefix
has full support-matching rank. Here is an infinite symbolic obstruction.
Let k>=2 be even, R=5k/2, and put

    B_k=Q[y,z]/(y^R,y^k*z^k,z^R),
    D_k={(i,j):0<=i,j<R, i<k or j<k},  |D_k|=4k^2.

Consider two tagged families in this actual cyclic quotient module:
the identity columns y^i*z^j and the convolution columns
A*y^i*z^j, for (i,j) in D_k, where A is the image of

    ((R-1)!)^2 exp(y+z)

truncated to D_k. All coefficients are integers. Assign input scores
tag*w_x+i*w_y+j*w_z with the positive integer weights

    (w_x,w_y,w_z)=(8k+1,8k+2,8k+4).

The first 4k^2 inputs consist exactly of the identity columns of total
degree <2k and the convolution columns of total degree <2k-1. Their
sizes are k(2k+1) and k(2k-1). The largest selected score is
16k^2-4, and the smallest unselected score is 16k^2+4k-1, so there is
a strict cutoff gap.

Define the nonzero polynomial

    P_k(y,z)=(y-z)^(2k-1)*(y+z).                            (HP1)

It is antisymmetric of total degree 2k, so its central coefficient
y^k*z^k is zero. All its other monomials belong to D_k. At the origin
it vanishes to order 2k, and at (1,1) to order 2k-1. If p_(i,j)
are its coefficients, the nonzero row functional

    ell_(i,j)=i!*j!*p_(i,j)                                (HP2)

annihilates every selected identity column. Its pairing with the
convolution column indexed by (a,b) is exactly

    ((R-1)!)^2 * (partial_y^a partial_z^b P_k)(1,1),

which is zero for a+b<2k-1. Thus these genuine geometric prefixes have
rank at most 4k^2-1 for every even k, without a genericity assumption.
This is a two-point Hermite interpolation defect, not a new theorem
about general fat-point postulation.

For k=2, after contracting the ten identity columns, the remaining
six-by-six matrix has a perfect support matching but rank exactly five.
The report gives the full matrix, matching, exact left null vector and
a nonzero five-by-five minor. The checker also verifies the uniform
null vectors at k=4,6. The next argument determines the exact defect
uniformly; a perfect matching is only asserted for the checked k=2 case.

These are two tagged generator families in B_k. A three-variable ideal
realizing the tag as multiplication has additional powers of that
variable, which have not been counted here. Consequently this is not
an initial-ideal count or an asymptotically stronger construction.
It does rule out the proposed automatic support-matching proof for
the coupled geometric-prefix gate.

#### Graded Jordan strings determine all exponential prefix ranks

There is a reusable calculation behind this defect. Let B be a graded
Artinian algebra over a characteristic-zero field, L in B_1, and

    B = direct_sum_j span(u_j,L*u_j,...,L^(ell_j-1)*u_j),
    deg(u_j)=nu_j,  L^(ell_j)*u_j=0,

be a homogeneous Jordan decomposition. These are the established
Jordan degree type data, not a new invariant. Write B_<a for the sum
of the graded components of degree less than a. For pairwise distinct
scalars tau_1,...,tau_s and integer cutoffs a_1,...,a_s, put

    c_(i,j)=min(ell_j,max(0,a_i-nu_j)).

Then the uniform rank formula is

    dim sum_i exp(tau_i*L)*B_<a_i
       = sum_j min(ell_j,sum_i c_(i,j)).                   (JD1)

All exponentials truncate because L is nilpotent. To prove (JD1),
decompose every prefix along the homogeneous strings. Exponentiation
preserves each string. On a string of length ell the input columns
are t^r*exp(tau_i*t), 0<=r<c_i, modulo t^ell. Their coefficient in
row h is tau_i^(h-r)/(h-r)! for h>=r. If q=sum_i c_i<=ell, the first
q rows, after scaling row h by h!, form the ordinary confluent
Vandermonde matrix on the distinct nodes tau_i. This is nonsingular:
its transpose is the Hermite evaluation map on polynomials of degree
less than q, whose kernel is divisible by a polynomial of degree q.
For q>ell shorten the column prefixes to total size ell. Thus each
string contributes exactly min(ell,q), proving (JD1).

Apply this to B_k above and L=y+z. The known two-variable strong
Lefschetz theorem says all powers of a general linear form have
maximal rank between graded pieces. Since the defining ideal is
monomial, diagonal coordinate rescaling transfers this property to
y+z: the nonempty open set of good linear forms meets the set where
both coefficients are nonzero. The relevant primary statement and
proof are [Harima--Migliore--Nagel--Watanabe, Proposition 4.4](https://arxiv.org/abs/math/0208201),
including its reduction from Borel-fixed ideals via generic initial
ideals. The theorem is credited as an input, not a new Lefschetz result.

The Hilbert function of B_k is

    h_d=d+1                         for 0<=d<2k,
    h_d=2k                          for 2k<=d<R,
    h_d=2(R+k-1-d)                 for R<=d<=R+k-2,

and zero afterward. Maximal rank for every power implies that any two
sets of Jordan strings active in two degrees are nested: their
intersection has the smaller cardinality. Thus the homogeneous strings
have starts and lengths, indexed by 1<=i<=2k,

    nu_i=i-1,
    ell_i=R+k+1-ceil(i/2)-i.                              (JD2)

Equivalently the end degree is R+k-1-ceil(i/2). A string starts at each
of the increasing Hilbert-function steps and the newest strings end
first. This derives (JD2) directly from the displayed Hilbert function
and the cited maximal-rank theorem.

Use (JD1) with nodes 0,1 and cuts 2k,2k-1. The sum of the two prefix
lengths in string i is 4k-2i+1. Its missing dimension is

    max(0,ell_i-(4k-2i+1))=max(0,floor(i/2)-k/2).

Summing over i gives the exact uniform conclusion

    selected-prefix rank = 4k^2-k^2/4.                     (JD3)

The defect is one sixteenth of the quotient length, rather than a
bounded accidental singular minor. Exact Gaussian ranks at k=2,4,6
are respectively 15,60,135 at lengths 16,64,144.

The homogeneous-string decomposition used in (JD1) is proved in
[Iarrobino--Macias Marques--McDaniel, Lemma 2.2](https://arxiv.org/abs/1802.07383v6).
The inspected scope includes its proof and Definition 2.28 of Jordan
degree type. This pass applies that established decomposition and
classical confluent Vandermonde interpolation; novelty of the combined
formula has not been established.

#### A consequence for all tag powers on diagonal weights

The same argument handles every power of an actual tag variable, not
only two selected generator families, on a specified weight slice.
Suppose B=K[y_1,...,y_r]/M with M homogeneous and of finite colength,
and define the actual finite algebra

    K[x,y_1,...,y_r]/(M,x-exp(L)) = B.                     (JD4)

Give x weight alpha>0 and all y variables weight beta>0. Within any
fixed x-power and total y-degree, replace a monomial spanning set by
the homogeneous Jordan basis. These are changes inside equal-cost
groups, so all ranks of cost sublevel spaces, and hence the minimum
basis cost, are preserved. The new configuration is a direct sum over
strings of t^j*exp(i*t), with costs alpha*i+beta*(nu+j).
For each string its first ell pairs (i,j) in N^2 form a downset, and
each fixed i gives a prefix in j. The distinct-node proof of (JD1)
therefore shows these first ell columns are independent. No cheaper
basis exists.

If c_m(alpha,beta) is the sum of the first m scores alpha*i+beta*j
in N^2, the exact restricted state support function is consequently

    min_basis_cost(alpha,beta,...,beta)
       = sum_strings [beta*nu*ell+c_ell(alpha,beta)].       (JD5)

This is a geometric consequence for the full quotient (JD4) in every
ambient dimension. It uses every x-power; x-exp(L) is a polynomial
relation because the exponential truncates. Scaling exp(L) by a
nonzero scalar only rescales columns and gives the same formula.
Equivalently the projected state support on this slice is the sum of
the known planar corner-cut support functions with a fixed translation.
No complexity estimate for computing Jordan data or for an enumeration
algorithm is asserted.

Importantly this is the diagonal slice of the weight space. Unequal
weights on the y variables break the equal-cost Jordan basis changes.
A face under a diagonal weight can contain many distinct full state
vertices. Thus (JD5) does not bound the whole Groebner fan, and cannot
be used as the unrestricted counting theorem.

The exact controls recover string intervals from all graded
multiplication-map ranks, independently of the Hilbert-function formula.
They check all total-degree cut pairs for k=2, selected pairs for k=4,6,
and both L=y and L=y+z on K[y,z]/(y^2,z^2). The last two cases have
the same Hilbert function but different Jordan data and different ranks
at cuts (1,2), so Hilbert data alone are not silently substituted for
the strings. Finally they construct all tag-power columns for k=2,
check their recursive multiplication, and compare actual minimum
coherent-basis costs with (JD5) at four diagonal positive weights.

#### Decision and remaining objective

Do not expand positive collision profiles to hunt for amplification:
(PC3) proves the whole fan is invariant. Do not use locality, maximal
multiplicity or intrinsic embedding dimension one to assert a bound
smaller than the generic corner-cut order: (ME4) forbids that shortcut.
For coupled boundary blocks, use actual Hermite rank conditions,
respecting (HP1), rather than independent support entries. An eventual
construction must still supply all monomial residues of an actual
quotient and a count proof exceeding the generic asymptotic lower order.
The concrete next distinction is between the solved total-degree
filtrations (JD1)--(JD5) and unequal-coordinate geometric prefixes.
Any use of Jordan data for the latter needs a separate argument; simply
applying the diagonal formula would miss the full-dimensional fan.

The complete strongest-primary-bound audit remains unfinished. The
working unrestricted upper remains the proof-checked Babson--Onn--Thomas
bound. None of this pass improves the exponent for G_3 or all-point F_3.
The unrestricted proof routes, further-consequence requirement and
honest-manuscript gate remain active. No new manuscript is started.

Reproduce these controls, together with the prior checks, using

```bash
python3 -B scripts/run_structural_complexity_checks.py
```

### Eleventh feasibility pass: coordinate flags survive the Jordan reduction

The preceding goal turn made progress in commit d58e0c2: it saved the
collision and monogenic reductions and proved the graded-Jordan rank and
diagonal state-support formulas. This pass tests the precise missing
unequal-weight transfer, rather than assuming it from those formulas.

#### The same graded algebra, operator and tag have different full fans

Work over Q and fix, once and for all,

    B=Q[y,z]/(y^2,z^3),  L=y+z,
    F=2*(exp(L)-1)=2y+2z+2yz+z^2+yz^2 in B.

The algebra has length six, Hilbert function (1,2,2,1), and homogeneous
Jordan strings for multiplication by L with (start,length) pairs
(0,4),(1,2). For lambda=0,1,2, choose generating coordinates

    (x,u,v) |-> (F,y+lambda*z,z).                           (CF1)

These are three coordinate presentations of the same fixed B, with
the same L, F and Jordan data. They are actual cyclic finite algebras:
y=u-lambda*v and z=v recover the original generators. Explicitly the
kernel of (CF1) is

    I_lambda=((u-lambda*v)^2,v^3,x-F_lambda(u,v)),
    F_lambda=2u+2(1-lambda)v+2uv+(1-2lambda)v^2+uv^2.        (CF2)

The exact full counts are

| lambda | independent six-element staircases | monomial initial ideals |
| --- | --- | --- |
| 0 | 16 | 11 |
| 1 | 21 | 9 |
| 2 | 32 | 15 |

Thus even the fixed graded algebra together with its fixed operator
and tag element does not determine the full fan without the generating
coordinates. In particular, graded Jordan degree type alone does not
determine the full fan count. This refutes exact determination and a
fan-preserving transfer; it does not rule out bounds in terms of Jordan
data with additional coordinate information.

The distinction is already visible in support functions, without relying
on the full counts. At weights (w_x,w_u,w_v)=(2,1,3), the three minimum
basis costs are 13,12,11 respectively. At diagonal weights (1,2,2) they
are all 10, and at (2,1,1) all eight, as required by (JD5). In fact the
entire diagonal support function agrees, since the change of generators
in u,v preserves each total-degree span. The affine rescaling from
exp(L) to F preserves monomial initial ideals and the same support costs.

For lambda=0 and 2, the ranks of cost sublevel spaces at thresholds
0,1,2,3,4 under (2,1,3) are respectively

    1,2,3,5,6   and   1,2,4,6,6.

The corresponding minimum basis costs are the sums of the six rank
birth times: 0+1+2+3+3+4=13 and 0+1+2+2+3+3=11. These are direct
residue-matrix checks. No extrapolation in lambda or length is made.

Every count in the table exceeds six. The saved witnesses explicitly
exhibit distinct initial ideals inside one strict ordering of the three
variable weights. Thus the full fan cannot be replaced by the six
variable permutations in these inhomogeneous examples.

#### Certificates cover both coherence and its failure

The bounded checker enumerates all 48 three-variable staircases of
cardinality six in each of the three presentations. It constructs the
commuting quotient matrices, verifies cyclicity and (CF2), and compares
independent determinant and inverse calculations for each candidate
basis D. For each independent D, let r_b be the unique residue normal
form of a border monomial b in the basis D. Then D is coherent exactly
when a positive weight satisfies

    w.(b-a)>0 whenever the coefficient of a in r_b is nonzero. (CF3)

Necessity follows from the standard-monomial normal form. For sufficiency,
the polynomials b-r_b lie in the kernel and have their border terms as
leading terms. Those terms generate the complement ideal of D, of
colength six. Since the kernel also has colength six, its initial ideal
must equal that complement. General global monomial orders are covered
by the previously proved finite-comparison positive-weight lemma; no
global scalar representation of an arbitrary order is assumed.

Feasible candidates have saved strict positive integer weight witnesses,
with variable-weight ties broken inside the strict region. Full-ground
Gaussian selection independently recovers their bases. Each independent
but incoherent candidate has a saved nonzero nonnegative dependence
among the normals in (CF3) and the three positivity normals e_i. Its
weighted sum is exactly zero. Were all inequalities strict, pairing
with w would give a strictly positive sum, a contradiction. Thus the
34 infeasibility certificates do not depend solely on polygon clipping.
Dependent candidates have exact zero determinants. All 144 candidates
are accounted for, including all 35 coherent bases.

#### A recent all-order theorem with different hypotheses

[Jonsson Kling--Lundqvist--Mohammadi--Orth, arXiv:2506.24028v1](https://arxiv.org/abs/2506.24028v1)
studies, in characteristic zero,

    (x_1^(m_1),...,x_d^(m_d),(x_1+...+x_d)^k), m_i>=2.

Theorem 4.11 makes the reduced basis depend only on the variable ordering;
Corollary 4.14 gives at most d! bases. The inspected scope is the setup,
Theorem 1.1, Proposition 4.9 and its proof, Theorem 4.11 and its proof,
Corollary 4.14 and its proof, and Proposition 4.15 with its adjacent
argument. The full lattice-path classification proof was not re-audited.

Here d denotes ambient dimension, not colength. This is a strong
restricted theorem and closes pure-power-plus-one-linear-form-power
constructions as an amplification route in fixed dimension. It does
not apply to (CF2), whose graph relation is inhomogeneous, or improve
the unrestricted G_d(n) bound. The existing result must be credited
if a related homogeneous family is used.

#### Later-source audit and the unresolved publication check

[Zhang--Hu--Liang--Dimitrova--Stigler, arXiv:2101.09384v1](https://arxiv.org/abs/2101.09384v1)
was inspected through its setup, Section 3.2 including the Corollary 3.17
count and adjacent arguments, and Section 4's database construction.
Its upper bound counts equivalence classes of finite-field data sets
under linear shifts. It does not count initial ideals of one ideal and
therefore does not replace the unrestricted colength bound.

The [2020 published Stigler--Zhang chapter](https://doi.org/10.1007/978-3-030-42687-3_9)
still needs its full theorem comparison. The publisher's primary
preview confirms the finite-field scope and identifies its bounds as
Theorems 1 and 8, whereas the inspected
[2019 preprint](https://arxiv.org/abs/1907.01080v1) numbers its modified
bound Theorem 4. The full NSF deposit could not be retrieved by either
the web reader or a direct timed download. The preview's references
include Babson--Onn--Thomas, an addition to the preprint's reference
list. Consequently the journal/chapter version has not been treated
as proof-checked or silently identified with the preprint. No purchase,
correspondence or request to an author was made.

The known finite-field cap does not by itself change a maximum over
all fields: field size may exceed colength, and characteristic zero
is included. The working unrestricted theorem remains the
proof-checked Babson--Onn--Thomas Corollary 2.5. The planar universal
order construction in Maruri-Aguilar Theorem 5.1 was rechecked: its
O(n^2) operation count is distinct from both the generic planar
corner-cut count and the number of initial ideals of one ideal.

This is an expanded, bounded primary comparison, not a certification
that no stronger theorem exists elsewhere. Unsuccessful searches and
the inaccessible full chapter do not establish novelty. The strongest
current-bound gate remains incomplete and is recorded accordingly.

#### Research decision

Retain the diagonal formulas with their proved scope. A full-dimensional
extension needs the monomial coordinate flags within each graded piece,
or an independent upper argument; the bare Jordan degree type cannot
determine the answer. Do not enlarge the finite lambda or length search
to guess asymptotics. The next lower proof mechanism must track those
flags while still accounting for all tag powers. The recent d! theorem
rules out falling back to its homogeneous family as a new lower result.

The unrestricted asymptotic goal, the further-consequence requirement,
the primary-priority audit and the manuscript gate remain unchanged.
No qualifying asymptotic advance or new-paper novelty claim is made.

### Twelfth feasibility pass: one-dimensional quotients of binomial algebras

The preceding goal turn made progress in commit b566928: it refuted
the Jordan-only full-fan transfer with exact certificates. The next
construction tested here replaces positive monomial-curve samples by
roots of unity. This has an exact rank model, and deleting an orbit
point really can increase the count. The uniform result below bounds
that increase without assuming monotonicity.

#### A uniform upper bound for the whole extension class

Let J be a binomial ideal of colength q>=2 in K[x_1,...,x_d], let
J be contained in I, and suppose dim_K(I/J)=1. Thus I has colength
n=q-1. Over every field, for d=2 and d=3 respectively,

    g(I)=O(M*alpha(M)),       M=|V_q^2|,                    (BE1)
    g(I)=O(M^2*alpha(M)),     M=|V_q^3|.                    (BE2)

Here alpha is the inverse Ackermann function. In particular (BE2) is

    O(q^2*(1+log q)^4*alpha(q*(1+log q)^2)).                (BE3)

This includes reduced and nonreduced algebras and quotients with a
unique relation involving many residue classes. It is a restricted
consequence of established normal-form and envelope theorems. Its
novelty has not been established. It does not improve the bound for
all colength-n ideals or all n-point configurations.

Here is the complete reduction. Use a monomial basis b_1,...,b_q of
K[x]/J. The binomial normal-form theorem used in (BQ1)--(BQ2) partitions
the nonzero monomial residues on V_q^d into q nonempty parallel classes
S_j, each mapping to the line K*b_j. Write the new one-dimensional
kernel as

    I/J=K*v,   v=sum_j c_j*b_j,   T={j:c_j is nonzero}.

The only dependence among these q reference columns in K[x]/I is v.
A basis consequently chooses one representative from every class
except one class in T. This statement also covers |T|=1, where that
class becomes a loop, and |T|=2, where its two classes merge.

For a generic positive weight w put

    m_j(w)=argmin_{a in S_j} w.a,   c_j(w)=w.m_j(w).

The notation c_j(w) for a cost is distinct from the fixed coefficient
c_j in v. The unique minimum-cost child basis is

    {m_1(w),...,m_q(w)} minus {m_k(w)},
    where k maximizes c_j(w) over j in T.                  (BE4)

To prove (BE4), first minimize the representative independently within
each selected class. Every possible omitted class has a basis cost
equal to sum_j c_j(w) minus its own cost, so the best omission maximizes
that cost. All parent and child standard monomials are covered: for
the same global order, in(J) is contained in in(I), so the child
staircase is contained in the parent staircase, which lies in V_q^d.
The finite-comparison positive-weight lemma covers arbitrary global
orders. No representation of their entire infinite order by one
scalar weight is needed.

For d=2, normalize w_1+w_2=1. Each c_j(w) is a concave polygonal
function with at most |S_j| affine intervals. There are O(M) intervals
in total. The upper envelope max_{j in T} c_j(w) of these line segments
has O(M*alpha(M)) pieces. Adding all O(M) individual breakpoints
determines every representative and the omitted class, proving (BE1).

For d=3, normalize to the closed weight triangle sum_i w_i=1, w_i>=0.
Each c_j is the lower support function of conv(S_j). A three-dimensional
polytope with m vertices has O(m) faces and edges; the same bound holds
for a polygon, segment or point. Its normal map on the weight triangle
therefore has O(|S_j|) edges and total polygon side incidences. Triangulate
each affine patch. Across all classes there are O(M) triangles, and
O(M) lines supporting their projected edges. Call this line set L.

The graph of max_{j in T} c_j is the upper envelope of O(M) triangles.
Its planar map has O(M^2*alpha(M)) total complexity by
[Pach--Sharir, *The upper envelope of piecewise linear functions and the
boundary of a region enclosed by convex plates: Combinatorial analysis*,
Theorem 1'](https://doi.org/10.1007/BF02187732).
Superimpose L. This step needs a count: multiplying the two map
complexities would lose the desired exponent.

The overlay argument in the proof of their Lemma 2.2 applies directly.
The lines have O(M^2) pair intersections. On any one line, each triangle
restricts to a segment, so the restricted upper envelope has
O(M*alpha(M)) pieces. Adding all O(M) lines thus creates at most
O(M^2*alpha(M)) further pieces and faces. On each open cell of this
refined map all class representatives and the maximizing class are
fixed. Equation (BE4) is then fixed, proving (BE2).

Degeneracies do not invalidate the bound. Choose finitely many strict
weights witnessing the different child bases, and preserve their
strict comparisons under an arbitrarily small generic perturbation of
the affine patches and their domains. Apply the envelope and overlay
bounds to the perturbed maps. Their cell count still bounds the number
of original strict bases. Boundary-only ties of a global order are
handled by its finite strict weight representation. The hyperbolic-cross
estimate M<=q*(1+log q)^(d-1) gives (BE3).

The geometric input was checked in the published primary full text:
the segment-envelope discussion, Section 2's Theorem 1 proof, Lemma 2.2
and its proof, Theorem 1', and the adjacent degeneracy argument.
No new envelope theorem is asserted. The combined application needs
its own priority comparison before it could be called novel.

#### Exact models refute monotonicity and cover a proper circuit

For integers h>=2 and e>=0, use the actual kernels of

    K[x,y,z] -> A=K[t]/(t^e*(t^h-1)),
    K[x,y,z] -> B=K[t]/(t^e*(1+t+...+t^(h-1))),
    (x,y,z) |-> (t,t^2,t^3).                              (BE5)

The parent modulus equals (t-1) times the child modulus. Both maps are
surjective since x maps to t, and their lengths differ by exactly one.
For e=0 over C, this deletes the point (1,1,1) from the full h-root
orbit. For e=1 there is also the point at zero. For e=2 both algebras
have a nonreduced component at zero.

The parent classes have labels 0,...,e+h-1. An exponent degree lambda
has label lambda for lambda<e, and otherwise

    e+((lambda-e) mod h).

The child has the single circuit on labels e,...,e+h-1, with every
coefficient equal to one. Labels 0,...,e-1 are coloops. This includes a
proper circuit, rather than relying only on a uniform truncation.

| h | e | parent length | parent initial ideals | child initial ideals |
| --- | --- | --- | --- | --- |
| 7 | 0 | 7 | 11 | 18 |
| 6 | 1 | 7 | 8 | 12 |
| 6 | 2 | 8 | 6 | 8 |

These are complete positive fan counts, not counts of orders sampled.
The checker enumerates all 514 relevant three-variable staircases.
It constructs commuting cyclic quotient matrices independently of the
class formula. Every candidate's determinant agrees with the circuit
criterion and is 0,1 or -1. Each coherent candidate has a strict integer
weight and independent Gaussian checks over Q,F_2,F_3,F_7; every
independent incoherent candidate has a nonnegative zero-dependence
certificate for its strict comparison normals. Each child witness
independently verifies (BE4). Two distinct child initials with the same
parent initial are saved, so the increase is not inferred only from
aggregate counts. The published bound does not assert monotonicity.

#### The prime connection is a known input, not an asymptotic advance

Let zeta be a primitive p-th root for a prime p, fix A_1=1 and other
integer exponents A_i, and take any n distinct parameters t=zeta^j.
For every finite monomial set E the evaluation rank over C is

    min(n, number of distinct A.a modulo p for a in E).    (BE6)

Indeed duplicate residues give identical columns; any square minor
on distinct residues and parameters is nonzero by the classical
Chebotarev theorem, proved in
[Tao, arXiv:math/0308286v6, Lemma 1.3](https://arxiv.org/pdf/math/0308286v6).
Its cyclotomic Lemma 1.2 and the full proof of Lemma 1.3 were inspected.
Thus, with p,n,A fixed, all n-element subsets have the same labeled
monomial matroid and positive fan. This is a direct consequence of
an existing theorem, not a new uncertainty principle.

For consecutive parameters j=0,...,n-1, ordinary Vandermonde already
proves (BE6) for every orbit size h, including composite h. In the
one-point deletion case, any h-1 Fourier columns restricted to the
other h-1 rows are independent for every h: the kernel on its h columns
is generated by the conjugate of the omitted row, whose entries are
all nonzero. Consequently primality is not needed for (BE5) or (BE3).

The full-orbit class-minimum description has an established routing
interpretation. [Sabariego--Santos, arXiv:0705.3631v1, Theorems 3.5--3.6](https://arxiv.org/pdf/0705.3631)
identifies the corresponding lattice ideal and its graded initial ideals
with coherent minimum distance diagrams. Its setup and these statements
were inspected. Their diagram counts refer to graded orders, whereas
our table uses all global orders; those counts must not be substituted
for each other. The periodic interpretation alone is not a new
algebraic correspondence.

#### Decision and the next proof gate

Deleting one orbit point genuinely increases some finite fans, but
(BE3) rules out an n^(2+epsilon) lower family by this operation for
every fixed epsilon>0. It does not rule out an improved logarithmic
factor. Stop enlarging one-point-deletion examples to fit an exponent.

A periodic lower construction must instead analyze a growing deletion
defect, or prove a logarithmic improvement. Before expanding its finite
search, the next concrete comparison is a defect-sensitive count for
the sum of the largest r class-minimum costs. For r=1 the envelope proof
above settles the exponent; for growing r the relevant level complexity
needs a separate bound. The general upper route still needs a
multiplication-compatible mechanism beyond the one-dimensional kernel.

The strongest-current-priority audit remains incomplete. No qualifying
unrestricted asymptotic advance, novelty certification or new manuscript
is claimed. The goal remains active and its acceptance criteria are
unchanged. Reproduce these controls together with the previous checks:

```bash
python3 -B scripts/run_structural_complexity_checks.py
```

### Thirteenth feasibility pass: deletion defect and the dual relation matroid

The preceding goal turn made progress in commit 4790dbb: it proved
(BE1)--(BE3) and saved exhaustive certificates. This pass proves a
defect-dependent extension for uniform truncations, states the exact
reduction for general quotients, and checks an obstruction to confusing
the two. The original unrestricted objective is unchanged.

#### The general reduction requires a matroid, not just a number of omissions

Let J be binomial of colength q, J contained in I, and dim_K(I/J)=r,
where 0<=r<q. Fix the parent monomial basis b_1,...,b_q and the
parallel classes S_j on V_q^d as in the previous pass. Put n=q-r.
Choose an r-by-q matrix H whose rows are a basis of I/J in this
reference basis. For positive generic w put

    m_j(w)=argmin_{a in S_j} w.a,   gamma_j(w)=w.m_j(w).

Then the child minimum basis cost is exactly

    sum_j gamma_j(w)
      - max_{T subset {1,...,q}, |T|=r, det(H_T)!=0}
          sum_{j in T} gamma_j(w).                        (DR1)

For r=0 the empty determinant is one. Here H_T is the square submatrix
on the omitted class labels. To prove the complementary basis criterion,
write R for the row space of H in K^q. The retained reference vectors
on B={1,...,q}\T form a basis of K^q/R precisely when
R intersects K^B trivially. Since dim R=|T|=r, this is equivalent to
the coordinate projection R->K^T being an isomorphism, or det(H_T)!=0.
Minimizing within each retained class now proves (DR1).

Thus the omitted classes form a maximum-weight basis of the matroid
represented by H, the dual of the projected reference-column matroid.
This is standard linear matroid duality, not a new correspondence.
Taking the r largest gamma_j is valid when every r-subset is a basis
of H. A quotient of a binomial algebra need not have this property.

#### A uniform bound when the relation matroid is uniform

Consider q nonempty classes S_j of points in R^d, with M input points
in total. Keep the n=q-r cheapest distinct classes, using their cheapest
representatives. For all 0<=r<q, the number of bases obtainable by
generic positive weights is bounded by

    O(M*(r+1)*alpha(M))       for d=2,
    O(M^2*(r+1)*alpha(M))     for d=3.                     (DR2)

The constants are independent of r. This extends the one-omission
envelope argument using an elementary sampling estimate. It is a
consequence of existing geometric tools; its novelty has not been
established. It is not asserted for a nonuniform relation matroid.

Here is the proof with the dependence on r explicit. Normalize weights
to the simplex. In d=3, triangulate the individual affine patches of
the concave functions gamma_j. As in (BE2), their total number t and
the number of their boundary support lines L are O(M), with t>=q.
For r=0, the arrangement of L already suffices and has O(M^2) cells.
Assume r>=1 and put p=1/(r+1).

Inside a cell of the arrangement of L, each class cost is affine.
An edge separating two different sets of r largest costs exchanges
classes tied at the cutoff. In general position an interior vertex
of this map is a triple equality. Two disjoint tied pairs at different
heights cannot both straddle the same cutoff. Consequently every such
vertex is an intersection of three active triangle interiors, with
at most r other triangles strictly above it.

Independently retain each triangle with probability p. A specified
such triple intersection becomes an upper-envelope vertex whenever
its three defining triangles are retained and all triangles strictly
above it are absent. This has probability at least

    p^3*(1-p)^r >= p^3/4.

If S is the number of sampled triangles, then

    E[S^2] <= (pt)^2+pt <= 2*(pt)^2,

since t>=r+1. The established triangle-envelope bound gives expected
upper-envelope complexity O(E[S^2]*alpha(t)). By summing the exposure
probabilities, the number of relevant interior triple vertices is
therefore O(t^2*(r+1)*alpha(t)). This argument counts events at the
cutoff, rather than the entire arrangement of pair-comparison lines.

It remains to count the overlay with L. Restrict the actual functions
gamma_j to a line in L. Each is the minimum of |S_j| affine functions
of one variable, so their restrictions have O(M) segments in total.
At a cutoff crossing two active segments are equal and at most r
segments are strictly above. Sampling segments with probability p
exposes that crossing on the upper envelope with probability at least
p^2/4. The expected segment-envelope complexity is O(p*M*alpha(M)).
Thus there are O(M*(r+1)*alpha(M)) cutoff crossings on that line.
Individual segment endpoints add O(M) events.

Summing over O(M) lines, and adding their O(M^2) pair intersections,
gives O(M^2*(r+1)*alpha(M)) overlay vertices. Between these vertices
the cutoff edges are straight. Every component meets an endpoint or
vertex, or is a polygonal cycle with vertices. Generic vertex degrees
are bounded, so Euler's formula bounds the refined map's cells by
the same order. Each cell fixes every m_j and the omitted set, proving
the spatial part of (DR2). The segment argument alone proves its
planar part, after adding the O(M) class breakpoints.

Degeneracies can be removed by perturbing the original input points
within their class labels. First choose finitely many strict weight
witnesses for the different bases and avoid every within-class and
between-class tie at those witnesses. Small perturbations preserve
these comparisons. Each perturbed class is still represented by a
convex polytope with O(|S_j|) normal-map complexity in dimension three.
The generic bound therefore also bounds the original strict bases.
The proof uses the segment and triangle bounds credited in the previous
pass; it does not claim a new general envelope or sampling theorem.

#### Arithmetic specialization and a necessary defect for amplification

For the prime cyclic point sets of (BE6), the reference-column matroid
is uniform after identifying equal residue classes. Thus (DR2) applies
to every n-subset of a fixed prime orbit of size q, with r=q-n.
It also applies to consecutive parameter indices 0,...,n-1 for every
orbit size q, by ordinary Vandermonde. With M=|V_q^3| this gives

    g(I) = O(q^2*(r+1)*(1+log q)^4
               *alpha(q*(1+log q)^2)).                    (DR3)

If q=O(n), obtaining g(I)>=c*n^(2+epsilon) from this uniform class
requires, for each fixed epsilon>0,

    r+1 = Omega(n^epsilon /
            ((1+log n)^4*alpha(n*(1+log n)^2))).           (DR4)

In particular any polylogarithmic deletion defect still rules out an
n^(2+epsilon) lower family. This does not exclude logarithmic gains,
prove that the bound is sharp, or bound arbitrary point configurations.

The bounded controls use q=7, A=(1,2,3), n=4,5. There are respectively
10 and 14 full positive coherent bases. All 37 candidate staircases
are classified, including positive-dependence certificates for each
independent incoherent candidate. Every coherent basis has a strict
integer weight and Gaussian checks on two different parameter subsets.

There is also an exact simultaneous certificate for the required
complex Fourier ranks. The prime ell=1000000009 and root z=440359144
satisfy Phi_7(z)=0 modulo ell. The products of all square Fourier minors
of sizes four and five are respectively 173122278 and 176714908 modulo
ell, both nonzero; 1225+441=1666 minors are checked. A complex-zero
determinant in Z[zeta_7] would reduce to zero under zeta_7->z, so these
products certify their nonvanishing over C as well. No inference from
complex nonvanishing to every finite characteristic is made.

#### A composite two-relation counterexample

Take the six-root orbit with t=1 and t=-1 removed, with coordinates
(x,y,z)=(t,t^2,t^3). The exact quotient over Q has modulus

    (t^6-1)/(t^2-1)=t^4+t^2+1.

In the parent reference classes 0,...,5, its relation matrix is

    H = [[1,0,1,0,1,0], [0,1,0,1,0,1]].

At positive weights (10,1,3), the parent class minima, in cost order,
have (class,cost) pairs

    (0,0),(2,1),(4,2),(3,3),(5,4),(1,5).

The unconstrained rule omits classes 1 and 5 and leaves the staircase
{1,y,y^2,z}. This is dependent because 1+y+y^2=0 in the child.
The correct dual basis omits classes 1 and 4 instead, leaving
{1,y,z,yz}, of cost eight. Independent quotient matrices and all 15
complementary minor decisions verify (DR1). This is a realizable
negative control, not an arbitrary matroid counterexample.

Consequently the next general dependency argument must retain the
dual relation matroid and its multiplication compatibility. Uniform
truncation cannot silently replace it, even at defect two.

#### A later parametric-matroid comparison

[Hausbrandt--Ruzika, arXiv:2503.08178v2](https://arxiv.org/pdf/2503.08178v2)
gives O(m^(2p)) optimal-value regions for m ground elements and p
parameters in Corollary 3.4. Its setup, the proof of Theorem 3.2 and
Corollary 3.4, and Theorem 3.6 with its proof were inspected. With two
weight parameters and m=|V_n^3| this gives O(m^4); it does not improve
the proof-checked unrestricted O(n^3) colength bound. Its general
matroid algorithm also does not turn arbitrary matroid examples into
commuting quotient algebras. No algorithmic novelty is claimed here.

The full published finite-field chapter remains unread after a fresh
repository and browser retrieval attempt failed. Its indexed text is
not a full proof check. The strongest-current-priority audit therefore
remains incomplete, rather than being certified by unsuccessful searches.

#### Application assessment: machine learning and transformers

The user asked whether the project's rank and factorization calculations
could apply to machine learning. There is a credible connection through
structured polynomial and exponential features. The established results
here concern matrix ranks and coordinate-ring factorizations; they are
not a new general theorem about CP rank, tensor-train rank, or compression
of arbitrary learned Transformer weights.

For example exp(y+z)=exp(y)*exp(z) has separable coefficients on a
rectangular truncation. But multiplication by exp(L) on any finite
nilpotent algebra is invertible, with inverse multiplication by exp(-L).
Its matrix rank is the full algebra dimension. Separable coefficients,
low tensor separation rank, and low matrix rank are different statements.
Sparse projection and prefix selection also need their own rank argument,
as (HP1)--(JD3) demonstrate.

The closest established application framework is
[Altschuler--Parrilo, arXiv:2106.02755v1, Section 3.1](https://arxiv.org/pdf/2106.02755v1).
Proposition 3.6 already expresses a polynomial kernel's rank after
restriction to an algebraic variety as rank(Phi*C*Phi^T); Corollary 3.8
gives optimal polynomial factorizations. The setup, Lemmas 3.5 and 3.7,
their proofs, Proposition 3.6 and proof, Corollary 3.8, and the
Proposition 3.2 proof were inspected. This prevents claiming the basic
coordinate-ring/kernel compression bridge as our discovery.

For attention, a feature factorization of the unnormalized kernel
A=Phi(Q)*Psi(K)^T permits computing A*V as Phi(Q)*(Psi(K)^T*V), with
the row denominator evaluated separately. This is established efficient
attention methodology. [TaylorShift, arXiv:2403.02920v1, Section 3.2](https://arxiv.org/html/2403.02920v1#S3.SS2)
explicitly uses tensor-product polynomial features; Sections 3.1--3.3
were inspected, including the distinction between its polynomial
kernel and exact exponential softmax. Our special Jordan-prefix rank
formula alone does not apply to arbitrary learned Q and K. Nonreduced
jet algebras likewise need a derivative-feature interpretation before
they can be used as ordinary data-evaluation models.

A classical Vandermonde example gives a useful limit. For N>=2 let scalar
queries be i=0,...,N-1 and scalar keys log(t_j), with distinct positive t_j.
The logit matrix has rank one, but its entrywise exponential is
A_ij=t_j^i, of determinant product_{j<k}(t_k-t_j) != 0. Row-softmax
normalization preserves its full rank N. The checker saves N=4,
t_j=1,2,3,4: det(A)=12 and det(row_softmax)=1/10000. This is not
a new softmax theorem and does not prohibit approximate low-rank methods.

A useful extension would need an explicit input model, a stable
approximate feature factorization, a bound on the resulting normalized
attention error, and a computation cost that includes finding the
representation. Any novelty comparison must include the existing
kernel-on-variety and polynomial-attention work. No ML performance,
new tensor-rank theorem, or change to the active mathematical objective
is asserted. This is an application lead, not a Transformer breakthrough.

#### Research decision

Retain (DR1)--(DR4) with their distinct hypotheses. Stop fixed-defect
uniform periodic searches for an exponent improvement; the uniform theorem
now covers a growing defect quantitatively. Before enlarging that
search, a symbolic family must explain how additive cyclic residue
classes realize many cutoff events. Arbitrary independently chosen
geometric patches do not supply that algebraic realization. For general
quotients the next lemma must instead control the dual matroid in
(DR1); the explicit defect-two example is a required negative control.

The unrestricted asymptotic theorem, further-consequence requirement,
complete priority audit and new manuscript remain unachieved. The goal
stays active. Reproduce the expanded exact checks with

```bash
python3 -B scripts/run_structural_complexity_checks.py
```

### Fourteenth feasibility pass: the exact boundary of residual duality

The preceding goal turn made progress in commit 76827ad: it saved
(DR1)--(DR4), their finite controls, and the ML application assessment.
This pass tests a possible general replacement for the uniform-matroid
assumption: recover a child ideal from its parent and its residual colon.
There is an exact reduction on box-initial cones. The proposed extension
to arbitrary complete-intersection parents is false. Neither conclusion
is the required unrestricted asymptotic advance.

#### An exact colon identity on every box-initial cone

Let S=K[x_1,...,x_d], let J be contained in I, and suppose that for a
global monomial order sigma,

    in_sigma(J)=C=(x_1^s_1,...,x_d^s_d),
    q=product_i s_i,   n=dim_K(S/I),   L=J:I.

Assume 0<n<q and put B={a:0<=a_i<s_i}, c=(s_1-1,...,s_d-1).
Then, over every field and without radicality or comaximality,

    dim_K(S/L)=q-n,
    in_sigma(L)=C:in_sigma(I),                            (RC1)
    D_sigma(L)={c-a:a in B\D_sigma(I)}.                    (RC2)

Here D_sigma denotes the standard exponent staircase. These statements
hold throughout every open cone where the parent has this box initial.

Here is a self-contained proof of the needed duality. In A=S/J define
lambda(a) as the coefficient of x^c in its reduced normal form. For
nonzero a take the leading monomial x^u of that normal form. Since
u belongs to B, multiply by x^(c-u). Its leading product is x^c, and
every other product is smaller under sigma. Reduction cannot increase
a monomial under sigma; x^c itself is standard. Consequently
lambda(a*x^(c-u)) is the nonzero leading coefficient of a. Thus the
bilinear form (a,b) -> lambda(ab) on A is nondegenerate.

Put R=I/J. Its orthogonal subspace is exactly Ann_A(R). One inclusion
is immediate. Conversely, if lambda(ar)=0 for every r in R but some
ar is nonzero, nondegeneracy gives b with lambda(bar)!=0; this is
impossible because br belongs to R. Therefore

    dim_K Ann_A(R)=q-dim_K R=n.

This gives the first identity in (RC1). The usual leading-product
argument always gives

    in_sigma(J:I) contained in in_sigma(J):in_sigma(I).

For the monomial box C, a monomial x^b in B belongs to C:in_sigma(I)
exactly when c-b belongs to D_sigma(I). Indeed, if c-b is outside that
staircase, its product with x^b is the surviving top monomial. In the
other direction, any initial monomial a with a+b<=c would force
c-b outside the staircase by upward closure. The right colon thus has
colength q-n, equal to the left colon. Inclusion and equal finite
colengths prove equality, and the reflected complement proves (RC2).
Applying the reflection twice recovers the original staircase.

This uses standard complete-intersection duality and the elementary
initial-colon inclusion, not a new general liaison theorem.

#### Grid parents: a general reduction, but no automatic smaller instance

For any monic nonconstant univariate polynomials f_i(x_i), the grid
parent J=(f_1(x_1),...,f_d(x_d)) has the same box initial under every
global order. Hence I and J:I have identical positive initial-ideal
fans, with pointwise reflected complementary staircases, for every
ideal I containing J. Repeated roots and local nonreduced quotients
are allowed. In particular, with r=q-n,

    g(I)=g(J:I)<=G_d(r).                                  (RC3)

Every finite-colength I admits a grid parent from its univariate
elimination polynomials. This makes (RC3) a general exact reduction,
but does not make it an asymptotic improvement. For n points with
pairwise distinct values in every coordinate, each of those univariate
polynomials has degree n. Every grid parent then has q>=n^d, and the
minimal one has residual length n^d-n. This includes generic point
configurations. A near-full-grid argument therefore cannot silently
replace the unrestricted problem by one of smaller length.

The established grid-complement comparison is
[Dimitrova--He--Robbiano--Stigler, arXiv:1704.03680v2, Section 4](https://arxiv.org/pdf/1704.03680v2).
Definitions 4.1 and 4.6, Lemmas 4.5--4.7 and their proofs, Theorem 4.8
and its proof, and Corollary 4.11 and its proof were inspected. Their
complementary-ideal definition splits primary components and supplies
comaximal factors; the radical grid case already gives the known
complementary-point fan correspondence. The present colon argument
does not require a Chinese-remainder decomposition. Its elementary
extension is not novelty-certified, and it supplies no new extremal count.

For a concrete noncomaximal case, take

    J=(x^2,y^2,z^2),   I=J+(x+y+z).

Over Q, I has length three and J:I has length five. Both ideals are
contained in (x,y,z), and xyz is a nonzero residue in both I/J and
(J:I)/J: (x+y+z)*yz=xyz modulo J, and xyz annihilates the maximal
ideal. Thus I+(J:I) is proper and I intersect (J:I) is strictly larger
than J. The arbitrary-nonreduced reading of Lemma 4.7(d) in the
inspected v2 cannot be used as written; the displayed witness refutes
its comaximality/Chinese-remainder claims in that scope. The published
version has not been checked here for a correction. The length and fan
duality used in (RC1)--(RC3) instead has the independent proof above.

#### A complete-intersection parent alone is insufficient

Use the already certified cyclic parent and child

    J=(y-x^2,z-x^3,x^7-1),
    I=(y-x^2,z-x^3,1+x+...+x^6).

The actual residual is L=J:I=(x-1,y-1,z-1), by the univariate
factorization x^7-1=(x-1)*(1+x+...+x^6). This remains valid in
characteristic seven. Thus L has the constant initial ideal (x,y,z).

At both weights (41,53,10) and (23,53,10), the parent staircase is

    D={1,z,z^2,z^3,z^4,x,xz},
    C=in(J)=(y,x^2,xz^2,z^5).

The first child staircase removes xz from D; the second removes z^4.
Both removed monomials are in the socle of S/C. Thus these two weights
give the same parent initial and the same residual initial but different
child initials. The parent J is a complete intersection; its displayed
initial C is not. Merely assuming that J is Gorenstein or a complete
intersection does not supply the missing box hypothesis.

Here even in(L)=C:in(I) holds for both children. What fails is invertible
recovery: C:(x,y,z)=C+(xz,z^4) adds both socle monomials, rather than
recovering either individual child initial. The proof route must retain
this information loss on nonbox parent cones. A common refinement of
the two ordinary parent/residual fans does not suffice.

#### Exact controls and an additional primary comparison

The checker constructs three eight-dimensional tensor-product grid
algebras and a principal relation ideal in each, over Q,F_2,F_3,F_5.
The grids include repeated roots. Multiplication matrices commute; the
relation and annihilator spaces are stable under all three actions.
The multiplication kernel and an independent top-pairing orthogonal
calculation give the same colon space. All twelve positive comparison
chambers are checked for every fixture. Descending relation pivots,
ascending quotient-column selection and the monomial-colon test agree:
12 fixtures, 144 chamber checks, 288 independent Gaussian selections.
The local cube has child/residual lengths 3/5 in characteristic zero,
three and five, and 4/4 in characteristic two; both fans have three
initials in each case. No cross-characteristic rank identity is assumed.

Two further checks use the nongrid parent
J=(x^2-y,y^2,z^2), its child J+(x), and residual J+(xy). The weights
(17,7,3) and (5,13,3) give the different parent boxes (2,2,2) and
(4,1,2), respectively; (RC2) holds in both cones. Independent cyclic
quotient matrices recheck both nonbox counterweights over
Q,F_2,F_3,F_7. These are bounded controls of the stated proof, not a
substitute for its uniform argument.

[Onn--Rozenblit, arXiv:1208.5639](https://arxiv.org/pdf/1208.5639)
was checked through the setup, Theorem 2.3 and proof, and Theorem 4.1
and proof. The small-coordinate projection estimate is
O((e(S)*||W||_infinity)^(d*(d-1))), with e(S)=2 for a matroid.
For the monomial exponent matrix on V_n^d the maximum entry can be
n-1, so its direct substitution is O_d(n^(d*(d-1))), or O(n^6) in
three variables. This does not improve the inspected O(n^3) bound.
The published locator is [Melamed--Onn, 2014](https://doi.org/10.1016/j.laa.2014.01.007);
the preprint's author name differs, and no version identity is assumed.
The full later finite-field chapter remains unread after the primary
repository again returned a retrieval error; that remains an incomplete
comparison rather than evidence of novelty or a reason to halt research.

#### Nonuniform fixed defect is a separate mechanism

The fixed-defect obstruction in (DR2)--(DR4) has the uniform-matroid
hypothesis. It does not rule out every nonuniform fixed-defect family.
In particular a general relation matroid can exchange two independent
pairs at different cost levels; there need not be one common cutoff.

Here is an actual quotient demonstrating the distinction. Let D be the
union of the divisor boxes of the six monomials

    x^3, y^3, xz, z^4, yz^3, x^2y,

and let J be its monomial complement. All six are in the parent socle,
|D|=17, and

    I=J+(x^3+y^3+xz, z^4+yz^3+x^2y)

has length fifteen over every field. The two relation rows have disjoint
three-element supports and are killed by every variable, so this is
an ideal rather than an arbitrary vector-space construction. Neither
row nor a nontrivial combination is a binomial modulo J; this child
does not reduce to the binomial quotient case.

At w=(1,1,1), the first circuit has two maximizers of cost three and
the second has two maximizers of cost four. The exchange lines are
w_x=w_y and w_z=w_y. Four adjacent distinct child initials are witnessed
by (101,100,102), (100,101,102), (101,100,99), (100,101,99).
They independently omit x^3 or y^3, and z^4 or yz^3. Thus the
single-cutoff vertex argument cannot be transferred to this matroid.

The checker verifies commuting cyclic quotient matrices over Z,
all 19 comparison chambers over Q,F_2,F_3,F_7, and exactly seven
initial staircases, with 76 independent Gaussian selections. This
example itself is within the full-socle class already bounded by
O(n^(3/2)) in (RSoc1); it is not an amplification family. A potential
nonuniform lower route must supply a mechanism beyond that class.
In particular, growing defect is necessary for the uniform route,
not an established necessity for every possible construction.

#### Decision for the unrestricted goal

Use residual duality only on proved box-initial cones. Do not replace
arbitrary quotient matroids by uniform matroids, arbitrary complete
intersections by box initials, or noncomaximal colon links by disjoint
point complements. Grid enlargement is not a general small-defect
construction: generic n-point inputs already force residual length
n^d-n. The upper route still needs a quantitative treatment of nonbox
initial socles and the compatible relation matroid. The uniform lower
route still needs a symbolic growing-defect family with an actual
additive-residue mechanism producing many distinct states. A nonuniform
route instead needs a compatible mechanism outside the already bounded
socle-support constructions; fixed defect there has not been ruled out.

The unrestricted asymptotic advance, further consequence, complete
priority audit and new manuscript remain pending. The goal is active;
these supporting deductions are not completion or a publication claim.

### Fifteenth feasibility pass: nonsocle exchanges and a global linear obstruction

The preceding turn made progress in commit 1c01f1f by proving the
box-initial colon identity and distinguishing uniform from nonuniform
exchange vertices. This pass asks whether independent exchanges can
persist outside the full-socle construction. They can, in an explicit
infinite family. A uniform global count nevertheless bounds the whole
family linearly, so it cannot supply the required asymptotic lower
advance. This is a route obstruction, not completion of the goal.

#### Primary input: a complicated basis is not a large fan

[Sturmfels--Weismantel--Ziegler, *Groebner bases of lattices, corner
polyhedra, and integer programming*, Section 4](https://www.mi.fu-berlin.de/math/groups/discgeom/ziegler/Preprintfiles/042PREPRINT.pdf)
was inspected through the setup, Theorem 4.2, Theorem 4.4, Lemma 4.5
and its proof, and Remark 4.6 and its proof. Its explicit trivariate
family has colength 2r(r-2) and a reduced basis of size 2r-2. The
exponential statement uses the bit length of r, not colength. This
does not improve a bound on the number of different initial ideals.
The example also prevents assuming bounded initial complexity solely
from a three-variable finite-group presentation.

The published family is a graph of a torus grid after writing x=y/z.
The construction below uses such graph grids with a divisibility
hypothesis; the resulting count is not asserted for every member of
the published family. Its group interpretation is established prior
art. No novelty claim follows from changing the grid periods.

#### The quotient and its multiplication-compatible relation space

Fix even V>=4 and a multiple U of V with U>=V. Over any field put

    J=(xz-y, y^U-1, z^V-1),     q=UV,
    h_0=sum_{a=0}^{U-1} sum_{0<=b<V, b even} y^a*z^b,
    h_1=z*h_0,
    I=J+(h_0),                  n=q-2.

The parent A=S/J is the group algebra on Z_U times Z_V: z is
invertible, x=y/z, and {y^a*z^b:0<=a<U,0<=b<V} is a basis. The
displayed h_0,h_1 have disjoint nonempty even/odd supports. On their
span y acts as identity and x,z interchange them. Thus they span an
ideal R=I/J of dimension two over every field, and I has length n.
All three variable actions on R are invertible. This does not fit the
variable-annihilated full-socle support construction of (RSoc1).

In characteristic zero, over a splitting field, I removes exactly
the points (1,1,1) and (-1,1,-1) from the q-point torus graph.
In characteristics dividing q the quotient may be nonreduced; its
length and the relation-space argument still hold.

Monomials have parent class labels

    (i,j)=(a+b mod U, c-a mod V) for x^a*y^b*z^c.

After choosing the cheapest representative of each class, the child
basis omits its most expensive even-j class and its most expensive
odd-j class. This is (DR1) for a rank-two relation matrix whose columns
are (1,0) or (0,1). It is not a uniform rank-two matroid. The rule and
all basis decisions here are independent of the characteristic.

#### A global bound for all parameters and all global orders

For this entire family,

    g(I)<=18*V^2+34*V+2=O(n).                              (NG1)

The constant in O(n) is absolute, since V^2<=UV=n+2. This is a
restricted family theorem, not an upper bound for G_3(n) or F_3(n).
Here is a uniform proof; write the positive weights as (X,Y,Z).
Finite strict weight witnesses let us avoid the finitely many boundary
equalities in the cases below.

If Y>X+Z, eliminate y=xz. Since V divides U, the remaining parent
is the grid (x^U-1,z^V-1), with fixed box 0<=a<U,0<=c<V.
Its top corner has even class parity. The odd maximum is one of its
two immediate predecessors, according as X<Z or Z<X. Thus this
entire region contributes at most two child staircases.

Suppose next that X>Y. Set A=X-Y>0 and S=A+Z. Reducing xz to y,
x^V to y^V, and the two periods shows that each class (i,c),
0<=i<U, 1<=c<V, has just two canonical candidates:

    x^(V-c)*y^((i-(V-c)) mod U),    y^i*z^c.

The c=0 candidate is y^i. Put a=V-c. If i>=a the difference
between the first candidate's cost and the second's is

    delta_0(c)=a*A-c*Z=V*A-c*S.

If i<a the difference is delta_1(c)=delta_0(c)+U*Y. Each is strictly
decreasing in c. Consequently two integer cutoffs 0<=k_0<=k_1<V
determine every parent representative: for c<=k_0 always choose yz;
for k_0<c<=k_1 choose yz only when i<a; for c>k_1 always choose xy.
There are at most V(V+1)/2 parent diagrams in this region.

Within any one of these regimes, each color maximum has at most six
affine candidates, independently of U,V. To see this, maximize first
over i. For an all-yz row use i=U-1; its maximum increases with c.
For an all-xy row use the largest y exponent U-1; its cost decreases
with c. In a hybrid row the xy maximum has cost
(U-1)Y+(V-c)A, decreasing in c. The yz maximum has cost
(V-1)Y+c(Z-Y), affine in c. For either parity, each applicable interval
therefore contributes an endpoint, or two endpoints for the last
affine expression. Include y^(U-1) separately for the even color.
The budget is at most 1+1+1+2+1=6 per color. Each parent diagram
therefore contributes at most 36 child staircases, one for each
possible pair of maxima.

It remains to consider X<Y<X+Z. Here y^V reduces to x^V and
x^U=1. The latter uses V|U. A canonical xy box is
0<=a<U,0<=b<V; the other candidate remains y^i*z^c. Except when
a+b=i+U, the xy candidate is cheaper. In the exceptional case
a=U-c and 0<=i<V-c, the cost difference is

    delta(c)=U*X-c*(X+Z-Y).

This decreases strictly in c. One integer cutoff 0<=k<V determines
the parent: for 1<=c<=k, replace the xy monomials
x^(U-c)*y^b with c<=b<V by y^i*z^c with 0<=i<V-c.
All other monomials of the xy box remain. There are at most V
parent diagrams.

Each color has at most four maximum candidates in such a diagram.
The full xy rows maximize at their largest a and b=V-1. A truncated
row has maximum cost U*X-Y+c*(Y-X), increasing in c, so take the
largest c of that parity. The yz row maximum has cost
(V-1)Y+c*(Z-Y), so take its two parity endpoints. Thus there are
at most 16 choices of the two maxima per parent diagram.

Adding the three bounds gives 36*V(V+1)/2+16*V+2, proving (NG1).
The proof counts all positive coherent bases; it does not infer a
global bound from weights sampled near the crossing below.

Nor does it assume bounded parent initial socle type. At weights
(U+1,U,U), the parent diagram is

    {x^a*y^b:0<=a<V,0<=b<U-a}
       union {y^b*z^c:1<=c<V,0<=b<V-c}.

It has 2V-1 socle monomials. Its unbounded type is compatible with
the endpoint argument, which controls maxima rather than listing
all corners separately.

#### An infinite family of independent nonsocle exchanges

Take U=V^2 and

    W=(2V-2, 4, 3V^2-7V+10).

This is strictly inside the parent regime k_0=0,k_1=2.
Indeed delta_0(1)=-V^2-V-4<0, while delta_1(1)=3V^2-V-4>0,
delta_1(2)=4V-8>0 and delta_1(3)=-3V^2+9V-12<0.
At W the two even maximizers are

    B=x^(V-2)*y^(U-V+1),    E=y^(V-3)*z^2,

both of cost 6V^2-10V+8. The two odd maximizers are

    A=x^(V-3)*y^(U-1),      C=x^(V-1)*y^(U-V),

both of cost 6V^2-8V+2. Other class maxima are strictly lower.
The two cost levels differ by 2V-6, and their equality planes
are independent.

For s,t in {-1,1}, explicit positive integer witnesses are

    W_(s,t)=(2V-2+s, 4,
             3V^2-7V+10+(V-2)*s/2+t).                    (NG2)

They preserve the parent diagram. Substitution gives C-A=2s and
E-B=2t in cost. The nearby parent inequalities remain strict:
delta_1(2) changes by -2t and the other margins are larger.
The next competing even full row stays at least 4+2s>=2 below B;
the other odd endpoint stays below A,C. Thus all four combinations
are realized for every even V>=4, over every field.

This is a symbolic mechanism outside the preceding full-socle model,
but not an asymptotic amplification: (NG1) bounds its entire fan.
An increasing number of local exchange sectors must not be confused
with a count exceeding n^2(log n)^2.

#### Independent finite controls and decision

The checker implements the threshold regimes and their endpoint lists,
clips their full positive weight regions exactly, and deduplicates the
child staircases. Cayley shortest paths independently check the parent
representatives at every resulting cell. The child multiplication
matrices are constructed over Z, commute, and have the stated cyclic
vector and two-dimensional relation ideal. Independent Gaussian checks
over Q,F_2,F_3,F_7 cover every initial basis in the two smallest fixtures.

| U | V | parent initials | child initials |
| --- | --- | --- | --- |
| 4 | 4 | 11 | 14 |
| 8 | 4 | 15 | 26 |
| 12 | 4 | 15 | 26 |
| 16 | 4 | 15 | 26 |
| 6 | 6 | 22 | 38 |

There are 166 regime/maximum cells and 160 independent Gaussian
checks. Formula (NG2) is checked at V=4,6,8 by independent shortest
paths, giving twelve exact sector witnesses; the uniform proof, not
these three parameters, establishes the infinite family.
The report encodes a staircase by triples (a,c,length), meaning exactly
x^a*y^b*z^c for 0<=b<length; this saves exact witnesses without repeating
every point of a long row. The checker verifies the encoding.

Park this divisible-grid two-point-deletion family as an exponent
amplification route. The nonuniform mechanism is real, but its global
threshold compression prevents the required growth. Retain that
compression argument as a possible ingredient for a broader dependency
bound, without assuming arbitrary quotients have two thresholds or
constant endpoint lists. The unrestricted upper route still needs a
quantitative multiplication-compatible statement beyond this special
group geometry. A lower route must escape the proved endpoint bound,
not enlarge these fixtures.

The strongest-current-primary audit remains incomplete; the full later
finite-field chapter was again unavailable from its public repository.
No qualifying unrestricted asymptotic advance, new-paper novelty claim,
further consequence for the target theorem, or manuscript is asserted.
The original goal remains active.

### Sixteenth feasibility pass: exact moment ranges and the remaining hard regime

The last recorded research pass, commit 0caf537, proved a uniform bound
that closed the divisible-grid amplification route. The intervening ML
assessment did not advance the unrestricted theorem. This pass resolves
a named baseline comparison at the level of its accessible repository
text and proves a parameter-sensitive bound for every finite quotient.
It does not improve the worst-case exponent or meet the goal's manuscript
gate. No further enumeration of the parked family is undertaken.

#### The later finite-field bound: relevant proof now inspected

The [NSF deposit](https://par.nsf.gov/servlets/purl/10177092) exposes the
statement and full displayed proof of Theorem 8 through indexed pages
15--16. The [publisher record](https://doi.org/10.1007/978-3-030-42687-3_9)
identifies the 2020 chapter and its Theorem 8. The deposit bounds ideals
of m rational points over F_p in d coordinates. With r=min(m,p^d-m),
its nonendpoint expression is

    O_d((p^2*floor(r/p)+(r mod p)^2)^(d*(d-1)/(d+1))).

The argument caps coordinate moments, applies a lattice-polytope vertex
bound, and uses complementary point sets. The cap and complement steps
are independently justified below and in (RC1)--(RC3). For p>m and d>=2,
r=m and the displayed base is m^2, recovering the old unrestricted
exponent. This resolves that named bound's scope; it does not supply an
improvement for a maximum over all fields. The complete chapter and a
line-by-line comparison with the publisher PDF remain unread. No identity
with the earlier arXiv version is assumed. Broader priority work remains
incomplete.

#### Exact coordinate ranges for all box staircases

Fix positive integer side lengths s_1,...,s_d, put Q=product_i s_i,
and let B={a:0<=a_i<s_i}. For 0<=n<=Q and coordinate i write
s=s_i and K=Q/s. Define

    U(n,s)=floor(n/s)*binom(s,2)+binom(n mod s,2),
    L(n,K)=K*binom(floor(n/K),2)+floor(n/K)*(n mod K),
    W_i(n)=U(n,s_i)-L(n,Q/s_i).                           (CM1)

For every n-element downset D in B, its moment satisfies

    L(n,Q/s_i) <= sum_{a in D} a_i <= U(n,s_i).            (CM2)

Both bounds are attained by downsets for every permitted n and box.
To prove the upper bound, write its i-parallel column heights as h_j.
Then 0<=h_j<=s, sum_j h_j=n, and the moment is sum_j binom(h_j,2).
If two heights are strictly between zero and s, transferring one unit
from the smaller to the larger does not decrease this convex sum.
Iteration leaves floor(n/s) full columns and at most one nonfull column,
giving U. Put those columns first in increasing lexicographic order on
the other coordinates, with the partial column next. This height
assignment is nonincreasing in the coordinatewise order, so it is a
downset and attains the bound.

For the lower bound, transferring one unit from a height exceeding
another by at least two strictly decreases the sum. The minimum thus
has heights floor(n/K) or floor(n/K)+1, with the longer columns first
in the same lexicographic order. This is again a downset and its moment
is L. These arguments optimize over all bounded column heights; the
explicit downsets show that the relaxation loses nothing. Individual
coordinate extrema need not be simultaneously attained, nor need their
extremizers be initial staircases of a given ideal.

The reflected complement E={c-a:a in B\D}, c=(s_i-1)_i, has length Q-n
and is a downset. Its moments obey the exact translation

    mu_i(E)=mu_i(D)+(Q-n)*(s_i-1)-Q*(s_i-1)/2.             (CM3)

All displayed half-integers combine to integers. Since reflection is
an involution, (CM3) also proves

    W_i(n)=W_i(Q-n).                                     (CM4)

In the colon setting of (RC1)--(RC3), this translates the convex hulls
of the corresponding initial-staircase moments, strengthening the
pointwise fan correspondence to an explicit moment translation.
The reduced grid-point case of this translation is already Lemma 2 of
Maruri-Aguilar--Saenz-de-Cabezon--Wynn (2013); the next pass records
the exact published comparison. It is not a new point-set consequence.

#### A parameter-sensitive bound over every field

Let I have positive colength n in K[x_1,...,x_d], d>=2. Choose monic
univariate polynomials f_i(x_i) in I of degrees s_i. Their grid parent
J=(f_1,...,f_d) is contained in I and has the same box initial for every
global order. Thus Q>=n and every initial staircase of I lies in B.
The established moment-state correspondence and (CM2) put its moment
polytope in the integral translate of a box with sides W_i(n). Adding
the integral unit cube makes it full-dimensional, refines its normal
fan, and increases its number of vertices. Its volume is at most
product_i(W_i(n)+1). The same lattice-polytope estimate used in the
proof of Babson--Onn--Thomas Corollary 2.5 therefore gives

    g(I) <= C_d*(product_i [W_i(n)+1])^((d-1)/(d+1)).      (CM5)

This covers all fields and nonreduced ideals, with no genericity,
comaximality or uniform-matroid assumption. It is a sharpening of the
coordinate cap in (MB1), obtained from established tools; novelty is
not certified. The unit cube also handles zero-dimensional or otherwise
lower-dimensional moment polytopes, where their ambient volume alone
would be zero.

Take f_i to be the minimal coordinate elimination polynomial. Its degree
s_i is at most n by Cayley--Hamilton for coordinate multiplication.
If Q=n, the inclusion J contained in I and equal lengths give I=J and
g(I)=1. Otherwise put r=min(n,Q-n)>=1. By (CM4) and

    W_i(r)+1 <= U(r,s_i)+1 <= r*s_i,

(CM5) implies

    g(I) <= C_d*(r^d*Q)^((d-1)/(d+1)).                   (CM6)

In three variables this is g(I)<=C*r^(3/2)*Q^(1/2).
Consequently, for each fixed delta>0, an elimination degree
s_i<=n^(1-delta) implies g(I)=O(n^(3-delta/2)). Conversely, an infinite
family with g(I)>=c*n^(3-epsilon) must have

    Q=Omega(n^(3-2*epsilon)),
    every s_i=Omega(n^(1-2*epsilon)).                    (CM7)

A lower family beating the generic exponent by a fixed epsilon>0,
g(I)>=c*n^(2+epsilon), requires Q=Omega(n^(1+2*epsilon)). These are
necessary conditions, not sufficient constructions. They make no claim
that a sublinear loss in degree is a linear approximation to cyclicity.

The hard regime is still unrestricted in an essential way. Generic
point configurations have s_i=n in every coordinate. So do the exact
seven-point monomial-curve controls (t,t^6,t^9), t=1,...,7, which have
37 initials versus 36 generic corner cuts. Large coordinate degrees
therefore do not force a generic fan. Arbitrary embeddings of local
curvilinear algebras can also have full coordinate degree. With Q=n^d
and r=n, (CM6) is precisely the old exponent. A universal improvement
still requires a new multiplication-sensitive argument in this regime.

For the finite-field comparison, with 1<=r and p>=2, set

    M=U(r,p),   H=p^2*floor(r/p)+(r mod p)^2.

Direct integer algebra gives H/4<=M+1<=H. Thus the positive cap used
here reproduces Theorem 8's asymptotic expression uniformly in p;
the added one has no effect on that comparison. Its ambient-dimension
convention is justified independently instead of treating a
lower-dimensional polytope's zero volume as an upper bound on vertices.
Rational point ideals over F_p contain x_i^p-x_i. General finite
F_p-algebras do not: F_p[x,y,z]/(x^n,y,z) for n>p is a simple negative
control. Field size cannot silently replace elimination degrees for G.

#### Independent exact controls and next decision

The checker independently optimizes the moment sum by dynamic programming
over all column heights. It checks every n from zero to box volume in
85 specified boxes: all side tuples in {1,2,3,4}^d for d=2,3, together
with (5,5), (5,5,5), (2,3,7), (2,2,2,2), and (2,2,3,3).
It reconstructs and checks 8,398 attaining downsets, their cardinalities,
moments and reflected complements. It also checks 432 coordinate
identities on all saved nonreduced/reduced grid-quotient chambers over
Q,F_2,F_3,F_5, and 192 finite-field bound comparisons for p=2,3,5 and
d=2,3. Full coordinate degrees and the nongeneric fan are checked
against the existing exact point fixture. Three nonreduced finite-algebra
controls verify the failure of the rational-point field equations. These
finite controls support
the uniform proof above; they do not prove an unrestricted exponent gain.

Use the full reproducible command:

```bash
python3 -B scripts/run_structural_complexity_checks.py
```

Retain (CM5)--(CM7) as a general screening lemma and the finite-field
comparison as a resolved named baseline. Do not enlarge small-coordinate
or near-full-grid searches to pursue the worst-case exponent: these
bounds already explain their limitation. The next upper argument must
handle large coordinate degrees and mixed dependencies. A next lower
family must escape all proved coordinate and endpoint caps and supply
a symbolic count before larger enumeration. The original unrestricted
objective, consequence requirement and manuscript gate remain active.

### Seventeenth feasibility pass: grading reduction and its discontinuity

The previous pass made progress in commit f5a9b3e: it proved exact box
moment ranges and a general degree-sensitive inequality. It did not
improve the unrestricted exponent. This pass tests a possible treatment
of the remaining large-coordinate-degree case: reduce by a grading, then
try to transfer the count to a nearby inhomogeneous embedding. The exact
grading bound follows from existing geometry. A uniform flat family shows
why coefficient closeness cannot justify the proposed transfer.

#### Primary comparisons: explicit prior art and an unsafe algorithm lead

[Maruri-Aguilar--Saenz-de-Cabezon--Wynn (2013), published full text](https://www.ism.ac.jp/editsec/aism/pdf/10463_2012_Article_390.pdf)
already proves the reflected-complement correspondence for fractions of
real Cartesian grids. Section 1.2.2 defines the reflection, Theorem 1
and Lemma 1 prove it, Corollary 1 identifies the fans, and Lemma 2 gives
exactly the moment translation in (CM3) in this reduced point-set scope.
The statements, proofs and Lemma 2's direct calculation in Sections
1.2.2--2.1 were inspected. Thus the translated point-state polytopes are
explicit prior art, earlier than the 2019 grid-complement paper. Our
nonreduced colon argument has an independent proof; no novelty is
claimed for its elementary extension.

[Onn, The Complexity of Vector Partition, inspected v1](https://arxiv.org/html/2109.06569v1)
Theorem 2.2 and its proof optimize attribute sums by a dynamic program.
Section 4 describes the established projected-edge zonotope method and
separates fixed parameters from unary and binary coefficient encodings.
These passages were inspected. Their partition-optimization bounds do
not give a sharper count for the multiplication-compatible matroids
in this project. Polynomial optimization and an improved number of
initial ideals remain separate claims.

[Heisel--Lauritzen, inspected v1](https://arxiv.org/html/2510.05103v1)
gives counterexamples to the original Hashemi--Kapur conversion algorithm.
The entire three-page note, including the specified processing order and
its remarks, was inspected. In its F_2 example, the returned polynomials
have leading monomials xz,z^2,y^2z, but an S-polynomial equals y^4 and
cannot reduce by those heads. This is a counterexample to that stated
algorithm, not a new extremal count, and its three-variable example has
positive dimension. A later corrigendum's primary text remains unchecked;
we make no claim that the repaired algorithm still fails. Do not use the
original procedure as a proved enumeration consequence.

The inspected Theorem 3, its proof and Theorem 4 in
[Berstein et al., arXiv:0808.3055](https://arxiv.org/pdf/0808.3055)
concern minimal linear aberration: optimization of weighted degree over
identifiable models. They are not bounds on the number of initial ideals.
The inspected Theorems 1.1 and 1.3 in
[Berstein et al., arXiv:0707.4618](https://arxiv.org/pdf/0707.4618)
likewise concern nonlinear matroid optimization, with explicit weight
encoding assumptions. Their full algorithm proofs were not audited here.
No stronger count is inferred from these optimization statements.
The broader strongest-current-bound audit remains incomplete.

#### Exact grading removes ambient moment dimensions

Let I be a colength-n ideal in d variables over any field. Suppose an
integer h-by-d matrix W of rational rank h defines a grading for which
I is homogeneous; weight columns need not be positive. Every monomial
residue has a degree in Z^h. Any monomial basis of the finite quotient
has the same number of elements in every degree, since the graded
components are a direct sum. Therefore, for every such basis D,

    W*mu(D)=sum_gamma gamma*dim_K(A_gamma),                (HG1)

an integer vector independent of D. In particular the convex hull P
of all initial-staircase moments has affine dimension at most r=d-h.
This uses integer degrees, not a trace identity reduced modulo the
characteristic.

Every moment coordinate lies between zero and binom(n,2). If P has
actual dimension k, its Euclidean k-volume is O_d(n^(2k)): it lies in
a section of that cube, of bounded diameter. Its affine integer lattice
has covolume at least one. Indeed the squared covolume is a sum of
squares of integer maximal minors of a lattice basis. Thus its relative
lattice volume is also O_d(n^(2k)). Applying the established lattice
vertex bound in dimension k gives

    g(I)=1                         if r=0,
    g(I)<=2                        if r=1,
    g(I)=O_d(n^(2r*(r-1)/(r+1)))    if r>=2.              (HG2)

If k<r, use the stronger bound in dimension k; its exponent is no
larger. Segments and points are handled directly. Constants are uniform
in W, including large grading coefficients. For d=3, any nontrivial
grading gives O(n^(4/3)). This is a dimension-sensitive consequence of
known state geometry, not a new unrestricted estimate.

The homogeneity space can be read from any reduced Groebner basis:
it is the real orthogonal complement of all differences between two
exponents in the support of the same basis polynomial. If I is
homogeneous, take the homogeneous component containing a basis
polynomial's marked leading term. Subtracting it leaves an element
of I supported only on standard monomials, so that difference is zero.
The converse follows from homogeneous generators. Using an arbitrary
nonreduced generating set instead would not prove this equality.
The homogeneity space and defining linear comparisons are established
Groebner-fan machinery; [Fukuda--Jensen--Thomas](https://arxiv.org/pdf/math/0509544),
Proposition 2.6 with its proof and the homogeneity-space discussion in
Section 6, were inspected. Their Section 5 does not provide a bound on
flip cost; its reverse-search operation count alone is not a new
colength-only runtime guarantee.

#### Full coordinate degrees and a flat perturbation obstruction

For n>=3 and fixed d>=2 work over Q and set B=Q[t]/(t^n). Put
A_i=n^(i-1) and define coordinate images, with parameter epsilon, by

    x_1=t,
    x_i=t+epsilon*((1+t)^A_i-1-t) mod t^n,  i>=2.          (HG3)

The kernel is generated by x_1^n and the d-1 displayed graph relations.
The quotient over Q[epsilon] is free with basis 1,t,...,t^(n-1), hence
flat; the abstract fiber algebra is always B. For epsilon>=0 every
coordinate has a nonzero t coefficient, so its minimal elimination
polynomial has degree n.

At epsilon=0 all coordinates equal t. A standard basis must select
one monomial of each total degree 0,...,n-1. A downset with that property
is a pure coordinate chain. Each of the d chains is initial by making
its coordinate uniquely cheapest, so g(I_0)=d. Its moment polytope
already has dimension d-1 despite its constant number of vertices.

For any n-element staircase D, its basis determinant P_D(epsilon) is
an integer polynomial. At epsilon=1 the tuple is (ME4), so every such
determinant is nonzero by the simultaneous Vandermonde proof. Thus
P_D is not the zero polynomial. There are finitely many staircases,
so every sufficiently small nonzero real epsilon makes all of them
bases at once. Under any global order the first n monomials are a
staircase; their independence forces them to be the standard basis.
Consequently

    g(I_epsilon)=C_d(n)
      for every sufficiently small nonzero epsilon,      (HG4)

while g(I_0)=d and every coordinate degree remains n. This does not
exceed the known generic lower bound. It prevents a continuity argument
from transferring (HG2) to merely near-homogeneous coefficients.

There is a uniform explicit rational choice, rather than only an
existence argument. Set A=n^(d-1), B_0=2+n*A^(n-1), and

    H=n!*B_0^(n*(n-1)/2).

The sum of absolute coefficients of each coordinate image, considered
as a polynomial in epsilon and t, is at most B_0. An n-element downset
has sum of total degrees at most n(n-1)/2: listing it in increasing
total degree, the j-th degree is at most j, starting with j=0. The
coefficient sum, and hence every coefficient's absolute value, of P_D
is at most H by the determinant expansion. After removing its power
of epsilon, its constant coefficient is a nonzero integer. For
0<|epsilon|<=1/(2H+2), the sum of absolute values of the other terms
is at most H*|epsilon|/(1-|epsilon|)<1. Hence P_D(epsilon)!=0 for all D.
One can choose epsilon=1/(2H+2), or any smaller positive rational;
its denominator has O_d(n^3 log n) bits. Choosing it smaller still
makes every graph coefficient as close to its epsilon=0 value as any
prescribed positive tolerance. No numerical rank decision is involved.

This argument uses (ME4) and standard open-minor reasoning. It does
not claim a new flatness theorem or established novelty. It supplies
a precise obstruction to the proposed upper-bound transfer in the
same full-coordinate-degree regime left open by (CM7).

#### Exact controls and research decision

For n=3,...,7 the checker tests all staircases against the actual
quotients t^n=0 with coordinate degrees (1,1,1), (1,2,3), (1,2,4).
The 531 independent determinant decisions agree with their degree-class
criterion. Every basis moment obeys (HG1), and the affine dimension is
checked independently. For the diagonal tuple it verifies that the
only bases are the three coordinate chains.

It also constructs (HG3) at five exact positive rational parameters,
chosen so every graph coefficient changes by less than 1/1000000.
All 177 staircase determinants are nonzero. The report saves their
integer-scaled determinants, coordinate numerators and denominators,
full coordinate elimination degrees, and the rank-three support
comparisons establishing trivial homogeneity space. At n=7 the count
changes from three to 36. These finite checks certify their fixtures;
the uniform determinant-height argument proves (HG4) for all n and d.
Reproduce them together with the earlier controls using

```bash
python3 -B scripts/run_structural_complexity_checks.py
```

A candidate improving the generic lower order in three variables must
have trivial homogeneity space, as well as escape the previous degree
and endpoint caps. Do not try to bound all nearby inhomogeneous ideals
by the fan count of a graded limit, even with fixed abstract algebra
and full coordinate degrees. Continue on genuinely mixed dependencies
with trivial homogeneity, rather than enlarging these perturbation
fixtures. The baseline audit, qualifying unrestricted advance, further
consequence and new manuscript are still incomplete. The goal remains
active; none of the known-tool consequences is substituted for success.

### Eighteenth feasibility pass: minimum product ranks and the information they lose

The intervening ML assessment did not change the theorem goal or advance
its unrestricted bounds. This pass tests a dependency mechanism in the
remaining full-coordinate-degree regime. Multiplication ranks are a
natural input to an upper argument, but a complete rectangular rank
profile loses information about the initial-ideal fan. The minimum-rank
case can nevertheless be classified uniformly, including nonreduced
algebras. This is a diagnostic endpoint, not the requested unrestricted
asymptotic advance or a novelty claim.

#### Primary comparison: code products and universal order lists

[Mirandola--Zemor, v2](https://arxiv.org/html/1501.06419v2#S3)
Theorem 3.3 proves the code-product Kneser inequality over every field,
using coordinatewise multiplication in K^n. Definitions and Lemmas
2.7--2.10, Lemma 3.6 and its proof, and the complete Theorem 3.3 proof
were inspected. Stabilizers give direct-sum decompositions of these
function spaces. Theorem 4.2 characterizes minimum-growth MDS pairs
under additional dimension hypotheses; its statement, rather than
its complete proof, was inspected. Thus code products and their inverse
problems are established tools. They do not directly count initial
ideals. A split point algebra is K^n, whereas an arbitrary nonreduced
quotient is not this algebra. The argument below proves its nonreduced
scope independently.

[Dimitrova--He--Stigler--Zhang, v2](https://arxiv.org/html/1811.01114v2#S3)
Theorem 3.12 and its proof explicitly give invariance of point initial
ideals under separate affine changes of coordinates. Corollaries
3.13--3.14 address equal counts and a sufficient condition for uniqueness.
These passages and the introduction's existing general bound were
inspected beyond the earlier abstract-only comparison. They supply no
stronger unrestricted extremal estimate. The affine invariance used
below is therefore also explicit prior art in the point setting.

[Maruri-Aguilar (2005)](https://rendiconti.dmi.units.it/volumi/37/05.pdf)
Definitions 4.1--4.5, Theorem 5.1 and its full proof, and the relevant
discussion and three-dimensional example in Sections 5--6 were now
read in the complete published PDF, beyond the earlier indexed excerpts.
The positive planar Hilbert-zonotope fan has primitive rays in the
integer square [0,n-1]^2. This yields a common O(n^2) order list for
all length-n planar inputs, not an O(n^2) improvement to the smaller
O(n^(4/3)) bound for one planar ideal. Its stated operation count does
not by itself establish bit complexity. The paper explicitly prevents
extending that ray description to three variables. The broader
strongest-current-bound audit remains incomplete.

#### A general additive rank lower bound

Let A=K[x_1,...,x_d]/I have dimension n and suppose each coordinate
alone generates A: [K[x_i]:K]=n. For positive integers a_i let C(a)
be the span of the monomials with 0<=u_i<a_i. Then

    dim C(a) >= min(n,1+sum_i(a_i-1)).                     (MR1)

Indeed, if a subspace U contains 1 and U+x_i U=U, it contains every
power of x_i, hence is A. Otherwise its dimension increases by at
least one. Starting from K*1 and successively applying U -> U+x_i U
a_i-1 times for each coordinate produces exactly C(a). This proves
(MR1) over every field and without radicality. The line presentation
x_i=t+c_i in K[t]/(f), deg(f)=n, attains equality for every box.
Thus full elimination degrees do not imply volume-proportional box
rank. No stronger bound follows merely by replacing the elementary
argument with Kneser's additive inequality.

#### Minimum bilinear rank forces fractional-linear coordinates

Suppose n>=3, all coordinates generate A, and

    dim span_K(1,x_1,x_i,x_1*x_i)=3,  i=2,...,d.           (MR2)

Then A=K[t]/(f), t=x_1, deg(f)=n, and every coordinate is a
nonconstant fractional-linear function of t whose denominator is a
unit modulo f. Here f may have repeated factors or be irreducible.

To prove this, write a nonzero bilinear relation as

    A_i+B_i*t+C_i*x_i+D_i*t*x_i=0.

Its two-by-two coefficient determinant A_i*D_i-B_i*C_i is nonzero.
If it were zero, the polynomial would factor into two nonzero affine
linear factors, one in t and the other in x_i. Multiplication by a
nonzero affine function of a cyclic coordinate has rank at least n-1:
in its univariate presentation, the gcd with its degree-n minimal
polynomial has degree at most one. Sylvester's rank inequality would
then give product rank at least n-2>0, contradicting that relation.
This also covers a constant factor, whose multiplication rank is n.

Put D(t)=C_i+D_i*t and N(t)=-A_i-B_i*t. These are coprime affine
polynomials. A Bezout relation alpha*D+beta*N=1, together with
N=D*x_i in A, makes D invertible. Consequently

    x_i=(-A_i-B_i*t)/(C_i+D_i*t) in A.                    (MR3)

The bilinear rank cannot be below three: either 1,t,x_i are already
independent, or x_i=a+b*t with b nonzero and t*x_i introduces t^2,
independent of 1,t since n>=3. Conversely (MR3) gives rank three.
The rational-function representation is unique: the numerator of the
difference of two representations has degree at most two<n, so it
cannot be a nonzero multiple of f. This makes its pole intrinsic.

Group the coordinates by their single simple pole on P^1(K), including
the pole at infinity of t. Let the group sizes be m_1,...,m_r. Two
coordinates have the same pole exactly when they are affine functions
of each other. This is also equivalent to
dim span(1,x_i,x_j)=2. The equivalence holds modulo f by the same
degree-two argument. These three-monomial flags retain information
that the four-monomial box rank in (MR2) omits.

The same conclusion follows from one larger minimum-growth box. Suppose
all a_i>=2 and

    dim C(a)=1+sum_i(a_i-1)<n.                            (MR2a)

Then (MR2) holds. Here is an inverse proof that does not use MDS or
reducedness. First take two coordinates x,y and sides a>=b>=2 with
rank a+b-1<n. Growing the x-prefix one power at a time increases rank
by at least one by (MR1). Equality at the final step forces the
two-by-b rectangle to have rank b+1. The relation space

    V={(P,Q):deg P,deg Q<b, P(y)+x*Q(y)=0}

has dimension b-1. For any two relations, eliminating x gives
P_1*Q_2-P_2*Q_1 evaluated at y equal to zero. Its degree is at most
2b-2<n, since n>=a+b>=2b and y has minimal polynomial degree n.
It is therefore the zero polynomial. Fix one relation and divide its
two polynomials by their gcd to get a primitive pair (P_0,Q_0).
Every relation is h*(P_0,Q_0) for a polynomial h. Put
m=max(deg P_0,deg Q_0). The multiplier space has dimension at most
b-m, so dim V=b-1 forces m<=1.

The case m=0 is impossible. For a nonzero h of degree at most b-1,
multiplication by h(y) has rank at least n-b+1. Multiplication by a
nonzero affine function P_0+x*Q_0 has rank at least n-1, so their
product has positive rank at least n-b. If Q_0=0, the second factor
is a nonzero constant and gives the same contradiction. Thus m=1.
The multiplier space now has dimension at most b-1, exactly dim V,
so all multipliers of degree at most b-2 occur, including 1. Hence
P_0(y)+x*Q_0(y)=0 itself is a bilinear relation, proving (MR2).

For a d-dimensional box in (MR2a), adding the other coordinate powers
to any two-coordinate subbox increases rank by at least one per step.
Equality of the full box forces equality in every such rectangle.
The two-coordinate proof therefore applies to each pair with x_1.
Consequently a non-fractional-linear cyclic presentation satisfies

    dim C(a)>=min(n,2+sum_i(a_i-1)) whenever all a_i>=2.   (MR2b)

This extra unit is an inverse boundary, not an asymptotic fan bound.
The strict inequality in (MR2a) is essential. In K[t]/(t^5), set
x=t and y=t+t^2. Both coordinates generate the algebra over every
field. The three-by-three rectangle has rank five, attaining its
additive lower bound after saturation, while span(1,x,y,xy) has
rank four. Thus a saturated minimum is not a certificate for (MR3).

#### Exact fan count for this entire class

Independent affine changes and nonzero scalings of coordinates preserve
every initial ideal under a global order: each monomial expands into
itself plus divisors. Normalize a representative at a finite pole
lambda to u=1/(t-lambda), and the infinity representative to t.
All finite-pole denominators are units modulo f.

Within one pole group, every variable except the cheapest one is a
leading term of an affine relation. Between different groups, every
mixed product is also a leading term. For distinct finite poles,

    u_i*u_j=(u_i-u_j)/(lambda_i-lambda_j),

and t*u_i=1+lambda_i*u_i. The product dominates all the terms on
the right in every global monomial order. Therefore any standard
staircase consists of 1 and pure coordinate chains, with at most one
active coordinate in each pole group.

Every such staircase of size n is independent. If its chain lengths
are ell_j, multiply a putative relation by the common denominator
product_(finite j)(t-lambda_j)^ell_j. The resulting polynomial has
degree at most sum_j ell_j=n-1. It vanishes modulo f, so is zero as
a polynomial. Uniqueness of partial fractions then makes every
coefficient zero. This proof includes all nonreduced cases.

Every allocation is coherent as well. Give its selected coordinate in
group j weight 2/(2*ell_j+1), and all other coordinates weight four.
The prescribed pure powers have cost below one, and the first omitted
pure power on each axis has cost above one. Mixed monomials reduce to
cheaper pure powers. Since the n prescribed residues are independent,
they form the initial staircase. Multiplying the weights by a common
denominator provides exact positive integer witnesses. Thus

    g(I)=[z^(n-1)] product_(j=1)^r (1+m_j*z/(1-z)).        (MR4)

An inactive group contributes one choice, not m_j choices. In three
variables the pole multiplicities (3), (2,1), (1,1,1) give respectively

    3,  2*n-1,  binom(n+1,2).                             (MR5)

For d distinct poles, (MR4) is binom(n+d-2,d-1). Its growth is below
the established generic n^(d-1)*(log n)^(d-1) order. This is an exact
restricted classification, not an improved unrestricted lower bound.

#### The complete box profile and all its stabilizers still lose the fan

In fact every presentation in (MR3) attains equality in (MR1) for
every anchored box, irrespective of its pole grouping. At pole j its
largest possible order is L_j=sum_(i in group j)(a_i-1). Partial
fractions show that the box spans precisely all rational functions
with these pole bounds. Multiplication by the common denominator
identifies this space with polynomials of degree at most
L=sum_i(a_i-1). The denominator is a unit modulo f, so

    dim C(a)=min(n,L+1).                                  (MR6)

If C(a) is proper, its stabilizer in A consists of scalars. Unit scaling
does not change a stabilizer in a commutative algebra, so reduce to
U=span(1,t,...,t^L), L<n-1. Any stabilizer element lies in U because
1 is in U. If its polynomial degree s were positive, multiplication
by t^(L-s+1), an element of U, would have degree L+1<n and could
not reduce back into U. Thus s=0. A full box span instead has
stabilizer A, of dimension n.

Over Q, choose f=t^n or f=product_(k=1)^n(t-k), and coordinates from
t,t+1,t+2,1/(1+t),1/(1+2*t) as appropriate. These give infinite
families with identical complete box ranks and stabilizer dimensions,
but fan counts ranging from three to binom(n+1,2). They also have
trivial homogeneity space: the displayed affine or inverse relations
are minimal monomial relations, so homogeneity would force all their
terms to have equal degree, and those equalities force every variable
degree to zero. This reasoning uses minimal supports, not an inference
from an arbitrary generating set. Hence the examples lie in the
full-coordinate-degree, trivial-homogeneity regime left by (CM7)
and (HG2).

This does not rule out an upper bound expressed through box ranks.
It rules out recovering a fan, or transferring a count within a bounded
factor, from those ranks and stabilizer dimensions alone. The minimum
case has the additional fractional-linear classification (MR3)--(MR4).
A useful general inverse argument must retain mixed monomial flags
and their multiplication compatibility when higher products have
larger rank. The critical-pair literature must be compared before
claiming any new inverse theorem.

#### Exact controls and decision

For n=3,...,7 the checker constructs 85 genuine quotients: reduced
point evaluations over Q and F_17, and nilpotent companion presentations
over Q,F_2,F_3,F_11. The third distinct-pole coordinate is deliberately
excluded in characteristic two, where 1/(1+2*t) is constant. No
characteristic-zero rank is silently transferred to a finite field.

It checks all 3,009 staircase decisions by determinants and independent
Gaussian ranks, and all 760 claimed initial ideals by full-ground
Gaussian selection at strict integer weights. It also verifies cyclicity
of every coordinate, the bilinear relations, their nonsingular coefficient
matrices, minimal-support grading constraints, 2,295 box ranks and 120
polynomial-prefix stabilizers. The report saves exact coordinates or
coordinate polynomials, determinants, weights and the complete staircase
decisions. These are finite controls for the uniform proof above.
Nine additional cyclic nilpotent controls over Q,F_2,F_3 check 81
rectangles for x=t,y=t+t^2 at lengths 5,7,9, including the saturated
three-by-three counterexample. They verify (MR2b) against actual
residue matrices; they do not substitute for the elimination proof.

Do not expand this minimum-rank family as a lower-bound search: (MR4)
already closes it below the generic order. Do not build a larger search
around box-rank or stabilizer dimensions alone. The next general upper
lemma must control compatible mixed flags beyond the minimum-growth
case; a lower family must escape this classification and the earlier
caps. The unrestricted asymptotic gap, complete priority audit, further
consequence and separate new manuscript remain open. Keep the original
goal active. Reproduce the controls with

```bash
python3 -B scripts/run_structural_complexity_checks.py
```

### Nineteenth feasibility pass: rank excess, exceptional ideals and extension bounds

The previous pass made progress in commit a694c84 by classifying minimum
product growth and identifying information lost by complete box-rank
profiles. This pass considers positive excess rather than expanding the
solved minimum family. In a cyclic algebra, a short relation module gives
an exact trade-off between rational degree and an exceptional ideal. An
elementary extension bound then controls the branch where the retained
coordinates are fractional-linear. Neither result closes the unrestricted
asymptotic gap.

#### Primary comparison: rational reconstruction is established machinery

[Neiger--Safey El Din--Tran, arXiv:2609.21051v1](https://arxiv.org/html/2609.21051v1)
defines the relation module R_M(a,b)={(p,q):ap+bq=0 mod M} in
Section 1.1. Its arbitrary-field setup, reduced/Popov definitions in
Section 2.2, and Section 2.3's Lemma 3 with its proof and Corollary 2
with its proof were inspected. These already provide relation bases and
their determinant characterization. The later fast algorithm proofs
were not audited. [Moreno Maza's primary course text](https://www.csd.uwo.ca/~mmorenom/CS424/Lectures/bugs/EuclideanMethods.html/node8.html)
distinguishes weak polynomial congruences from rational reconstruction
with a unit denominator; Theorem 5 and its proof were inspected. Thus
the reconstruction mechanism below belongs to established interpolation
and relation-module mathematics. No new reconstruction algorithm or
complexity bound is claimed. These sources do not supply an inspected
improvement to the unrestricted initial-ideal count. The strongest-current
baseline audit remains incomplete.

#### An exact trade-off for two-by-b products

Let A=K[t]/(f), deg(f)=n>=3, over any field, and let x also generate
A. For b>=2 assume 2b-2<n. Set

    R=dim span(1,t,...,t^(b-1),x,x*t,...,x*t^(b-1)).

If R<2b, put e=R-b-1. Then e>=0, and there are coprime polynomials
P,Q, with Q nonzero, and an ideal E of A such that

    m=max(deg P,deg Q)>=1,
    delta=dim_K E,
    R=b+m+delta,  hence m+delta=e+1,                      (EX1)

and, on the quotient A/E of length n-delta, Q is a unit and
x=-P(t)/Q(t). This includes repeated or irreducible f; delta is an
ideal dimension, not a count of distinct points unless A is split and
reduced.

Here is a complete module proof. The space

    V={(U,W):deg U,deg W<b, U(t)+x*W(t)=0 in A}

has dimension 2b-R>0. Eliminating x from any two relations produces
U_1*W_2-U_2*W_1 of degree at most 2b-2<n, zero modulo f, hence
zero as a polynomial. Divide one nonzero pair by its polynomial gcd
to obtain the primitive pair (P,Q). Every pair in V is h*(P,Q),
with h a polynomial. This is polynomial proportionality, using
coprimality; it does not cancel h(t) in A.

Put r=P(t)+x*Q(t) in A and E=r*A. Represent x by its polynomial
of degree below n. Write f_0=gcd(f,r) and h_0=f/f_0, with monic
gcds. Then

    dim E=n-deg(f_0)=deg(h_0)=delta.

A multiplier h gives a relation exactly when f divides h*r, equivalently
h_0 divides h. The truncated relation space therefore has the basis

    h_0*t^j*(P,Q),  0<=j<b-m-delta.

Its dimension is b-m-delta, proving (EX1). If m=0, r is a nonzero
affine function of the cyclic coordinate x, or a nonzero constant.
Its multiplication rank delta is at least n-1. But (EX1) and R<2b
give delta<=b-1<n-1, a contradiction. Thus m>=1 and delta<=e.

On A/E, the primitive relation P+xQ is zero. A polynomial Bezout
identity U*P+W*Q=1 gives Q*(W-U*x)=1 there, so Q is a unit. This
proves the rational representation without assuming that the common
factor removed from a weak reconstruction was already a unit in A.

For a larger a-by-b rectangle with a>=b>=e+2 and rank
a+b-1+e<n, the cyclic growth argument in (MR1) implies that its
two-by-b subrectangle has rank at most b+1+e. Also n>=a+b+e>2b-2.
Applying (EX1), with its actual smaller excess, yields

    rational degree + exceptional dimension <= e+1.      (EX2)

The same statement can be applied to each coordinate against one
common cyclic parameter t. Taking the sum of the exceptional ideals
retains one common quotient; its lost dimension is at most the sum
of their individual dimensions. This is coordinate-compatible quotient
construction, not a comparison between unrelated abstract presentations.

#### Why the alternatives and degree threshold matter

Excess one has two possibilities: m=2,delta=0, or m=1,delta=1.
Both occur over every field in nilpotent examples of sufficient length.
For x=t+t^2, the relation is already quadratic and exact. For
x=t+t^(n-1), the primitive relation x-t fails in A but holds after
quotienting by the one-dimensional ideal (t^(n-1)). The weak relation
t*(x-t)=0 must not be divided by t inside A.

At excess two, the examples x=t+t^3, x=t+t^2+t^(n-1), and
x=t+t^(n-2) realize (m,delta)=(3,0),(2,1),(1,2). Unit rational
denominators and reduced one-outlier examples are checked separately.
Small rank excess alone cannot identify which branch occurs.

The degree threshold is real. In K[t]/(t^6), take x=t+t^2+t^3 and
b=4. The rank is six and there are two independent short relations.
Their eliminant can be a nonzero multiple of t^6, of degree 2b-2=n.
They need not have the same rational ratio as polynomial pairs. The
argument above therefore does not extend across the degree threshold
by treating equality modulo f as equality in K[t].

#### A general bound for coordinate-compatible extensions

Let I be contained in J in K[x_1,...,x_d], of respective finite
colengths n and n-e. Suppose every monomial initial ideal of J has
at most q minimal monomial generators. Then, over every field,

    g(I)<=g(J)*product_(j=0)^(e-1)(q+(d-1)*j).            (EX3)

For a fixed global order, in(I) is contained in in(J). Its staircase
is obtained by adding exactly e boxes to the staircase of in(J).
Any finite downset extension admits an ordering of its added boxes
in which every intermediate set is a downset. Its first addable boxes
are precisely the minimal generators of in(J). Adding one box removes
that candidate and creates at most d new candidates, its immediate
successors. After j additions there are at most q+(d-1)*j choices.
Their product counts all ordered addition paths and hence bounds the
number of resulting staircases. Sum over the at most g(J) starting
staircases to prove (EX3). Multiple paths or starting staircases leading
to the same result only overcount. No reducedness or point deletion
interpretation is used.

This direction is important: it bounds the larger algebra from a
smaller quotient with controlled borders. It is not a monotonicity
claim for fan counts under arbitrary quotienting or deformation.
The containment is in the same polynomial ring with the same coordinate
images. An arbitrary isomorphism between abstract algebras is insufficient.

For a fractional-linear quotient in (MR3) of length n-e>=3, each
initial staircase is a union of coordinate chains. If r axes are active,
its minimal generators are the d inactive-or-terminal axis monomials
and the binom(r,2) mixed quadratic monomials. Thus q<=d*(d+1)/2,
and (MR4) is at most binom(n-e+d-2,d-1). Consequently

    g(I)<=binom(n-e+d-2,d-1)
          * product_(j=0)^(e-1)(d*(d+1)/2+(d-1)*j).       (EX4)

The factor is exp(O_d(e*log(e+2))). Fixed e gives O_(d,e)(n^(d-1)),
below the existing generic logarithmic order. If
e=o(log n/log log(n+2)), the bound is n^(d-1+o(1)), so this regime
cannot produce n^(d-1+epsilon) for any fixed positive epsilon. The
latter statement does not exclude logarithmic gains when e grows.

Combining (EX1)--(EX2) with (EX4) controls the branch where every
retained primitive rational map has degree one. If some map has degree
two or more, (EX4) has not been justified. It is not legitimate to
call every small-excess presentation a bounded-error fractional-linear
one. Likewise, (EX3) alone gives no unrestricted exponent improvement:
general smaller quotients need not have uniformly bounded borders or
a known smaller fan count.

#### A sharp extension of the bilinear test to large coordinate degrees

There is also an all-ideal screening statement that removes the initial
assumption of full coordinate degrees. Write s_i=dim_K K[x_i] inside
an arbitrary length-n quotient A generated by its d coordinates. If

    s_1+s_i>n+2 and dim span(1,x_1,x_i,x_1*x_i)<=3
      for every i>=2,                                    (EX5)

then A=K[x_1], every s_i=n, and the fractional-linear classification
applies. Indeed multiplication by a nonzero affine function of x_i
has rank at least s_i-1 by restriction to K[x_i]. A degenerate bilinear
relation would factor into affine functions and have positive product
rank at least s_1+s_i-n-2, impossible for a zero product. Any remaining
bilinear relation is nondegenerate, and Bezout gives x_i=N(x_1)/D(x_1)
with D a unit in A. Its inverse already belongs to the finite subalgebra
K[x_1]: its minimal polynomial has nonzero constant, so expresses that
inverse polynomially. Thus all coordinates lie in K[x_1], which must
be all of A. Each nonconstant fractional-linear coordinate also generates
it. For n>=3, its bilinear rank is then exactly three by (MR2).

Strictness cannot be dropped. For m>=3 take

    A=K[x,y]/(x*y,x^m,y^m,x^(m-1)-y^(m-1)), z=y+1.

This has length n=2m-2, all three coordinate degrees equal m, and
both bilinear ranks against x equal three. But 2m=n+2 and none of
these coordinates generates A. The connected-sum presentation is a
standard small algebraic counterexample, not a new lower family.
The larger-box inverse proof for full cyclic degrees is not transferred
automatically to partial degrees: its primitive constant component can
be supported in an affine-factor annihilator.

#### Exact verification and research decision

The checker saves 98 reconstruction fixtures over Q,F_2,F_3,F_37,
including nilpotent and split point algebras. It computes full relation
kernels independently, removes their polynomial gcds, checks exact
polynomial proportionality, and compares the exceptional multiplication
rank with the gcd quotient dimension. The annihilator multiples reproduce
the complete truncated kernel. Three degree-threshold counters retain
their nonzero eliminants, and nine connected-sum controls check the
sharp high-degree boundary with actual commuting cyclic presentations.

There are also 24 finite extension controls over Q,F_3,F_11,F_37.
The smaller fractional-linear algebras have lengths n-e with n=5,6,7
and e=1,2; the larger algebras are nilpotent extensions or reduced
one-outlier configurations. All 1,264 staircase rank decisions are
saved, with the staircase tables shared across fixtures. All independent
staircases in these particular controls lie in the enumerated extension
candidates. This stronger fixture fact is not asserted for general
extensions. At 323 strict integer weights the checker independently
selects both quotient bases, verifies containment, and checks that
exactly e monomials were added. It does not report an exact full fan
count for the larger algebra from these selected weights.

These results close bounded-exception fractional-linear amplification
and identify rational degree versus exceptional dimension as separate
parameters. They do not settle the degree-two-and-higher branch, a
growing-exception family, or the general fast-growth alternative. A
uniform dependency-sensitive upper bound still has to control these
branches; no unrestricted asymptotic advance, completed strongest-priority
audit, further consequence at that target, or new manuscript is claimed.
Keep the full goal active. Reproduce all controls with

```bash
python3 -B scripts/run_structural_complexity_checks.py
```

### Twentieth feasibility pass: a characteristic-free multiple-pole bound

The last published mathematical pass, commit 95dcbc3, separated rational
degree from exceptional dimension. The intervening ML assessment changed
no theorem or authoritative artifact and is classified as no progress
toward the unrestricted goal. This pass returns to the degree-two branch.
It proves a bound for a whole multiple-pole family, including arbitrary
nonreduced divisors over every field. The proof removes an apparent small
characteristic exception by retaining the least admissible gap, rather
than assuming gap one is independent. This is a restricted structural
component, not the unrestricted advance required by the goal.

#### A later counting comparison and the remaining baseline audit

[Jonsson Kling--Lundqvist--Mohammadi--Orth--Saenz-de-Cabezon,
arXiv:2411.10209v1](https://arxiv.org/html/2411.10209v1#S2.SS4)
studies the characteristic-zero ideals

    (x_1^2,...,x_r^2,(x_1+...+x_r)^k).

The setup, Definition 2.26, Lemma 2.27 with its proof, Proposition 2.29
with its proof, and Remark 2.30 were inspected. Their number of distinct
reduced bases is the multinomial coefficient with blocks k,2,...,2 and
a final singleton when needed. Variable rankings determine these bases.
Here r is the number of variables, not colength. With fixed d this family
has bounded length (at most 2^d), so its growing-r count supplies neither
an unrestricted fixed-d upper bound nor an infinite fixed-d lower family.
No full audit of the paper's main basis construction or its later
Lefschetz proofs is claimed. This is a precise comparison with a named
later counting result, not evidence that all relevant literature has been
exhausted.

A fresh NSF-deposit retrieval again failed. The already inspected
finite-field Theorem 8 and proof remain the comparison recorded in the
sixteenth pass: when p>n that expression recovers the old exponent.
The failed retrieval is not a new mathematical blocker and is not used
as evidence of novelty. The broader strongest-current-bound audit remains
incomplete; the inspected unrestricted benchmark is still Babson--Onn--
Thomas's O_d(n^(2d(d-1)/(d+1))).

#### The family, its length, and its complete parent normal domains

Let K be any field, d>=2, and choose distinct nonzero elements
lambda_3,...,lambda_d in K when d>2. In K(t) set

    x=t,   y=t^(-2),   z_j=(t-lambda_j)^(-1), j=3,...,d.

Let f be monic of degree n>=1 with
f(0)*product_j f(lambda_j) nonzero, and let I_f be the kernel of these
coordinate images in A=K[t]/(f). The quotient has length n since x=t.
It can be nonreduced or nonsplit. In fact these are all finite quotients
of this coordinate ring: x*y=t^(-1), so the ring is the localization
K[t,t^(-1),(t-lambda_3)^(-1),...,(t-lambda_d)^(-1)], a localization of
the principal ideal domain K[t]. Its finite-colength ideals are given
by the stated divisors f after removing inverted factors.

Call M={x,z_3,...,z_d} the simple-pole axes. At a generic positive weight,
choose its unique cheapest member u. Relations between distinct members
of M have their mixed products leading, by partial fractions. Also

    x^2*y=1,
    lambda_j^2*y*z_j+x*y+lambda_j*y-z_j=0.                (QP1)

If u=x, (QP1) makes y*z_j leading for every j. If u=z_s, its equation
makes x*y leading; subtracting its equation from the others makes every
y*z_j with j!=s leading. Finally

    y*(1+lambda_s*z_s)^2=z_s^2

makes u^2*y leading, including in characteristic two. The complete
parent normal domain therefore consists of pure chains on M, together
with y^k and u*y^k for k>=1. These monomials are independent in K(t):
at t=0 the pairs y^k,u*y^k separate even and odd pole orders, since
u=t or u has a nonzero first derivative at zero. After these coefficients
vanish, uniqueness of partial fractions separates the pure chains.
The displayed reductions span the ring and decrease the chosen order,
so independence proves completeness of the parent domains, not merely
an inclusion in them. Every global initial ideal of I_f has a positive
weight witness avoiding the finitely many ties between axes in M.

Consequently each length-n initial staircase has parameters

    a>=1, b0>=b1>=0, c_v>=0 for v in M minus {u},
    n=a+b0+b1+sum_v c_v,
    a>=2 if b1>0.                                        (QP2)

It contains 1,u,...,u^(a-1), y,...,y^b0,
u*y,...,u*y^b1 and v,...,v^c_v for the other axes v.
Put L=b0-b1 and h=(b0+b1)/2. Its three relevant moments are

    mu_u=a*(a-1)/2+h-L/2,
    mu_y=h^2+h+L^2/4,
    mu_v=c_v*(c_v+1)/2.                                  (QP3)

The cases a=1,b1=0 contribute only O_d(n^(d-2)) patterns and can be
counted separately. All arguments below concern a>=2.

#### Rank flags which retain the missing dependency information

If u=x, hold the other chain lengths c_j fixed and set
C(t)=product_j(t-lambda_j)^c_j. The paired y and x*y chains, the pure
x chain and the other poles span all rational functions with respective
pole bounds 2b1,a-1,c_j. This space has dimension n-L. Multiplication
by the unit t^(2b0)*C(t) identifies the entire candidate span with

    t^(2L)*K[t]_(<n-L)
       + span(t^(2k)*C(t):0<=k<L), modulo f.             (QP4)

This follows either by partial fractions or by multiplying those full
pole-bound spaces by their common denominator. In particular its rank
depends only on L and the c_j, not on h or on a separately.

If u=z_s, write c0 for the length on x and c_j for the other finite
simple poles. Set v=1/t and

    g(v)=v^n*f(1/v)/f(0),
    C_s(v)=(1-lambda_s*v)^(a-1)
                  *product_(j!=s)(1-lambda_j*v)^c_j.

Multiplication by the unit v^c0*C_s(v) identifies the full pole-bound
core with K[v]_(<n-L). The L remaining columns are

    v^(n-a-sum_(j!=s)c_j-L+2k)*C_s(v), 1<=k<=L,
    modulo g.                                           (QP5)

Thus the rank depends only on (a,(c_j)_(j!=s),L), not on c0 and h
separately. Projecting the columns in (QP5) onto the top L coefficients
after reduction modulo g gives an L-by-L rank test; total rank is
n-L plus that rank. Both (QP4) and (QP5) are identities of subspaces
after multiplication by units, so they hold over every field, with
arbitrary multiplicities in f. A collection of independent partial
fractions before reduction alone would not prove them modulo f.

For each fixed flag (a,(c_j)_(j!=s)) and each parity, let L_* be its
smallest admissible L giving rank n, if one exists. Even parity always
allows L_*=0. Odd parity need not allow L=1. For instance over F_2,

    f(t)=t^5+t^4+t^3+t+1, lambda_s=1, a=2

has no basis at odd gap one, whereas gap three is a basis. This exact
countercontrol prevents replacing an arbitrary admissibility set by
the condition on its first minor. No bound on L_* is assumed below.

#### A separable monotonicity count in the x sector

Fix c=(c_j) and one parity of L. Put N=n-sum c_j, so a=N-2h. By
(QP4), feasibility is independent of h, except for the inequalities

    L/2<=h<=(N-2)/2.

The pairs (h,L) are closed under coordinatewise minimum and maximum,
even if admissible L values have arbitrary gaps. Normalize w_y=1 and
put r=w_x. The cost in (QP3) is a sum of a function of h and a
function of L; mu_x strictly decreases when either h or L increases
on this parity lattice. If p minimizes at r1 and q at r2>r1, write
m=p meet q and j=p join q. Separability gives

    F_(r2)(j)-F_(r2)(q)
      = F_(r1)(p)-F_(r1)(m)
             +(r2-r1)*(mu_x(p)-mu_x(m)).                (QP6)

If p is not coordinatewise at most q, the last term is strictly
negative and the first is nonpositive, contradicting optimality of q.
Strict coherent minimizers therefore form a coordinatewise monotone
chain as r increases. Each coordinate has O(N) possible values, so
there are O(N) such pairs for fixed c and parity. Summing over the
O_d(n^(d-2)) other chain vectors proves O_d(n^(d-1)) initials in this
sector. Costs of the fixed other axes do not affect this argument.

#### The least-admissible-gap count in every finite-pole sector

Fix u=z_s, a>=2, the c_j at j!=s, and a parity with some admissible
gap. Normalize w_y=1, write s=w_u and r=w_x. The sector inequality
is r>s>0. Take any coherent pattern with gap L and put Delta=L-L_*.
Replacing L by L_* at the same h,c0 gives another basis by (QP5),
with integral b1'=h-L_*/2>=0. If Delta>0, strict cost optimality and
(QP3) imply

    0>F(L)-F(L_*)
      =(L-L_*)*((L+L_*)/4-s/2),
    hence L+L_*<2s.                                     (QP7)

If c0>=2, increasing both b0,b1 by one and decreasing c0 by two
also preserves the rank flag and gives another basis. Its cost minus
the original cost is 2h+2+s-(2c0-1)r>0. Since s<r,

    (c0-1)*r<h+1.                                       (QP8)

For c0>=2 and Delta>=2, combine (QP7)--(QP8) to get

    (c0-1)*(Delta-1)<2*(h+1)
       =n-a-sum_(j!=s)c_j-c0+2<n.                       (QP9)

The same product with each factor replaced by its nonnegative part
is zero in the other cases. For fixed a,(c_j) and parity, c0 and
Delta determine L,h,b0,b1. The number of positive integer pairs
R,S with R*S<n is at most

    sum_(R=1)^(n-1) floor((n-1)/R)=O(n*log(n+1)).

The zero-factor boundary contributes O(n) more pairs. There are only
O_d(n^(d-2)) choices of a,(c_j); parity and the choice of u cost
constant factors. Combining this with the x sector proves

    g(I_f)=O_d(n^(d-1)*log(n+1)), d>=3,
    g(I_f)=O(n), d=2.                                   (QP10)

The constants are independent of f, the poles and the characteristic.
In particular there is no extra n^d/p term from persistent small-field
dependencies. The useful mechanism is comparison to the least feasible
flag, with its potentially large offset removed before counting integer
pairs. The proof does not require a generic divisor, a nonsingular rank
minor at gap one, or an equality of fans for different divisors.

#### A linear rank-delay counterexample even at rational degree two

The rank flags cannot be replaced by a claim that every parent-normal
prefix completes rank after n+O(1) entries. Take d=2, n=2m+1,
f=t^n-1, x=t and y=t^(-2), over any field. Both coordinates generate
the quotient: x=y^m modulo f, and the y powers permute the t-power
basis. With w_y=1 and w_x=m+epsilon, 0<epsilon<1, the parent-normal
monomials of cost below 2m are exactly

    y^j, 0<=j<2m;       x*y^j, 0<=j<m.

There are 3m of them. The latter m columns repeat y^m,...,y^(2m-1),
so their combined rank is n-1. The next monomial y^(2m) completes
rank at position 3m+1=n+m. The delay is linear in n despite fixed
rational degree and full coordinate elimination degrees. This is a
rank-completion obstruction, not a large fan or an unrestricted lower
bound. The matrices are standard unit columns in every characteristic,
so the argument does not rely on a numerical rank tolerance.

#### Verification, scope and next mechanism

The checker constructs twelve actual cyclic quotients over Q,F_2,F_3
in dimensions two through four. It compares direct residue ranks with
the independent flag formulas, and independently enumerates every
downset in each parent domain. Eight three-dimensional fixtures have
complete positive fans: every candidate has a saved rank, and every
independent candidate has exact border comparisons and a full polygon
clip; nonempty regions have strict integer weights and independent
cross-wide greedy checks. Four other fixtures have selected strict
weights, explicitly not a claimed full fan. Checks also cover the
monotone-chain and minimum-gap counting inequalities, and the three
rank-delay controls m=2,3,5. These are bounded proof controls; the
uniform asymptotic statement is established by (QP1)--(QP10).

Shared pattern tables store the chain parameters in (QP2), which reconstruct
each exact basis. Border normals and polygon vertices use comma-separated
exact coordinate strings. This avoids repeating full bases across models;
the checker validates the expanded bases before encoding them.

For this whole degree-two multiple-pole family, (QP10) lies below the
existing generic n^(d-1)*log^(d-1)(n) lower order when d>=3. Park it
as an exponent-amplification source, including in small characteristic.
No novelty for this restricted bound has been certified. The retained
counting mechanism suggests a next proof gate: for higher rational
degree, determine whether independent common shifts leave a bounded
collection of rank flags and whether least-admissible comparisons still
bound all coherent patterns. An arbitrary rational parametrization may
have several interacting pole gaps, and general finite algebras need
not have such a parametrization. A dependence on rational degree and a
complementary bound for the remaining algebras are still required for
an unrestricted result. Do not expand enumeration merely by increasing
n in the now-bounded family.

The original goal, further-consequence requirement, primary-priority
audit and honest new-manuscript gate remain open. Reproduce with

```bash
./build.sh
python3 tests/test_smoke.py
python3 -B scripts/run_structural_complexity_checks.py
```

### Twenty-first feasibility pass: the staircase-moment relaxation is sharp

The previous turn made progress in commit 13f37ec: it proved the all-field
quadratic-pole bound, saved exact controls and recorded the user's broader
ML architecture direction. It did not improve the unrestricted initial-ideal
bound. This pass tests whether individual staircase geometry can improve
that bound without additional algebraic information. An explicit moment
cube shows that this relaxation already has the full old exponent, in
every fixed dimension. The remaining compatibility gap is made precise;
no family of ideals attaining that exponent is constructed.

#### The relaxation and the exact scope of the result

Let O_d(n) be all n-element downsets in N_0^d and write
mu(D)=sum_(a in D) a. For a family F contained in O_d(n), call D in F
positively exposed in F if some w in R_(>0)^d satisfies

    w.mu(D)<w.mu(E) for every E in F minus {D}.

Let H_d(n) be the largest possible number of such D over all choices of
F. Different members with the same moment cannot both be strictly exposed.
If F is the initial-staircase family of one colength-n ideal, every member
is positively exposed. Hence

    F_d(n)<=G_d(n)<=H_d(n).

The converse is not assumed: F in the definition of H need not be the
basis family, or initial family, of any common quotient or even a common
linear matroid. The result proved below is

    H_d(n)=Theta_d(n^(2d(d-1)/(d+1))), d>=2.              (MCB1)

This determines a geometric relaxation of the target, not G_d or F_d.
In three dimensions it gives H_3(n)=Theta(n^3), while the established
gap for actual ideals remains Omega(n^2*log^2 n)<=G_3(n)<=O(n^3).
Thus a uniformly smaller bound cannot follow from downward closure,
cardinality, positive exposure and moment geometry alone.

#### An arithmetic interval of moments at fixed planar cardinality

Fix an integer R>=1. Put M=4R+1 and B=12M+2R+13=50R+25. For every
s in {0,...,R^2}, use Lagrange's classical four-square theorem to choose
nonnegative integers z_0,...,z_3 with sum z_j^2=s. Each z_j is at most R.
Make a vector v on rows i=0,...,12 by setting v_0=0 and, for j=0,...,3,

    (v_(3j+1),v_(3j+2),v_(3j+3))=(z_j,-2z_j,z_j).

Its four supports are disjoint, and direct summation gives

    sum_i v_i=0, sum_i i*v_i=0, sum_i v_i^2=6s.          (MCB2)

Set ell_i=B-M*i+v_i. Since |v_i|<=2R, these are strictly decreasing
positive row lengths with ell_i>=13. The first row ell_0=B is fixed.
Let A_s be the planar downset with row i containing x=0,...,ell_i-1.
Then every A_s has the same size

    S=sum_i ell_i=13B-78M=338R+247.

Its second-coordinate moment is independent of s by (MCB2). Its first
moment is sum_i binom(ell_i,2). Expanding this expression, the linear
terms in v vanish because the baseline lengths are affine in i, while
one half of sum v_i^2 is 3s. Thus, with constants U,V depending on R,

    |A_s|=S, mu(A_s)=(U+3s,V).                          (MCB3)

This realizes every third integer in an interval of length 3R^2 while
preserving size and the other moment. The four-square input is established
number theory, not a new arithmetic consequence here. One accessible
statement is [Steven V Sam, Theorem 1](https://mathweb.ucsd.edu/~ssam/notes/lagrange-4sq.pdf);
only that classical theorem is needed for this step.

#### Independent arms give a full moment cube in fixed dimension

For d>=3 place one copy of A_(s_i) in each cyclic coordinate plane
(i,i+1 mod d), taking i as its long coordinate, and let D_s be their
union. Different arms meet only on pure coordinate axes or at the origin.
On each pure axis the overlap is exactly the thirteen points 0,...,12,
independent of s, because every row length is at least thirteen. Counting
nonorigin pure-axis points and the disjoint mixed-plane points gives

    |D_s|=d*(S-13)+1=338dR+234d+1.

Coordinate i gets the long moment U+3s_i from its own arm and the short
moment V from the preceding arm. The duplicated axis contributes
sum_(j=0)^12 j=78. Therefore

    mu(D_s)=(U+V-78)*1+3s, s in {0,...,R^2}^d.          (MCB4)

For d=2 the two oppositely oriented planar arms overlap in the fixed
thirteen-by-thirteen square. The same calculation gives

    |D_s|=2S-169=676R+325,
    mu(D_s)=(U+V-1014)*1+3s.                            (MCB5)

These unions are genuine downsets. They are not being called initial
staircases of a common ideal.

The fixed first row makes exact-cardinality padding possible. Extend the
first pure axis by delta points starting at B. This adds delta to size
and delta*B+binom(delta,2) to the first moment, independently of s. It
preserves downward closure and (MCB4)--(MCB5) up to a fixed translation.
Define (a_d,b_d)=(676,325) for d=2 and
(a_d,b_d)=(338d,234d+1) for d>=3. For any n>=a_d+b_d, choose

    R=floor((n-b_d)/a_d)>=1,
    delta=n-a_d*R-b_d, 0<=delta<a_d.

It follows that for every sufficiently large n, moments of actual
n-element downsets contain a translated grid

    c_n+3*{0,...,R^2}^d, R=Theta_d(n).                  (MCB6)

All points of the grid are realized, not just its corners, and the padding
handles every large n rather than an infinite subsequence alone.

#### The classical lattice bound is already sharp on these moment sets

[Barany--Larman (1998), Theorem 1](https://www.renyi.hu/~barany/cikkek/73.pdf)
states that the integer hull P_r=conv(Z^d intersect r*B^d) has
Theta_d(r^(d(d-1)/(d+1))) vertices. The theorem, the k=0 deduction in
Section 3 from Theorem 4, and Section 5's full proof of that volume
estimate were inspected. Their upper theorem and the usual volume estimate
are established inputs, not a new lattice-polytope theorem in this project.
Proofs for the other face numbers and the later smooth-body extension were
not audited. This does not complete the strongest-current-initial-ideal
priority comparison.

Take r=floor(R^2/2). After translation by r*1, all lattice points of this
ball fit in the grid in (MCB6). Its vertices therefore give distinct actual
n-staircases with that moment configuration, after scaling by three.
At least a 2^(-d) fraction of these vertices have a strict positive
minimizing normal. Indeed every vertex has a strict minimizing normal
with all coordinates nonzero. Coordinate reflections preserve P_r and
can make that normal positive; any resulting vertex has at most 2^d
preimages under these reflections. Choosing the ball-indexed family of
staircases gives

    H_d(n)>=c_d*R^(2d(d-1)/(d+1))
           >=c'_d*n^(2d(d-1)/(d+1)).                   (MCB7)

For the upper bound, every coordinate of every n-staircase moment is
between zero and binom(n,2). The established lattice-polytope vertex bound,
with the unit-cube device of (MB2) for lower-dimensional hulls, gives the
matching upper estimate for any subfamily. This proves (MCB1).

The lower estimate is inherited from the classical integer-ball theorem
after an explicit staircase encoding. It does not establish a new
initial-ideal asymptotic. The missing common-algebra realization is exactly
the distinction between H and G; it cannot be discharged by convexity or
by the fact that each individual downset defines a monomial ideal.

There is an established stronger necessary condition even before one uses
commutation. [Babson--Onn--Thomas, Proposition 2.2 and its proof](https://arxiv.org/pdf/math/0207135)
show that the monomials of an initial basis, ordered by its monomial order,
are no later at each position than the ordered monomials of any other
basis. For a positive weight, this implies componentwise dominance of
the sorted lists of costs. Rechecking the full displayed proof gives a
simple usable filter: below any cost threshold, an initial basis must
contain at least as many entries as any independent competing prefix.
Total moment minimization alone does not supply this condition for an
arbitrary family. It follows from the actual residue-rank structure.

For example, take d=3,R=2,delta=7, so n=2738, and index the encoded
family by lattice points p in the radius-two ball, using s=p+(2,2,2).
At weight (4,1,1), p=(-2,0,0) is a strict total-moment minimizer in
that family. But compared with q=(-1,-1,-1), its 1262nd monomial has
cost 61 while the other's has cost 60. If both encoded staircases are
bases of one quotient, the first cannot be initial at that weight or
under any order refining it. The checker supplies analogous witnesses
for all four positive exposed vertices of this small ball. These are
weight-specific disqualifications, not a proof that the full family
cannot be realized at some other collection of weights.

#### Exact controls and the next algebraic question

The checker supplies four-square witnesses and directly expands the planar
partitions for R=1,...,6 and every s<=R^2. It verifies the row constraints,
cardinalities and both moments without using the closed moment formula to
compute the answer. Cyclic unions in dimensions two through five are
expanded as point sets; checks include downward closure, exact moments,
and padding at zero, one and the end of each cardinality period.

For d=3, R=2,3, it also enumerates all lattice points of balls of radii
two and four. Exact polygon clipping finds every positive exposed ball
vertex, then gives a strict integer normal and checks the corresponding
expanded staircase against the entire ball-indexed moment family. These
certificates concern exposure in a relaxation, with no assertion of a
common quotient or residue-rank certificate. The general asymptotic
argument uses (MCB2)--(MCB7), not extrapolation from the finite ball counts.
Four exact sorted-prefix countercontrols separately check the stronger
necessary comparison, including its strict total-cost and prefix-count
inequalities.

Park moment-only improvements to the unrestricted exponent, even if they
retain the complete set of possible individual staircase moments. A viable
upper argument must use compatibility among bases of one multiplication
representation. A viable lower argument must construct one common algebra
and prove its rank and coherence conditions; the independently encoded
moment cube is insufficient. The next discriminating gate is a quantitative
restriction on rank-feasible staircase families from commuting cyclic
multiplication matrices, retaining the ordered-prefix dominance required
of a greedy basis, beyond the properties used in H. Arbitrary matroids
must not be substituted for such quotient representations. The newly located
Hilbert-function curve-carrier results are only a possible comparison to
inspect, not a theorem or reduction being invoked here.

The original unrestricted goal, further consequence, primary-priority audit
and separate manuscript requirements remain open. No publication novelty is
certified for this relaxation result. Reproduce all controls with

```bash
./build.sh
python3 tests/test_smoke.py
python3 -B scripts/run_structural_complexity_checks.py
```

### ML follow-up: an exact reference for a compressed attention hot path

The user asked specifically whether the factorizations can improve existing
ML hot paths. The most concrete candidate is the feature expansion and
streaming-state update in polynomial attention. This is an architectural
constraint with an exact implementation benefit. It is not evidence that
our initial-ideal counting results accelerate ordinary softmax, and its
underlying kernel factorization is already established mathematics.

#### Several latent coordinates instead of a single moment curve

Let s>=1, m>=1, A={0,...,m}^s minus {0}, and D=|A|=(m+1)^s-1.
Use learned latent query and key encoders q,k in R^s and positive parameters
theta_a. Design the score and the positive attention kernel as

    S(q,k)=sum_(a in A) theta_a*q^a*k^a,
    K(q,k)=1+S(q,k)+S(q,k)^2/2.

Equivalently the D-dimensional query/key embeddings have entries
sqrt(theta_a)*q^a and sqrt(theta_a)*k^a. The product relations are built
into these embeddings. They do not hold for unrestricted D-dimensional
queries and keys. Collect equal monomials in the quadratic expansion:

    K(q,k)=sum_(b in B) c_b*q^b*k^b,
    B={0} union A union (A+A)={0,...,2m}^s,
    c_b=1[b=0]+1[b in A]*theta_b
        +(1/2)*sum_(a,a' in A, a+a'=b) theta_a*theta_(a').       (MLH1)

All coefficients are positive. The identity is elementary collection of
monomials, an instance of the already audited
[Altschuler--Parrilo polynomial-kernel framework](https://arxiv.org/pdf/2106.02755v1).
The M=(2m+1)^s features are also the minimal exact separation width for
this particular kernel on nonempty open latent domains: the M monomials
are linearly independent on each domain and the diagonal coefficient
matrix is nonsingular. Choosing unisolvent evaluation sets on both sides
gives a rank-M kernel matrix. This standard rank argument is not being
claimed as a new tensor-rank result.

For value vectors v_j of width v, maintain for each b a denominator
statistic Z_b=sum_j k_j^b and a value statistic H_b=sum_j k_j^b*v_j.
At query q the numerator is sum_b c_b*q^b*H_b and the denominator is
sum_b c_b*q^b*Z_b. Restrict the sums to j<=i for causal attention. This
computes exactly the designed quadratic kernel, including normalization.
Since K=(S+1)^2/2+1/2, every nonempty denominator is at least T_visible/2.
That lower bound does not establish good conditioning of the monomial
basis or of its gradients.

Counting the latent encoders separately, straightforward monomial tables
give O(T*M*(s+v+1)) forward arithmetic and O(M*(v+1)) streaming state.
Direct coefficient preparation costs O(D^2*s), charged whenever theta
changes. Backpropagation storage/recomputation needs a separate analysis;
the forward state bound is not a constant-memory-training claim.
For fixed s, M=Theta_s(D), compared with 1+D+D(D+1)/2 generic symmetric
quadratic features. Concrete feature/state-row counts are

| Latent coordinates s | Box degree m | Embedding D | Generic symmetric features | Collected features |
| --- | --- | --- | --- | --- |
| 1 | 64 | 64 | 2145 | 129 |
| 2 | 7 | 63 | 2080 | 225 |
| 3 | 3 | 63 | 2080 | 343 |

These are exact representation and streaming-state counts for the same
constrained embeddings. They are not measured wall-clock speedups, and
they are not savings against a method already using the collected kernel.
Larger s costs more: (2m+1)^s is exponential in latent dimension. Multiple
heads and nonlinear encoders can be tested, but cannot be assumed to
remove this interaction-class restriction. Arbitrary learned mixing of
the expanded coordinates generally destroys the diagonal form in (MLH1).

#### Implementation and prior-art boundary

The existing research checker now includes a dependency-free rational
reference for five one-, two- and three-dimensional fixtures. Its direct
oracle forms S and K from the uncollected score, independently of (MLH1).
The streaming implementation collects exponents and maintains the finite
state. Exact controls cover pair kernels; query, key and tied-weight
kernel derivatives; causal and noncausal normalized outputs; and normalized
query/key derivatives, including a key masked out of early causal prefixes.
This validates the algebraic reference, not an optimized training kernel.

[TaylorShift, Sections 3.1--3.3 and 4](https://arxiv.org/html/2403.02920v1#S3)
already expands quadratic attention features, changes contraction order,
and analyses normalization and crossover costs. Its efficient expression
uses the full tensor product; the symmetric feature count in the table
already improves that elementary redundancy before applying (MLH1).
[Linear attention, Sections 3.2--3.3](https://proceedings.mlr.press/v119/katharopoulos20a/katharopoulos20a.pdf)
already supplies the causal sufficient-statistic construction. Thus neither
linear context scaling nor a finite recurrent state is a novelty claim.

Other credible hot paths are structured token/channel mixing and recurrent
state transitions. The project's separable exponential coefficients explain
how a full-matrix-rank operator can still have a cheap structured application;
they do not imply low matrix rank. Existing
[Monarch matrices](https://arxiv.org/abs/2204.00595),
[structured state-space duality](https://arxiv.org/abs/2405.21060), and
[Tensor Product Attention](https://arxiv.org/abs/2501.06425) are mandatory
comparisons before presenting a new factorized block. The prior inspection
scopes in the goal manifest still apply; this follow-up rechecked their
primary abstracts and TaylorShift's displayed feature formulas, not all
their proofs or experiments.

#### The experiment that would establish an ML contribution

Prioritize (MLH1) as a small architectural experiment within an existing
polynomial-attention block. Compare its actual implementation with the
same constrained kernel through generic symmetric features to measure the
implementation benefit. Separately compare learned models against ordinary
polynomial/linear attention at matched state, parameters and training cost;
include associative recall to expose losses in selective retrieval. Charge
encoders, coefficient updates, forward/backward computation and state traffic.
Use a numerically suitable basis or bounded latents and measure both output
and gradient error against a high-precision reference. No GPU or model-quality
experiment has yet been run. Do not report the feature-count ratios as
throughput or quality improvements.

A worthwhile new theorem would have to add a stable approximation or
representation guarantee for a useful interaction class at a lower complete
cost than the strongest applicable prior method. A practical paper could
instead demonstrate a useful accuracy/memory/throughput tradeoff with a
proper controlled comparison. The current exact factorization alone does
not establish either result. This ML application work preserves the active
unrestricted mathematics goal rather than silently substituting for it.

### Twenty-second feasibility pass: basis exchange completes the missing moments

The previous turn made progress in commit 45e9061: it established the sharp
staircase-moment relaxation and saved weight-specific prefix failures. That
result did not improve the unrestricted bound. This pass replaces those
individual failures by an all-weight obstruction for the radius-two-ball
template, valid uniformly as its row lengths grow. The mechanism is a
general basis-exchange certificate. It uses established matroid theory and
does not constitute the asymptotic advance required by the goal.

#### A rank-two exchange and a coordinate-cycle certificate

Let D be a monomial basis of any finite quotient, over any field. Suppose
u,u+h*e_i are in D, h>0, and a,b are not in D, with

    E=(D minus {u,u+h*e_i}) union {a,b}

also a basis. Contract the span of the n-2 common monomials. The images
of the old pair and of the new pair are both bases of a two-dimensional
space. In particular a has a nonzero image, so at least one of

    D minus {u} union {a},
    D minus {u+h*e_i} union {a}                           (BC1)

is a basis. If a-u=-alpha*e_i+beta*e_j with alpha,beta>0, the moment
change of that forced basis is at most

    -alpha*e_i+beta*e_j                                  (BC2)

coordinatewise; removing the second old monomial subtracts h*e_i too.
Consequently, initiality of D at a positive weight requires

    beta*w_j>=alpha*w_i.                                (BC3)

The argument is independent of all multiplication coefficients. It uses
the linear matroid of one common quotient, not unrelated rank-two spaces.

Suppose such exchanges are available around a directed cycle of distinct
coordinate indices i=0,...,ell-1, with j=i+1 modulo ell. If

    product_i alpha_i > product_i beta_i,                (BC4)

then the necessary inequalities (BC3) are inconsistent. More strongly,
there is a positive rational convex combination of the forced actual
bases whose moment is strictly below mu(D) on every cycle coordinate
and equal to it on the others. This proves noninitiality without choosing
a weight. Here is an explicit certificate, avoiding logarithms. Put
L_i=alpha_i/beta_i and choose a rational t>1 with

    t^(ell-1)<product_i L_i.

Set q_0=1 and q_i=t*q_(i-1)/L_i for i=1,...,ell-1. Take convex
coefficients lambda_i proportional to q_i/beta_i. The upper bound on
coordinate i of the combined moment change, before normalization, is

    q_(i-1)-L_i*q_i=(1-t)*q_(i-1)<0, i>=1,
    q_(ell-1)-L_0*q_0<0, i=0.                           (BC5)

Extra old-entry removals in (BC1) only decrease the corresponding
coordinate further. Each term is a genuine basis if the hypotheses
hold, so this convex combination belongs to the common basis polytope.
It has lower value than mu(D) under every strictly positive weight.
Every initial staircase has a realizing open positive weight cone; D
therefore cannot be initial under any global monomial order.

This is an explicit use of established basis exchange and state-polyhedron
theory. [Babson--Onn--Thomas](https://arxiv.org/pdf/math/0207135) already
uses basis exchange in Proposition 3.4 and represents the basis polytope
as a projected matroid polytope in Theorem 3.5. Those displayed proofs,
together with Proposition 2.2, were rechecked. No novelty is claimed for
basis exchange, gain inequalities, or the general polyhedral framework.
The certificate is a sufficient obstruction, not a characterization of
initiality or an upper bound on the number of initial ideals.

#### The exact interval forced by a balanced row patch

Consider three consecutive rows r,r+1,r+2 of an affine baseline of
slope M. A patch (z,-2z,z) changes to (z+1,-2z-2,z+1). Write L for
the unperturbed length of the first row, so their old lengths are

    L+z, L-M-2z, L-2M+z.

The two old monomials removed in the exchange have exponents

    u=(L-M-2z-2,r+1), u+e_i,

while the added exponents are

    a=(L-2M+z,r+2), b=(L+z,r).

Here i is the long direction and j the row direction. Thus

    a-u=-(M-3z-2)*e_i+e_j,
    b-u=(M+3z+2)*e_i-e_j.                               (BC6)

Common monomials cancel from every cumulative cost comparison. The
two old costs are w.u and w.u+w_i. For their sorted list to dominate
the two added costs, its smaller cost must be at most both added
costs. This is exactly

    M-3z-2 <= w_j/w_i <= M+3z+2.                        (BC7)

The larger old cost imposes no additional condition: the two added
costs, relative to w.u and divided by w_i, sum to 6z+4, so their
maximum is at least 3z+2>=2>1. This proves the equivalence in (BC7)
for prefix dominance against this one competing basis. It does not
make that condition sufficient against other bases. In particular,
competing bases with a larger moment can still force a smaller
one-for-one basis through (BC1); moment comparisons among the listed
candidates alone miss that completion.

#### Uniform exclusion of the radius-two-ball template

Return to the thirteen-row construction, with any R>=2, M=4R+1 and
any padding delta>=0. Only the profiles A_0,A_1,A_2 are needed here.
The canonical four-square witnesses give A_1 a unit patch on rows
10,11,12, and A_2 adds a second unit patch on rows 7,8,9. Each of
A_0 to A_1 and A_1 to A_2 is therefore the z=0 exchange above.

Index all 33 points p in the integer ball ||p||_2<=2 by s=p+(2,2,2),
and use the earlier canonical profile for each s-coordinate. Suppose
all these staircases are bases of one common quotient. For positive
weights, the only possible strict moment minimizers in this family are

    s=(0,2,2),(2,0,2),(2,2,0), and (1,1,1).              (BC8)

Indeed a lattice point of the radius-two ball either has one coordinate
of absolute value two and the others zero, or has all coordinates of
absolute value at most one. Under a positive minimizing normal, the
latter are dominated by p=(-1,-1,-1). Strictness excludes ties and
nonvertices. An initial member must be one of (BC8), by the established
strict total-moment condition for a realizing generic weight.

For the first axial candidate, compare s=(0,2,2) with its inward neighbor
(1,2,2), also in the ball family. The exchange requires
w_2>=(M-2)*w_1. But strict axial moment minimality requires
w_1>w_2+w_3. Since M-2=4R-1>=7, these are incompatible. The other
two axial candidates fail cyclically.

For s=(1,1,1), all three competing bases obtained by increasing one
s-coordinate to two belong to the family. Their exchanges have
alpha=M-2 and beta=1. One actual basis is forced for each axis, with
moment change -((M-2)+epsilon_i)*e_i+e_(i+1), where epsilon_i is zero
or one. Averaging the three gives a point of the actual basis polytope
at least

    (M-3)/3=(4R-2)/3                                    (BC9)

below mu(D) in every coordinate. This rules out the central candidate
under every positive weight. Hence none of the 33 listed staircases
can be initial if they are all common quotient bases. The conclusion
holds uniformly for n=1014R+703+delta, not just at one saved weight.
It neither excludes other encodings nor proves failure of the larger
balls used in the asymptotic relaxation result.

The quotient can have other initial bases outside the listed family.
There is also no assertion that common-basis realization itself fails:
generic point evaluation can make all these downsets bases. What fails
is their proposed simultaneous role as initial bases. The exchange
completion produces missing basic sets with smaller moments, explaining
the gap left by independently exposed staircase moments.

#### Verification, literature scope and the remaining main problem

The checker expands the two planar partitions and subtracts them before
testing (BC7), including the interval endpoints and both sides. It checks
every invertible rank-two new pair over F_2,F_3,F_5,F_7,F_11 after
normalizing the old pair to a basis. Rational unequal-gain cycles check
(BC5) and every possible forced old-entry removal. Product-one and
product-less-than-one controls retain feasible weak weight ratios, without
claiming initiality. Finally the exact arm unions and padding are rebuilt
for R=2,...,6 and delta=0,7; all six required whole-staircase exchanges
and the uniform average drop are independently checked. Contraction
matrices are local implication controls, not a purported construction of
one algebra from independently chosen matrices.

The published finite-field chapter is now explicitly located at
[Stigler--Zhang, 2020, pp. 139--158](https://link.springer.com/chapter/10.1007/978-3-030-42687-3_9).
The publisher page and the existing indexed Theorem 8 comparison do not
replace the full chapter: the public deposit again failed direct retrieval.
No full publisher/preprint identity or exhaustive baseline audit is claimed.
[Migliore's primary survey](https://academicweb.nd.edu/~jmiglior/Hilb.pdf),
Theorems 5.3--5.4, was also inspected at the statement level: it works in
characteristic zero with reduced projective point schemes. Its maximal-growth
and reduction-number hypotheses must not be imposed on arbitrary finite
algebras. The original proofs of those carrier theorems have not
been audited here and no curve-carrier reduction is used.

This pass supplies a general exchange-closure screen and closes the small
ball template as a lower-bound construction. It supplies no new bound for
G_d or F_d. The next unrestricted mechanism must quantify which proposed
exchange alternatives have to be dependent, and use their common
multiplication representation to control all coherent bases. Merely
adding the cycle test to larger independent-staircase searches would not
satisfy that gate. The original asymptotic theorem, further consequence,
strongest-current baseline audit and separate manuscript remain open.
The goal remains active. Reproduce the exact controls with

```bash
./build.sh
python3 tests/test_smoke.py
python3 -B scripts/run_structural_complexity_checks.py
```
### ML implementation gate: learned mixing, state closure and honest comparisons

The current ML application priority is the feature expansion and recurrent
state update in a deliberately structured polynomial-attention block. The
existing exact example has 63 expanded query/key coordinates, from two
latent coordinates of box degree seven, and 225 collected features instead
of 2080 generic symmetric quadratic features. This is a 2080/225 ratio in
feature/state rows for the same constrained kernel. It is not a throughput
measurement or equivalence to unrestricted 63-dimensional queries and keys.

#### Dense learned mixing retains the state reduction but adds computation

The earlier diagonal restriction is sufficient for the cheapest algorithm,
but is not necessary for an M-row recurrent state. Let
A={0,...,m}^s minus {0}, D=|A|, B={0,...,2m}^s and M=|B|.
For any real D-by-D matrix H, including a nonsymmetric matrix obtained from
different learned query/key projections, put

    S(q,k)=sum_(a,b in A) H_(a,b)*q^a*k^b,
    K(q,k)=1+S(q,k)+S(q,k)^2/2.

Writing u(q)=(q^r)_(r in B), collection of the two sets of exponents gives

    K(q,k)=u(q)^T*C*u(k),
    C_(r,t)=1[r=t=0]+1[r,t in A]*H_(r,t)
      +(1/2)*sum_(a+c=r,b+d=t) H_(a,b)*H_(c,d).              (MLH2)

Thus store Z=sum_j u(k_j) and V=sum_j u(k_j)*v_j^T, using M*(v+1)
scalar entries. At each query first form u(q)^T*C, then contract it with
Z and V. The same formulas hold for causal prefixes. They compute exactly
the designed kernel, including normalization, since K>=1/2 pointwise.
That positivity does not imply that C is positive semidefinite.

The proof of (MLH2) is direct multiplication: linear terms have exponent
pair (a,b), while quadratic terms have pair (a+c,b+d). Every pair lies in
B-by-B. It is an instance of the already established coordinate-ring
kernel principle, not a new factorization theorem.
[Altschuler--Parrilo, Proposition 3.6 and Corollary 3.8](https://arxiv.org/html/2106.02755v1#S3.SS1.SSS3)

For a general dense H the straightforward costs are O(M*v) state
update/read, O(M^2) additional query transformation, O(M^2) coefficient
storage and O(D^4*s) naive coefficient preparation whenever H changes.
These are explicit upper bounds, not optimal preprocessing claims. Latent
encoders, backward computation and training storage remain additional.
The 225-row state survives mixing, but the dense C contains 50625 scalar
entries. For value width 64, the recurrent state alone shrinks from
135200 entries to 14625; storing a dense C as well gives 65250 entries
before other parameters/workspace. The 9.24-fold state-row reduction must
not be advertised as a 9.24-fold total-memory or speed improvement.
Positive diagonal H retains the earlier diagonal C and cheaper contraction.

The public checker independently forms the direct bilinear score and the
collected coefficient matrix for four small dense nonsymmetric examples.
It verifies 100 pair kernels, 400 query/key kernel derivatives and 80
normalized causal/noncausal output coordinates using rational arithmetic.
These are forward and kernel-derivative checks, not full backward or
learned-model validation.

#### Factorized features do not ensure a factorized accumulated state

For two latent coordinates, u(k)=u_1(k_1) tensor u_2(k_2) is separable.
However the state is a sum of these rank-one tensors. Take N distinct real
numbers t_j and keys k_j=(t_j,t_j). With powers 0,...,L-1 its denominator
state is

    Z_(r,c)=sum_(j=1)^N t_j^(r+c)=V*V^T,
    V_(r,j)=t_j^r.

The Vandermonde matrix V has rank min(N,L), and over the reals V*V^T
has the same rank. Thus the state can reach full matrix rank L, although
every token update has rank one. The checker saves exact positive leading
minors for ranks one through five. This classical witness rules out an
unjustified fixed-rank-state inference from the project's separable
coefficients. It does not preclude approximate state compression under
additional assumptions or nonlinear finite-precision memory schemes.

#### Practical decision

Begin with the diagonal structured block, because its representation
benefit has the clearest complete forward cost. Compare it with the same
kernel evaluated through generic symmetric features, then with a learned
225-feature linear/polynomial baseline and the relevant existing methods
at matched memory, parameters and training compute. Measure associative
recall as context grows before investing in a custom GPU kernel. A smaller
state that loses necessary retrieval capacity does not establish an
architectural advantage.

[TaylorShift, Sections 3--4](https://arxiv.org/html/2403.02920v1#S3)
already treats polynomial attention, normalization and crossover costs.
[Power Attention, Section 4.1](https://arxiv.org/html/2507.04239v1#S4.SS1)
specifically fuses feature expansion with matrix multiplication and retains
some redundant features to obtain regular GPU tiles. It demonstrates why
minimum mathematical feature width is insufficient as a hardware claim.
For structured MLP/channel mixers, also compare with
[Block Tensor-Train, Sections 2--3](https://arxiv.org/html/2406.06248v1#S2),
including its structure-dependent initialization and learning rates.
No new advantage over those architecture families has been established by
the project's current special factorizations.

An additional local NumPy diagnostic compared the diagonal generic and
collected implementations against an independent direct float64 kernel
and analytic query derivatives, using float32 feature/state arithmetic.
Its four bounded random fixtures each used 96 tokens, value width 16,
causal and noncausal evaluation, and both representations: 16 comparisons
in total. Maximum absolute output error was below 6e-6 and maximum absolute
query-derivative error below 1e-5. Local diagnostic sources and raw values
remain in ignored storage; these observations are not the public exact
certificate, a uniform conditioning theorem, a mixed-kernel numerical
test, a timing result, a GPU experiment, a full gradient check or a
trained-quality result. The public reproducible assertions are (MLH2)
and the rational checks above.

The actionable lead is an exact, explicitly restricted attention block
with a smaller recurrent state. A new mathematical or ML contribution
still requires a retained interaction/recall guarantee or a useful measured
quality--cost tradeoff against the strongest relevant method. The active
unrestricted initial-ideal theorem goal is preserved.
### Twenty-third pass: a uniform bound from lattice facets and planar shadows

The preceding turn committed the ML application checks in `afdccf3`. That
was progress on the application assessment, not an advance on G_d(n).
The present pass gives a uniform proof of a smaller upper bound than the
proof-checked Babson--Onn--Thomas bound. Priority remains a separate gate:
this derivation is not yet a certification of the strongest current bound
or a completed publication manuscript.

For every fixed integer d>=3 define

    D_d = d^2-d-1,
    E_d = (d-1)*(6d(d-2)+4)/(3D_d),
    L_d = (d-1)^2/D_d.

The resulting all-fields, all-ideals bound is

    F_d(n) <= G_d(n) = O_d(n^E_d * log(n+1)^L_d).             (LS1)

In particular,

    G_3(n) = O(n^(44/15) * log(n+1)^(4/5)),
    G_4(n) = O(n^(52/11) * log(n+1)^(9/11)).

The power gain over 2d(d-1)/(d+1) is exactly

    2(d-1)(d-2)/(3(d+1)(d^2-d-1)) > 0.                       (LS2)

This covers nonradical ideals and every characteristic. The logarithm
does not remove the strict asymptotic power gain. No sharp growth claim,
prime-distribution consequence or improved running-time claim is made.

#### A geometric lemma with an explicit proof mechanism

Let P be a full-dimensional lattice polytope in R^d, d>=3. Write V for its
number of vertices, S for its Euclidean boundary (d-1)-measure, and s(P)
for the expected number of vertices of its orthogonal projection onto a
uniform random two-dimensional linear subspace. Put

    a = d(d-2)/(d-1).

Then, with constants depending only on d,

    V <= C_d*(s(P)*H^a + S/H), for every H>=1,               (LS3)
    V <= C_d*S^(a/(a+1))*s(P)^(1/(a+1)).                    (LS4)

The maximum shadow size may replace the expectation. In dimension three
(LS4) reads V <= C*S^(3/5)*s(P)^(2/5), equivalently
s(P) >= c*V^(5/2)/S^(3/2). The proof uses integer facet normals and the
boundaries of vertex normal cones. It does not use an assumed relation
between total primal edge length and expected shadow size.

**Facet incidence budget.** Every facet F has a primitive integer outward
normal z_F. Triangulate F using all its k_F vertices. Such a triangulation
has at least k_F-d+1 top-dimensional simplices: its dual adjacency graph is
connected, the first simplex has d vertices, and each subsequent simplex
along a spanning tree can introduce at most one new vertex. For an integer
(d-1)-simplex in that hyperplane, the vector of signed maximal minors of
its edge matrix is a nonzero integer multiple of z_F. Its volume is
therefore at least ||z_F||/(d-1)!. Since k_F>=d, this gives

    k_F <= d!*area(F)/||z_F||.                               (LS5)

Call a vertex bad when it belongs to a facet with ||z_F||>H. By (LS5) the
number of bad vertices is at most d!*S/H. All remaining vertices are good.

**Size of a good normal cone.** At any vertex, the incident facet normals
span R^d and generate its full-dimensional pointed normal cone. At a good
vertex choose d independent incident normals z_1,...,z_d, all of length
at most H. Their determinant is a nonzero integer. For u_i=z_i/||z_i||,

    |det(u_1,...,u_d)| >= H^(-d).

The simplex conv(0,u_1,...,u_d) lies in the normal cone intersected with
the unit ball. Radial integration shows that the spherical (d-1)-area A
of the cone section satisfies

    A >= 1/((d-1)!*H^d).                                    (LS6)

This argument does not assume that the vertex is simple.

**Boundary from area.** A pointed cone section is contained in an open
hemisphere. Stereographically project that hemisphere from its opposite
pole to the unit ball in R^(d-1). In these coordinates the spherical metric
has conformal factor 2/(1+|x|^2), between one and two. If E is the image,
then A <= 2^(d-1)*vol(E), while the spherical boundary measure p is at
least the Euclidean perimeter of E. The classical Euclidean isoperimetric
inequality, applied to this piecewise smooth domain, consequently gives

    p >= c_d*A^((d-2)/(d-1)) >= c'_d*H^(-a).                (LS7)

Thus no spherical isoperimetric or integral-geometric extension is being
silently assumed. Euclidean isoperimetry is a classical input, not a new
claim of this project.

**Boundary from shadows.** For an edge e of P let A_e be the spherical
(d-2)-measure of its normal cone. Its cone is a wall shared by precisely
the two normal cones of its endpoints. Lower-dimensional intersections
have zero relevant measure, so

    sum_(vertices v) p_v = 2*sum_(edges e) A_e.

Let E_2 be a uniform random two-plane and let H_e be the (d-1)-dimensional
linear span of the edge normal cone. Almost surely E_2 intersects H_e in
a line. By rotational invariance within H_e that line is uniform among
unoriented lines. Since the normal cone is pointed, its section and its
negative are disjoint up to measure zero. The probability that the line
meets the cone is exactly 2*A_e/area(S^(d-2)).

For almost every E_2, its great circle avoids every normal cone of
codimension at least two. Its crossings of the walls correspond exactly
to edges of the shadow polygon, hence also to vertices of that polygon.
Each wall contributes at most one crossing because its cone is pointed.
Linearity of expectation now proves

    sum_v p_v = area(S^(d-2))*s(P).                          (LS8)

Together (LS7)--(LS8) bound the good vertices by C_d*s(P)*H^a; adding the
bad vertices proves (LS3).

For completeness the additive form can be removed. Each vertex belongs
to at least d facets, and ||z_F||>=1; summing (LS5) gives
S>=V/(d-1)!. Choose H=2d!*S/V, which is at least one. At most V/2 vertices
are bad. Inserting this choice in the good-vertex bound yields
V^(a+1)<=C_d*s(P)*S^a and proves (LS4).

#### Matroid shadows and the universal monomial ground set

An established theorem of Dey, combined with Eppstein's token-chain
reduction, bounds the number of successive optimal bases for an m-element
rank-r matroid with affine one-parameter element costs by O(m*r^(1/3)),
after generic perturbation. The inspected proof sources are
[Dey's six-page conference manuscript, Section 4.1, Theorem 4.1 and its proof](https://www.researchgate.net/profile/Tamal-Dey-2/publication/3728283_Improved_bounds_on_planar_k-sets_and_k-levels/links/0912f50aa92321dfc8000000/Improved-bounds-on-planar-k-sets-and-k-levels.pdf)
and [Eppstein's author-hosted manuscript, Section 4, Theorem 4 and its proof](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf).
Dey bounds r edge-nonoverlapping convex polygons on m+r lines by
O((m+r)*r^(1/3)+(m+r)^(2/3)*r^(2/3)); since r<=m this is O(m*r^(1/3)).
Eppstein traces r tokens through the optimal bases, obtaining convex
chains with one decreasing-slope exchange per change. Both relevant
proofs were inspected. The conference version is not asserted to be
identical to Dey's longer 1998 journal version.

Consequently any planar linear projection of a rank-r matroid basis
polytope, including a projection already mapped into R^d, has at most
O(m*r^(1/3)) vertices. Indeed directions (1,t) and (-1,t), t in R, cover
the two open half-planes of objective directions up to positive scaling.
Every polygon vertex is strictly exposed in one of them. Ties among
different bases with the same projected point cause no problem: choose
one strict exposing parameter for each of the finitely many distinct
points. At these parameters there is a positive gap to every basis
outside the exposed group. A sufficiently small generic perturbation of
all element cost lines preserves those groups. The resulting sequence
therefore has at least as many distinct bases as the original projected
polygon has vertices on that half-plane. This justifies using the
generic theorem without counting exponentially many tied bases.

For a colength-n ideal I over any field use the established ground set

    V_n^d = {u in N_0^d : product_i(u_i+1)<=n},
    m = |V_n^d| = O_d(n*log(n+1)^(d-1)).

The residues of these monomials represent a rank-n matroid. Every initial
staircase is contained in V_n^d because it contains all divisors of each
of its monomials. The set spans the quotient since it contains at least
one such staircase. Define B(I) as the convex hull of the sums of exponent
vectors over all matroid bases in this ground set. The established
state-polyhedron correspondence injects the initial staircases into the
vertices of B(I). These are the basis polytope and correspondence of
[Babson--Onn--Thomas, Theorem 2.4 and Section 3](https://arxiv.org/pdf/math/0207135),
not new definitions or correspondences claimed here.

Each exponent coordinate is at most n-1, so B(I) is a lattice polytope in
[0,n(n-1)]^d. It need not be full dimensional. Put Q=B(I)+[0,1]^d.
The normal fan of this Minkowski sum refines that of B(I), so Q has at
least as many vertices. More explicitly, an objective uniquely exposing
a vertex of B(I), chosen also generic for the cube, exposes its sum with
a cube vertex; uniqueness of an exposed Minkowski decomposition makes
this assignment injective. Each planar projection of the cube is a
zonogon with at most 2d vertices. Adding it introduces at most 2d new
polygon edge directions, including when the projected B(I) is a segment
or a point. Therefore

    s(Q) <= O_d(m*n^(1/3))+2d
          = O_d(n^(4/3)*log(n+1)^(d-1)).                    (LS9)

Q is contained in a cube of side R=n(n-1)+1. Its surface area is at most
2d*R^(d-1). An elementary proof sums the areas of its facets multiplied
by the absolute coordinate components of their unit normals: each
coordinate projection is covered twice, and the sum of the absolute
components is at least one. Each projected area is at most R^(d-1).
Thus S(Q)=O_d(n^(2(d-1))). Substitution in (LS4) proves (LS1). The
arithmetic simplification of the exponent and logarithm gives (LS2).
All constants are independent of the field and the presentation of I.

The same argument provides a reusable bound for any rank-r matroid on
m elements assigned integer vectors in [0,p]^d. If B is the convex hull
of their basis sums, then

    f_0(B) <= O_d((rp+1)^(d(d-1)(d-2)/D_d)
                     *(m*r^(1/3)+1)^((d-1)/D_d)).          (LS10)

This is a cardinality bound; it does not by itself provide a faster
enumeration algorithm or independence oracle.

#### Further algebraic consequence

Every n-element staircase has O_d(n^((d-1)/d)) minimal outside monomials.
Here is a direct proof of the estimate used in the earlier literature.
Set h=ceil(n^(1/d)). Minimal outside exponents form an antichain. Those
in [0,h-1]^d number at most h^(d-1), by projecting along one coordinate.
For each coordinate i, a minimal outside exponent u with u_i>=h gives
the h inside points u-e_i,...,u-h*e_i. These chains are disjoint for
different u at fixed i, again by the antichain property. There are at
most n/h such u per coordinate. Hence the total is at most
h^(d-1)+dn/h=O_d(n^((d-1)/d)).

Each reduced Groebner basis has one polynomial per minimal outside
monomial. Taking the union over initial ideals and using (LS1) gives

    |U(I)| <= O_d(n^(E_d+(d-1)/d)*log(n+1)^L_d).            (LS11)

For d=3 this is O(n^(18/5)*log(n+1)^(4/5)), improving the proof-checked
historical exponent 11/3. The multiplication argument is the known
argument of Babson--Onn--Thomas Corollary 2.5(2); the changed quantitative
input is (LS1). This counts distinct elements of the universal basis,
not their total coefficient encoding or computation time.

#### Certificates, primary comparisons and remaining publication gate

The dependency-free checker now records eight integer polytope fixtures,
including nonsimple vertices, a unimodularly skewed cube, rational and
finite-field point-quotient basis polytopes, and a lower-dimensional basis
polytope regularized by the cube. It saves 63 facet triangulation
certificates, 68 independent incident-normal determinants and 64 exact
planar projections. Projected convex hulls and intersections with the
primal edge normal cones are calculated independently and agree edge by
edge. Three point-matroid fixtures supply 66 rational parameter samples
and ten single-exchange token transitions, independently checked against
all bases. Exact exponent arithmetic covers dimensions 3--32; the proof,
not this finite range, establishes the identities for every d>=3.
No sampled integral is used as a proof of (LS8) or of isoperimetry.

The directly matched historical comparison is
[Babson--Onn--Thomas Corollary 2.5 and its complete proof](https://arxiv.org/pdf/math/0207135).
Its second assertion gives the historical universal-basis exponent just
compared. The more recent
[multi-parametric matroid preprint, Corollary 3.4 and Theorem 3.6](https://arxiv.org/html/2503.08178v1#S3)
was checked through the definitions, arrangement argument and complete
Theorem 3.6 proof. For p parameters it gives O(m^(2p)) regions and an
algorithm with that factor; p=2 does not imply (LS1). Its input has
arbitrary cost coefficients, whereas the additional lattice/surface
budget is essential here.

[Black--Criado's inspected manuscript, Theorem 1.1 and Corollary 1.4](https://arxiv.org/html/2406.06936v1)
bounds expected shadows from diameter and edge lengths, giving an upper
bound O_d(R) for a polytope in [0,R]^d. It does not supply (LS4), which
bounds vertices using surface area and a shadow upper bound. The served
HTML is labeled arXiv v1 but bears an August 2026 date; no unverified
identity with other versions is assumed. The inspected principal
statements of [Small Shadows of Lattice Polytopes](https://arxiv.org/pdf/2204.09129)
concern selected short pivot paths and monotone diameter, not a bound on
all vertices from a maximum shadow. The earlier Moser construction
explains why a shadow bound without a lattice-size budget is insufficient.
These comparisons have the reading scopes stated here; a search failure
is not evidence that (LS4) or (LS1) is new.

The next publication task is a focused priority audit of (LS4), (LS10)
and the resulting unrestricted initial-ideal exponent, followed by a
separate manuscript if that audit supports novelty. This is now an
actual uniform upper-bound proof and further algebraic consequence,
not merely another family obstruction. The strongest-current-priority
gate, independently assessed proof and publication artifact remain
unfinished, so the goal stays active. Existing manuscripts and kernels
remain unchanged; no outbound correspondence or submission is authorized.
### Twenty-fourth pass: manuscript, classical averaging and consequence checks

The preceding turn is substantive progress: commit `b846187` saves the
uniform all-fields bound (LS1), its geometric proof and exact certificates.
This pass prepares a separate complete working manuscript and makes the
attribution and remaining priority question more precise. It does not
replace the unrestricted objective with manuscript preparation alone.

#### Classical projection identities and the finite-field comparison

[Grunbaum, Grassmann angles of convex polytopes (1968)](https://faculty.washington.edu/moishe/branko/BG64%20Grassmann%20angles.pdf)
provides the classical framework behind (LS8). Section 2 defines the
rotation-invariant Grassmannian measure. Equation (2.5), with its proof,
reduces intersection probabilities to the span of a cone; Section 3's
equation (3.1) and averaging argument express projected-face counts as
sums of local probabilities. These relevant definitions and proofs were
inspected, together with the adjacent cone-polarity and averaging results.
Our direct wall-count proof is retained for precision, with explicit
credit to this existing principle. It is not a new Grassmann-angle or
Crofton identity.

The precise later Stigler--Zhang comparison can be made without silently
equating its 2019 preprint with the 2020 chapter. The indexed primary
[later deposit, Theorem 8](https://par.nsf.gov/servlets/purl/10177092)
bounds the number for m points in d-dimensional prime-field space using,
in its first range, the factor

    (p^2*floor(m/p)+(m mod p)^2)^(d(d-1)/(d+1)).

For p>m this is exactly m^(2d(d-1)/(d+1)), recovering the earlier exponent.
The relevant theorem and displayed proof were inspected through indexed
primary pages 15--16; complete PDF retrieval remains unconfirmed. This
supports the named quantitative comparison, not an assumption about every
statement in the chapter. Its fixed-field or high-density bounds can be
stronger in their own parameter regimes. Likewise the Euler-zonotope
bound can be stronger for a fixed integer coordinate range. Neither
comparison is presented as blanket dominance by (LS1) or (LS10).

#### Complete working manuscript

A separate anonymous working manuscript is saved in the ignored project
`results/local/lattice_shadow_paper/`, with `main.tex` and the standalone
11-page PDF at `output/pdf/main.pdf` within that project. It contains the
following complete arguments:

1. The integer facet incidence budget, with a triangulation argument.
2. The determinant lower bound for spherical normal-cone area, including
   nonsimple vertices, and the stereographic Euclidean-isoperimetry step.
3. The classical normal-fan averaging identity, derived explicitly.
4. The Dey--Eppstein planar matroid bound, with a separate perturbation
   argument for tied bases and projected vertices.
5. Cube regularization and the coordinate-projection surface-area bound,
   covering lower-dimensional basis polytopes.
6. The general integer-weight matroid estimate and the all-colength
   application, including an elementary hyperbolic-cross cardinality proof.
7. The universal Groebner-basis cardinality consequence, with an explicit
   proof of the established minimal-generator packing estimate.

The paper states exactly which inputs are existing results, includes a
version-aware priority appendix and an AI-assistance disclosure, and
avoids claims of sharp growth, a new enumeration algorithm, prime
distribution or an ML performance consequence. It is a complete working
paper for internal human review. Independent human refereeing and external
submission have not occurred; neither is introduced as an extra condition
on the user's research goal. The unresolved gate is the defensible novelty
assessment for the unrestricted quantitative result.

The editor compiler and terminal LaTeX compiler both confirm compilation.
All 11 pages were rendered for layout inspection, with the proof and
certificate pages additionally inspected at full size. The final build has
no unresolved references or overflowing boxes, and the PDF author metadata
is empty. The goal manifest records exact source/PDF hashes; the paper and
its private reading/QA records are not included in the public commit.

#### Independent certificates for the algebraic consequence

The new `staircase_generator_controls` checks 325 staircases: sizes one
through seven in dimensions two and three, and sizes one through five
in dimension four. It finds their 1462 minimal outside generators in two
independent ways: a one-step minimal border and the first missing point
on coordinate fibers, followed by removal of dominated caps. It verifies
the antichain property, the injective small-coordinate projection and
the disjoint inside chains used in the proof of (LS11). Each fixture
saves its staircase, generators and exact rational packing bound.
The estimate is established mathematics, and finite agreement is not
the uniform proof or a new generator theorem.

Run the existing complete research command:

```bash
python3 -B scripts/run_structural_complexity_checks.py
```

The public report retains all earlier controls and adds the separate
`staircase_generator_control` block. The historical completed manuscript,
original kernels, CLI and ML assessment remain intact.

#### Publication judgment at this point

The all-colength power improvement is mathematically substantial enough
to justify the completed working draft: its scope is every field and
every colength-n ideal in each fixed d>=3. This is a contribution
assessment, not a prediction of journal acceptance or broad impact.
The precise comparison with the proof-checked Babson--Onn--Thomas bound
is established. The remaining priority audit must determine whether an
equivalent mixed lattice/angle estimate or an existing unrestricted
counting result already implies the claimed advance. Failure to find a
keyword match is not evidence that it is new. The goal remains active
with this question explicit; the draft alone does not complete it.

### ML architecture extension: exact power features and delta-state closure

The user asked whether the efficient factorizations can improve existing ML
hot paths. The actionable connection is polynomial feature expansion and
recurrent memory reads/writes. The prime-specific and initial-ideal counting
theorems do not themselves accelerate arbitrary Transformer weights. This
extension records a more useful architectural option than merely shrinking
the quadratic example: retain higher-order interactions within a modest
state budget, and check that compression preserves the entire memory update.

#### Higher powers with an explicit representation constraint

Let A={0,...,m}^s, theta_0=1 and theta_a>0, shared across the sequence.
Set z(x)=(sqrt(theta_a)*x^a)_(a in A), D=|A|, r(x)=||z(x)||_2>0.
For an even positive integer p use the normalized power kernel

    K_p(q,k)=(z(q).z(k)/(r(q)*r(k)))^p.

The encoders producing q,k in R^s may be learned and distinct. Their
expanded coordinates are constrained; this is not unrestricted R^D.
Write P(t)=sum_(a in A) theta_a*t^a and let c_b be the coefficients of
P(t)^p. Every c_b is positive and its support is B={0,...,pm}^s.
Indeed every integer coordinate between zero and pm is a sum of p
integers between zero and m; these coordinatewise decompositions combine
to p vectors in A. Consequently

    phi_b(x)=sqrt(c_b)*x^b/r(x)^p, b in B,
    K_p(q,k)=phi(q).phi(k), ||phi(x)||_2=1.              (MLP1)

For fixed theta the M=(pm+1)^s functions in (MLP1) are linearly
independent on any nonempty open latent domain: multiply a proposed
linear relation by r(x)^p to obtain a relation among distinct monomials.
Positive coefficients then give kernel separation rank M, as in the
established coordinate-ring framework. This is not a new rank theorem.

For s=2,m=7, the augmented embedding has D=64 coordinates:

| Power p | Generic symmetric power width | Collected width M |
| --- | --- | --- |
| 2 | 2080 | 225 |
| 4 | 766480 | 841 |
| 6 | 119877472 | 1849 |

These compare two exact representations of this same constrained kernel.
The generic column describes an algorithm that does not collect the
additional monomial identities. It is not an assertion that our two
latent coordinates retain arbitrary 64-dimensional interaction capacity.
The power kernel also differs from the earlier quadratic Taylor kernel;
their coinciding width counts do not make their outputs equal.
[Power Attention, Section 4 and Lemma 4.2](https://arxiv.org/html/2507.04239v1#S4)
already supplies the symmetric-power kernel, width and hardware-aware
implementation. [Altschuler--Parrilo, Proposition 3.6 and Corollary 3.8](https://arxiv.org/html/2106.02755v1#S3.SS1)
already supplies the restriction/rank principle. No novelty for their
combination has been established.

Coefficient preparation need not enumerate D^p tensor entries. Pad the
coefficient box to length L>=pm+1 in every coordinate, evaluate P by an
s-dimensional discrete Fourier transform, raise each value to p, and
inverse-transform. No wraparound occurs because the true degree in each
coordinate is at most pm. With power-of-two padding this costs
O_s(M*(log(M+1)+log(p+1))) unit-cost arithmetic for fixed s, including
exponentiation by squaring. It is a standard convolution algorithm, not
a bit-complexity or numerical-error theorem. Charge it whenever theta
changes; differentiating it is part of training cost. Forward feature
evaluation plus a recurrent update/read costs O(M*(s+v)), with value
width v, in addition to encoders and gates. Training memory, recomputation
and actual GPU utilization still require measurement.

#### The compressed subspace is closed under scalar-gated delta updates

This point is stronger than agreement of pairwise kernels. Index generic
symmetric features by multiplicity vectors nu with sum(nu)=p. Their
coefficient is w_nu=p!*product(theta_a^nu_a)/product(nu_a!), and their
monomial exponent is b(nu)=sum(nu_a*a). The generic feature is
sqrt(w_nu)*x^b(nu)/r(x)^p. Since c_b=sum_(b(nu)=b) w_nu, the matrix

    U_(nu,b)=1[b(nu)=b]*sqrt(w_nu/c_b)

satisfies U^T*U=I and f_generic(x)=U*phi(x). Initialize both memory
states at zero. For scalar gates alpha_t,beta_t and value row v_t^T,
use the ordinary gated delta recurrence

    S_t=alpha_t*S_(t-1)
        +beta_t*f(k_t)*(v_t^T-alpha_t*f(k_t)^T*S_(t-1)). (MLP2)

Substitute S_(t-1)=U*T_(t-1) and f=U*phi into (MLP2). The identity
U^T*U=I gives S_t=U*T_t, where T follows the same recurrence in M
coordinates. Every read f(q)^T*S_t equals phi(q)^T*T_t. This proves
exact equivalence for every sequence, not just an additive attention sum.
It also covers gates determined by already identical outputs. Arbitrary
feature-coordinate gates, feature-wise nonlinearities or a learned initial
state outside this subspace require a separate closure argument.

For unit features, beta in [0,2] makes I-beta*phi*phi^T nonexpansive.
With |alpha|<=1 the homogeneous state transition remains nonexpansive.
This does not bound the sum of forced writes. At beta=1, the just-written
key reads back v_t exactly, before subsequent writes can interfere.
These are established delta-rule properties; exact latest-write recall is
not a guarantee for all previous keys. The full mathematical identity also
implies equality of derivatives where both implementations differentiate
all theta-dependent operations correctly; no full backward implementation
or optimizer-cost advantage is claimed here.

[DeltaNet, Sections 2.2--3.2](https://arxiv.org/html/2406.06484v1#S2.SS2)
already gives the correction, Householder-product representation and
chunkwise training. [Gated DeltaNet, Section 3.1, equation (8)](https://arxiv.org/html/2412.06464v1#S3.SS1)
already gives exactly the scalar-forgetting recurrence above. Those
definitions and displayed derivations were inspected. No novelty is claimed
for delta memory or parallelizing it. Its existing chunk algorithm is a
required implementation comparison, rather than a sequential Python loop.

#### Verification and the practical gate

The dependency-free public checker now builds six exact examples in one
and two latent dimensions at powers two, four and six. It compares direct
pair kernels against convolution coefficients (150 pairs), independently
expands generic symmetric multinomials, verifies their grouped coefficients,
and compares generic and collected scalar-gated delta states and outputs
at every prefix (150 query reads). It checks unit-feature exact writes at
beta=1. These finite controls support (MLP1)--(MLP2); the derivation above
establishes the general identity.

An additional ignored local NumPy diagnostic used the 64-coordinate
embedding, 128 random latent keys in [-1,1]^2, value width 16 and fixed
positive random weights. At powers two, four and six, FFT coefficients
agreed with independent direct convolution to relative L2 error below
5e-16, and collected float32 attention outputs agreed with direct float64
pair kernels to maximum absolute error below 5e-7. This is bounded forward
evidence, not conditioning, training, runtime or GPU validation.

The same diagnostic deliberately queried stored keys. Although every key
was its own top-scoring match, the degree-six normalized attention assigned
the correct key only about 0.0784 weight on average. For independent
mean-zero unit-variance scalar values, its analytic average retrieval MSE
was about 0.886, computed as ||A-I||_F^2/128. One pass of ungated beta=1
delta writes had MSE about 1.437, despite exact newest-key writes. The zero
predictor has MSE one in this diagnostic. Thus neither this fixed embedding
nor a delta correction by itself establishes useful recall. These are
untrained, single-fixture diagnostics, not comparative model-quality results.

The next useful experiment is learned multi-query associative recall with
matched total state, parameters and training compute: ordinary learned
features in DeltaNet, polynomial Power Attention and this collected feature
family. Include encoder, coefficient-update, backward and state-traffic
costs. Test several latent dimensions/heads and out-of-distribution context
lengths before custom GPU work. Continue only for a reproducible useful
quality/memory/throughput tradeoff or an additional representation theorem.
Exact compression alone does not establish a fundamental advantage over
existing architectures. This application assessment leaves the active
unrestricted mathematical theorem and priority audit unchanged.

### Twenty-fifth pass: the primary-source boundary of the lattice-shadow claim

The last mathematical goal work was progress: `b846187` proved the uniform
bound and `4aecba3` saved the manuscript record and generator certificates.
The intervening user-directed ML work in `e1f169c` changed the application
assessment, not the unrestricted asymptotic bound. This pass returns to the
remaining theorem-level priority requirement. It supplies positive primary
comparisons and a clearer proof boundary, without declaring novelty from
an unsuccessful search or a missing source.

#### What the classical area budget gives, and what the shadow adds

Retain the notation of (LS3)--(LS8). Vertex normal cones partition the unit
sphere up to null boundaries, so their spherical areas sum to a constant.
By (LS6), a good vertex contributes at least c_d/H^d to that area. Hence
the classical area budget alone gives

    V <= C_d*(H^d+S/H), H>=1,
    V = O_d(S^(d/(d+1))).                               (PA1)

For bounded S, fixed-dimensional constants absorb the cutoff H>=1;
the earlier incidence budget gives S>=V/(d-1)! in any case. Substituting
S=O_d(n^(2(d-1))) in (PA1) recovers precisely the historical all-colength
exponent 2d(d-1)/(d+1). This is a reconstruction of a classical surface
bound, not a new theorem.

The proposed improvement uses the boundary budget instead: (LS7) gives
p_v>=c_d/H^(d(d-2)/(d-1)), and the classical projection identity (LS8)
sums those boundaries using the mean shadow size. This produces (LS3)
with s(P)*H^a in place of H^d. The distinction identifies the additional
quantitative step; it does not prove that no earlier theorem implies it.

#### Nonlinear angle sums and discrete intrinsic volumes

A public university archive contains the original coauthored chapter
[McMullen--Schneider, Valuations on convex bodies (1983)](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/McMullen/McMullen9.pdf).
Section 8, especially Theorems (8.5)--(8.9), gives cone dissections and
internal/external angle identities. The displayed derivations were read,
and (8.9) was checked against the scanned page. The section explicitly
paraphrases McMullen's 1975 argument. Its identities are classical
background; the inspected section does not state our joint cardinality
bound or the primitive-facet-normal cutoff. This is not a claim that the
entire 1975 article, including its lattice application, has been read.

[Dospolova, arXiv:2107.06549v1](https://arxiv.org/html/2107.06549v1)
provides a full primary preprint behind the later 2024 article identified
in the previous audit. Theorems 1.1--1.4, the cone definitions, Theorems
3.1--3.2 with their displayed angle arguments, and Sections 4.2--4.4 were
inspected. These define lattice-point sums weighted by tangent-cone
angles and prove valuation/dilation/reciprocity properties. In particular
the preprint has A_0(P)=1 and G_0(P)=|P intersect Z^d|-1. Those quantities
are not f_0(P) or s(P): on P=[0,t]^d they are 1 and (t+1)^d-1, while
f_0(P)=2^d and the generic planar shadow has 2d vertices. The cited
statements do not themselves state (LS4). No equivalence to the full
2024 translation or exclusion of every possible corollary is asserted.

#### A bounded citation trace selects further actual theorem comparisons

A recorded public OpenAlex lookup of the Babson--Onn--Thomas DOI returned
35 citing work records. The inventory is metadata, not evidence that all
later research has been indexed or read. It selected the following
positive comparisons; the scientific judgments use the primary texts.

[Wagner's 2003 thesis](https://ti.inf.ethz.ch/ew/research/phdtheses/wagner03.pdf),
Section 3.1, explicitly restricts the corner-cut/state-polyhedron
correspondence to generic point ideals and then cites Babson--Onn--Thomas
for finite-colength ideals. Theorem 3.15 gives the corner-cut upper bound
O((k log k)^(d-1)); Section 3.5 extends it to flags of that polyhedron.
The definitions, final count argument and flag extension were inspected.
They do not replace the multiplication-compatible matroid of a general
quotient by a uniform matroid, so they do not supply an unrestricted
G_d(n) bound of that order.

[Lemmens, arXiv:1812.09527v3](https://arxiv.org/html/1812.09527v3),
Theorem 4 and Corollary 5, concerns lattice-point saturation of distinct
subset sums of planar lattice polygons and of the two-dimensional corner
cut polyhedron. The full Corollary 5 proof in Section 3 was inspected.
The application is absent from v1, so the revised version was checked
separately. Lattice-point saturation of a sum set does not enforce
residue-matroid independence or count initial ideals of arbitrary
finite algebras. No identity of the preprint and 2022 article is assumed.

[Berstein et al., arXiv:0707.4618v1](https://arxiv.org/pdf/0707.4618v1),
Theorems 1.1 and 1.3, solves nonlinear matroid optimization. Its
Proposition 4.2 and Lemmas 4.3--4.4 encode attainable basis-weight profiles
by det(A*Y*A^T) and interpolate on the box {0,...,rp}^d. Those statements
and the complete determinant/interpolation derivation were read. They
provide an algorithm for integer matrix representations and unary-bounded
weights, not the smaller vertex cardinality in (LS10). Our comparison is
to this version; a failed published-copy link does not establish its
identity with the 2008 journal article.

[Onn, arXiv:2109.06569v1](https://arxiv.org/html/2109.06569v1),
Theorem 2.2 and its full dynamic-programming proof, optimizes vector
partitions with fixed part count and attribute dimension. Its feasible
objects are partitions satisfying size constraints, without an arbitrary
residue-matroid independence constraint. No unstated reduction or
improved enumeration count is imported from that algorithmic result.

[Onn's author-hosted 2010 monograph, archived in 2013](https://web.archive.org/web/20131015172108/http://ie.technion.ac.il/~onn/Book/NDO.pdf),
is now retrieved and inspected in Section 6.3.3, including Definition 6.29,
Corollary 6.30 with its full proof, and the attribution notes on page 128.
The copy bears the publisher's 2010 imprint, matching ISBN and pagination.
It credits the general-colength framework and universal-basis algorithm
to Babson--Onn--Thomas. Its edge-direction route uses at most binom(N,2)
differences for N ground elements; Lemma 2.13 gives the usual
O(N^(2(d-1))) zonotope bound. These inspected results do not state the
smaller count (LS1). This closes the previously unread book-section
comparison; no identity with another archived copy is assumed.

[Deza--Pournin--Rakotonarivo (2021)](https://www.cas.mcmaster.ca/~deza/cm2021.pdf),
Theorems 1.2--1.4, estimates primitive-zonotope vertices with fixed
coordinate bounds and growing dimension. The primary offprint and its
definitions/principal statements were inspected. These statements do not
establish the optional fixed-dimension, growing-generator-bound sharpness
construction mentioned in the private continuation note. That construction
remains unproved here; no sharpness claim is added to the manuscript.

#### Artifact and remaining gate

The private manuscript now includes these source/version comparisons and
the area-versus-boundary explanation. Its previous source and PDF are
preserved in ignored version storage. This is a substantive priority
audit update, not a new numerical improvement to (LS1).

The original [McMullen 1975 article](https://doi.org/10.1017/S0305004100051665),
including its lattice application, remains incompletely inspected.
Publisher metadata identifies this target; neither an abstract nor a
failed public download settles its relevant content. The author-archived
Onn book comparison above is now complete at its stated section scope.
The later finite-field deposit still has only the exact indexed theorem
comparison already recorded. These are explicit limits, not reasons to
assume a theorem is new or to contact anyone. The strongest-current-bound
and novelty flags remain false, and the full mathematical goal stays
active. Kernel and historical-manuscript bytes are preserved.

### ML hot-path decision: preserve selective forgetting in the compressed layer

The strongest present application candidate is a structured recurrent
attention layer. The saving concerns redundant polynomial features and
memory reads/writes. Neither prime factorization nor the lattice-shadow
initial-ideal bound supplies a faster algorithm for arbitrary learned
weights. A useful architecture needs both economical computation and
retained learning/recall quality; the exact identity alone settles only
the former's representation count.

#### Which memory transitions survive compression?

Use the fixed orthonormal-column grouping map U and unit feature phi from
(MLP1), with generic feature f=U*phi. For any linear state transition A,
the subspace im(U) is invariant exactly when

    (I-U*U^T)*A*U=0.

Indeed invariance implies A*U=U*B with B=U^T*A*U, and conversely this
identity implies invariance. This is elementary linear algebra, not a new
compression theorem. For a diagonal generic gate diag(g_nu), the grouping
map's positive entries give the sharper criterion

    diag(g)*U=U*diag(gamma)
    iff g_nu=gamma_b whenever b(nu)=b.                    (MLG1)

Every nonzero entry U_(nu,b) forces the corresponding equality; it also
suffices. Thus scalar forgetting is not required. Choose arbitrary learned
gates gamma_(t,b) in [0,1] on the M collected coordinates and use

    T_t=(I-beta_t*phi(k_t)*phi(k_t)^T)*diag(gamma_t)*T_(t-1)
        +beta_t*phi(k_t)*v_t^T.                           (MLG2)

The generic model with tied gates (MLG1) satisfies S_t=U*T_t at every
prefix, including the delta correction and all query reads. The proof is
substitution using U^T*U=I. At beta=1 the current key reads its new value
exactly. For beta in [0,2] the homogeneous transition is nonexpansive,
since both factors are contractions. None of this guarantees that earlier
values remain retrievable. Shared theta must remain fixed within the
sequence, as in (MLP1); arbitrary time-dependent changes of U are excluded.

Equation (MLG2) is the established channel-gated delta recurrence applied
to our deliberately restricted features.
[Kimi Linear, Section 3, equation (1)](https://arxiv.org/html/2510.26692v1#S3)
already gives that recurrence and its chunkwise implementation; Sections
3.1--3.2 discuss parallelization and precision. We inspected those displayed
formulas and adjacent explanation, not all experiments or kernel code.
[GLA, Section 4.1](https://arxiv.org/html/2312.06635v6#S4.SS1)
already allows channel-dependent forgetting and discusses gate-parameter
costs. No new gate, parallel algorithm or improvement over these methods
is claimed. The architecture can use their recurrence interfaces with
collected width M, although this integration has not been implemented or
benchmarked against their kernels.

The dependency-free checker adds six exact two-value fixtures: independent
collected-channel gates and cheaper exponent-separable gates
gamma_b=g_0*product_j g_j^(b_j), with g_j in [0,1]. It verifies 30 invariant
state steps, 300 query-output coordinates and 24 exact latest-write
coordinates. The separable family needs only s+1 gate outputs but imposes
an additional restriction; independent M-channel gates do not.

A negative control uses P(x)=1+x+x^2, power two, and a unit-feature write
at x=1 with value one. Exponent two has two generic coordinates of weights
two and one. Two gate-only steps retain disjoint members of this pair.
The full state then reads zero. Replacing each gate by its weighted group
average first matches the read 2/9, but after the second step incorrectly
reads 2/27. Projection after each transition therefore does not justify
compressing an arbitrary existing gated model. This counterexample is
saved exactly alongside the positive controls.

#### What an improvement would have to beat

For s=2,m=7,p=2, the same constrained polynomial model has 2080 generic
symmetric coordinates or 225 collected coordinates. With value width 64,
their delta states contain 133120 or 14400 scalars: 520 or 56.25 KiB in
float32. These are state counts, not total-model memory or a measured
9.24-fold speedup. Ordinary KDA with 225 features already has the same
state size; our candidate must earn its extra feature-encoder cost through
a better quality/cost tradeoff. Raising polynomial degree while retaining
a small latent dimension does not restore unrestricted query/key capacity.

Forward recurrence arithmetic remains O(M*v), plus feature evaluation,
gate/latent encoders and coefficient preparation. Generating M gates
with a dense model-width-to-M projection can be expensive. Chunkwise
training, all derivatives, activations, padding and state traffic must be
charged. [Power Attention, Section 4.1](https://arxiv.org/html/2507.04239v1#S4.SS1)
explains why a mathematically smaller feature representation need not
produce a faster GPU kernel. Its polynomial features are also a mandatory
baseline, not evidence that combining existing ingredients is novel.

For dense channel/FFN mixing, [Monarch](https://arxiv.org/abs/2204.00595)
already designs fast structured matrices. For attention-cache factorization,
[Tensor Product Attention](https://arxiv.org/abs/2501.06425) and
[DeepSeek-V2's MLA](https://arxiv.org/abs/2405.04434) are existing architecture
families. Their primary abstracts were rechecked here; this is not a
proof or experiment audit. No project-specific advantage in those hot paths
has been established. Rank-one token updates also need not leave a low-rank
accumulated memory, as the earlier exact Vandermonde control demonstrates.

The next implementation gate is one trained multi-query recall comparison:
native learned-feature KDA, scalar-gated structured features, and (MLG2),
at matched total state, parameters and training compute. Include Power
Attention in the subsequent sequence-model comparison. Report held-out
keys, longer contexts, interference and overwrite behavior across seeds;
then measure complete forward/backward and decoding on the same GPU.
Only a useful quality/memory/latency frontier justifies a custom fused
kernel or wider architecture program. No trained-model or GPU result was
produced in this pass. This application assessment preserves the active
unrestricted mathematical goal and its open priority audit.

### Twenty-sixth pass: sharpness of the lattice-shadow mechanism

The preceding mathematical pass made progress by comparing the proposed
bound with eight named primary sources and extending the manuscript;
`217265a` saved that audit alongside user-directed ML checks. It did not
change (LS1) or close its priority gate. This pass proves that the two
powers in the geometric bridge (LS4) are independently optimal. It keeps
the original all-colength objective and its unresolved literature gate.

#### A matching family in every fixed dimension

Let G_d(p) contain one representative of each primitive integer direction
in [-p,p]^d, chosen by making its first nonzero coordinate positive. Set

    Z_d(p) = sum_(g in G_d(p)) [0,g].

These are established primitive zonotopes, already defined and used for
integer-weight matroid projections by
[Deza--Manoussakis--Onn, Section 3.1](https://optimization-online.org/wp-content/uploads/2015/12/5266.pdf).
Neither the family nor its pair-matroid realization is claimed new. The
fixed-coordinate, growing-dimension statements compared in pass 25 do
not supply the following fixed-dimension argument; no novelty for the
zonotope vertex asymptotic itself is inferred from that distinction.

For every fixed d>=3, as the integer p tends to infinity, we prove

    f_0(Z_d(p)) = Theta_d(p^(d(d-1))),
    S(Z_d(p))   = Theta_d(p^(d^2-1)),
    s(Z_d(p))   = Theta_d(p^d).                              (SH1)

Here S and s have exactly the meanings used in (LS3)--(LS8). With
alpha=d(d-2)/(d^2-d-1) and beta=(d-1)/(d^2-d-1),

    alpha*(d^2-1)+beta*d = d(d-1).                          (SH2)

Thus (SH1) saturates (LS4) up to constants. Decreasing either alpha or
beta by a fixed positive amount, while holding the other fixed, makes
the proposed upper bound grow strictly slower than this family's vertex
count. This is a sharpness theorem for the geometric bridge, not a
matching lower bound for G_d(n) or for its Dey-based specialization.

#### Primitive tuples without an unproved density assumption

Put delta=1/(100d). Let C_j(h), j=1,...,d-1, be the integer box near
h*e_j with diagonal interval [h,h+floor(delta*h)] and all other intervals
[0,floor(delta*h)]. For all sufficiently large h, at least a fraction
8^(-(d-1)) of ordered column tuples Q=[q_1 ... q_(d-1)] have gcd one of
the maximal minors. In particular there are Omega_d(h^(d(d-1))) such
tuples. The following finite-prime argument establishes the density for
these moving boxes; a global natural-density statement alone would not.

At stage r assume the previous r-1 columns, in their first r rows, form
a primitive r-by-(r-1) matrix U. The empty first stage is allowed. Its
cofactor normal is primitive. For the new first-r-coordinate column u,
D=det[U u] is nonzero: the first-r square matrix divided by h differs
from I in maximum row-sum norm by at most r*delta<1. The coarse bound
abs(D)<=2^(d^2)*h^(2(d-1)) holds uniformly.

For each prime ell, D=0 modulo ell is a nontrivial equation in fresh
coordinates, of probability at most 1/ell+a/h, a=100d. Conditional on
this equation, u=U*lambda modulo ell has a unique lambda because U has
full column rank. Failure of the full next-row matrix is one congruence
on the fresh last entry, again with conditional probability at most
1/ell+a/h. Thus failure at a specified prime has probability at most
(1/ell+a/h)^2.

Let t=floor(log_2(h)/2). The union over primes ell<=2^t is bounded by

    3/4 + 2a(t+1)*2^(-2t) + a^2*2^(-3t).

This uses sum_(j>=2) j^(-2)<3/4 and H_(2^t)<=1+t. For larger primes,
first condition on U,u. Each bad prime divides the same nonzero D.
Once t>=d^2+4(d-1), at most b=4d-3 such primes occur. The last-entry
congruence bounds their total probability by

    b*(2^(-t)+a*2^(-2t)).

All error terms tend to zero for fixed d. Hence each conditional stage
succeeds with probability at least 1/8 for every sufficiently large h.
Earlier stage events use only earlier projected coordinates, so the new
row and column retain their independent interval distributions. Multiplying
over the d-1 stages proves the claim. This does not exchange an infinite
prime sum with a limit.

Gcd one of the maximal minors also gives an integer left inverse LQ=I:
combine the adjugate identities of the maximal row submatrices by integer
Bezout coefficients. Thus Q*R^(d-1) intersect Z^d equals Q*Z^(d-1), and
each column is individually primitive.

#### From saturated spans to vertices, surface and shadows

A hyperplane W spanned by a good tuple contains at most 5^(d-1) integer
points in [-2h,2h]^d. Indeed write such a point as Q*c; saturation makes
c integral. The top square block divided by h is within (d-1)*delta<1/100
of I, so its inverse has infinity norm at most 1/(1-(d-1)*delta). Hence
norm_inf(c)<=2/(1-(d-1)*delta)<3, giving at most five choices per entry.
At most 5^((d-1)^2) ordered tuples therefore span the same W. There are
Omega_d(h^(d(d-1))) distinct W.

Take p>=2h. The central arrangement of g-perpendicular hyperplanes,
g in G_d(p), has a different line flat W-perpendicular for each W just
counted. A generic affine hyperplane avoiding zero meets these lines in
distinct vertices of its restricted arrangement. Any affine arrangement
has at least as many regions as vertices: deletion--restriction adds
exactly the restricted region count, and all new vertices are vertices
of the restriction; induction on dimension and hyperplane count proves
the assertion. Restricted regions inject into central chambers by
convexity. Central chambers correspond bijectively to zonotope vertices
by strictly maximizing each segment endpoint. This gives the lower count
in (SH1), taking h=floor(p/2). The standard deletion--restriction upper
bound is O_d(|G_d(p)|^(d-1))=O_d(p^(d(d-1))).

The tuple count yields Omega_d(h^d) distinct primitive first columns,
because each has only O_d(h^(d(d-2))) completions. Permuting coordinates
gives the same count near every positive coordinate axis. Consequently
|G_d(p)|=Theta_d(p^d). For any unit objective w, choose i with
abs(w_i)>=1/sqrt(d). Each generator near h*e_i has
abs(w.g)>=h*(1/sqrt(d)-delta*sqrt(d))>=h/(2sqrt(d)). The centered
zonotope's support function, one half of sum_g abs(w.g), is therefore
at least c_d*p^(d+1) in every direction. It contains a ball of that
radius and lies in a cube of side C_d*p^(d+1). Coordinate projection
bounds give both surface estimates in (SH1): surface area is at most
the cube bound and at least twice any projected ball area.

For almost every planar projection all projected generators are nonzero
and pairwise nonparallel. Each forbidden condition has measure zero,
as witnessed for a pair by the plane they span. The projected zonogon
has exactly 2*|G_d(p)| vertices. This proves the final part of (SH1).
All arguments hold in every fixed dimension, not only tested dimensions.

Finally, write g=g^+-g^- with both parts in {0,...,p}^d. The known
partition-matroid construction chooses one element of weight g^+ or g^-
for every generator. Its basis-sum polytope is a translate of Z_d(p).
Assigning the same unit vector to the two members of a pair and distinct
unit vectors to different pairs represents it over every field. Thus
representability alone cannot remove this sharpness obstruction. A
residue matroid of a polynomial quotient is a further constraint; no
such realization is asserted for this family.

#### Reproducible controls and effect on the main goal

The new `lattice_shadow_sharpness_control` in the existing checker uses
three-dimensional primitive zonotopes at p=1,2,3. Their 13,49,145 generators
give 96,1248,10944 vertices. An independent great-circle Euler count agrees
with vertices built from strict integer support objectives, checking 12288
vertices in total. Generic shadows have 26,98,290 vertices. The report
saves generator lists, flat multiplicities, vertex-set hashes and exact
support witnesses. It also checks 4160 small axis-box tuples in dimensions
three and four (26 and 1560 primitive cases), four integer left inverses,
and rational exponent balances in dimensions 3--32. These finite cases
support implementation correctness, not the uniform density or asymptotic
proof. The small boxes are below the sufficient asymptotic range.

Run `python3 -B scripts/run_structural_complexity_checks.py` as before.
The separate working manuscript now includes the complete sharpness
argument and certificate scope. The original smooth-number kernels and
completed corner-cut manuscript are unchanged. The unrestricted bound
(LS1), its universal-basis consequence and the unresolved priority flags
remain intact. A stronger unrestricted estimate must exploit information
beyond these two geometric budgets or the existing shadow estimate; this
result does not rule out such an advance. Complete the focused primary
comparison before claiming the resulting paper is new or publishable.

### Twenty-seventh pass: coordinate covers, residual colength and a precise crossover

The preceding goal turn was progress: `a80e26f` added a uniform sharpness
proof for the geometric powers, independent exact controls and the updated
manuscript. This pass retains the unrestricted bound (LS1) and derives a
parameter-sensitive theorem valid over every field, including nonreduced
quotients. It also makes the finite-field comparison quantitative. It does
not replace the full objective with a grid-only problem.

#### A positive lattice-valuation comparison

The available eight-page
[McMullen, Lattice invariant valuations on rational polytopes (1978)](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/McMullen/McMullen11.pdf)
was read through the simple-valuation construction, Theorems 6--8 and the
Section 4 application. Theorem 6, the angle-inversion formula and the
application were checked on the page images. These concern polynomial
expansions and reciprocity for translated lattice-point enumeration.
The later coauthored account, Section 7, printed page 203, explicitly
attributes the lattice versions of simple-valuation polynomiality and
reflection parity to the 1975 and 1977 papers. This is a comparison of
named statements, not a claim to have read the full 1975 article or to
have classified all consequences. No joint vertex/surface/shadow count
is a stated result of these inspected passages. The original article's
remaining scope and the strongest-current-priority assessment stay open.

#### A uniform bound with a separated univariate cover

Let f_i in K[x_i] be monic of degree b_i>=1 over any field, and suppose

    J=(f_1,...,f_d) is contained in I,
    Q=product_i b_i,  n=dim_K K[x]/I,  0<n<Q.

Set Omega=product_i {0,...,b_i-1}, and define

    h=min(n,Q-n),  p_i=min(b_i-1,h-1),
    M=#{u in Omega : product_i(u_i+1)<=h},
    A=sum_i product_(j!=i)(h*p_j+1),
    alpha=d(d-2)/(d^2-d-1),  beta=(d-1)/(d^2-d-1),
    B(h,b)=A^alpha*(M*h^(1/3)+1)^beta.

For each fixed d>=3 we prove

    g(I) = O_d(B(h,b)),
    M <= min(Q,h*(1+log h)^(d-1)).                         (CV1)

The constants do not depend on the field, coefficients or multiplicities.
At n=0 or Q, the initial-ideal count is one. For a cover with all b_i>=n
and h=n, (CV1) recovers the unrestricted estimate (LS1), so the full
scope is retained.

Here is a complete justification of the residual step. The algebra
T=K[x]/J has basis Omega. Extracting the coefficient of x^(b-1) defines
a nondegenerate multiplication pairing: for nonzero a, choose a
componentwise maximal exponent u in its support and multiply by
x^(b-1-u). Every other support exponent has a coordinate below u and
cannot contribute to the top coefficient, even after univariate reductions.
For the ideal R=I/J, its annihilator equals its orthogonal complement,
because R is closed under multiplication. Thus J:I has colength Q-n,
and double annihilation recovers I.

For any global order, in(J)=(x_i^b_i), and the standard initial-colon
inclusion is

    in(J:I) subset in(J):in(I).

If D is the standard staircase of I in Omega, the right side has
standard set {b-1-u : u in Omega minus D}. To check this directly,
v is in the monomial colon exactly when b-1-v is in D: the existence
of a nonstandard u with u+v inside Omega is equivalent, by upward
closure, to b-1-v being nonstandard. Both sides have colength Q-n,
so the inclusion is equality. Reflection is invertible, hence I and
J:I have identical initial-ideal counts over every field. For point
grids the reflected-complement correspondence is established prior art:
[Maruri-Aguilar--Saenz-de-Cabezon--Wynn, Theorem 1, Lemma 1 and Corollary 1 (2013)](https://www.ism.ac.jp/editsec/aism/pdf/10463_2012_Article_390.pdf).
The general box-colon proof was already established in this project;
no new duality is claimed.

Choose I or J:I to have colength h. Its initial staircases are all in
Omega intersect V_h^d. This capped set represents a rank-h residue matroid,
contains every such staircase and has M elements. Restricting the usual
strict greedy exposure preserves the injection of initials into vertices.
All basis-sum coordinates are at most h*p_i. Adding the unit cube gives
surface at most 2*A by coordinate projections, and planar shadows of size
O_d(M*h^(1/3)+1) by the established matroid bound. Substituting these in
(LS4) proves (CV1). The ground-size estimate is the hyperbolic-cross bound.

Every reduced basis of I also has at most Q-n+d polynomials: minimal
outside monomials inside Omega are among its Q-n nonstandard monomials;
those outside Omega must be among the d pure powers x_i^b_i. Combining
with the classical generator estimate gives the further consequence

    |U(I)| = O_d(B(h,b)*min(n^((d-1)/d),Q-n+d)).             (CV2)

This counts elements, not coefficient sizes or enumeration work.

#### Exact comparison with the named prime-field expression

For a set in F_q^d use f_i=x_i^q-x_i. Then Q=q^d and h is the smaller
of the set and complement sizes. Formula (CV1) applies to every prime
power q, with p_i=min(q-1,h-1). To compare only with the scope actually
inspected in the later Stigler--Zhang Theorem 8, take q=p prime. Its
first-range upper-bound expression is

    (p^2*floor(h/p)+(h mod p)^2)^(d(d-1)/(d+1)).

For h=Theta(p^gamma), fixed 1<=gamma<d, its power in p is
(gamma+1)*d(d-1)/(d+1). The power in (CV1) is

    (gamma+1)*(d-1)*alpha+(4*gamma/3)*beta,

with logarithmic factor (1+log p)^L_d. The difference is exactly

    (d-1)*(3d-gamma*(d+4))/(3*(d+1)*(d^2-d-1)).             (CV3)

There is a strict asymptotic power improvement when
1<=gamma<3d/(d+4). At the threshold the powers agree and the logarithm
prevents claiming an improvement; beyond it the earlier expression has
the better power. This is a comparison of upper-bound expressions, not
of achieved maximum counts or unknown finite-size constants. For d=3
and gamma=6/5 the powers are 82/25 and 33/10, a gain of 1/50. The
finite-field comparison is not blanket dominance or a prime-distribution
consequence. Use the smaller applicable proved bound outside this range.

#### Verification and remaining scope

The new `bounded_grid_shadow_control` checks all 258 staircases in four
small boxes. A direct monomial-product colon oracle independently checks
reflection; capped grounds and the residual generator budgets are also
verified. It checks 144 chambers of the existing reduced, repeated-root
and local quotient fixtures in characteristic zero and finite fields.
Exact rational calculations verify (CV3) and its sign in dimensions 3--32.
Finite tests support the implementation; the preceding arguments prove
the uniform statements.

The working manuscript includes (CV1)--(CV3), credits the existing duality
and records the actual 1978 reading scope. Its stale sentence disclaiming
geometric sharpness was corrected to distinguish the proved geometric
sharpness from the unresolved initial-ideal lower bound. The original
kernels and completed corner-cut manuscript are preserved. The main
unrestricted exponent is unchanged and the goal remains active: the
strongest-current-priority audit is still incomplete. No failed retrieval,
new parameter specialization or finite check is used as novelty evidence.

### Twenty-eighth pass: completion audit and the remaining primary-source gate

The preceding turn was progress: `e70ca14` saved the all-field separated-cover
bound, residual generator consequence, finite-field crossover and exact
controls. This pass checks the original objective against the current proof
package. It does not declare the goal complete from a working manuscript,
passing tests or the comparison with one historical bound.

| Original requirement | Current evidence | Completion assessment |
| --- | --- | --- |
| Broad uniform theorem | (LS1), proved for every colength-n ideal over every field and every fixed d>=3; the separate manuscript's main theorem and geometric/matroid lemmas supply the argument. | Proof supplied, with established inputs explicitly credited. Finite checks are supporting evidence. |
| Advance over the strongest relevant primary baseline | The exponent is strictly smaller than the inspected Babson--Onn--Thomas exponent; later named bounds and their parameter regimes have scoped comparisons. Originality and the strongest-current comparison are still explicitly open. | Not established. This prevents completion of the full objective. |
| Reusable dependency or geometric mechanism | (LS3)--(LS8), the primitive-facet incidence budget, cone boundary estimate and planar matroid reduction; (SH1)--(SH2) delimit geometric sharpness. | Uniform structural argument supplied; priority of its quantitative application remains open. |
| Initial upper and infinite lower attempts | (RS1)--(RS3) and their rank-delay counterexamples; the symbolic nonsocle family (NG1)--(NG2), including its global O(n) obstruction. | The staged proof attempts and their limitations are recorded. The upper-bound route does not require pretending the tested lower family is extremal. |
| Further consequence | (LS11), the integer-weight projection estimate, geometric sharpness, and (CV1)--(CV3) with the residual generator budget. | Precise derivations supplied. No enumeration speedup or prime-distribution claim follows. |
| Reproducible certificates | Existing dependency-free research command and report, with source hashes, independent finite oracles and exact witnesses. | Sources and all prior blocks checked; computations support the proofs and do not prove novelty. |
| Honest separate manuscript and preservation | Current 18-page working draft, matching recorded source/PDF hashes, prior compilation and page review; the seven protected kernel/manuscript files match their saved baseline. | Artifact present and internally checked. Priority remains a stated limit; no independent referee or journal decision is imposed as an additional requirement. |

The source access gap has recurred through the three consecutive mathematical
passes ending in `a80e26f`, `e70ca14` and this audit. The precise remaining
comparisons are:

1. [McMullen, Non-linear angle-sum relations for polyhedral cones and
   polytopes (1975), pp. 247--261](https://doi.org/10.1017/S0305004100051665),
   particularly its lattice application. The 1978 and 1983 readings establish
   the named angle/valuation comparisons already recorded; they are not
   asserted to reproduce every statement of this original article.
2. [Stigler--Zhang, The Number of Groebner Bases in Finite Fields (2020),
   pp. 139--158](https://doi.org/10.1007/978-3-030-42687-3_9), with the later
   [public deposit](https://par.nsf.gov/servlets/purl/10177092). Its indexed
   Theorem 8 and displayed proof support the numerical comparisons above,
   but the full later source has not been retrieved. The arXiv history was
   rechecked and lists only v1; a later arXiv version cannot be assumed.

Current public requests still did not retrieve either complete text. This is
an access observation, not evidence of novelty or mathematical non-equivalence.
There is no live download or background research process to wait on.

The proof package and all exact artifacts are retained. Adding another
restricted consequence would not settle this remaining requirement. The
appropriate status is blocked on the primary-source priority comparison,
with the original full objective preserved. A complete copy of the named
sources, or a primary replacement that actually resolves their relevant
content, permits this audit to resume. No purchase, author contact, external
message or submission is authorized or has been made.
### ML follow-up: preserve decoupled erase and write in the smaller state

The next comparison should include newer recurrent architectures, rather
than stopping at the scalar-gated delta rule. In particular,
[Gated DeltaNet-2, Section 3.1, equations (8)--(10)](https://arxiv.org/html/2605.22791v1#S3.SS1)
separates an erase gate on the key axis from a write gate on the value
axis. We inspected these definitions and the adjacent derivation; this is
not a replication of its training or GPU results. Its update is

    Sbar=D_alpha*S,
    Snew=Sbar+f*((w*v)^T-(b*f)^T*Sbar).                  (MLE1)

Here * between vectors denotes elementwise multiplication, D_alpha is
diagonal, and f is the key feature. The recurrence is published work.
The following observation extends our compression compatibility check;
it does not establish a new architecture or a novelty claim.

Fix the positive grouping weights from (MLP1)--(MLP2), so f=U*phi and
U^T*U=I. Suppose the decay is constant within each equal-exponent group,
as required in (MLG1), and S=U*T. Set

    bbar_a=sum_(nu: b(nu)=a) (w_nu/c_a)*b_nu.

The supports of different columns of U are disjoint. Therefore
U^T*diag(b)*U=diag(bbar), even when b varies inside a group. Substitution
in (MLE1) gives the exact smaller recurrence

    Tbar=diag(alpha_collected)*T,
    Tnew=Tbar+phi*((w*v)^T-(bbar*phi)^T*Tbar).             (MLE2)

It follows inductively from zero initialization that every full state
equals U times its smaller state, and every query read agrees. Values
and their write gates can be arbitrary. The feature weights and grouping
map must remain fixed within the sequence. Gates may vary with time; if
they depend on past reads, their evaluation must receive the same reads.

This distinguishes two operations. Untied diagonal forgetting can move
the stored state outside im(U), as the earlier two-step counterexample
shows. Untied erasing only changes the scalar read before the update
writes along f, which is in im(U). Consequently arbitrary generic erase
gates have an exact weighted-average reduction. For gates in [0,1], the
averaged gates remain in [0,1]; every averaged gate vector is attainable
by tying its generic entries. This is a statement about the recurrence's
gate values, not equivalence of arbitrary gate-encoder parameterizations.

The independent rational control `algebraic_erase_write_controls` checks
four one- and two-dimensional polynomial fixtures at powers two and four,
including deliberately different erase gates within duplicate groups.
It compares 24 state steps and 288 prefix-query coordinates exactly.
It also retains a stability limit: a unit key (3/5,4/5) with erase gate
(1,0) maps (-1/2,1) to (-8/25,31/25). Squared norm increases from 5/4
to 41/25. Thus the previous nonexpansiveness guarantee for a symmetric
scalar delta correction cannot be transferred to arbitrary channel-wise
erase gates. Decoupled erasing and writing likewise does not automatically
give exact current-key writes.

The implementation opportunity remains a small feature encoder attached
to an existing recurrent update. Generate M collected erase/decay gates
directly and keep the state at M*v scalars. Computing F generic gates
and averaging them still costs at least O(F); an arbitrary nonlinear
gate encoder cannot be folded by simply averaging its weights. The
M=225 versus F=2080 example retains the previous 9.24-fold state-count
ratio for the same constrained feature family, with no measured speedup.
An ordinary 225-feature Gated DeltaNet-2 has the same state size and is
the stronger direct learned-feature comparison.

Another recent comparison is
[Complex KDA, abstract and Section 3](https://arxiv.org/html/2609.24797v1#S3).
Its signed gates and delta coefficients in [0,2] already provide stable
rotations and state-tracking results. The inspected statements establish
a relevant baseline, not an audited proof comparison. Our scalar-gated
grouping identity also permits signed collected gates with absolute value
at most one; this elementary compatibility alone adds no new rotation,
state-tracking or learning guarantee over that work.

Proceed to a trained, matched-state associative-recall comparison before
custom GPU development. Include native learned-feature KDA/Gated
DeltaNet-2, the structured scalar and channel-gated variants, and the
decoupled variant (MLE2). Use held-out keys, overwrite/interference and
longer contexts across seeds; charge encoders, gates, coefficient updates,
backward and state traffic. Compare with Power Attention in a subsequent
sequence-model test. Select based on measured quality at equal memory
and time, not feature count alone. No trained-model or GPU benchmark was
run here. The prior untrained poor-recall diagnostic remains a reason to
test this, not evidence of an improvement. Prime-specific theorems and
initial-ideal counts still supply no direct acceleration of dense weights.

### Learned associative recall: the structured block loses the matched-budget test

The user requested a learned comparison against existing recurrent memory
updates. We ran the experiment rather than treating feature-count savings
as a model-quality result. The tested two-latent structured block did not
show an advantage. Preserve its exact factorization, but park this version
as a proposed general ML replacement; these results do not justify a
custom GPU kernel for it.

#### Protocol and the comparison class

`scripts/run_learned_recall.py` implements the published KDA and Gated
DeltaNet-2 recurrences with a common small encoder interface, each with
either unconstrained learned features or the collected quadratic features
of (MLP1), s=2 and m=7. These are one-layer memory-block references, not
complete Kimi Linear/Gated DeltaNet-2 models or their fused kernels. No
official Zoology MQAR or language-model benchmark was replicated.

All four blocks have exactly 225*16=3600 float32 recurrent-state scalars
(14400 bytes). They share a cap of 32768 active trainable parameters,
with actual counts 32331--32686, and have no persistent basis/index
buffers. Hidden widths use the largest size fitting that cap; no unused
padding parameters are counted. Training activations and peak process RAM
are not matched. Thus memory matching concerns recurrent state and the
closely matched parameter storage, not total training-memory consumption.

Each architecture receives the same two-learning-rate tuning budget on a
separate seed, then 30 process CPU seconds per main run at seeds 11,23,37.
The budget includes batch generation, feature/coefficient construction,
gate encoders, forward/backward, clipping and AdamW. Evaluation and
checkpoint I/O are excluded. Actual main budgets are 30.000--30.053
seconds; the excess is the final step. One CPU thread is used throughout.
The local host has 16 logical x86_64 CPUs and 16 GiB RAM; the environment
is Python 3.12.14, PyTorch 2.2.2 and NumPy 1.26.4. The PyTorch wheel is
prebuilt; no custom compiler flags or GPU kernels are involved. These
optional dependencies do not change the core toolkit's dependencies.

Episodes contain fresh normalized Gaussian keys and independent uniform
16-way values. Models receive explicit key/value writes, then shuffled
read-only queries for those keys. Query/key encoders are shared. Training
mixes 4/8/16 target associations, overwrites and unrelated writes. Test
episodes are separate from tuning and training; they include 32/64-key
contexts, overwrites and 32 unrelated writes. Key dimensions two and eight
are tested separately. Chance classification accuracy is 6.25%.

A preflight linear feature head would have limited the native query span
to hidden_width+1. We discarded that tuning before any main test run and
restarted with a final SiLU activation. An independent rank control checks
that the native feature span exceeds that artificial bound. No final test
data were used for this correction. All subsequent choices and the
superiority criterion were saved in `protocol.json` before the main tests.

#### Results, including equal-exposure controls

The primary score is the mean of 16-key recall and 16-key overwrite
accuracy. Percentages below are means over the three training seeds.

| Memory block | Two-dimensional keys | Eight-dimensional keys | 64-key recall, eight-dimensional keys |
| --- | --- | --- | --- |
| Learned-feature KDA reference | 75.98% | 98.14% | 84.65% |
| Structured KDA | 67.77% | 70.64% | 37.01% |
| Learned-feature Gated DeltaNet-2 reference | 70.43% | 97.24% | 82.19% |
| Structured Gated DeltaNet-2 | 68.58% | 71.12% | 37.31% |

For eight-dimensional keys, structured-minus-native primary-score gaps
are -27.50 percentage points for KDA and -26.12 for Gated DeltaNet-2.
Their descriptive three-seed 95% paired t intervals are respectively
[-28.69,-26.32] and [-28.64,-23.61] points. The two-dimensional gaps
are -8.21 and -1.85 points; the latter interval includes zero. These are
small-pilot descriptive intervals, not multiple-comparison-adjusted
confirmatory significance claims. None of the four comparisons passes
the predeclared structured-superiority criterion.

The equal-example checkpoint control also favors native learned features
in every seed and recurrence family. For eight-dimensional keys all
models are evaluated after exactly 1000 updates; the mean structured
gaps are -28.58 and -25.62 points. In the two-dimensional condition the
largest predeclared common milestones are 1000,1000,500 by seed, with
mean paired gaps -7.03 and -2.27 points. Thus fewer training examples
within the CPU budget do not explain away the negative result.

Complete reference forwards at batch one, with 16 writes and 16 reads,
also cost more for the structured variants on this CPU: for eight-dimensional
keys, median timings are 1.173 versus 1.606 ms for native/structured KDA
and 1.483 versus 1.610 ms for native/structured Gated DeltaNet-2. These
include the encoders and feature preparation. They are not production
decoding or GPU-throughput comparisons.

#### Verification, artifacts and decision

The saved controls check all four recurrences against independent dense
transition matrices, full-network directional gradients, overwrite labels
against a latest-match oracle, and an exact orthogonal-address positive
control. Additional checks of the actual 225-feature two-dimensional
encoder compare against the direct 64-coordinate squared kernel, compare
parameter gradients through those independent expressions, and check
full-network directional gradients. Maximum pair-kernel discrepancy in
double precision is below 1.6e-15 in these fixtures. These checks support
implementation correctness, not universal numerical conditioning.

`scripts/analyze_learned_recall.py` independently verifies all 24 result
cells, memory caps, CPU budgets, training exposure counts, test query
counts, aggregate scores and paired intervals. The protocol, every tuning
trial and main run, equal-exposure scores, raw stdout, source hashes and
machine-readable tables are under `results/benchmarks/learned_recall_v2/`.
Private checkpoints and the optional environment remain ignored locally.
Exact CLI arguments are in `protocol.json`; the audit command is

    python3 scripts/analyze_learned_recall.py results/benchmarks/learned_recall_v2

Add `--check-model` inside the optional PyTorch environment to repeat the
actual-configuration kernel and gradient checks. Budgeted reruns can
complete different numbers of updates on different hosts; raw rows save
the exact count, and the exposed `train(..., fixed_steps=count)` function
permits deterministic exposure replay in the recorded environment.

This experiment finds no practical or fundamental advantage for this
particular structured block under its tested budgets. It does not rule
out every structured feature family, training regime or hardware target.
The exact 2080-to-225 compression remains correct for its constrained
model, while a learned 225-feature model already uses the same recurrent
state and learns better in this comparison. A new ML direction needs a
stronger representation or task-specific mechanism before another kernel
development program. The separate mathematics goal and its unresolved
primary-source priority comparison remain unchanged.

### Full-dimensional factorization: a precise barrier and a higher proof target

This follow-up answers whether the two-latent restriction can be removed
before any more model training. No new neural experiments were run. The
previous learned-recall protocol, implementations and raw results are
unchanged. The conclusions below are established-literature consequences
and elementary controls, not a new architecture or a novelty claim.

#### What exact compression does, and what it does not imply

Let A be a finite set of nonnegative integer exponent vectors, containing
0, with D=|A| and affine dimension r>=1. Write
z_a(x)=sqrt(theta_a)*x^a, where every theta_a>0. On an open subset of the
positive orthant, consider the normalized kernel

    K_p(x,y) = (z(x).z(y)/(||z(x)|| ||z(y)||))^p,  p>=1 integer.

Expanding the numerator and collecting equal exponents gives one feature
for each element of the sumset pA. All collected coefficients are
positive. Distinct monomials are linearly independent on an open set;
division by the same nonzero normalization preserves that independence.
An invertible evaluation matrix at |pA| suitable points then makes the
kernel Gram matrix have rank |pA|. Thus its minimal exact separation rank
is M=|pA|, not merely an upper bound from one implementation. This is a
special case of the kernel-rank framework in
[Altschuler--Parrilo, Proposition 3.6 and Corollary 3.8](https://arxiv.org/html/2106.02755v1#S3.SS1).
It compares exact representations of this same kernel. It does not show
superiority to an unrelated learned M-feature kernel.

The classical Freiman dimension bound gives

    |2A| >= (r+1)*D - r*(r+1)/2.

For general p the sharp bound is

    |pA| >= binom(p+r-1,r)*D - (p-1)*binom(p+r-1,r-1).

See [Lopez Bruch--Jing--Mudgal, Theorem 1.5 and Section 7](https://arxiv.org/html/2409.05638v1).
The statement, sharp construction and Section 7 argument were inspected;
the underlying cited compression lemmas were not independently re-proved.
These are published bounds, not results originating in this project.

For D=64, the quadratic consequences are:

| Independent exponent directions r | Minimum exact quadratic width |
| --- | --- |
| 2 | 189 |
| 3 | 250 |
| 8 | 540 |
| 32 | 1584 |
| 63 | 2080 |

Since this lower bound increases for 1<=r<=D-1, D=64 and M=225 force
r<=2. The tested support {0,...,7}^2 has |2A|=225. Within this family we
cannot keep those two widths and simply add independent continuous
directions. This is not a lower bound for arbitrary neural architectures,
approximate features, discrete input sets or other kernel families.

There is a concrete escape if we change the base width. The sharp support

    A = {0,e_2,...,e_r} union {e_1,2e_1,...,(D-r)e_1}

has D elements and attains the bound. For r=8,D=29 it gives exactly 225
quadratic features. It includes every coordinate monomial: the ratio of
the collected e_i feature to the constant feature is a known positive
constant times x_i. Hence the normalized feature map retains all eight
input coordinates on the positive orthant. But degrees up to 21 are
concentrated on one axis. Injectivity alone guarantees neither good
conditioning nor balanced interactions or successful training. Keeping
D=64 instead allows eight directions at width 540, versus 2080 for the
generic symmetric square of 64 independent coordinates. That is a
conditional exact feature-count saving, not an advantage over a freely
learned 540-feature memory.

#### Exact rank can coexist with failed recall

A two-write example makes the missing implication explicit. Take
A={0,1}, theta=1 and p=2, so the unit-norm features have kernel

    K(x,y) = (1+xy)^2 / ((1+x^2)*(1+y^2)).

Write the distinct keys x_1=-1/10 and x_2=1/10 with scalar values +1
and -1. Their correlation is mu=9801/10201. The two-key Gram matrix has
determinant 1-mu^2>0 and rank two, so the features have enough linear
capacity to interpolate these two associations.

Nevertheless, the ordinary ungated, unit-step delta rule

    S <- S + f(x)*(v - f(x)^T*S)

starting at zero gives, after the two writes,

    S = f(x_1) - (1+mu)*f(x_2),
    read(x_1) = 1-mu-mu^2 = -91979201/104060401,
    read(x_2) = -1.

The first read has the wrong sign, about -0.884 instead of +1, in exact
arithmetic. As the two distinct keys approach zero, mu approaches one
and the first read approaches -1. No dimensional information was lost
in this example. Exact compression faithfully preserves a kernel whose
geometry permits destructive interference; full rank is not a stability
or single-pass recall guarantee. This supplies a failure mechanism,
not a causal diagnosis of the trained models' score gaps.

Two apparent remedies also need precise limits. A sum of kernels on
disjoint coordinate blocks retains coordinates but gives additive read
functions across blocks. On the four points (x,y) in {-1,1}^2, even
arbitrary one-coordinate functions span only {1,x,y}; they cannot equal
the binding xy. The corresponding evaluation ranks are three and four.
Taking products restores interactions but generally increases the
expanded feature count. Rank-one individual writes also need not give a
low-rank sum of many writes.

Nor can an arbitrary exact state sketch automatically retain independent
channel gates. If a linear compression E represents all states and every
diagonal gate can be updated from ES alone, ker(E) must be invariant under
each coordinate projector. If E*e_i is nonzero for every original
coordinate, any nonzero vector in ker(E) would imply E*e_i=0 for one of
its nonzero coordinates, a contradiction. Therefore E is injective and
cannot reduce dimension. For example E=(1,1) identifies states (1,0)
and (0,1), whereas gate diag(1,0) separates them. The scope is exact
linear compression with all those gates and states, not all recurrent
models. Tied gates, restricted reachable states or approximation change
the question.

#### Existing architectures that avoid two scalar latents

The clearest exact alternative is
[Large Memory Layers with Product Keys](https://arxiv.org/html/1907.05242v2),
Sections 3.1--3.2. Its two query parts are vectors, not scalar latents.
For query (q_1,q_2), two codebooks of n vectors give N=n^2 addresses with
score s_ij=q_1.a_i+q_2.b_j. A globally top-k address must use a top-k
entry in each codebook: otherwise keeping the other index fixed already
provides k better addresses. Compatible deterministic tie-breaking gives
the same conclusion for ties. Thus full-dimensional queries can have
exact factorized access to this structured address family. Computing
subscores costs O(d*sqrt(N)); at most k^2 candidate sums remain, plus
selection and O(k*v) costs to read v-dimensional values. The N*v value
table still occupies memory. This advantage is over exhaustive scoring
of the same Cartesian-product codebook, not arbitrary learned addressing
or equal-memory neural models. The defining construction, top-k proof
and cost argument were inspected; experiments were not replicated.

[Sparse Delta Memory](https://arxiv.org/html/2607.07386v1), Sections 2--3,
already combines product-key indexing with sparse delta updates to an
explicit memory table. Its selected addresses permit reads and writes
that leave other slots untouched. Its equal-FLOP framing allows increased
state memory; it is not an equal-memory superiority theorem. The
definition and update rule were inspected, not its code or experimental
claims. Consequently, combining high-dimensional factorized access with
delta updates is already occupied prior art.

Another route keeps rich polynomial kernels but permits approximation.
[PolySketchFormer](https://arxiv.org/html/2310.01655v3), Theorem 1.1 and
Theorems 2.2/2.4, gives polynomial-kernel sketches with matrix-error
guarantees and a positivity construction. These do not require a
two-dimensional encoder. The displayed statements and proof arguments
were inspected. A matrix approximation guarantee alone does not prove
uniform error control through long selective-overwrite trajectories or a
trained-model advantage. The original
[high-degree polynomial sketch paper](https://arxiv.org/abs/1909.01410)
was checked at abstract level, not independently proof-audited.
[Incremental sketching for online kernel learning](https://proceedings.mlr.press/v97/zhang19h.html)
already studies online approximation and regret; only its primary abstract
was inspected here. These prevent treating a generic sketch-plus-online-
update proposal as a demonstrated gap.

#### Higher-impact target and decision gate

The useful objective is a theorem relating memory size, active access
cost and retention for a model that preserves high-dimensional queries
and cross-coordinate binding. An exact factorization theorem for a
redundant implementation is insufficient. Before another model run:

1. Specify the query/address family and allowed overwrite and forgetting
   operations. Count all stored values, indices, feature parameters,
   metadata and required precision. Distinguish an exact codebook search
   saving from compression of the value table.
2. Prove either a new access/retention guarantee beyond product-key and
   sparse-delta memory, or a sharp approximation/capacity tradeoff for
   full recurrent trajectories beyond existing kernel-sketch/online
   learning results. The guarantee must support cross-coordinate
   associations and state its dependence on horizon, overlap and gates.
3. Establish the primary-source boundary for that exact theorem. Compare
   against the strongest applicable architecture, not just its dense
   implementation. If it follows from existing work, park that proposal.
4. Only then use experiments to test the theorem's practical regime.
   Explicit exact-key recall also admits ordinary finite dictionary
   baselines when keys fit in memory; it alone cannot establish a
   fundamental architectural advantage. Meaningful extensions must
   specify noisy/semantic queries or other binding demands and charge
   the corresponding memory.

No project-specific architectural advantage or publishable novelty has
been established by this review. The dimension obstruction and the
interference control explain why the earlier exact mathematics was
insufficient to imply one. The eight-direction construction is an
algebraic feasibility witness, not a recommended new training run.
The strongest current direction is full-dimensional structured access
with a retention guarantee; its first deliverable must be a theorem and
an audited comparison, not another benchmark sweep. This does not replace
or resolve the separate blocked mathematical priority audit.

Reproduce the finite certificates with

    python3 scripts/run_feature_dimension_checks.py

The output is `results/research/full_dimensional_factorization_checks.json`.
It records 45 sharp higher-sumset fixtures, the five dimension-frontier
cases, the eight-direction witness, rational two-write interference,
additive-binding and channel-gate counterexamples, and 57 exact product-
key top-k comparisons including ties. These computations check concrete
certificates; the uniform claims rest on the arguments and cited theorems
above. They do not test model quality or establish novelty.

### Associative-memory goal: first proof iteration

The user has explicitly requested a new active goal: obtain a substantive
provable advantage in total storage, active read/write cost or retention,
with high-dimensional queries and cross-coordinate binding. The complete
criteria and current milestones are in
`results/research/associative_memory_goal.json`. This supersedes the old
blocked algebraic-complexity objective as the active task; it does not
claim that the old objective was achieved. Previous manuscripts and
results remain intact. No new neural training has been run.

The first iteration establishes the model and three rigorous baselines.
They are useful proof tools and feasibility tests, not a completed novel
theorem or a claim of superiority over existing methods.

#### Model and information accounting

Fix n distinguishable keys. A scalar value at each key can be replaced
independently; an exact overwrite changes one requested value and leaves
every other stored-key read unchanged. For an alphabet of size a, a
zero-error memory needs at least n*log2(a) mutable bits: the a^n possible
latest-value vectors must produce different distinguishable states.
This elementary bound includes any mutable parameter/index storage. It
does not apply unchanged to distributional, lossy or correlated values.

For the initial algebraic model, stored state s has m scalar cells and
fixed-key reads are Fs, where F is an n-by-m real feature matrix. The
feature encoder and its parameters are fixed during an overwrite. Cost
t means the number of changed scalar cells, a lower bound on physical
write access, rather than a claim about total running time. Values with
several independent components require corresponding storage and work.
Computing an address, constructing a feature, reading the old value and
updating metadata must also be charged. Adaptive/nonlinear read schemes
or changes to the encoder need separate analysis.

#### A sharp locality/interference identity

Choose a target key i and a set T of cells that an overwrite may change.
Let C=F[:,T], let P_T be the orthogonal projector onto range(C), and set
ell_i(T)=(P_T)[i,i]. Among all such changes producing a unit change at
key i, the minimum total squared change at the other keys is exactly

    (1-ell_i(T))/ell_i(T),          if ell_i(T)>0.

If ell_i(T)=0, a unit target change is impossible. To prove the formula,
write the read change as r in range(C), with r_i=1. Then

    1 = <P_T e_i,r> <= sqrt(ell_i(T))*||r||_2.

Thus sum_{j!=i} r_j^2=||r||_2^2-1 is at least 1/ell_i(T)-1. Equality
is attained by r=P_T e_i/ell_i(T). A target change delta multiplies the
optimal interference energy by delta^2. Minimizing over allowed sets
|T|<=t equivalently maximizes this leverage score. This is an elementary
projection argument using classical statistical leverage, not a newly
introduced invariant. See the projection-diagonal definition in
[Drineas--Mahoney, Section 3](https://arxiv.org/html/1005.3097v1#S3), and
the related underconstrained-regression discussion in
[Drineas et al., Section 6.2](https://arxiv.org/html/1109.3843v2#S6.SS2).
The definition and relevant displayed identities were inspected; a full
audit of their randomized algorithms is not claimed.

One exact consequence is particularly restrictive. Suppose deleting
any row of F leaves every set of at most n-1 columns independent. Then
every nonzero exact isolated overwrite changes at least n cells,
irrespective of how many redundant columns m>=n are available. Indeed,
for a smaller support T the unchanged reads require

    F[all rows except i,T]*delta_s[T] = 0.

Column independence forces delta_s[T]=0, contradicting the requested
change. The update algorithm itself can be nonlinear or depend on the
current state; its before/after difference still obeys this equation.
The hypothesis holds for generic evaluations of distinct monomials on
an open domain: each relevant minor is a nonzero polynomial in the
evaluation points, and finitely many such nonvanishing conditions can
hold simultaneously. It is stronger than merely having full row rank.
Special sparse address matrices need not satisfy it.

An exact approximate-retention frontier is available for a Hadamard
read matrix H of order n. Its columns are orthogonal with norm sqrt(n),
so every t-column projector has ell_i(T)=t/n for every i. Hence

    minimum off-target squared interference = n/t - 1.

This is attained by delta_s[j]=H[i,j]/t on the selected columns.
Orthogonal rows alone therefore do not make a sparse coordinate update
safe. Changing to a private-address basis can eliminate this cost, but
the computation/storage needed for the resulting query representation
must be included. Neither statement lower-bounds all memory algorithms.

In particular, storing a linear encoding By of latest values and reading
with F gives exact stored-key retention precisely when FB=I_n. More
general specified linear queries W require W=RB. Column sparsity of B
is update locality and row sparsity of R is query locality. This
correspondence already belongs to dynamic linear-data-structure and
linear-circuit research; the definition and proof in
[Brody--Larsen, Section 3](https://arxiv.org/html/1208.2846v2#S3)
were inspected. That preprint uses different lemma numbering from its
journal version. It would be incorrect to claim a new theorem merely
by translating this correspondence into neural-memory terminology.

#### A constructive baseline with horizon-independent retention

Fixed sparse features also give a positive control. Let G=FF^T and
choose alpha>0 such that

    A=I-alpha*G,       ||A||_infinity = rho < 1.

For integer L>=0 define

    P_L=alpha*sum_{j=0}^L A^j,       B_L=F^T*P_L.

Store the current latest-value matrix V and s=B_L*V. An overwrite
V_i<-v changes s by B_L[:,i]*(v-V_i). This preserves the encoding
identity exactly. Since G*P_L=I-A^(L+1), the stored-key reads obey

    F*s - V = -A^(L+1)*V,
    maximum absolute entry error <= rho^(L+1)*max_abs_entry(V).

The bound depends on the current value range, not the number or order
of preceding overwrites. This distinction requires storing the latest
V_i separately: subtracting an approximate read instead would not
maintain the stated identity. Existence of alpha is a genuine hypothesis;
the result does not assert contraction for arbitrary learned features.

This construction also specifies a richer query function. For any query
feature z, compare z^T*s with the exact minimum-norm kernel interpolant
z^T*F^T*G^(-1)*V. The difference is at most

    ||z^T*F^T||_1 * alpha*rho^(L+1)/(1-rho) * max_abs_entry(V)

in each value component. No reduction of z to two coordinates is used,
but the bound alone does not construct a suitable high-dimensional
encoder, establish semantic robustness, or support arbitrary binding
tasks. Those requirements remain part of the goal.

Here is the cost ledger. Suppose every stored-key feature has at most k
nonzeros, and the off-diagonal support graph of G has maximum degree
Delta. Let b_L bound the number of keys within graph distance L of any
key. A column of P_L is supported in this ball, so a column of B_L has
at most c=min(m,k*b_L) nonzeros. With an explicitly stored sparse B_L:

- a replacement touches at most c state rows plus the latest-value row;
- a query with k_q nonzero features reads k_q state rows, after paying
  the encoder and routing cost;
- state and latest-value storage together use (m+n)*value_width scalars;
- coefficient/index storage can use n*c entries, in addition to keys,
  encoder parameters, F/graph data and any query-routing tables;
- given the sparse Gram graph, a direct column construction takes at
  most O(n*L*(Delta+1)*b_L+n*k*b_L) scalar operations, with an additional
  O(n*k) term for L=0. Forming/validating the Gram graph is separate.

For a general bounded-degree graph, b_L can grow exponentially with L.
The cost does not disappear because the coefficients are fixed. Exact
matrix operations over arbitrary reals are also not a bit-complexity
or hardware guarantee.

Finite precision can be made explicit. Round nonzero entries of B_L to
dyadic coefficients with absolute error at most eta, retaining zero
entries. Let h<=n bound the number of stored coefficients in a state
row. Then the retention bound becomes

    [rho^(L+1) + ||F||_infinity*h*eta] * max_abs_entry(V).

Represent bounded dyadic values by integers, and maintain the resulting
integer state sums exactly. If their integer magnitudes are at most R
and coefficient numerators at most C, each accumulator is bounded by
n*C*R at every time. Signed word length O(log(n*C*R+1)) therefore
suffices independently of the overwrite horizon. Every coefficient,
value and accumulator bit is charged. Repeated floating-point additions
without this or another error-control mechanism do not inherit the
horizon-independent bound automatically.

This baseline follows geometric-series inversion and fixed linear
encoding. Sparse inverse locality and sparse matrix-function
approximation are established subjects, including
[Demko--Moss--Smith](https://www.ams.org/mcom/1984-43-168/S0025-5718-1984-0758197-9/S0025-5718-1984-0758197-9.pdf)
and [Benzi--Razouk](https://etna.ricam.oeaw.ac.at/vol.28.2007-2008/pp16-39.dir/pp16-39.pdf).
Benzi--Razouk's Theorem 3.4 and its proof were inspected in the full
paper, including graph-distance support of matrix powers and the
eigenvector-conditioning factor; Section 4.1 gives polynomial-approximation
algorithms. This establishes a direct classical baseline for the locality
mechanism. Demko--Moss--Smith was checked at searchable-statement level;
its full proof remains to audit. Our elementary proof is self-contained,
but establishes neither priority nor an advantage over those methods.

#### Literature boundary and next proof test

Regression-based retention is not new. The model and online inverse
updates in [Flexible Kernel Memory](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0010955),
Sections 2.1/2.4, already use kernel interpolation. Formula images were
not independently transcribed, so this is not a full proof audit.
[Kalman Delta Networks](https://arxiv.org/html/2609.07816v1) and
[Memory by Design](https://arxiv.org/html/2605.31163v1) explicitly track
memory uncertainty or approximate it. Their relevant model definitions
and covariance recurrences were inspected; experimental claims are not
adopted. Adding a covariance state without a new cost/retention guarantee
does not establish a contribution.

Two further primary sources identify strong comparison targets.
[Adaptive exact graphical inference](https://jmlr.org/papers/v12/sumer11a.html)
was inspected beyond its abstract: Section 1.2 and Theorems 4/5 give
preprocessing O(n*d^(3w)), marginal queries O(d^(2w)*log n), and a
single factor/model update O(d^(3w)*log n), for finite domain size d
and the specified elimination width w. These are not dimension-free
bounds, and a continuous Gaussian specialization is not established by
that statement alone. The relevant query/update proof arguments were
inspected; the entire paper was not independently re-proved.
[Online structured kernel interpolation](https://proceedings.mlr.press/v130/stanton21a.html)
was checked at primary-abstract level. Its inducing-grid, precision and
approximation costs remain to audit before claiming a gap.

A [preliminary sparse-factorization note](https://lucablanchi.it/presentation-theory/local-linear-presentations-sparse-factorization)
also explicitly states the workload factorization correspondence. It
labels itself AI-assisted and not fully verified. Its elementary overlap
is a priority flag; it is not an authoritative theorem audit or a reliable
basis for claiming an open problem.

The next constructive test is insertion/deletion or motion of keys,
because the fixed-key baseline leaves a significant issue unresolved:
changing F changes G and the entire encoding B_L. Reusing the old
coefficients silently invalidates the bound. We will derive the affected
supports and residual through the inverse/resolvent identity, charge
updates to both values and coefficient metadata, and compare with dynamic
inference and online kernel-regression guarantees. A new result must
improve an applicable published bound for a clearly specified richer
query workload; it cannot simply outperform a dense implementation or
ordinary delta writes. Exact-key dictionary lookup remains a baseline
where it satisfies the whole query specification.

The goal remains active. None of the bounds in this first iteration is
being presented as a publishable novelty or as completing it. No neural
experiment should be launched from these controls alone.

Run the exact certificates with

    python3 scripts/run_retention_tradeoff_checks.py

The report is `results/research/retention_tradeoff_checks.json`. It checks
row-deleted monomial evaluation minors, sharp constrained-leverage
optimizers, Hadamard support frontiers and local polynomial encodings.
The fixed-feature examples use F=I+(1/16)*adjacency(cycle), without private
feature coordinates. They independently compare sparse replacement
updates with the full current-value encoding and include dyadic integer
accumulator controls. Finite streams and sampled assignments verify the
implementation; the horizon-uniform statement follows from the encoding
identity, not from extrapolating a tested stream length.

#### Second proof pass: local repair when keys change

The fixed-key construction extends algebraically to local key edits.
This is a candidate comparison baseline, with novelty unresolved; it is
not an established improvement over dynamic kernel/inference algorithms.
Keep alpha and L fixed in advance, and suppose every permitted feature
matrix F_t satisfies the same contraction bound rho<1. A local edit
changes one row i of F. Then G'=F'F'^T differs from G only in row/column
i. Let H be the union of the old and new off-diagonal Gram support
graphs, and C the keys within H-distance L of i.

For A=I-alpha*G and A'=I-alpha*G', the noncommutative identity

    A'^j-A^j = sum_{r=0}^{j-1} A'^r*(A'-A)*A^(j-1-r)

shows that every changed column of P_L is indexed in C. In each term,
the changed entry lies at i or an incident Gram edge, and no product
can reach a column farther than j<=L from i. Moreover,

    B'_L-B_L = (F'-F)^T*P'_L + F^T*(P'_L-P_L).

The first term also has columns only in C, since the changed feature row
is i. Hence B'_L[:,j]=B_L[:,j] whenever j is outside C. Recompute the
columns in C using their new sparse polynomial expansions. With old
values V and a concurrent value change delta at slot i, update

    s <- s + sum_{j in C}(B'_L[:,j]-B_L[:,j])*V_j
             + B'_L[:,i]*delta.

This gives s=B'_L*V' exactly after the edit. The current retention bound
therefore remains rho^(L+1)*max_abs_entry(V') for arbitrarily many such
edits. Quantizing each current encoding column to a canonical dyadic
value and updating its integer contribution gives the same fixed
coefficient-error bound as before, without history-dependent drift.
One must change the coefficients as well as the value state: a local
residual update that silently leaves B_L unchanged does not prove this.

To include insertions and deletions, reserve a maximum of N slots and
mark unused slots inactive. Set their feature rows and latest values to
zero, and use the auxiliary matrix

    G_hat=F*F^T+alpha^(-1)*diag(inactive slots).

Inactive rows are isolated, their A entries are zero, and their encoding
columns F^T*P_L are zero. They require occupancy metadata rather than
physical feature coordinates. On active rows, G_hat equals the active
Gram matrix and the same retention identity applies. Activating or
deactivating a slot is still a single star-supported matrix edit. This
is a bounded-capacity construction; resizing is not free.

If the union overlap graph has L-balls of at most b_L vertices and every
feature row has at most k entries, a direct repair recomputes at most
b_L encoding columns. A conservative scalar-work bound is

    O(b_L^2*[L*(Delta+1)+k+k*value_width])

plus finding and updating the affected Gram entries, feature incidence
lists and routing metadata. Here Delta bounds the union graph degree.
An inverted feature index with at most h keys per coordinate can find
the changed inner products in O(k*h) scalar operations; maintaining
that index, the Gram graph and all sparse encoding columns uses
additional storage. Large h is not excluded just by a small Gram degree,
because inner-product cancellations can hide shared coordinates.
All these quantities must be charged. The construction visits old/new
graph neighborhoods and compares coefficient changes; it does not prove
that encoder evaluation or graph discovery is cheap for arbitrary keys.

The fixed alpha is essential. Retuning alpha after each key change can
alter all diagonal entries of A and invalidate this locality argument.
Insertion requests violating the uniform contraction or support bounds
are outside the guarantee. Global changes to encoder parameters are
also outside the one-row edit model.

The exact script now checks this extension against an independently
recomputed full encoding, including inserts, deletes, feature replacements,
concurrent value changes and fixed-point state updates. It uses sparse
selected-column polynomial expansion for the local correction, with a
dense Gram/full-encoding oracle for verification. It is not a runtime
implementation or performance benchmark.

The literature gate has also become more specific.
[Dynamic Kernel Graph Sparsifiers](https://arxiv.org/html/2211.14825v2)
already studies moving points, dynamic kernel graphs and sketches of
Laplacian solutions. Its introductory statements and the displayed
Theorems 2.1--2.4/G.1/I.5 were inspected. Its Laplacian/sketch output,
error norm, kernel regularity, dimension/aspect-ratio assumptions and
adversary model differ from the explicit per-key maximum-error guarantee
above; a full proof comparison is still required. A low-dimensional
sketch of a solution must not silently be treated as individually readable
memory values. This review does not assert a flaw or an improvement over
that work.

The next decisive question is whether a useful high-dimensional binding
encoder can satisfy the locality/separation assumptions with a total
storage and precision cost that improves an applicable prior bound.
The polynomial identity alone is elementary and dynamic graph-filter
maintenance is an existing subject. The example cycle features are
controls, not evidence of a novel architecture. The unresolved work is
the strongest-prior-art comparison and a construction/separation with
meaningful query robustness, not a larger sweep of these examples.

#### Generality check: a concrete binding workload and its strong baselines

One sufficiently explicit workload is a matrix of independently
replaceable bindings V_ij with bounded values, addressed by two vector
queries u and v. Asking for u^T*V*v preserves interactions between the
two index sets while allowing arbitrary, full-rank value tables. The
two query parts can each have many coordinates; they are not two scalar
latents. All table entries, query construction and index storage count.

Requiring fast answers for every *unnormalized* Boolean u,v at absolute
error strictly below 1/2 already contains Boolean vector-matrix-vector
multiplication: on Boolean V, threshold the integer-valued arithmetic
answer at 1/2 to decide whether any requested binding is present.
[Henzinger et al., Theorem 2.7/Corollary 2.8](https://arxiv.org/html/1511.06773v1#S2.SS2)
place this problem behind the OMv conjecture even with polynomial
preprocessing. The statements and reduction argument were inspected.
Thus an unrestricted truly subquadratic guarantee here would be a major
conditional-complexity breakthrough. This is a known conditional
barrier, not an unconditional impossibility theorem or a new result.

Normalization and approximation change the baseline substantially.
For nonnegative probability vectors u,v and scalar values in [-1,1],
draw independent I~u,J~v and average s samples of V_IJ. This estimates
u^T*V*v with error at most epsilon except with probability delta when

    s >= 2*epsilon^(-2)*log(2/delta).

This is the standard bounded-variable concentration argument. A current
table supports exact latest-value replacement in one cell; one-hot
queries also read their value exactly. General binding queries use s
table probes after paying to construct/sample the two distributions.
For dense input vectors, constructing the samplers generally costs at
least linear work in their lengths. The independent table still costs
all its stored bits. Fresh randomness gives a conditional per-query
guarantee under adaptive past updates; a simultaneous guarantee for H
queries needs its own failure accounting, for example delta/H per query.
It is not an unlimited-horizon deterministic retention theorem.

This is an existing Monte Carlo baseline, not a proposed contribution.
A constant-cost approximate binding query by itself would not beat it.
Our next comparison must distinguish exact stored associations, robust
mixture queries, deterministic versus probabilistic error, and the costs
of encoding/routing. We should seek a useful improvement to that complete
tradeoff while avoiding both an OMv-hard unrestricted promise and a
repackaging of ordinary lookup or sampling.

### Associative-memory goal: address geometry and signed writes

The previous goal turn made progress: it established supporting bounds,
defined the binding workload and supplied exact certificates. It did not
establish the new architectural advantage required by the active goal.
This iteration tests whether a useful encoder can meet those bounds.
The objective and acceptance criteria remain unchanged. No additional
neural training is justified by the results below.

#### A positive-address restriction on the local-correction route

Let F have n nonnegative rows f_i in R^m, each of Euclidean norm one.
Write G=F*F^T, s_i=||f_i||_1^2 and
R=max_i sum_{j != i} G_ij. Here s_i measures effective width in this
particular memory basis; it is **not** intrinsic input dimension.
Then

    (sum_i sqrt(s_i))^2 <= m*n*(1+R).

Indeed, Cauchy--Schwarz applied to the m coordinates of sum_i f_i gives

    (sum_i ||f_i||_1)^2 <= m*||sum_i f_i||_2^2
                         = m*1^T*G*1 <= m*n*(1+R).

Positivity is used in replacing the sum of row norms by the norm of the
sum. If every s_i >= s, this implies m >= n*s/(1+R).
For alpha>0, the exact infinity-norm contraction constant is

    ||I-alpha*G||_infinity = |1-alpha| + alpha*R.

Some alpha gives a strict contraction if and only if R<1. Consequently
the first iteration's guarantee, at m approximately n cells, requires
effective positive address width below approximately two. Broad equal
positive weights on k cells require m > n*k/2 under this particular
guarantee. Tiny nonzero tails can still encode many input coordinates;
this is not a dimension or general memory impossibility theorem.

There is also a spectral variant. If lambda_max(G)<=K, the same proof
gives m >= (sum_i sqrt(s_i))^2/(n*K). In particular a positive-definite
G with condition number at most K has lambda_max(G)<=K, since its
eigenvalues average one. Thus replacing infinity-norm contraction by
bounded spectral conditioning does not automatically remove the issue.

The bound is sharp. For a cyclic permutation matrix P on n>=3 cells,
F=c*I+s*P with c,s>=0 and c^2+s^2=1 has effective width
(c+s)^2=1+2*c*s and off-diagonal Gram row sum R=2*c*s.
Signed Hadamard features give the opposite control: m=n, G=I and
effective width n. They are a dense change of basis of private cells,
with dense writes; their existence proves no architectural advantage.

These are elementary frame inequalities. Their use here is to reject
an unjustified positive-encoder assumption, not to claim a paper.
Strictly positive softmax addresses also have a complete overlap graph;
the earlier sparse graph repair guarantee needs genuine sparsity or a
separately bounded truncation error.

#### A capacity bound for normalized additive squared-kernel recall

Let u_1,...,u_n be unit vectors in R^p or C^p. Put P_i=u_i*u_i^*,
G_ij=Tr(P_i*P_j)=|u_i^*u_j|^2 and R_i=sum_{j != i}G_ij.
Consider exactly the canonical additive memory with *latest* scalar
values v_j in [-1,1], queried at stored key i as

    y_i = sum_j G_ij*v_j / sum_j G_ij.

This grants correct replacement of old contributions, so it does not
confuse stale duplicates with interference. Maintaining that canonical
state may require a latest-value/key cache, which must still be charged.

For each i the largest possible |y_i-v_i| over these values is exactly
2*R_i/(1+R_i), attained by v_i=1 and all other v_j=-1. On the other hand

    sum_ij G_ij = ||sum_i P_i||_F^2 >= n^2/p,

by Cauchy--Schwarz on the p eigenvalues, whose sum is n. Therefore

    max_{i,v} |y_i-v_i| >= 2*max(0, 1-p/n).

A uniform error guarantee <=epsilon<2 requires
n <= p/(1-epsilon/2), despite a quadratic feature state of order p^2.
This is a Welch/frame consequence for this additive normalized reader.
It does not bound delta corrections, dual writes, nonlinear reads,
stochastic or distributional recall, or all architectures at equal memory.
For a degree-2r dot-product kernel the same argument applies to normalized
symmetric r-th tensor powers in dimension D=binomial(p+r-1,r), replacing
p by D. A larger exact factorization rank alone is not a uniform
arbitrary-value retention guarantee.

The directly relevant published baseline is
[KATA, arXiv:2607.17419v1](https://arxiv.org/html/2607.17419v1).
Its PSD features, MUB dictionary, idealized retrieval statement and delta
recurrences were inspected (Theorems 5/7, Sections C.5/C.7 and D.1/D.3).
Theorem 7 uses isotropic independent values and idealized distractor
correlations; it is not the uniform adversarial-value guarantee above.
Equation 28 corrects the current key using that key's PSD feature; it
does not assert preservation of all earlier keys. These distinctions
are consistent with the paper's scope, not counterexamples to its stated
theorems. We have not reproduced its experiments or audited every proof
and implementation detail. Its larger feature state must be counted in
any comparison. Neither polynomial lifting nor MUB addressing is new.

#### Exact signed dual writes, and why they do not yet pass our goal

Suppose r mutually unbiased orthonormal bases in R^p or C^p are available:
their rank-one projectors satisfy

    Tr(P_ai*P_bj) = delta_ij  if a=b, and 1/p otherwise.

Choose one reference projector P_a0 in each basis and retain the other
p-1 projectors as query addresses. For each retained address use the
signed write matrix W_ai=P_ai-P_a0. Directly from these overlaps,

    Tr(P_bj*W_ai) = delta_ab*delta_ij

for every retained query. Thus n=r*(p-1) independently replaceable
values can be kept in S=sum_ai v_ai*W_ai with exact reads Tr(P_ai*S).
The update

    S <- S + (v_new-Tr(P_ai*S))*W_ai

replaces one value and preserves every other retained read. It uses no
latest-value cache in exact arithmetic. Each W has trace zero and
matrix rank two. The state is generally indefinite: this changes the
write/normalization architecture, rather than preserving a convex
attention average. The ordinary unit-strength delta direction P_ai
instead changes a cross-basis read by the correction divided by p.

A complete complex MUB set has r=p+1 at prime-power p, yielding p^2-1
independent labels, the dimension of traceless Hermitian matrices.
These are full p-coordinate quadratic queries, including cross-coordinate
terms. The classical MUB decomposition already explains this construction.
[Filippov--Man'ko, Section 3, equations 6--10](https://arxiv.org/html/1008.2675v1)
explicitly uses omitted reference projectors, the centered Gram matrix
delta_ab*(delta_ij-1/p), and tomographic inversion. Those definitions,
the displayed reconstruction and its proof were inspected. Our rank-two
dual identity is an elementary consequence, not a defensible novelty claim.

The arithmetic and access accounting are restrictive. A generic dense
matrix read/write takes O(p^2*d_v) scalar operations for value width d_v;
rank-two writes do **not** touch only two memory cells. A Hermitian state
uses p^2-1 real coefficients per value after omitting its trace; real
symmetric states use p*(p+1)/2-1. An explicit codebook can cost O(r*p^2)
real coordinates, and explicit dual matrices O(n*p^2). An implicit known
construction can save that metadata only after counting its evaluation,
index decoding and precision costs. Addressing arbitrary continuous new
keys is not supplied by this fixed-codebook result. A dynamic dictionary
already solves discrete-address recall, so that is not a sufficient
separation; any claimed benefit must include the richer query semantics.

The general identity above is algebraic. Our exact finite-precision
controls use real MUBs at p=2 and p=4, with dyadic projectors. If D*W_ai
is integral and |v_ai|<=R are integers, store T=D*S as integers. Then
|T_jk|<=n*D*R, and exact current reads followed by integral replacement
maintain this bound for any overwrite horizon. This is O(log(n*D*R))
bits per state entry, plus signs/indices and read-arithmetic costs.
It is not a claim that arbitrary prime-power complex formulas can be
used at fixed precision without a further error and number-field audit.

Noise is a material limitation even in exact arithmetic. Give every
nonreference key of one basis value +1 and all retained keys of other
bases value zero. Then S=I-p*P_a0. For a stored unit vector u orthogonal
to the reference vector w, a nearby unit query
q=sqrt(1-t^2)*u+t*w reads

    q^*S*q = 1-p*t^2.

Thus this construction's worst-case constant-error neighborhood can
shrink as p^(-1/2). Exact retention at codewords is not automatically
robust semantic retrieval. No robustness advantage over existing memory
methods follows from this construction.

#### Bounded continuous reads impose another capacity gate

There is a broader supporting bound even for signed writes. Suppose
each assignment v in {-1,+1}^n has some Hermitian state S_v, with
|u_i^*S_v*u_i-v_i|<=epsilon<1 at the fixed unit keys. Suppose also
|q^*S_v*q|<=B for every unit query q and every assignment. Then

    n <= p*B^2/(1-epsilon)^2.

The second condition is equivalent to ||S_v||_op<=B, so
||S_v||_F<=sqrt(p)*B. For independent signs sigma_i let
X=sum_i sigma_i*P_i. The matched-value guarantee gives
Tr(S_sigma*X)>=(1-epsilon)*n for every assignment, while

    E Tr(S_sigma*X) <= sqrt(p)*B*E||X||_F
                   <= sqrt(p)*B*sqrt(E||X||_F^2)
                    = B*sqrt(p*n).

The last equality uses independence and ||P_i||_F=1. Rearranging proves
the bound. This does not assume a linear update or encoder for values.
It does assume a single quadratic reader with bounded outputs everywhere.
Clipping, nonlinear or routed reads escape that reader model and must be
examined separately. In particular, n of order p^2 requires worst-case
output range B at least of order sqrt(p) at constant matched error.
The signed construction's large outputs are not just an omitted
normalizer one can replace by a constant without checking retention.

This is a classical Hilbert-space margin/shattering argument.
[Gyurik--van Vreumingen--Dunjko, Proposition 5 and Appendix A.3](https://arxiv.org/html/2105.05566v3)
already bounds fat-shattering by observable Frobenius norm squared divided
by margin squared. Those statements and the displayed reduction/proof
were inspected. Taking norm bound sqrt(p)*B and margin 1-epsilon gives
the same order; the constant-one proof above does not supply substantive
novelty. This is a known-corollary gate, not our paper's main theorem.

#### Literature decisions and reproducible checks

Other tempting shortcuts already have substantial prior art:

- [SANTA/S2ANTA](https://arxiv.org/html/2605.01910v2) analyzes sampling
  and equal-mass stratification, including unbiasedness and covariance
  domination (Appendix A). Those arguments were inspected. Value-fetch
  savings do not eliminate score formation or sampling setup costs.
- [Three Tokens Force Exponential Feature Rank](https://arxiv.org/html/2608.11427v1)
  already gives a short-context lower bound for a nonnegative kernel
  head. Its task, reader model, rank argument and stated limitations
  were inspected; its entire proof chain was not independently verified.
  Signed readers and recurrence are outside that rank theorem. A generic
  positivity/rank separation alone is not an available novelty claim.
- [Learnable Low-Rank Polynomial Sketch](https://arxiv.org/html/2610.14219v1)
  already studies learned tensor factorizations. Theorem 3.1 and its
  truncated-SVD proof were inspected. That fixed-data existence statement
  is not a uniform dynamic-key retention theorem; experiments and the
  remaining proof chain were not audited here.
- [ZeroS](https://arxiv.org/html/2602.05230v1) already studies signed,
  zero-sum attention and stability. Proposition 3.5 and its displayed
  Appendix A.1.5 proof were inspected: its Lipschitz guarantee assumes
  bounded logits, residual-weight scaling and a decaying output factor.
  That stability statement alone does not assert exact latest-value
  retention. A signed reader proposal must compare both guarantees;
  removing a positive background is not itself new.

Commands:

```bash
python3 -B scripts/run_retention_geometry_checks.py
python3 scripts/run_retention_tradeoff_checks.py
python3 scripts/run_feature_dimension_checks.py
./build.sh
python3 tests/test_smoke.py
python3 -B scripts/run_structural_complexity_checks.py
```

The new report is `results/research/retention_geometry_checks.json`.
It checks positive-width equality cases, signed-basis escapes, additive
PSD worst-value witnesses and dyadic dual-write streams against direct
reconstruction of all latest values. Exact arithmetic supports the
displayed proofs; small dimensions do not establish novelty or an ML win.

The decision is to park a naive positive local-correction encoder and
standalone MUB dual-writing as proposed contributions. Retain their
proofs and controls. Continue toward a genuine **robustness/access/retention
tradeoff** for a richer reader: specify the allowed routing or nonlinearity,
prove a useful uniform guarantee, count its whole state and precision,
and compare against dictionaries, sampling, kernel regression and signed
attention on the identical workload. No new publishable advantage has
yet passed that gate. The goal stays active; its success criteria do not
weaken to accommodate these supporting results.

### Associative-memory goal: local correction with shared features

Previous turn classification: progress. The address-geometry proofs,
literature decisions and exact certificates changed the next construction
to investigate. They did not achieve the active goal. This iteration first
tests routed retrieval, then removes an unnecessary limitation of the
earlier local-correction model. The full novelty, usefulness and publication
criteria remain in force; no additional neural training is authorized by
these supporting results.

#### Routed sparse retrieval is already a strong baseline

For the Euclidean simplex projection sparsemax, a score gap of at least
one makes the output exactly one-hot. Perturbing all scores by at most
eta preserves this conclusion if the original gap is at least 1+2*eta.
Keeping latest values in a table then gives exact, history-independent
retention in that query neighborhood, with ordinary point replacements.
This is an established margin/projection consequence, not our result.

Factorized access does not make this a new architecture theorem:
[Correia et al., Section 4.1](https://arxiv.org/html/2007.01919v2)
already evaluates structured sparsemax through a top-k oracle, and
certifies the full solution when the restricted solution has support
strictly below k. Those definitions, the support certificate, Section 4.2
and Appendix A's displayed forward/backward formulas were inspected.
Combining that primitive with known product-key top-k search is a known
construction route. Its active support can be the whole table in the
worst case. A sparse output is not a proof that finding its support is
cheap, nor an improvement over an appropriately indexed dictionary.

[Nonparametric Modern Hopfield Models, Section 3.2](https://arxiv.org/html/2404.03900v2)
also studies masked/sparse retrieval. Its displayed construction, three
mask examples and the target-in-mask premise were inspected, not its
entire proof chain. Searching for the mask must be counted; assuming
the target is included does not supply an arbitrary-query routing theorem.
[U-Hop+, Section 3.1 and Algorithm 1](https://arxiv.org/html/2410.23126v2)
was checked as well: its displayed O(1/T) rate concerns feature-map
optimization iterations. We must not read this as a sublinear online
lookup-time bound. No claim about these papers' empirical results is
adopted here. The simple routed-lookup proposal is parked as a sufficient
contribution; richer binding semantics and total-cost comparisons remain.

#### A conditional model that separates local and shared overlap

The previous positive-address bound does not justify excluding every
dense Gram matrix. Dense overlap may come from a small shared component
that can be handled separately. Let an encoded key have features
f_i=(b_i,u_i), with B in R^(n*m), U in R^(n*r). The rows b_i have at most
k nonzeros; u_i can use all r shared coordinates. The value matrix V has
width v, entries bounded by R, and independently replaceable rows. There
is no low-rank assumption on V or on a binding table formed by its rows.

Use n fixed-capacity slots. Inactive slots have b_i=u_i=V_i=0. Set

    G0 = B*B^T + diag(inactive),  A = I-G0,
    P = sum_{j=0}^L A^j.

Assume ||A||_infinity<=rho<1. Let C be P rounded entrywise to a fixed
dyadic grid, with symmetric rounding and error at most eta per supported
entry. Graph distance is measured in the off-diagonal support of G0.
If each radius-L ball has at most b vertices, then

    ||G0*C-I||_infinity <= e_L,
    e_L = rho^(L+1) + (1+rho)*b*eta.

This follows from G0*P-I=-A^(L+1), the graph support of the powers, and
||G0||_infinity<=1+rho. Assume C is positive definite and ||C||_infinity<=K.
These are genuine conditions: for instance, symmetry and
b*eta < (1-rho^(L+1))/(1+rho) ensure positivity, and one may take
K=1/(1-rho)+b*eta. An inactive slot retains its isolated unit diagonal
exactly; its coefficient, value and shared-feature contributions vanish.

Define the following current quantities:

    Z=C*V,  X=C*U,  H=I+U^T*X,  t=U^T*Z,  h=H^(-1)*t,
    S0=B^T*Z,  Y=B^T*X.

A query encoded as (b_q,u_q), with b_q supported on k_q coordinates,
reads

    read(q) = b_q^T*S0 + (u_q^T-b_q^T*Y)*h.

The state is represented by these cached quantities; it is not required
to materialize the dense vector of interpolation coefficients
d=Z-X*h. In fact

    B^T*d = S0-Y*h,  U^T*d = h.

The second identity follows from U^T*X=H-I. Thus the displayed reader is
exactly f_q^T*[B^T*d;U^T*d]. Its effective query coefficients may be signed.
It does not preserve the convex-average interpretation of attention.

#### Uniform retention and a precisely scoped robustness bound

Suppose, uniformly over admissible key edits,

    ||u_i||_2 <= M,   H >= c*n*I  for a fixed c>0.

The inequality is in the positive-semidefinite order; it is stronger than
mere invertibility. Degenerate shared directions or arbitrary insertions
are not silently admitted. Since

    ||Z||_max <= K*R,   ||t_column||_2 <= n*M*K*R,
    ||h_column||_2 <= M*K*R/c,

we have ||V-U*h||_max <= Gamma*R, with

    Gamma = 1 + M^2*K/c.

Here ||.||_max is the maximum absolute entry, not a matrix row-sum norm.
For the active keys, use G=G0+U*U^T. Direct multiplication gives the
residual identity

    G*d - V = (G0*C-I)*(V-U*h).

For inactive slots d_i=0, so the ghost diagonal has no effect on the
actual feature reads. Consequently every active key has maximum value
error at most e_L*Gamma*R, independently of overwrite horizon. This bound
depends on the current local geometry and shared conditioning, not on
the contraction constant of the full, generally dense Gram matrix.

It also yields a conditional query-neighborhood statement. Suppose each
private feature coordinate occurs in at most h0 key rows and |B_ij|<=J.
Then ||S0_column||_infinity<=h0*J*K*R and each row of Y has Euclidean norm
at most h0*J*K*M. If the encoded query differs from an active key by
delta_b in l1 norm and delta_u in Euclidean norm, its error relative to
that key's latest value is at most

    R*[e_L*Gamma + h0*J*K*Gamma*delta_b + M*K*delta_u/c].

This follows by subtracting the two displayed read expressions and
applying the bounds on S0, Y and h. It is a guarantee in **encoded feature
space**. No raw-input semantic encoder or dimension-independent routing
margin is supplied by this lemma. Those costs and hypotheses cannot be
removed when proposing a model. The parameters r, c, M and K also count;
holding them fixed is not a theorem for arbitrary key sets.

#### Why key edits need not invalidate every cached row

An edit to B_i changes G0 only on row/column i, including the ghost
diagonal when occupancy changes. The matrix-power difference identity
from the first iteration shows that C changes only inside the square
indexed by the radius-L ball S of i in the union of the old/new local
overlap graphs. Rounding preserves that support statement. Editing U_i
or V_i in the same operation then changes X=C*U and Z=C*V only in S.

For each j in S, replace its old row contribution by its new contribution
in H, t, Y and S0:

    H <- H + sum_j (u'_j^T*x'_j - u_j^T*x_j),
    t <- t + sum_j (u'_j^T*z'_j - u_j^T*z_j),
    Y <- Y + sum_j (b'_j^T*x'_j - b_j^T*x_j),
    S0 <- S0 + sum_j (b'_j^T*z'_j - b_j^T*z_j).

Row vectors are used in these outer products. Recompute h from the new
r-by-r H and t. The sums include unchanged u_j or b_j whenever x_j or
z_j changes. Updating only the edited key's outer product is incorrect.
All caches remain the exact functions of the *current* rounded C, keys
and latest values above. Retuning the global polynomial or changing all
encoder parameters is outside this one-row-edit guarantee.

In scalar operations, with degree at most Delta and union-ball bound b,
a conservative key-edit bound is

    O(b^2*[L*(Delta+1)+r+v] + k*b*(r+v)
      + b*r*(r+v) + r^3 + r^2*v),

plus incidence discovery and encoder evaluation. It allows recomputing
at most b polynomial columns and their supported row products. Value-only
replacement is cheaper: X, H and Y stay fixed; symmetry of C gives
delta_t=X_i^T*delta_V_i. Its cost is O(k*b*v+r^2*v), plus the value cache
write. The reader costs O(k_q*(r+v)+r*v), plus obtaining its encoded query.
These are algebraic operation bounds, not machine-word or measured time.

All caches must be charged: latest V and Z use 2*n*v scalars; S0 uses
m*v; U and X use 2*n*r; Y uses m*r; H and its optional inverse use
O(r^2); t and h use O(r*v). The sparse C needs O(n*b) coefficients and
indices; B and incidence lists need O(n*k), with slot masks and any raw
keys/codebooks/encoder parameters added separately. Preprocessing C,
all row products and the small solve is charged, too. In particular this
does not compress independent values below the dictionary information
bound. A bounded degree alone can make b exponential in L; a small
update cost needs an actual ball-growth or depth bound.

#### Finite precision without accumulation over the stream

Assume B, U and V lie on fixed dyadic grids, with bounded coefficient
ranges. Use a fixed dyadic scale for C and keep Z, X, H, t, Y and S0 as
exact integer numerators at their known product scales. The row-difference
updates above replace current canonical contributions. Their bit lengths
are therefore bounded by the current sums, not by the number of edits.
For example, if the scales of C,U,V are D,E,Q respectively, then

    |Z_numerator| <= D*Q*K*R,
    |X_numerator| <= D*E*K*M,
    |H_numerator| <= D*E^2*(1+n*M^2*K),
    |t_numerator| <= D*E*Q*n*M*K*R.

Use ceiling integer bounds where needed. The h0,J bounds similarly
bound numerators for Y and S0. Invert or solve the small current H with
exact rational arithmetic, rather than accumulating rounded inverse
updates. Clearing denominators and the determinant bound give bit
length O(r*[log(n)+log(r)+input_scale_bits+coefficient_bits+range_bits])
for this solve. Its bit-operation cost must multiply the scalar-operation
accounting; it is not an unlimited-precision unit-cost claim. Rounded
floating-point cache updates do not inherit the history-independent
guarantee without a separate numerical analysis.
Forming P exactly before rounding also uses temporary integers with bit
length growing linearly in L and the input grid precision (plus degree
and range terms). Their arithmetic must be charged when repairing columns.

#### What this advances, and what remains unproved

The static identities are classical Woodbury/Schur-complement algebra.
[WISKI, Section 4 and Appendix A.2--A.3](https://arxiv.org/html/2103.01454v1)
already combines sparse interpolation, cached sufficient quantities and
Woodbury updates. Its displayed formulas and update arguments were
inspected in this iteration, beyond the earlier abstract check. Its
constant-time-in-data-count statement retains inducing-dimension costs
and O(m^2) caches; its accumulating-observation posterior is not a
latest-value overwrite specification. This difference is not evidence
that our conditional extension is new: dynamic regression, sparse inverse
approximation, latent-variable inference and local graph-filter methods
still require an exact theorem comparison. We have not established a
strongest-baseline separation. In particular B=I reduces to ordinary
table storage plus small regression summaries.

There is also close numerical-linear-algebra prior art:
[Zheng--Xi--Saad, Sections 3.1--3.2](https://arxiv.org/html/2002.00917v1)
already combines a power expansion with a low-rank Schur correction.
Its displayed inverse formula (22), Theorem 3 and telescoping proof were
inspected. The static power-plus-low-rank construction is therefore
occupied. Whether the fully counted dynamic repair and maximum-coordinate
retention guarantee above is a useful extension still requires a complete
comparison; a different application name or error statement alone is
insufficient evidence of a substantive gap.

The concrete advance is an explicit, precision-accounted retention and
encoded-query bound allowing shared dense overlap, with local repair
under key insert/delete/replace. This invalidates an overly broad rejection
of dense-address models based only on full-Gram contraction. It does not
turn the supporting algebra into a new paper or establish raw-input
robustness. An encoder realizing the promises for a useful binding
workload, together with a better full tradeoff than existing methods,
remains necessary for the active goal.

Reproduce the exact controls with

```bash
python3 -B scripts/run_shared_feature_retention_checks.py
```

The report `results/research/shared_feature_retention_checks.json` checks
n=9,27 slots, r=1,3 shared dimensions and L=0,2,4, including 108 key edits
with concurrent value changes and 36 value-only overwrites. Sparse-column
recomputation and local integer cache repairs
are compared with full current reconstruction. A small independent dense
inverse checks the Woodbury solve; direct feature reads check the implicit
state. The depth-four controls certify a worst-case value error below
one tenth under their stated promises and test encoded perturbations.
These are exact proof controls with dense oracles, not runtime or neural
benchmarks. Small r controls do not establish a high-dimensional end-to-end
architecture; the displayed proof quantifies the dependence on arbitrary r.

### Associative-memory goal: correlated routing and a stronger-baseline check

Previous turn classification: progress. It established the local/shared
retention lemma with precision and repair costs, without establishing a
useful input encoder or a novel strongest-method advantage. This iteration
tests an explicit higher-dimensional binding workload and proves a
head/support tradeoff. Stronger-baseline controls then show why that
tradeoff is insufficient to complete the goal.

#### A uniform signed-head bound, with its exact resource model

Store arbitrary independent scalar values V_ij in an n-by-n Cartesian
table. A head returns sum_ij P_h(i,j)*V_ij, where P_h is nonnegative.
The final linear output may use arbitrary signed coefficients lambda_h.
All heads use the same table coordinates; separate row/column reorderings
are allowed. No pair-specific table remapping, value-dependent query
weights or additional value summaries are included in this model.

Suppose each head support is a chain/Ferrers bipartite graph: whenever
(i,j) and (i',j') are present with distinct rows/columns, at least one of
(i,j') and (i',j) is present. Additive scores a_i+b_j with a threshold
or top-k cutoff have this property. Ties require a compatible rule, such
as lexicographic row/column infinitesimals; an arbitrary pair-specific
tie rule could itself encode a correlated mask and is excluded.

[Original PKM, Section 3.3](https://arxiv.org/html/1907.05242v2) was checked
for this comparison: its heads have separate query networks/subkeys,
share the value table and sum their outputs. Our signed combination is
more permissive than that reader. Applying the bound to independently
replaceable literal bindings still requires the stated value-storage
model; it does not constrain every learned encoding of those bindings.

Let H_+ be the number of heads with positive lambda_h, and let

    C_support = sum_{h: lambda_h != 0} |support(P_h)|.

This counts **unfused intermediate head-support entries**, not unique
physical value probes, memory size, or an implementation-independent
runtime. Suppose the resulting coefficient matrix W equals the uniform
matching D, with D_ii=1/n and all off-diagonal entries zero. Then

    H_+ * C_support >= n^2.

The same statement holds for any permutation matching by renaming columns.
Proof: each positive target diagonal must have a positive contribution.
Assign every diagonal to one such head, obtaining disjoint groups of
sizes a_h summing to n. Within group h, every pair of selected diagonals
forces an off-diagonal support entry in that head. Choose one such entry
per pair. These entries are distinct, including across groups, because
the groups have disjoint row/column labels. There are
F=sum_h binomial(a_h,2) forced entries. Positive heads use at least n+F
support entries. Every forced off-diagonal contribution must be cancelled
by at least one negative head, using at least F additional entries.
Thus C_support>=n+2*F=sum_h a_h^2>=n^2/H_+.

In the unsigned case there can be no cancellation, so every head contains
at most one target edge and H_+>=n. This is the classical induced-matching
obstruction for chain-graph covers, not a new head-rank principle.
[Alecu et al.](https://arxiv.org/html/2104.04471v1) was inspected for the
chain/quasi-chain definitions and decomposition context. The explicit
induced-matching argument also appears as Lemma 4.1 in
[Edge partitions into induced-2K2-free bipartite graphs](https://combinatorialpress.com/ars-articles/volume-166/edge-partitions-into-induced-2k2-free-bipartite-graphs/),
whose displayed definition, lemma and proof were inspected. Its unsigned
edge-partition result is not a claim about signed attention.

The signed bound is tight in this support-count model. Partition the n
indices into g equal blocks of width m=n/g. In each block use separable
positive raw weights 2^(j-i), normalized on the upper triangle including
the diagonal, and on the strict upper triangle. Very small row/column
weights outside the block keep both supports within it. Each is a top-k
head for the same additive ordering. If their normalizers are Z_+,Z_-,
combine them with coefficients Z_+/n and -Z_-/n. Off-diagonal weights
cancel and each diagonal has weight 1/n. The two supports together use
m^2 entries. Summing blocks gives H_+=g and C_support=n^2/g, with at most
2*g heads. Singleton blocks need only the positive head.
These are exact rational-weight witnesses; their coefficient bit lengths
and cancellation gain grow with block size. This is not a fixed-precision
or numerically stable model implementation claimed for free.

#### An approximate version requires a weight condition

For standard monotone normalizations of additive scores, two positive
diagonal weights also force a cross weight at least their minimum.
The checker verifies this for normalized separable weights with top-k
and compatible ties. Assume this condition in addition to the support
condition. Write epsilon=||W-D||_1, summing absolute values of all entries.
For 0<=epsilon<=1/2,

    C_support >= n^2*(1-2*epsilon)^2/H_+ - 2*n*H_+*epsilon.

Here is the proof, including its dependence on H_+. Let S be the diagonal
indices with W_ii>=1/(2*n). There are s>=n*(1-2*epsilon) such indices.
Assign each to its largest positive head contribution, which is at least
1/(2*n*H_+). A pair assigned to the same head now forces a cross
contribution at least this large. Among the F>= (s^2/H_+-s)/2 distinct
forced entries, at most 2*n*H_+*epsilon can lack a negative-head support:
each uncovered entry would contribute at least 1/(2*n*H_+) to the error.
Hence C_support>=s+2*F-2*n*H_+*epsilon, proving the formula.
It can be vacuous when H_+ is large; do not omit the negative term.

For fixed value-independent linear coefficients, this norm is exactly
the worst output error over V_ij in [-1,1], attained by sign values.
Thus it is a finite-alphabet value witness as well as a coefficient bound.
It does not cover value-dependent routing, arbitrary nonlinear decoding,
precomputed value summaries or multiple recurrent stages. Nearby head
lower bounds use different models: the scalar normal form and Theorem 1
statement in [The Head Complexity of Boolean Functions](https://arxiv.org/html/2609.04046v1)
were inspected, without a full proof audit. Its Boolean attention-only
model is not automatically our truncated Cartesian-support model.

No novelty is established for the signed counting statement. Its relation
to signed chain decompositions and architecture resource bounds needs
an exact prior-art comparison. A search that did not find the formula
would not establish a mathematical gap.

#### Explicit correlated queries, and the dictionary that matches them

An extra bilinear score can produce a matching support in one head.
For general score arrays

    s_ij = a_i+b_j+x_i*y_j,

sorting x and y gives the anti-Monge identity

    s_ij+s_i'j' - s_ij'-s_i'j = (x_i-x_i')*(y_j-y_j') >= 0.

Row-max search in totally monotone matrices is classical
[Aggarwal et al.](https://research.ibm.com/publications/geometric-applications-of-a-matrix-searching-algorithm);
its primary abstract and stated resource were inspected. This alone
does not establish an output-sensitive global sparse-normalization
algorithm, a new selection theorem, or a novel architecture advantage.

A concrete raw-query family is enough to test the proposed comparison.
Let x=(1,...,n), let y be an arbitrary permutation of those tags, and set

    s_ij = -(4/n)*(x_i-y_j)^2.

The query consists of two n-dimensional vectors, not two scalar latents.
Sparsemax over all n^2 scores has weight 1/n exactly at the matching
pairs and zero elsewhere. Arbitrary values in all n^2 table positions
are still independently replaceable; no low-rank binding assumption is
being made. Perturb each tag by at most eta<=1/8. Matched squared gaps
are at most 4*eta^2; unmatched gaps are at least (1-2*eta)^2. The simplex
threshold lies between -(1+16*eta^2)/n and -1/n. All matched weights stay
positive, all unmatched weights vanish, and the coefficient l1 error
relative to uniform matching is at most 8*eta^2. The last bound uses
mean absolute deviation at most half the range of the matched squared
gaps. All these statements hold uniformly in n.

However, an ordinary indexed lookup is stronger on this very workload.
Round the noisy tags to integers, build an inverse-tag array from y,
then read V_i,j at the unique matching j and average. This takes O(n)
word operations/probes and O(n) temporary index words on the stated
integer-tag family, with exact answers for tag noise below 1/2. Latest
point replacement takes one value-word write for any horizon. It uses
the full n^2 value table, plus the query and index storage, with no hidden
continuous state. The tags need O(log n) integer bits and their specified
fractional bits; the sum needs the value precision plus O(log n) bits.
Sorting gives a standard deterministic alternative for more general
separated comparable tags. The bilinear example therefore does not beat
the strongest applicable dictionary/join baseline.

Approximate probabilistic queries have an even cheaper familiar baseline
on this particular family. Since the first tags are the canonical order,
sample a column j uniformly, round y_j to get its row i, and read V_ij.
The sample mean is unbiased for the matching average. For values in
[-1,1], O(epsilon^(-2)*log(1/delta)) samples suffice for additive error
epsilon with per-query failure probability delta. This avoids building
the inverse-tag array. It retains the whole latest-value table; fresh
randomness gives a per-query guarantee, while simultaneous guarantees
over a horizon need separate failure accounting. The deterministic
coefficient-norm lower bound does not exclude such per-value randomized
guarantees. These standard sampling quantifiers must remain explicit.

#### Two controls that prevent a false fundamental separation

First, arbitrary table layout is outside the common-grid bound. The
bijection (i,j)->(i,(j-i) mod n) moves the logical diagonal to one physical
column. A single additive top-n head can select that column and return
the uniform logical diagonal average with n physical accesses. Individual
logical replacements still change one cell. This example concerns one
fixed matching, rather than solving every permutation through one fixed
layout, but already refutes an unrestricted reading of the head bound.

Second, signed-head fusion changes what must be read. The two triangular
heads above have n^2 intermediate support entries but their combined
matrix has only n nonzeros. An evaluator that recognizes the cancellation
can read only the n target values. Hence C_support is not a lower bound
on physical probes for the fused computation. Computing normalization,
discovering cancellation and handling numerical error have their own
costs; we have not lower-bounded those for all implementations. Claiming
a value-access win solely against unfused evaluation would violate our
strongest-baseline criterion.

The decision is to preserve the common-grid signed tradeoff as a scoped
supporting lemma, and park the matching-query construction as a claimed
architectural advantage. This iteration supplies an explicit raw-query
family and a sharp signed support count, but the honest decision remains
that no new publishable strongest-method advantage has been established.
The next candidate must survive arbitrary permitted indexing, algebraic
fusion and useful cached summaries on the *same* dynamic workload, not
merely a head-by-head implementation comparison. The earlier local/shared
retention lemma remains a conditional tool, not a completed contribution.

Reproduce the controls with:

```bash
python3 -B scripts/run_correlated_routing_checks.py
```

The report `results/research/correlated_routing_checks.json` records
6,885 additive-head cases including ties, 19 tight signed decompositions,
19 approximate mixtures, 456 raw-query cases, three layout/fusion
counterexamples and 192 dictionary replacements/matching reads. All
arithmetic is rational. The checker uses dense oracles and is not a
runtime benchmark, a novelty proof or a new neural training run.

### Associative-memory goal: regression equivalence and the dimension audit

Previous turn classification: progress. The signed-support lemma survived
its exact checks, but indexing, fusion and sampling defeated its use as
a general access-cost separation. This iteration tests the remaining
local/shared construction against an explicitly specified regression
baseline. The original research goal and publication gates remain intact.

#### The reader is an ordinary weighted regression and residual decoder

Keep the notation and all hypotheses of the local/shared section. For
each current value column v, consider the strictly convex objective

    J(h) = (v-U*h)^T*C*(v-U*h)/2 + h^T*h/2.

Here C is the current symmetric positive-definite rounded local inverse
approximation; it need not be diagonal or entrywise nonnegative. The
normal equations give exactly

    h = (I+U^T*C*U)^(-1)*U^T*C*v,
    d = C*(v-U*h),
    read(b_q,u_q) = b_q^T*B^T*d + u_q^T*h.

The Hessian is I+U^T*C*U>0, so these equations identify the unique
minimizer, not merely a stationary point. Substituting d gives precisely
the earlier reader. Equivalently, Woodbury gives

    d = (C^(-1)+U*U^T)^(-1)*v.

These equalities hold for every value assignment and encoded query.
The local/shared retention and encoded-perturbation bounds therefore
apply without change to this regression baseline. Current values are
replaceable responses, not accumulating independent observations. There
is no statistical assumption on v in this equivalence.

For fixed C, replacing v_i changes U^T*C*v by
(C*U)_i^T*delta_v_i. This explicitly resolves the apparent distinction
between accumulating-observation caches and latest-response overwrites:
subtracting the old contribution and adding the new one implements the
same current normal equations. When a key edit changes C, every affected
row contribution must also be repaired, exactly as in the prior lemma.
The equivalence does not excuse omitting those changes.

#### A baseline family that includes the candidate's full cache choice

One baseline maintains the same C, current B,U,V and the four aggregates

    H=I+U^T*C*U,  t=U^T*C*V,
    Y=B^T*C*U,    S0=B^T*C*V.

It solves the small current system for h, and uses the displayed reader.
It may store X=C*U and Z=C*V, giving exactly the earlier algorithm and
its costs. Thus the current architecture has no separation from this
baseline family on any of the specified inputs, errors or horizons.
This is a constructive equivalence, not a lower bound against every
future extension or proof that every dynamic guarantee is published.

A second baseline omits persistent X and Z. For every affected row j in
the union repair ball S, reconstruct

    x_j = sum_l C_jl*u_l,  z_j = sum_l C_jl*V_l

from its sparse C row. Subtract the old row's outer-product contributions
to H,t,Y,S0 and add the new contributions. Outside S these products and
the input rows stay unchanged, so all aggregates remain canonical
functions of the current state. Reconstructing both versions costs
O(b^2*(r+v)) scalar operations when |S|<=b and every filter row has at
most b nonzeros. This fits within the earlier conservative key-edit
bound; it does not turn an exponential-size ball into a small one.

The second variant saves exactly n*(r+v) persistent scalar entries.
Local old products and input rows require O(b*(k+r+v)) temporary entries,
in addition to local polynomial-repair storage and the small solve. It
retains all keys, values, C coefficients/indices, incidence lists and
the four aggregates. Value-only updates reconstruct x_i at O(b*r) cost;
their total bound becomes

    O(b*r + k*b*v + r^2*v),

with H's inverse cached while keys stay fixed. The original cached
variant avoids the added b*r term. Query costs are identical. These are
two explicit storage/time choices, not uniform dominance in every
parameter regime. Preprocessing, encoder work, integer bit lengths and
small-solve costs remain charged as before. Exact canonical integer
sums give the same horizon-independent precision guarantee to either
baseline; renaming the aggregates does not change numerical behavior.

#### Shared-dimension dependence cannot be held constant for free

There is a further constraint on the earlier retention theorem's
parameters. For symmetric C>0, ||C||_infinity<=K implies C<=K*I. If
||u_i||_2<=M and H>=c*n*I, taking traces yields

    c*n*r <= tr(H) = r+tr(U^T*C*U) <= r+K*n*M^2,
    r*(c*n-1) <= K*n*M^2.

When c*n>1, r<=K*n*M^2/(c*n-1). Consequently M,K,c cannot all stay
fixed as both capacity and shared dimension grow without bound. The
conditional theorem remains valid for arbitrary finite r with its
stated parameters. It does not give a dimension-free guarantee by
suppressing their dependence on r. In particular, if c scales as 1/r
with M,K fixed, its stated Gamma=1+M^2*K/c scales as r. Holding the
stated error bound fixed can then require a larger polynomial depth
and repair ball. This concerns that sufficient bound; it is not a
universal capacity or error lower bound for other readers.

The constraint is attained by C=I and rows of U consisting of m copies
of each of the r standard basis vectors. Then n=m*r, M=K=1 and
H=(1+m)*I, so c=(1+m)/(m*r) and both trace inequalities are equalities.
This elementary control prevents small-r examples being reported as
dimension-independent architecture results.

#### Additional primary comparisons and their limits

[Malioutov--Johnson--Willsky, Section 3.2 and Appendix A](https://www.jmlr.org/papers/volume7/malioutov06a/malioutov06a.pdf)
already interprets inverse powers as weighted graph walks and proves
equivalent walk-summability criteria. The Neumann expansion, Proposition
1 and its displayed proof were inspected. Our norm contraction implies
convergence of the absolute power series. Local inverse expansion is
therefore established inference machinery. This comparison does not
assert their paper contains our complete dynamic cache/bit-cost bound.

[Liu--Chandrasekaran--Anandkumar--Willsky, Sections III--IV](https://arxiv.org/html/1105.1853v1)
already combines inference on a remaining graph with a small feedback
solve. Figure 4, Theorems 1/2 and their displayed block-matrix proofs
were inspected. The tree case is exact; approximate feedback inference
has exact means when its loopy inference converges. This occupies the
broad local-plus-global correction idea. It is not a direct reduction
of our model: a positive Gaussian block precision eliminates to a
*subtracted* positive-semidefinite term, while our covariance correction
is C^(-1)+U*U^T. Its batch inference bounds also do not alone certify
our dynamic costs or rounded finite-precision bound.

[van den Brand--Nanongkai--Saranurak, Theorems 4.1--4.2](https://arxiv.org/html/1905.05067v1)
provides general dynamic inverse tradeoffs for column/row and
entry/entry updates and queries, with field-operation accounting. The
statements and transformation-maintenance reduction were inspected;
the entire fast-algorithm proof was not audited. Its general inputs,
exact inverse queries and algebraic model differ from our local,
approximate reader. It cannot establish our finite-word costs, and
its general-case exponents are not a lower bound on our structured
inputs. Comparing only against those exponents would not establish
superiority over the strongest applicable structured method.

The earlier WISKI and Schur-preconditioner comparisons remain in force.
These readings do not establish that every part of the fully costed
dynamic statement is already published. They do establish that an
architectural claim based on regression, local inverse expansion and
a small global solve has not identified a substantive mathematical
gap. The explicit regression baseline above closes the claimed
architecture distinction on the present workload, regardless of the
remaining bibliographic question about the combined dynamic bound.

#### Decision and the next admissible problem

Park the current local/shared construction as a proposed main
contribution, while retaining its conditional theorem and certificates
as tools. Neither this baseline equivalence nor the trace inequality
is offered as a publishable theorem. No strongest-method advantage,
useful semantic encoder or defensible publication claim has emerged.

Before another architecture construction, specify an independent
mathematical query target and an admissible dynamic key/query class.
Exact stored-key recall already has a dictionary baseline; the prior
matching averages already have indexed and sampling baselines; the
present interpolation reader already has the regression baseline above.
A new candidate must improve a quantified tradeoff for the *same*
target and assumptions, allowing the baseline useful indexing, fused
evaluation and auxiliary summaries. Finite-precision, randomized-query
and horizon quantifiers must agree. Improving only a dense or
head-by-head implementation does not pass this gate. A raw-input
encoder theorem or a genuinely stronger dynamic solver theorem remains
a possible route, not an established open problem or promised paper.
Do not add training runs or enlarge a search framework on the strength
of the present identities.

Reproduce the new baseline controls with:

```bash
python3 -B scripts/run_regression_baseline_checks.py
```

The report `results/research/regression_baseline_checks.json` checks six
instances, 108 key edits, 48 value replacements, 12 independent dense
surrogate solves and eight trace-equality controls. A baseline retaining
neither full X nor Z is compared with the earlier full reconstruction,
the regression normal equations and the direct feature decoder. The
checker uses dense rational oracles; it makes no runtime, novelty,
architecture-superiority or new neural-training claim.

### Associative-memory goal: robust recall beyond a bounded-output assumption

Previous turn classification: progress. The exact regression baseline and
dimension audit changed which constructions remain worth investigating.
This iteration checks deletion prior art and proves a fixed-key robust
capacity bound that also applies after output clipping. It does not
establish the strongest-method advantage required by the active goal.

#### Deletion promises need an explicit counterfactual

[Exact Record Omission in Delta Attention, Theorem 1](https://arxiv.org/html/2609.06872v1)
already writes the difference between two affine recurrent trajectories
as a transported boundary difference plus accumulated forcing from
changed later updates. The statement, proof and write-rule corollary
were inspected. Subtracting the transported receipt is exact precisely
when the forcing sum vanishes. The paper explicitly does not give a
lower bound on every deletion algorithm. We must distinguish frozen
external writes from a model whose later write inputs depend on the
deleted history. Our local current-state constructions do not prove
counterfactual deletion for an entire contextual model.

[Forgetful Attention, Sections 3.1--3.4 and 4](https://arxiv.org/html/2607.12204v1)
already uses incremental/decremental support-vector optimization for
context-memory deletion. Its displayed objective, fixed-C deletion
definition and maintained-versus-batched distinction were inspected.
Returning to the retained-data optimizer is its exact-deletion target;
it is not a promise that every retained query prediction stays unchanged.
Its exact maintained path and approximate training path have different
claims. A new decay, erase or small regression mechanism alone would
not establish a new deletion theorem. We have not replicated either
paper's model experiments or audited all implementation details.

#### A homogeneous quadratic reader cannot robustly realize too many labels

Fix unit keys u_1,...,u_n in real p-dimensional Euclidean space and
0<delta<1. Suppose that for **every** sign assignment sigma in {+1,-1}^n
there is a symmetric matrix S_sigma such that

    sigma_i * x^T*S_sigma*x > 0
    whenever ||x-u_i||_2 <= delta, for every i.

Then

    n <= p*(2/delta^2-1)^2.

No bound on S_sigma or its precision is assumed. The matrix can be
obtained by an arbitrary nonlinear encoding of all labels. This is a
bound on robustly realizing independent labels at the fixed keys, not
merely on a particular update rule. For fixed delta its order is linear
in p, despite the quadratic number of matrix parameters.

Proof. Write r=delta^2. For a mixed-sign assignment, S has a negative
eigenvalue -a and a positive eigenvalue b, where a,b>0 are the extreme
magnitudes. For a positive key u, perturb by delta times a unit
eigenvector of -a, choosing its sign so the cross term is nonpositive.
The quadratic value at that point is

    u^T*S*u - 2*a*delta*|<u,w_negative>| - a*r > 0.

Consequently u^T*S*u>a*r. For a negative key the corresponding
perturbation toward a b-eigenvector gives u^T*S*u<-b*r. Since unit-key
values lie between -a and b, existence of both signs also implies
a>b*r and b>a*r.

Normalize and shift on the unit sphere:

    T = S/(a+b) - [(a-b)*r/(2*(a+b))]*I,
    gamma=r/2,  B=1-r/2.

At every key, sigma_i*u_i^T*T*u_i>gamma. Also

    ||T||_op = max(a,b)*(1+r)/(a+b)-r/2 <= B,

because max(a,b)/(a+b)<1/(1+r). For either constant sign assignment
use T=+gamma*I or -gamma*I instead. Thus for every assignment we have
a matrix T of operator norm at most B and margin at least gamma.

Set P_i=u_i*u_i^T and X_sigma=sum_i sigma_i*P_i. For independent
uniform signs, E||X_sigma||_F^2=n. The margin and trace Cauchy--Schwarz
give, for each assignment,

    gamma*n <= tr(T_sigma*X_sigma)
             <= sqrt(p)*B*||X_sigma||_F.

Take expectations and use E||X||_F<=sqrt(n); then
n<=p*B^2/gamma^2, which is the displayed bound. The argument also
works for Hermitian matrices and complex unit keys, using phases to
choose the eigenvector perturbation. Exact rational controls below
cover the real cases; the complex extension follows from the same
spectral and trace identities.

For clipped recall within error epsilon<1 of every binary label,
the signs must be correct throughout each ball, so the same bound
applies. It also covers a quadratic numerator divided by a denominator
that is strictly positive on the balls. A subsequent monotone scalar
decoder with a fixed decision threshold preserves that classification
reduction when its threshold can be absorbed into a quadratic form.

The domain and reader restrictions matter. The result is for a fixed
homogeneous quadratic form in the prescribed query coordinates. It
does not cover arbitrary nonlinear feature maps, query-dependent
routing, sums of separate quadratic ratios, affine-quadratic terms,
iterative attractor dynamics or complete networks containing such
operations. In particular it is not a lower bound on an entire
published transformer architecture. It bounds fixed-center robust
realizability, not general adversarial-PAC sample complexity.

The last trace/margin step is established machinery:
[Gyurik--van Vreumingen--Dunjko, Proposition 5 and Appendix A.3](https://arxiv.org/html/2105.05566v3)
already bounds fat-shattering by observable Frobenius norm and margin;
those displayed arguments were inspected previously. The present
geometric step obtains a bounded-norm margin representative from
full input-ball robustness, rather than assuming an output-norm bound.
Its precise priority remains to be established. An unsuccessful search
for this formulation is not evidence that it is new.

#### An explicit nonlinear reader escapes the quadratic restriction

For each unordered coordinate pair i<j, keep an independently replaceable
value V_ij in [-1,1]. There are N=p*(p-1)/2 such bindings. Define

    Q_V(q) = [sum_{i<j} V_ij*q_i^2*q_j^2]
             /[sum_{i<j} q_i^2*q_j^2].

Its numerator and denominator are quartic; equivalently, square the
p input coordinates and apply a quadratic pair reader. Use the fixed
unit key u_ij=(3*e_i+4*e_j)/5. At this key Q_V=V_ij for every assignment
of values. These are arbitrary independent bindings, not a low-rank
assumption on V. General queries can use all p coordinates.
However, the Jacobian of coordinate squaring has rank two at these
two-sparse keys, although it has rank p at dense queries. This code
does not supply a uniformly conditioned p-direction semantic encoder
or a theorem for arbitrary learned keys. Its special local geometry
must not be hidden behind the dimension of the ambient query vector.

For every real q with ||q-u_ij||_2<=delta=1/10, write t for the squared
mass outside coordinates i,j. Then t<=delta^2, the target pair mass
is at least (3/5-delta)^2*(4/5-delta)^2, and all other pair masses sum
to at most (1+delta)^2*delta^2+delta^4/2. The latter follows from
||q||_2<=1+delta and sum_{a<b outside}q_a^2*q_b^2<=t^2/2. Therefore

    |Q_V(q)-V_ij| <= 243/1225 < 1/5.

The denominator is positive in every such ball. This proves a fixed
noise-radius, all-assignment guarantee for N quadratic-in-p bindings.
At sufficiently large p it exceeds the preceding linear-in-p upper
bound for homogeneous quadratic readers. This is an explicit reader-
class separation, not a statement that higher-order memories or the
benefits of nonlinear features are new.

The representation is just an indexed binding table. It needs N value
words, plus an implicit pair-index formula; no N-by-p explicit codebook
is required. Indices need O(log(N)) bits; word-operation claims require
words large enough for those indices and the chosen value precision.
An overwrite or a clear-to-zero replaces one value word. Clearing
a value retains its address in the denominator: it is not removal of a
key followed by renormalization. Dense exact reads cost O(p^2*v) scalar
operations for value width v, with O(p) squaring work. On fixed bounded
rational grids, exact products and sums have current-state bit lengths
O(value_bits+4*query_bits+log(N)); overwrite history does not enter.
The N independent value words are still charged. This is not a storage
compression below the independent-value information bound.

The strongest baseline prevents an architectural-superiority claim.
On each key neighborhood, the largest two coordinate magnitudes identify
the pair: both target magnitudes are at least 1/2, all other magnitudes
are at most 1/10. A linear scan and one indexed value read return the
latest label **exactly** at O(p+v) cost. On arbitrary queries, the same
table directly evaluates Q_V with exactly the same formula and costs
as this reader. Thus the complete specified query function has no
separation from the indexed-table baseline. The neighborhood shortcut
alone does not compute Q_V exactly on general queries.

#### Prior-art and continuation decision

Higher interaction order increasing associative-memory capacity is
established. [Krotov--Hopfield, Section 2](https://arxiv.org/html/1606.01164v1)
already studies higher-order interactions, their capacity scaling and
iterative binary-state updates. The displayed model, update and capacity
formulas were inspected. Its pattern-recovery criterion differs from
independent scalar labels on Euclidean balls; we have not established
that this distinction makes the present bound new. Its architecture and
the prior product-key/dictionary baselines remain part of the comparison.

[Arriaga--Vempala, *Robust Concepts and Random Projection*](https://doi.org/10.1007/s10994-006-6265-7)
uses Euclidean input-ball robustness and a unit-ball data support in
Definition 1. After the author's PostScript failed to render, the
author-uploaded journal text was retrieved; Definition 1 and Section
4.3's Theorem 9 with its proof were inspected. The latter learns robust
balls of radius at most one by random projection. It does not state
the arbitrary indefinite-quadratic fixed-key capacity bound above.
This checks an important model comparison; it does not settle priority
across that paper's other results or the broader robust-learning literature.

The candidate worth auditing is the precise robust-capacity theorem,
including clipping and no a priori weight bound, and whether a sharper
noise-radius dependence or broader reader theorem is known or feasible.
The general message that nonlinear features can increase capacity is
already occupied. Do not market the pair table as a new architecture,
infer a practical advantage from an asymptotic class separation, or
claim a paper from these elementary constructions. In particular the
active goal's strongest-applicable-method criterion remains unmet.
Continue the exact prior-art and structural proof work before any
additional architecture experiment, training run or manuscript.

Reproduce finite proof controls with:

```bash
python3 -B scripts/run_robust_reader_checks.py
```

The report `results/research/robust_reader_checks.json` records 720
whole-ball normalization checks using exact distance to quadratic zero
lines, including very large and small reader rescalings; 320 exact sign
assignments for the trace identity; 74 exhaustive pair-label assignments
and 4,454 rational noisy queries; and 96 replacements/zero-clears. The
whole-ball theorem is proved above; finite probes alone do not prove it.
The checker also verifies the indexed top-two escape on every tested
key neighborhood. No neural or runtime experiment was run.

### Associative-memory goal: affine quadratics and normalized queries

Previous turn classification: progress. The robust homogeneous capacity
bound, quartic countercontrol and indexed baseline were saved. This
iteration establishes a direct prior-art boundary and extends the proof
to linear terms and biases, including a separate proof on the unit sphere.
The full architecture/novelty/publication goal is still unachieved.

#### The homogeneous geometric step is published prior art

[Fawzi--Fawzi--Frossard, Lemma 6.1](https://arxiv.org/pdf/1502.02590)
already perturbs toward an opposite-sign extreme eigenvector and bounds
distance to a homogeneous quadratic decision boundary by the square
root of its value divided by that eigenvalue magnitude. The lemma and
its displayed proof, Section 6's assumptions and Theorem 6.2 were
inspected. The latter additionally uses an extreme-eigenvalue ratio
condition and a second-moment distinguishability measure. Our earlier
geometric step is the same mechanism as the lemma; combined with the
established norm/margin argument, the homogeneous capacity bound is
an elementary consequence, not a standalone novelty claim. This does
not establish that the affine or spherical extensions below are new.

#### Two domains, with explicit affine-quadratic bounds

Let p>=2, let u_1,...,u_n be fixed real unit keys, and let 0<delta<1.
For every binary assignment suppose there is an arbitrary real
affine-quadratic function

    f(x)=x^T*S*x+b^T*x+c,  S symmetric,

with the prescribed strict sign on the whole specified neighborhood
of every key. Coefficients and the algorithm that obtains them are
unrestricted. There are two different neighborhood specifications:

* Full balls: all real x with ||x-u_i||_2<=delta. Then

      n <= 300*(p+8)/delta^4.

* Spherical caps: all unit x with ||x-u_i||_2<=delta. Put
  C0=64*660^2=27,878,400. Then

      n <= 3*(p/4+2)*(2*C0)^2/delta^24.

The usual quadratic feature-dimension bound
n<=(p+1)*(p+2)/2 can be taken in the minimum in both cases.
For p=1 there are only two distinct unit keys, so n<=2.
Both bounds are linear in p at any fixed positive noise radius.
The spherical constant and radius dependence are extremely loose;
they do not establish a useful practical memory advantage. No optimal
dependence on delta or priority is claimed.

Here is a uniform proof of both bounds. For a mixed assignment, let
u be a positive key with smallest f-value and v a negative key with
largest f-value. Write

    g=f(u)-f(v)>0, w=lambda_max(S)-lambda_min(S),
    B=||b||_2, D=w+B,
    m=(lambda_max(S)+lambda_min(S))/2, T=S-m*I.

D>0, since otherwise f is constant on the unit sphere. On unit inputs,
f(x)=x^T*T*x+b^T*x+c+m and ||T||_op=w/2.

#### Full-ball gap: spherical averaging handles the linear term

Let r=delta^2 and alpha=r/4. If B>0 put beta=b/B. For any unit u the
vector mu=(1-alpha)*u-alpha*beta belongs to the convex hull of its
spherical cap: ||mu||<=1 and <mu,u> >= 1-r/2. The convex hull equals
the unit ball cut by that half-space; its cross-sections prove this
identity for p>=2. Hence a random unit Q in the cap can have mean mu.
With h=Q-u,

    E h=-alpha*(u+beta), E||h||^2=2*alpha*(1+<u,beta>)<=4*alpha,
    |E[Q^T*T*Q]-u^T*T*u|<=4*alpha*w.

The last inequality follows by expanding into 2*u^T*T*E h+E[h^T*T*h].
Positivity throughout the cap gives

    f(u)>alpha*(B+b^T*u)-4*alpha*w.

Use mean (1-alpha)*v+alpha*beta at the negative key to get

    f(v)<-alpha*(B-b^T*v)+4*alpha*w.

Since g=b^T*(u-v)+u^T*T*u-v^T*T*v, and the last two terms differ
by at most w, subtraction yields

    (1-alpha)*g>alpha*(2*B-9*w).

Separately, averaging f at u plus/minus delta times a minimum
eigenvector gives f(u)>-delta^2*lambda_min(S). Averaging at v with
a maximum eigenvector gives f(v)<-delta^2*lambda_max(S). Thus
g>delta^2*w. If B>=9*w the preceding spherical estimate gives
g>delta^2*B/4 and D<=10*B/9. Otherwise D<10*w and the eigenvector
estimate suffices. In either case

    g>delta^2*D/10.

For B=0 the eigenvector estimate alone proves the same conclusion.
This proves the full-ball gap without assuming bounded weights,
well-conditioned S, or a bound on the original constant term.

#### Spherical gap: approximate extrema force an almost constant form

The full-ball eigenvector step is unavailable for a sphere-only reader.
For example f(x)=||x||^2-1+x_1/100 has the correct opposite signs
around +e_1 and -e_1 on small spherical caps. At x=(9/10)*e_1 its
value is -181/1000, so the same form fails the positive full-ball
promise. An off-sphere proof cannot simply be reused after normalization.

We instead show g>=D*delta^12/C0. If g>=D*delta^4 this is immediate.
Otherwise put eta=sqrt(g*D). Along any unit-speed great circle the
second and fourth derivatives of f have magnitude at most 2*D and
8*D, respectively; the radial term m*I is constant there. Taylor's
formula at u, with steps sqrt(g/D), bounds the tangent gradient by
2*eta. A symmetric Taylor average with steps (g/D)^(1/4) bounds the
tangent Hessian below by -8*eta/3. Both steps are less than delta,
and their chord distances are at most their arc lengths, so they stay
inside the cap. The same argument for -f at v gives the opposite bounds.

In detail, f(u)<=g and f(v)>=-g. The one-sided Taylor estimate gives
|f'(0)|<=g/t+D*t=2*eta at t=sqrt(g/D). The symmetric estimate gives
f''(0)>=-2*g/t^2-(2*D/3)*t^2=-8*eta/3 at t=(g/D)^(1/4).
These estimates hold in every tangent direction.

Define the sphere Lagrange quantities

    lambda_u=u^T*S*u+b^T*u/2, r_u=S*u+b/2-lambda_u*u,
    lambda_v=v^T*S*v+b^T*v/2, r_v=S*v+b/2-lambda_v*v.

Then ||r_u||,||r_v||<=eta. On u-perpendicular directions S-lambda_u*I
is bounded below by -4*eta/3. On v-perpendicular directions
S-lambda_v*I is bounded above by 4*eta/3. Set d=1-<u,v> and z=u-v.
The negative key cannot belong to the positive cap, so d>delta^2/2.
Direct subtraction gives

    d*(lambda_u-lambda_v)=g-<u,r_v>+<v,r_u>,
    (S-lambda_u*I)*z=r_u-r_v+(lambda_u-lambda_v)*v.

Since g<D and g<=eta, these imply

    |lambda_u-lambda_v|<=3*eta/d,
    ||(S-lambda_u*I)*z||<=7*eta/d.

For completeness, the following conservative constants make the global
step explicit. Put A=S-lambda_u*I. For any unit x, subtract
(<x,u>/d)*z to make x perpendicular to u, or add (<x,v>/d)*z to
make it perpendicular to v. Each adjusted vector has norm at most
4/d; the scalar adjustment is at most 1/d. Using ||z||<=2, the
quadratic-form error from either adjustment is at most

    (10/d^2)*||A*z||<=70*eta/d^3.

The tangent lower and upper bounds therefore put every Rayleigh
quotient of A between -338*eta/(3*d^3) and 482*eta/(3*d^3). Hence

    ||A||_op<=161*eta/d^3,
    w<=322*eta/d^3,
    B<=2*||r_u||+2*||A||_op<=338*eta/d^3,
    D<=660*eta/d^3.

Squaring and using d>delta^2/2 gives
g>=D*d^6/660^2>D*delta^12/(64*660^2), as required.
This argument is also a quantitative version of a simple rigidity fact:
a quadratic function on the sphere cannot have a local minimum and
a distinct local maximum at the same value unless it is constant on
the sphere. At exact extrema the residual identities give equal
Lagrange multipliers and A*z=0; the two tangent semidefinite conditions
then force A=0 and b=0. This fact is not claimed as new optimization theory.

#### Converting either gap into a capacity bound

Choose tau=(f(u)+f(v))/2 and divide the on-sphere coefficients by D.
Every key then has the required sign with margin at least gamma, where

    gamma=delta^2/20                 for full balls,
    gamma=delta^12/(2*C0)            for spherical caps.

The shifted constant obeys |c+m-tau|<=w/2+B, because tau lies between
two values on the unit sphere. In the Hilbert feature representation

    Phi(x)=(x*x^T,x,1),
    W=(T,b,c+m-tau)/D,

we have ||Phi(u_i)||^2=3 and ||W||^2<=p/4+2. Constant sign assignments
use just the final scalar coordinate. The same independent-sign trace
argument as before, now in the direct-sum Hilbert space, gives

    n<=3*(p/4+2)/gamma^2.

Substitution proves both displayed theorems. The feature-dimension
alternative follows from linear dependence of more feature vectors
than their ambient dimension: assigning signs to agree with the
nonzero coefficients of a dependence contradicts strict separation.

#### Scope, primary comparison and publication gate

The spherical result covers normalized query coordinates, arbitrary
linear terms and biases. Clipping does not evade a sign-realizability
bound. A positive quadratic denominator and a fixed scalar decision
threshold reduce to an affine-quadratic numerator; sums with a common
denominator can be fused before applying the result. Different ratio
denominators, nonlinear feature maps, routing and iterative dynamics
remain outside it. This is a restriction on the effective reader in
specified coordinates, not a lower bound on all complete ML models.

The final bounded-margin step is established learning theory.
[Ivanov--Jung--Naszodi, Theorem 3.1 and its proof](https://arxiv.org/html/2607.10751v1)
gives the relevant general expanded-half-space bound through balanced
signed sums and Rademacher type. Those statements and proof were
inspected; the paper's expanded-ball domain is not an arbitrary
indefinite-quadratic reader. The normalization/gap argument still needs
its own precise priority audit. Different domains alone do not establish
a new theorem of publishable significance.

There is also substantial prior optimization theory:
[Martinez, *Local Minimizers of Quadratic Functions on Euclidean Balls and Spheres*](https://doi.org/10.1137/0804009).
Its primary abstract was inspected, including its characterization of
local nonglobal minimizers; the full proof was not accessed or audited.
Do not claim the rigidity fact as new or infer from the abstract that
our quantitative robust-capacity bounds are absent from the literature.

The advance is a complete affine-quadratic extension on both domains,
not a new efficient architecture. Its constants, especially on the
sphere, are far too weak to justify a practical capacity claim. The
earlier quartic table and indexed-dictionary comparisons still apply;
no strongest-method access, storage or retention advantage has been
proved. The full goal remains active. Audit the exact geometry and
sharper bounds before training or preparing a contribution manuscript.

Reproduce the new controls with:

```bash
python3 -B scripts/run_affine_robust_checks.py
```

They check rational cap measures and their moments, derivative bounds
with arbitrarily scaled radial terms, residual-transfer identities,
the full-ball gap algebra, all mixed labels on small orthogonal key
sets with linear terms/biases, the affine lifted sign-sum identity and
the explicit spherical/full-ball counterexample. These finite controls
support the proofs; they do not establish priority or verify the uniform
theorem merely by sampling query neighborhoods.

### Associative-memory goal: the capacity frontier and uniform constants

Previous turn classification: progress. The direct homogeneous prior art
and the affine/full-ball and sphere-only proofs are preserved. This pass
checks whether optimization error bounds supply the missing uniform
guarantee, and establishes a constructive capacity benchmark. The full
architecture, novelty and publication objective remains unachieved.

#### Fixed-problem exponents do not give uniform memory guarantees

[Gao--Liu--Chen--Yuan, Theorem 2.1](https://arxiv.org/html/1611.08781v2)
states a sphere-constrained quadratic gradient inequality with exponent
3/4. Its preliminary equations (14)--(15), theorem, displayed proof
cases and Example 2.1 were inspected. The constants and neighborhood
are for the specified problem; the proof uses nonzero singular values.
Example 2.1 already exhibits quartic flatness on a circle. This is
relevant prior art, not a coefficient-uniform capacity bound. We have
not independently verified every proof case.

[Jiang--Li, Definition 1.1 and Theorem 3.13](https://arxiv.org/html/1911.11955v2)
give a global error bound to the optimizer set of a fixed quadratic
problem on the unit ball. The exponent is 1/4 in the specified ill case
and 1/2 otherwise; the multiplicative constant is problem dependent.
The definitions, main statement and displayed proof structure were
inspected, without independently auditing every supporting lemma.
Distance to global optimizers differs from our opposed classification
caps. An exponent alone does not resolve the comparison.

The absence of uniform constants is visible without relying on either
paper's proof. For 0<epsilon<=1/2 consider on the two-dimensional unit
ball

    f_epsilon(x)=epsilon*x_1^2+x_2^2.

The coefficient matrix has operator norm 1 and spectral width at least
1/2. The optimizer set is {0}. At x=e_1 the distance to that set is 1
and the objective gap is epsilon. Consequently any global bound
dist(x,{0})<=tau*f_epsilon(x)^(1/2) requires
tau>=epsilon^(-1/2). Even weakening the exponent to 1/4 requires
tau>=epsilon^(-1/4). These examples are not in the paper's ill case;
the latter calculation simply shows that weakening an exponent does
not by itself make its constant uniform.

Similarly, on the unit sphere in three dimensions put

    F_epsilon(x)=epsilon*x_2^2+x_3^2, x_star=e_1.

Its coefficient operator norm and spectral width are both 1. For a
unit q=(a,b,0), a>0, b>0, elementary tangent projection gives

    F_epsilon(q)=epsilon*b^2,
    ||grad_sphere F_epsilon(q)||^2=4*epsilon^2*a^2*b^2,
    (F_epsilon(q)^(3/4)/||grad_sphere F_epsilon(q)||)^4
        =b^2/(16*epsilon*a^4).

In any fixed positive neighborhood of e_1 choose such a q with fixed
b>0. The displayed ratio diverges as epsilon tends to zero. Hence a
tangent-gradient inequality of this form cannot have both a common
constant and a common positive neighborhood over this coefficient-
normalized family. Fixed-instance gradient inequalities remain
compatible with the example. Neither comparison proves novelty of
our earlier geometric estimates.

#### A spherical counterfamily prevents a quadratic center-gap bound

The previous sphere proof used a very weak twelfth power of the noise
radius. Improving that proof directly to a second power is impossible
for an arbitrary given classifier. Here is a uniform, explicit example;
it is related to the established quartic flatness phenomenon, not a
standalone novelty claim.

For p>=2 and 0<t<=1 define

    a=(1-t^2)/(1+t^2), s=2*t/(1+t^2), d=1-a,
    delta=t/4, u=e_1, v=(a,s,0,...,0),
    f_t(q)=(q_1-a)^2-d^2/4.

Its quadratic spectral width is 1 and its linear coefficient norm is
2*a, so D=1+2*a. The center values are 3*d^2/4 and -d^2/4,
respectively, with gap g=d^2.

For any unit q with ||q-u||<=delta, q_1>=1-delta^2/2 and
delta^2/2<=d/32. Therefore

    f_t(q)>=(31*d/32)^2-d^2/4=705*d^2/1024>0.

For a unit q with ||q-v||<=delta, decompose
e_1=a*v+(e_1-a*v), where the second summand has norm s.
The identity <v,q-v>=-||q-v||^2/2 gives

    |q_1-a|<=delta^2/2+s*delta<=9*d/32,
    f_t(q)<=-175*d^2/1024<0.

Thus the strict opposite signs hold on both entire spherical caps.
Nevertheless

    g/(D*delta^4)=1024/[(1+t^2)*(3-t^2)]
                lies between 256 and 1024/3.

In particular g/D is of order delta^4. No uniform positive lower
bound g>=c*D*delta^r for r<4 can hold for every such classifier.
This does not prove a universal fourth-power lower bound, nor does
it prove a fourth-power capacity obstruction. The same labeling may
have another, much better classifier: for these two keys a separating
linear form suffices. A poor witness is not a bound on the best
witness for each labeling.

#### A standard quadratic table retains order min(p^2,p/delta^2) labels

Before pursuing sharper upper bounds, record a constructive lower
bound with all labels independent. This is an ordinary matrix-unit
binding table; no new architecture or efficient factorization is
being claimed.

Partition p>=2 coordinates into consecutive blocks of size at most
r, using as many full blocks as possible, where 2<=r<=p. For every
unordered pair i<j in one block store an arbitrary V_ij in [-1,1].
There are N such pairs. Set

    u_ij=(3*e_i+4*e_j)/5,
    S_ii=0, S_ij=S_ji=(25/24)*V_ij within a block,
    S_ij=0 otherwise, R(q)=q^T*S*q.

Then R(u_ij)=V_ij exactly. The encoder is the identity on R^p and
the values have no low-rank restriction. The prescribed keys are
special two-coordinate codewords; this is not a claim about arbitrary
learned keys or a conditioned nonlinear encoder.

The maximum absolute row sum gives

    ||S||_op<=(25/24)*(r-1).

At each stored key, the two corresponding coordinates of S*u_ij
have total squared magnitude at most (25/24)^2. Each of the other
at most r-2 coordinates has magnitude at most (25/24)*(7/5).
It follows that

    ||S*u_ij||_2<=(35/24)*sqrt(r-1).

For every full Euclidean perturbation ||h||_2<=delta, expansion now
proves, uniformly over all bounded value assignments,

    |R(u_ij+h)-V_ij|
        <=(35/12)*delta*sqrt(r-1)+(25/24)*delta^2*(r-1).

If delta^2*(r-1)<=1/64, the right side is at most

    35/96+25/1536=195/512<1/2.

Clipping R to [-1,1] cannot increase its error. For vector values the
same guarantee holds coordinatewise, charging every value coordinate.
The statement covers the full ambient ball and hence also its unit-
sphere restriction.

For 0<delta<=1/8 choose r=min(p,1+floor(1/(64*delta^2))). Counting
the full blocks gives

    N>=p*(r-1)/4
      =(p/4)*min(p-1,floor(1/(64*delta^2))),
    N<=p*(r-1)/2.

Indeed floor(p/r)>=p/(2*r), and each full block contributes
r*(r-1)/2 pairs. The remainder can only increase the lower bound.
This establishes capacity Omega(min(p^2,p/delta^2)) throughout the
stated regime, with explicit constants and a uniform proof rather
than extrapolation from finite examples.

#### Storage, overwrites, precision and a matching indexed baseline

Store the N independent V entries once, with the common scale 25/24
implicit. The block layout and pair addresses require p,r and a fixed
indexing rule, not a learned dense codebook. A word-address model must
pay at least log_2(N) address bits. Initializing arbitrary contents
costs Omega(N) value writes. Given an address, replacing its scalar
value takes one value-word write; zero-clear retains the address and
does not delete a key or renormalize a query. Canonical replacement
has no accumulation with the stream horizon.

The complete quadratic query costs O(p+N) scalar work, or O(p+N*v)
for value width v. With bounded dyadic value and query coordinates,
exact evaluation uses O(value_bits+2*query_bits+log(N)) arithmetic
bits, plus fixed rational constants. These arithmetic costs are not
free real-number operations. Other finite-precision implementations
must charge their own rounding error within the slack below 1/2.

A table with the same layout evaluates exactly the same quadratic
function on every fresh query for exactly the same asymptotic cost.
For local codeword recall it does better: since delta<=1/8, the two
intended coordinate magnitudes are at least 3/5-delta>delta, while
every other coordinate magnitude is at most delta. Scanning q for
its two largest magnitudes identifies the pair, then a direct lookup
returns the latest value exactly in O(p+v) work. Both tables use the
same independent values and implicit metadata. This construction
therefore supplies a capacity benchmark, not a strongest-baseline
memory, access or retention advantage.

#### Consequence for the active research direction

For the all-binary-label robust affine-quadratic model, the currently
proved bounds leave the following window (p>=2, 0<delta<=1/8):

| Neighborhood | Constructive lower capacity | Current upper capacity |
| --- | --- | --- |
| Full Euclidean ball | Omega(min(p^2,p/delta^2)) | O(min(p^2,p/delta^4)) |
| Unit-sphere cap | Omega(min(p^2,p/delta^2)) | O(min(p^2,p/delta^24)) |

The constants in the upper bounds, especially the sphere bound, remain
very large. This table is not a claim that the displayed powers are
sharp or absent from the literature. The general norm/margin step is
already known. A gap estimate for every possible witness and optimal
robust shattering capacity are different targets; the cap counterfamily
cannot substitute for a capacity lower bound.

The next mathematical target is a direct, substantially sharper
capacity theorem, or a construction disproving the natural candidate
order min(p^2,p/delta^2), with a primary-source priority audit. Merely
improving an intermediate constant or rediscovering the block table
does not complete the goal. Any proposed architecture must additionally
beat the strongest indexed, routed, kernel or regression method on the
same complete workload, including fresh queries, precision and updates.
No such advantage, publishable novelty or new neural result is claimed
in this pass. The full goal remains active; no training was launched.

Reproduce the new exact controls with:

```bash
python3 -B scripts/run_robust_frontier_checks.py
```

The report records 148 value assignments, 1,656 exact key/gradient
checks, 6,224 noisy-query and indexed-baseline controls, 28 capacity
counts, 96 point replacements/zero-clears, eight opposed-cap scales
and twelve fixed-problem constant controls. Rational neighborhood
probes check the formulas independently of the whole-neighborhood
proofs; they do not establish those proofs by sampling or settle
priority.

### Associative-memory goal: sharp order for robust numerical quadratic recall

Previous turn classification: progress. The standard block lower bound,
opposed-cap counterfamily and fixed-problem constant audit changed the
next proof target. Here we close the capacity-order gap for numerical
recall in two reader classes. This is a supporting theorem with unresolved
priority, not completion of the architecture/strongest-baseline goal.

#### Specify the value guarantee before bounding capacity

Fix p>=2 keys u_1,...,u_n of Euclidean norm at most 1, independently
of their stored values. Fix a>epsilon>=0 and delta>0. For every
y in {-1,+1}^n suppose some affine quadratic

    f_y(q)=q^T*S_y*q+b_y^T*q+c_y, S_y symmetric,

satisfies |f_y(q)-a*y_i|<=epsilon on the entire radius-delta ball
around u_i. Encoding the values into the coefficients may be arbitrary,
nonlinear and computationally expensive. There is no coefficient norm
or precision assumption. Then, writing r for the affine dimension of
the keys, we prove

    (n-1)^2/n <= (r+1)*epsilon^2/[delta^2*(a-epsilon)^2],
    n <= 2+(p+1)*epsilon^2/[delta^2*(a-epsilon)^2].

The first displayed bound also holds with the smaller geometry quantity
Lambda defined below in place of r+1. These statements require two-sided
numerical accuracy. Merely assigning the correct sign is a weaker
promise, and the earlier sign-only bounds are not replaced by this one.

Clipping to [-1,1] is covered when a+epsilon<1: a clipped value within
epsilon of either +/-a lies strictly inside the clipping interval,
so the raw value must equal it. A memory supporting every value in
[-1,1] can be tested on these interior values. Clipping only binary
endpoints has a different behavior. For example f(q)=8*q_1, at keys
e_1,-e_1 and radius 1/4, clips to the correct +/-1 exactly, although
its gradient is nonzero. The gradient bound below cannot be inferred
from that endpoint-only promise.

#### Full-ball proof: local accuracy bounds an empirical gradient norm

For any unit direction w, exact quadratic central differences give

    f_y(u_i+delta*w)-f_y(u_i-delta*w)
        =2*delta*<2*S_y*u_i+b_y,w>.

Both values are within epsilon of the same target. Therefore

    ||2*S_y*u_i+b_y||_2 <= K, K=epsilon/delta.

This step uses numerical accuracy, not a global Lipschitz or weight
bound. Center the keys and define

    mu=(1/n)*sum_i u_i, x_i=u_i-mu,
    G=sum_i x_i*x_i^T, r=rank(G),
    t_y=2*S_y*mu+b_y, ell_i=x_i^T*G^dagger*x_i.

The pseudoinverse acts only on the span of the centered keys; no
conditioning assumption is needed. We have sum ell_i=r,
||x_i||<=2 and Tr(G)<=n. Since sum x_i=0, the gradient bounds imply

    E_y := 4*Tr(S_y^2*G)+n*||t_y||^2
          =sum_i ||2*S_y*u_i+b_y||^2 <= n*K^2.

Choose independent uniform signs y_i, put bar_y=(1/n)*sum_i y_i,
z_i=y_i-bar_y, and set

    M_y=sum_i z_i*x_i*x_i^T, v_y=sum_i z_i*x_i.

The constant coefficient vanishes in this centered sum. Weighted
Frobenius Cauchy--Schwarz gives the deterministic identity and bound

    sum_i z_i*f_y(u_i)=Tr(S_y*M_y)+<t_y,v_y>,
    |sum_i z_i*f_y(u_i)|^2 <= E_y*D_y,
    D_y=(1/4)*Tr(G^dagger*M_y^2)+||v_y||^2/n.

For clarity, use the direct-sum pairs
(2*S_y*G^(1/2),sqrt(n)*t_y) and
((1/2)*M_y*(G^dagger)^(1/2),v_y/sqrt(n)). Their inner product is
the displayed sum because M_y is supported on the range of G.

The independent-sign second moment can be computed exactly:

    E D_y = Lambda
          =(1/4)*sum_i ||x_i||^2*ell_i+3*Tr(G)/(4*n)
          <= r+1.

Indeed M_y=sum_i y_i*(x_i*x_i^T-G/n), so its squared-matrix moment
has trace against G^dagger equal to
sum_i ||x_i||^2*ell_i-Tr(G)/n. Also E||v_y||^2=Tr(G).
This explains the positive 3/4 coefficient; dropping the centering
correction would give a different expression.

Numerical accuracy gives a lower bound for each assignment, even
though f_y may depend arbitrarily on all signs:

    sum_i z_i*f_y(u_i)
        >=(a-epsilon)*n*(1-bar_y^2).

Here z_i has the same sign as y_i and sum |z_i|=n*(1-bar_y^2).
Taking expectations, using E[n*(1-bar_y^2)]=n-1, and applying
Cauchy--Schwarz once more yields

    (a-epsilon)*(n-1) <= sqrt(n)*K*sqrt(Lambda).

Squaring proves (n-1)^2/n<=Lambda*K^2/(a-epsilon)^2. Finally
(n-1)^2/n=n-2+1/n>n-2, which proves the simpler stated upper
bound. Singular, nearly singular and translated key sets are included.

#### Sphere-only proof for a quadratic plus bias

Now require unit keys and only unit queries at chord distance at most
delta, where 0<delta<=1. Consider readers q^T*S_y*q+c_y, with no
linear term in these coordinates. Absorb the bias into
T_y=S_y+c_y*I. On all queried points the same reader is q^T*T_y*q.
For every binary assignment assume the same numerical accuracy
|q^T*T_y*q-a*y_i|<=epsilon. If r_0 is the linear rank of the unit
keys, then

    n <= r_0*[(a+epsilon)^2+4*epsilon^2/(3*delta^2)]
             /(a-epsilon)^2.

To prove it, fix a unit tangent direction w perpendicular to u_i.
Choose sin(t)=delta/2, cos(t)=sqrt(1-delta^2/4), and use the two
unit queries q_+=cos(t)*u_i+sin(t)*w and
q_-=cos(t)*u_i-sin(t)*w. Their chord distances to u_i are at most
delta. The exact odd-part identity gives

    f_y(q_+)-f_y(q_-)=4*cos(t)*sin(t)*u_i^T*T_y*w,
    ||(I-u_i*u_i^T)*T_y*u_i||
        <=epsilon/[delta*sqrt(1-delta^2/4)]
        <=2*epsilon/(sqrt(3)*delta).

The radial component is f_y(u_i), whose magnitude is at most
a+epsilon. Hence ||T_y*u_i||^2<=L^2, where
L^2=(a+epsilon)^2+4*epsilon^2/(3*delta^2).

Put G_0=sum_i u_i*u_i^T and M_0=sum_i y_i*u_i*u_i^T.
The same weighted Frobenius inequality, without centering, gives

    n*(a-epsilon) <= |Tr(T_y*M_0)|
        <=sqrt(n)*L*sqrt(Tr(G_0^dagger*M_0^2)).

Since every key has norm 1, independence gives
E Tr(G_0^dagger*M_0^2)=sum_i u_i^T*G_0^dagger*u_i=r_0.
Taking expectations and squaring proves the result.

This argument does not include a separate linear term on the sphere.
Nor does it include a state-dependent rational denominator, routing,
arbitrary nonlinear encoders or iterative retrieval. These cases
remain distinct mathematical targets. The full-ball theorem above
does include all linear terms and biases.

#### Matching the order at a fixed numerical accuracy

For arbitrary values in [-1,1], fix error tolerance 1/4 and take
0<delta<=1/16. In the previous block construction choose

    r_block=min(p,1+floor(1/(256*delta^2))).

The same whole-ball expansion now bounds the error by

    35/192+25/6144=1145/6144<1/4.

It retains at least
p*min(p-1,floor(1/(256*delta^2)))/4 independent values, which is
Omega(min(p^2,p/delta^2)) with absolute constants.
The upper theorems use interior test values +/-1/2 and epsilon=1/4,
so they cover both raw readers and clipping to [-1,1]. Together with
the elementary feature-dimension bound, they give:

| Numerical reader and query domain | Upper bound | Capacity order |
| --- | --- | --- |
| Affine quadratic, full balls | min((p+1)*(p+2)/2, 2+(p+1)/delta^2) | Theta(min(p^2,p/delta^2)) |
| Quadratic plus bias, unit-sphere caps | min(p*(p+1)/2, p*(9+4/(3*delta^2))) | Theta(min(p^2,p/delta^2)) |

The capacity order is for the stated fixed accuracy and noise range,
maximized over fixed key sets, with independent bounded numerical
values. It is not a theorem about the weaker sign-only model, a
precise optimal constant, arbitrary transformer architectures or a
general upper bound for affine-quadratic sphere readers.

The lower construction keeps the identity encoder in R^p and actual
cross-coordinate pair values. As previously accounted, it stores one
finite-precision word per independent scalar value, plus implicit
block/pair metadata; a given-address overwrite changes one value
word. Complete quadratic queries cost O(p+N*value_width), with
precision and initialization charged. Its indexed-table baseline
evaluates the same full function and gives cheaper exact local
recall. The theorem therefore does not by itself yield an algorithm
beating that baseline. Upper bounds allow arbitrary real states,
under the continuous-neighborhood promise, so restricting state
precision alone cannot evade them. A promise only on a finite query
grid is different and needs its own interpolation/rounding audit;
we do not infer a continuous guarantee from such a grid. Lower bounds
use ordinary bounded dyadic values and a common fixed rational scale.
Bounded evaluation roundoff must be added to the numerical tolerance.
Any overwrite procedure promising all independent current values on
the fixed keys inherits the snapshot bound at every horizon. Adapting
the keys or encoder to the stored labels changes the model.

#### Primary comparison, remaining novelty audit and next action

The weighted norm and independent-sign method are classical.
[Bartlett--Mendelson, Lemma 22 and its proof](https://jmlr.csail.mit.edu/papers/volume3/bartlett02a/bartlett02a.pdf)
give the Hilbert feature trace bound underlying the final step.
Those pages and the surrounding kernel-margin theorem were inspected.
Our gradient-to-empirical-norm reduction and exact centered moment
identify the relevant geometry; they should not be presented as a
new general Rademacher inequality.

Gradient-based complexity control is also established. The publisher
abstract of [Vieira et al., the Q-norm/minimum-gradient paper](https://doi.org/10.1109/TNN.2008.2000442)
was inspected; the full journal proof was not accessed. In the author's
[2006 thesis, Chapter 4, Sections 4.1--4.4](https://repositorio.ufmg.br/bitstream/1843/BUOS-8CTH6W/1/douglas_alexandre_gomes_vieira.pdf),
the quadratic-form complexity and gradient regularizer, including its
linear SVM comparison, were inspected. These are close methodological
prior art. Exact equivalence to, or coverage of, the two capacity
statements above remains to be established, rather than inferred from
titles or a keyword search.

[Falk--Karlov's robust quadratic separation paper](https://doi.org/10.1016/S0305-0548(99)00134-3)
defines robustness by distance to the zero set. The primary abstract
and available definition/algorithm excerpts were inspected; its full
proof was not audited. That sign-separation objective differs from
two-sided numerical recall, and the present theorem does not settle
its all-label robust capacity.

A current transformer comparison is
[MLPs are Hebbians, Appendix B.10](https://arxiv.org/html/2607.10034v1).
The noisy-margin definitions, Theorem B.62's displayed argument and
Corollaries B.63--B.65 were inspected. Their isotropic-value decoding
guarantee and sufficient noise ceiling in equation (78) depend on
fact count, dimension and a kernel Lipschitz constant. They are not
an arbitrary-value, fixed-radius numerical-retention guarantee.
The concentration argument needs its conditional-independence scope
retained; no uniform adaptive-query extension or independent audit
of the entire paper is being claimed here.

This pass proves a structurally sharper numerical-capacity statement
and explains the role of the extra value promise. Publication novelty
and an advantage over the strongest applicable architecture remain
unestablished. The next work is to audit the precise empirical-gradient
capacity result against this prior art, and determine whether it can
support a new constructive retention/access tradeoff on complete
queries that an indexed or routed baseline cannot already match.
Keep the original goal active. No training or outreach was performed.

Reproduce the exact checks with:

```bash
python3 -B scripts/run_numerical_capacity_checks.py
```

They independently verify rational pseudoinverses through all Penrose
identities, centered sign moments and weighted Cauchy inequalities on
singular and arbitrarily thin fixtures, ambient and spherical central
differences, radial-bias cancellation, explicit whole-ball numerical
families, block counts and the clipping distinction. Finite controls
support the uniform proofs; novelty and significance require the
remaining literature and architecture comparisons.

### Associative-memory goal: affine numerical recall on the sphere

Previous turn classification: progress. Matching numerical-capacity
bounds on full balls and for quadratic-plus-bias spherical readers were
proved and saved. This pass closes the linear-term gap for normalized
queries and audits the projection method against relevant prior art.
The original architecture-advantage objective remains unachieved.

#### The full affine spherical theorem

Fix p>=2 unit keys u_1,...,u_n, independent of their stored values,
and 0<delta<=1. For every y in {-1,+1}^n suppose an arbitrary
affine quadratic f_y(q)=q^T*S_y*q+b_y^T*q+c_y has numerical error
at most epsilon<a relative to a*y_i on every unit q with
||q-u_i||<=delta. There are no coefficient or conditioning bounds.
Then, with r_0=rank(sum_i u_i*u_i^T),

    n <= (r_0+p)*[(a+epsilon)^2+4*epsilon^2/delta^2]
                  /(a-epsilon)^2.

In particular r_0<=p. The feature-dimension bound
n<=p*(p+3)/2 can also be taken in the minimum: on unit queries
the constant is absorbed by S_y+c_y*I, leaving a symmetric quadratic
and a linear term. Clipping to [-1,1] remains covered for interior
test values with a+epsilon<1, as in the preceding proof.

For arbitrary bounded numerical values, accuracy 1/4 and
0<delta<=1/16, choose the interior test values +/-1/2. Combining
the upper theorem with the already verified block construction gives

    n <= min(p*(p+3)/2, 2*p*(9+4/delta^2)),
    maximum capacity = Theta(min(p^2,p/delta^2)).

Thus the matching order now covers full affine quadratics on both
full balls and spherical caps. This is still numerical recall on
continuous neighborhoods; the weaker sign-only problem, grid-only
promises and arbitrary nonlinear/routed readers remain different.

#### Numerical accuracy controls the tangent derivative

For a unit tangent w perpendicular to a key u, write
A=u^T*S*w, B=b^T*w. Put s=delta/2 and c=sqrt(1-s^2).
Use the symmetric unit queries at angles t and 2*t, where sin(t)=s:

    q_small,+/- = c*u +/- s*w,
    q_large,+/- = (2*c^2-1)*u +/- 2*c*s*w.

Both pairs are within the cap; the larger pair has chord distance
exactly 2*s=delta. Their normalized odd differences are

    d_small=[f(q_small,+)-f(q_small,-)]/(2*s)=2*c*A+B,
    d_large=[f(q_large,+)-f(q_large,-)]/(4*c*s)
            =2*(2*c^2-1)*A+B.

Numerical accuracy implies |d_small|<=2*epsilon/delta and
|d_large|<=epsilon/(c*delta). With
lambda=2*(1+c)/(1+2*c), exact extrapolation gives

    2*A+B=lambda*d_small-(lambda-1)*d_large,
    |2*A+B| <= [2*lambda+(lambda-1)/c]*epsilon/delta
             =(2+1/c)*epsilon/delta <= 4*epsilon/delta.

Here c>=sqrt(3)/2>1/2. Since w was arbitrary,
||grad_sphere f(u)||<=4*epsilon/delta. This proof allows an
unbounded linear term and an arbitrarily large radial coefficient
canceling a bias; neither is bounded separately.

#### A value-and-derivative projection counts the extra directions

We use elementary finite-dimensional projection facts, rather than
assuming separate norm bounds for canceling components. For any fixed
linear function class form the sampling map

    J(f)=(f(u_i))_i direct_sum ((1/2)*grad_sphere f(u_i))_i.

Let C be its range, P_C its orthogonal projector, and R the selector
of the first n (value) coordinates. Define d_C=Tr(R*P_C*R^T).
If each sign assignment has a realizing f_y with
||J(f_y)||^2<=n*L^2, then

    n <= d_C*L^2/(a-epsilon)^2.

Indeed let e_y=(y,0) in the sampling space. Accuracy gives
n*(a-epsilon)<=<e_y,J(f_y)>. Membership in C and Cauchy--Schwarz
bound this by sqrt(n)*L*||P_C*e_y||. Independent signs satisfy
E||P_C*e_y||^2=d_C. Taking expectations and squaring proves the
claim, including singular sampling maps and nonlinear value encoding.

For the homogeneous quadratic class C_0, with f(q)=q^T*S*q,

    ||J(f)||^2=sum_i [(u_i^T*S*u_i)^2
                    +||(I-u_i*u_i^T)*S*u_i||^2]
               =sum_i ||S*u_i||^2=Tr(S^2*G_0),
    G_0=sum_i u_i*u_i^T.

The i-th value coordinate obeys

    |u_i^T*S*u_i|^2
        <=(u_i^T*G_0^dagger*u_i)*Tr(S^2*G_0).

This is the same weighted Frobenius inequality used previously. The
squared coordinate norm on a subspace is the corresponding projector
diagonal, so summing gives

    d_C0 <= sum_i u_i^T*G_0^dagger*u_i = r_0.

Now add the p linear functions. The constant jet is already the jet
of q^T*I*q on the sphere. Write C=C_0+E and
W=(I-P_C0)*E. Then C=C_0 orthogonal_direct_sum W, dim(W)<=p,
and P_C=P_C0+P_W. Consequently

    d_C <= d_C0+Tr(P_W) <= r_0+p.

The orthogonal residual W need not lie in E. This distinction matters
when the linear and quadratic parts nearly cancel. No conditioning
or angle between these two spaces was assumed.

Finally the center values have magnitude at most a+epsilon, and the
tangent derivative bound gives

    ||J(f_y)||^2 <= n*[(a+epsilon)^2+4*epsilon^2/delta^2].

Insert this L^2 and the projection trace bound into the general
claim to finish the theorem. All bounds concern the total realized
function; large individual coefficients do not evade them.

#### Primary comparison and an exact prior-art counterexample

The sampling and leverage machinery is established. In
[Knight, Sections 1--3](https://www.utstat.toronto.edu/keith/papers/leverage.pdf),
the augmented least-squares matrix, projection formula (3), and
subset-trace identities were inspected. Our projector on values plus
half-gradients has this standard penalized-least-squares interpretation.
The projection/rank facts alone are not a new architecture or theorem
of publication significance.

[CS4ML v2](https://arxiv.org/html/2306.00945v2), Definition 2.6,
Appendices A.2 and B.1--B.3, and Lemma E.1 with its proof were
inspected. They cover generalized leverage measurements and stacked
value/derivative regression. However, Lemma E.1(iii)'s unqualified
subspace-sum estimate fails in a two-dimensional example. Its exact
statement and proof were also visually checked on PDF page 42.

Take X=R^2 with its Euclidean norm, scalar measurement L(x)=x_2,
F_1=span((1,0)), F_2=span((1,1/2)). Their squared measurement
norms on unit vectors are 0 and 1/5, respectively; their sum is R^2
and has measurement norm squared 1. Thus the claimed bound's two
sides are 1 and 2/5. If aggregate nondegeneracy is desired, use both
coordinate measurements with counting measure; their energy is
exactly ||x||^2. The counterexample persists at the second coordinate.

For these spaces F_2 intersect F_1-perpendicular is zero and cannot
supply the complement used in that proof. The projected residual
(I-P_F1)*F_2 is instead span((0,1)). Our rank argument uses that
residual. We do not use the invalid sum estimate, infer novelty from
this counterexample, or audit its consequences for the other results.

#### Implication for the architecture search and full goal

Adding an unrestricted linear branch or bias cannot restore quadratic
numerical capacity at fixed positive noise and accuracy. The complete
affine reader still has capacity only linear in p. The more general
projection argument also limits a finite-dimensional fixed feature
augmentation when the total value-and-derivative energy is bounded:
adding m columns raises the value trace by at most m. It does not
infer such an energy bound for arbitrary oscillatory or learned
features; those require their own assumptions.

This finishes the obvious affine extension of the numerical reader
restriction. Do not make further constant improvements or ordinary
Sobolev least-squares reformulations the main project contribution.
The existing block-table construction and its indexed full-query
baseline retain exactly the same storage, overwrite and precision
accounting; no architecture advantage was obtained here.

The exact priority of the capacity application remains unresolved.
A publishable contribution still needs a substantive new consequence
or constructive tradeoff on a complete query/update workload beyond
the strongest applicable indexed, routed, kernel and regression
baseline. The next architecture candidate must state such a workload
before proposing an experiment. The full original goal stays active,
with no additional training, messages or submissions.

The existing command now includes the new independent controls:

```bash
python3 -B scripts/run_numerical_capacity_checks.py
```

The extension checks exact value/derivative projectors on seven key
sets, including duplicate, embedded and thin configurations; 296
independent-sign moments; 592 projected Cauchy controls; 162 two-angle
identities with large linear/radial coefficients; 372 whole-cap affine
states covering 1,920 key caps; and three subspace-sum counterexamples.
Previously saved certificate blocks are retained. These exact finite
controls support the uniform proof and its stated limitations.

### Associative-memory goal: conditional sampling meets the weighted-join baseline

The preceding turn was **progress**: the full affine numerical-capacity
case was closed and checked. The next concrete candidate changes the
reader to a normalized distribution on allowed coordinate pairs. It
retains two full query vectors and independently replaceable pair values.
It avoids the two-scalar bottleneck and permits inexpensive approximate
reads. **It fails the strongest-baseline gate:** a weighted indexed join
implements the same distribution and all the same memory/update guarantees.
This candidate is therefore parked as a main contribution. The original
architecture/publication goal remains active and unachieved; no training
or new contribution manuscript follows from this construction.

**Complete workload and model.** Fix m rows, n columns and prefix lengths
0<=l_i<=n, with allowed pairs E={(i,j):1<=j<=l_i}. Monotone lengths give a
Ferrers mask, although the following formula needs only the row prefixes.
Let N=sum_i l_i>0. Store an independent current value
V_ij in [-1,1]^v at every allowed pair. A value update receives (i,j,V_new)
and replaces that payload. The mask and pair identities stay fixed in this
workload; insertion/deletion of rows, changes to the mask and learned-key
updates are separate workloads with unproved additional costs.

A fresh query supplies nonnegative vectors a in R^m and b in R^n. Define

```text
B_k = sum_{j<=k} b_j,  B_0=0,
Z = sum_i a_i B_{l_i} > 0,
P_ij = a_i b_j / Z  on E,
R(a,b,V) = sum_{(i,j) in E} P_ij V_ij.
```

The encoder can be the identity on the m+n query coordinates: none are
compressed to two scalars. Any learned map producing these vectors must
be charged in both methods. The query has compositional product weights
conditioned on a fixed mask; it does not represent every joint weighting
of the N pairs or arbitrary dense learned pairwise scores. Zero weights
are permitted, but Z=0 is an invalid query rather than an arbitrary
normalization convention. All value assignments are allowed independently.

**Uniform sampling identity and full costs.** Build the column cumulative
weights B_k and the row cumulative weights of w_i=a_i B_{l_i}. Choose i
with probability w_i/Z and then j<=l_i with probability b_j/B_{l_i}.
Rows with w_i=0 are never chosen. Their product is exactly P_ij, including
zero-weight and empty-row cases. This proves the conditional sampler
without enumerating the N pairs. Independent row and column choices do
not implement it: for l=(1,2) and a=b=(1,1), the correct distribution is
uniform on three pairs, while the product of its marginals puts mass
1/9 on the forbidden pair (1,2).

For integer weights with B bits per coordinate, both cumulative arrays
are constructed with O(m+n) additions/multiplications on integers of at
most L=2B+ceil(log2(mn))+1 bits. Inverse cumulative lookup takes
O(log(m+1)+log(n+1)) comparisons per sample. The second ticket is drawn
from a range that depends on the selected row. An exact uniform ticket
uses rejection on random bits with fewer than two trials in expectation
and O(L) expected random bits per ticket. There is no rejection of
completed join pairs, and no unit-cost arbitrary-real random oracle.
Bit-operation costs of arithmetic/comparisons on L-bit operands remain
charged; these are word-operation bounds only with adequate word width.

Use row offsets to store all payloads in one packed table. If each value
coordinate has S dyadic fractional bits, persistent state costs
N*v*(S+O(1)) payload bits plus O(m*(log(n+1)+log(N+1))) bits for prefix
lengths/offsets. Query input and cumulative workspace cost
O((m+n)*(B+L)) bits. Initialization writes all N*v coordinates. A canonical
point overwrite costs O(v) value-word writes and has no accumulation of
old update errors. K samples cost K*v value-word probes, O(K*v) additions,
the cumulative lookup costs above, and the charged O(m+n) query setup.
The accumulator requires O(v*(S+log(K+1))) bits; exact averaging retains
the integer numerator and denominator K, or output rounding adds its
own fixed error. Workspace and preprocessing are not amortized away
over unspecified future queries. An exact O(N*v) read is an available
fallback when K is too large, for both implementations.

At a valid point code a=e_i,b=e_j, the output is exactly V_ij under every
sample. Parsing those vectors still costs O(m+n) if the input format is
dense. Explicit pair-address input can bypass that scan for both methods.
All N independently overwritten values are retained, because all N are
stored. There is no compression below the information required by them.

**Error, adaptivity and retention.** Conditional on the complete history,
choose the current query and payload state, then draw K fresh independent
pairs. The sample average is unbiased for the current R. For epsilon>0,

```text
Pr(||R_hat-R||_infinity > epsilon | history)
    <= 2*v*exp(-K*epsilon^2/2).
```

This is the classical bounded-variable estimate, not a new concentration
result. One direct proof writes h(lambda)=log E exp(lambda*(X-E X))
for a coordinate X in [-1,1]. The second derivative is the variance under
an exponential tilt and is at most 1; h(0)=h'(0)=0. Thus
h(lambda)<=lambda^2/2. Independence, exponential Markov and optimizing
lambda give the two scalar tails, followed by a union over coordinates.
See [Hoeffding (1963)](https://doi.org/10.1080/01621459.1963.10500830)
for the established bounded-sum method.

For a specified finite horizon of Q adaptive reads, the conditional bound
and a union bound give total failure probability at most
2*v*Q*exp(-K*epsilon^2/2). It suffices to take
K>=2*epsilon^-2*log(2*v*Q/eta) for failure probability eta. The current
payload/query can depend on previous samples; it cannot be selected after
seeing the fresh samples of its own read. No fixed-K, all-reads-success
claim for an infinite horizon is made. No analogous adaptivity guarantee
for an old public sketch is inferred.

For example, on a full 4-by-4 uniform query, retain four diagonal samples
and later set their values to -1 and the remaining twelve values to +1.
The retained-sample answer is -1, the true mean is 1/2, and the error is
3/2. Fresh sampling of this current state remains unbiased. This refutes
only the proposed reuse-of-these-samples guarantee; it does not rule out
other adaptive sketches with different hypotheses/resources. A canonical
payload table gives indefinitely stable stored values; its statistical
read-error probability still needs the separate horizon accounting.

**Noise promises are also different.** Normalize a,b to probability
vectors when discussing geometric stability. Around a valid point code,
if their respective total-variation distances from e_i,e_j are at most
alpha,beta, then P_ij>=(1-alpha)*(1-beta), because the unconditioned target
mass has that lower bound and Z<=1. Therefore every independent payload
assignment satisfies

```text
||R-V_ij||_infinity <= 2*(alpha+beta-alpha*beta).
```

This bound is attained on a full mask by putting the target value at one
endpoint and all other values at the opposite endpoint. It is a simplex
noise guarantee, not the previous continuous ambient-ball or sphere-cap
promise. In particular l1 and l2 errors must not be interchanged without
dimension factors. The normalized bilinear reader has a denominator and
a restricted nonnegative domain; it does not contradict the quadratic
capacity theorem.

For general simplex queries write mu_ij=a_i*b_j*1_E and define mu' and
Z' similarly. If D=||a-a'||_1+||b-b'||_1, then ||mu-mu'||_1<=D and

```text
||mu/Z-mu'/Z'||_1 <= min(2, 2*D/max(Z,Z')).
```

Indeed, adding/subtracting mu'/Z bounds this by
(||mu-mu'||_1+|Z-Z'|)/Z<=2D/Z; exchanging the two queries gives
the other denominator. If Z,Z'>=gamma>0 this is at most 2D/gamma,
and it also bounds the output error for every bounded payload. Without
a mass promise there is no uniform Lipschitz bound. For l=(1,2), compare

```text
query A: a=(1,0),     b=(t,1-t),
query B: a=(1-t,t),   b=(0,1),        0<t<1.
```

Their combined l1 distance is 4t, and both mask masses are t. Their
conditional distributions are concentrated on (1,1) and (2,2),
respectively. Set V_11=-1,V_22=+1 to obtain output gap 2 for arbitrarily
small t. Finite-precision or encoding noise therefore cannot be declared
harmless merely because the cumulative sampler is exact on its supplied
integer weights. The stability promise is shared by both implementations.

**Primary baseline comparison.** Shekelyan, Cormode, Triantafillou,
Shanghooshabad and Ma,
[Weighted Random Sampling over Joins, v1](https://arxiv.org/html/2201.02670v1),
Definitions 2.1--2.2, Section 3.2/Algorithm 1 and Section 3.3, explicitly
use products of row weights, downstream completion sums and conditional
continuation sampling. Section 3.2 also treats inequalities by cumulative
weights. These definitions and the displayed algorithm were inspected,
rather than inferring equivalence from the title. Our prefix sampler is
the already ordered two-table specialization with predicate j<=l_i.
The paper's more general streaming implementation need not have the
same constant factors as this indexed specialization. Its results are
not being cited as a theorem about dynamic payload overwrites: those
costs follow directly here because changing V does not change the query
distribution. We have not audited every algorithm or experiment in that
paper.

The independent reference implementation can use the identical integer
cumulative arrays, packed value table and random tickets. Coupled on
those tickets, it returns the same pairs and answers after every payload
overwrite, for every legal query/history. Thus its distribution, failure
bound, retention, metadata, precision, setup and active value costs all
match. Comparing only against materializing the dense masked matrix,
against independent-coordinate sampling, or against an unindexed
product-proposal sampler would miss this baseline. For the last example,
l=(1,2),a=(M,1),b=(1,M) has proposal acceptance
(2M+1)/(M+1)^2 and at least M/3 expected trials. Conditional sampling
avoids those rejections, but that advantage is already available in
weighted joins. No optimality claim for binary-search CDF lookup is made.

Two additional primary comparisons constrain the next candidate:

- [Kuang et al., Sections 2.1--3.3, v1](https://arxiv.org/html/2509.07963v1)
  explicitly replace low-rank attention scoring by full/high-rank BTT and
  MLR bilinear forms, and combine them in MLBTC. Those definitions and
  their scoring equations were inspected. Full rank plus inexpensive
  structured evaluation is therefore already an architecture direction,
  not a new claim by this project. We have not independently reproduced
  their learning results or established uniform retention from them.
- [Song--Xie--Zhou, v1](https://arxiv.org/html/2510.03678v1), Definition 1.1,
  Theorem 4.4, Algorithm 2 and Theorem 7.6, address streamed attention
  sampling. These statements and their surrounding arguments were
  inspected, not independently validated in full. The tensor construction
  studies squared amplitudes of (A1 tensor A2)x, with fixed x and updated
  factors; Theorem 7.6 is conditional on auxiliary estimators. A score
  sampler alone is not a theorem retaining N independent payloads or
  handling adaptively chosen fresh queries and changing factors for free.
  Its formal guarantee must be separated from a whole memory architecture
  before adopting or claiming an improvement to it.

**Decision and next proof obligation.** Park the masked product-mixture
candidate as a contribution: its complete guarantee is matched by the
published-method specialization above. Do not run a learned comparison
of two implementations of that same mathematical algorithm. Preserve it
as a baseline for any genuinely richer construction. A next proposal must
identify a new complete workload, such as full-dimensional scores with
changing bindings, and explicitly beat the applicable structured-scoring,
dynamic-index or tensor-sampling method under that workload. It must
include access to arbitrary current payloads and charge factor/key edits,
normalizers, precision, query setup and adaptive failure probability.
Another known sampler, rank calculation or parameter-count advantage
does not fulfill the unchanged goal. There is no established new
architecture advantage or publication-ready contribution from this turn.

The existing checker was extended in place:

```bash
python3 -B scripts/run_correlated_routing_checks.py
```

It independently enumerates 12,580 positive-mass queries and rejects 4,934
zero-mass queries; checks every integer row ticket and conditional child
ticket in those cases (73,881 and 49,254 respectively); verifies 640
canonical vector overwrites and point-code reads; covers 26 adaptive
histories with fresh sampling; checks 52 simplex-noise cases, four
vanishing-mask-mass witnesses through 128-bit exponents, four rejection
controls, the fixed-sample counterexample and 72 exact Bernoulli-tail
controls. The slow direct completion enumeration is an independent proof
oracle, not a runtime benchmark for the production indexed baseline.
These finite certificates support the uniform identities and expose
overstated guarantees; they do not establish novelty. All prior
certificate blocks and the smooth-number kernels remain preserved.
