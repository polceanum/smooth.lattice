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
Proposition 17 is the first candidate under the active research goal: an
integer representation of an entire ordered prefix with a cutoff score
linear in its length for fixed dimension. Its novelty and usefulness are
unresolved; it is not yet a paper result or a faster unranking algorithm.
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

The active quantitative goal is recorded in
`results/research/quantitative_goal.json`. The final section below records a
multiscale proof candidate for cumulative corner cuts in every fixed dimension,
with stronger exact-size lower bounds in dimensions four and higher. Candidate
22 gives a uniform finite sieve bound and the general argument. Its internal
proof audit now includes closed constants and a step-by-step ledger; priority
and publication significance remain unresolved. The structural note and its
elementary consequences do not satisfy this goal.

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
