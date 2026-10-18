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
  Section 5, supplies a finite-field specialization. Its parameters and the
  later published theorem must be compared without identifying fixed field
  size with the maximum over all fields. The earlier pass found that the
  v1 expression recovers the old exponent when p>n.
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
