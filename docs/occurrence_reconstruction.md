# Returning to the occurrence and prime-exponent observations

## Scope and current assessment

This track returns to the observations that started the project. Order the
products of a fixed set of distinct primes and record their exponent vectors.
Deleting entries with a zero in one coordinate, then subtracting one in that
coordinate, reproduces the vector sequence. Counting earlier occurrences of
one exponent stream twice also returns that stream.

These facts, their established signature-sequence connections and the earlier
results are preserved in [the original research note](interleaving_research.md).
The ML investigation is paused and preserved separately in
`associative_memory_timing.md`. None of the results below uses its Hamming
kernel or memory architecture.

The question here is precise: **how much lower-dimensional ordering information
determines a finite higher-dimensional interleaving?** We prove a logarithmic
upper bound and a logarithmic worst-case lower bound. The uniform proofs are
given below; finite checks support their delicate steps. Priority and
publication significance remain unresolved. This is a candidate structural
addition, not a completed publication claim or a faster unranking algorithm.

## The occurrence partner already identifies the missing object

For a word, write

    T(a)_n = number of m<n with a_m=a_n.

Let s_n be the n-th P-smooth number, E_i(n) its p_i exponent and C_Q(x) the
number of Q-smooth numbers at most x, including 1. Proposition 8 of the
original note proves

    T(E_i)(n) = C_(P without p_i)(s_n / p_i^E_i(n)) - 1,
    T(T(E_i)) = E_i.

For two primes, the first expression is just the other exponent. For three
or more, it is the rank of an entire remaining exponent vector. This explains
why the two-stream picture survives as a duality but stops exchanging two
individual coordinates. Ordinary occurrence duality is established background.
The remaining question concerns reconstructing the ordering hidden in those
face ranks. The results below address that finite information problem.

## Definitions

Let w be a positive real d-vector with rationally independent coordinates.
Write W_n(w) for the ordered first n vectors of the nonnegative integer
lattice under x -> w.x, starting with 0 at rank 1. Prime logarithms satisfy
the independence assumption by unique factorization.

For a nonempty coordinate subset S, the *face word* W_n(w_S) orders the
lattice points supported on S. The input retains the coordinate labels.
It consists of n exponent vectors, not numerical smooth values and not the
first n entries of a projection of the full word. Those are different data.
Coordinate order means label order; increasing numerical weights are not
assumed. The higher-dimensional lower examples can also be relabelled by a
common permutation so both prime tuples increase numerically.

For an integer vector z put

    V(z) = product_i (1+abs(z_i)).

In particular, 2^|support(z)| <= V(z). Mixed signs are necessary for the
hyperplane w.z=0 to meet the strictly positive weight cone.

## A dimension-independent difference bound

**Lemma 1.** If a and b are distinct vectors among the first n points, with
b the later one at rank j, then

    V(a-b) <= 2(j-1) <= 2(n-1).

Consequently their difference involves at most floor(log_2(2n-2)) coordinates.
This concerns coordinates in which they differ, not the size of the prime set.

**Proof.** Cancel their common coordinatewise minimum c. Put u=a-c, v=b-c
and m=u+v=abs(a-b); the supports of u and v are disjoint. In the integer box
Q={x:0<=x_i<=m_i}, pair x with m-x. The two scores sum to

    w.u+w.v < 2 w.v.

Every two-element pair therefore has a member strictly below w.v. A possible
fixed point has score (w.u+w.v)/2 and is also strictly below it. This gives
at least ceil(|Q|/2) box points below w.v. The pair {u,v} supplies both
members at or below w.v, adding one more when counting at or below. Thus
the rank of v is at least ceil(|Q|/2)+1. It is at most the rank j of b,
because cancellation decreases its score or leaves it unchanged. Hence
|Q|=V(a-b)<=2(j-1). The support bound follows by taking logarithms. QED.

The same involution explains the comparison-wall threshold. It is useful to
keep the strict-prefix and tied-wall conventions distinct.

**Lemma 2.** Let z be primitive, integral and of mixed sign. Write u=z_+,
v=z_-. At a positive weight on w.z=0 with no other independent rational
relation in their box, the first rank occupied by the tied pair is at least

    V(z)/2.

**Proof.** The box involution now pairs scores about the common height h.
Its only height-h points are u and v. Indeed a further such point differs
from u by an integer multiple of the primitive z; the box and nonnegativity
permit only u and v. Exactly (V(z)-2)/2 box points lie below h. At least
that many full-lattice points lie below h, so the first tied rank is at
least V(z)/2. The product is even: primitivity gives an odd coordinate of z.
QED.

The constant 2 in this universal wall bound cannot be decreased. In dimension
two, for coprime a,b>=1 and weights (b,a), the tied points (a,0),(0,b) first
occur at rank (a+1)(b+1)/2. This planar threshold is already recorded in
[Adams-Watters's OEIS A163823](https://oeis.org/A163823); it is not new here.

## A finite comparison certificate and face reconstruction

**Theorem 3.** For n>=2, the signs of w.z for all primitive mixed-sign integer
vectors with V(z)<=2n determine W_n(w).

**Proof.** Two weight vectors having those same signs lie in the same open
convex cone cut out by the corresponding strict inequalities and positivity.
There are finitely many such z, since every |z_i| is at most 2n-1.

Join the two weights within this cone by a generic polygonal path. Only a
finite arrangement needs consideration: every point that can occur among
the first n, or tie for a rank at most n, has each coordinate at most n.
Having a coordinate greater than n would already give more than n strictly
smaller pure powers. Perturb the path within the open cone to avoid
intersections of distinct comparison hyperplanes of this finite box.

At a crossing, let z be its primitive normal. Any pair changing order there
has difference an integer multiple of z. The lowest such tied pair is u=z_+
and v=z_-: other pairs are positive translates of positive integer multiples
of these, and have strictly higher common height. If the crossing changes
W_n, this lowest pair has first tied rank at most n. The generic crossing
has no second independent relation in its box, so Lemma 2 gives V(z)<=2n.
Such a hyperplane cannot be crossed inside the fixed-sign cone. Thus the
word never changes along the path. QED.

**Theorem 4.** Let n>=2 and r=min(d,1+floor(log_2 n)). The collection

    { W_n(w_S) : nonempty S subset {1,...,d}, |S|<=r }

determines W_n(w). Face prefixes of length n suffice; extra terms are not
assumed. The one-term word is always (0) and requires no face data.

**Proof.** The set of positive weights realizing a specified face word is
an open convex cone. One finite description imposes its consecutive strict
comparisons and compares its last point with every coordinatewise minimal
excluded point in that face. The minimal border is finite and every excluded
point dominates a border point. These inequalities therefore specify exactly
that prefix. Intersect the cones for the supplied face words.

If two realizing weights gave different full words, a generic path inside
their common open cone would cross a comparison wall that changes W_n.
As in Theorem 3, its primitive normal z has V(z)<=2n. Thus its support S has
size at most floor(log_2(2n))=1+floor(log_2 n). The lowest tied pair is
supported on S. Its first rank in that face is no greater than its full
rank, and hence at most n. Crossing this wall changes that face's n-term
word, contradicting the common cone. QED.

This is a uniqueness theorem under an existence promise. Arbitrarily supplied
face words need not be jointly realizable. It does not prove that pairwise
compatibility is sufficient. Nor is it an efficient algorithm merely to
request all these faces: their number can be sum_(j<=r) binomial(d,j).
An LP using all the strict prefix constraints can find a positive rational
realization when one exists; a generic perturbation and enumeration then
recover the full word. That elementary procedure gives no new cost advantage
over the established rank/unrank methods.

## Logarithmically many coordinates are necessary in the worst case

The earlier three-prime obstruction says that arbitrarily long pair prefixes
can agree while the full words differ at rank 5. The following extends the
same separating-wall mechanism to arbitrarily many coordinates.

**Theorem 5.** For every k>=1 and M>=1, there are two ordered sets of 2k
distinct primes, sharing their first 2k-1 prime bases, such that all
corresponding proper coordinate faces have identical M-term exponent words,
but the full words differ at a rank at most binomial(3k,k).

**Proof.** Put d=2k and

    u=(1,...,1,0,...,0), v=(0,...,0,1,...,1), z=u-v,

with k ones on each side. Choose a positive weight w on w.z=0 within eta
of the all-one vector, where eta<1/(2k+1). Choose it generically on that
hyperplane: its only rational relations are multiples of z. This is possible
by avoiding the countably many other proper hyperplanes. In particular
every proper coordinate subvector is rationally independent.

The common height h=w.u=w.v is at most k(1+eta). Every lattice point of
total degree at least k+1 has score at least (k+1)(1-eta)>h. Thus at most
binomial(d+k,k)=binomial(3k,k) points lie at or below h. The only tied points
at h are u,v: a point with that score differs from u by an integer multiple
of z, and nonnegativity allows only the two endpoints.

Perturb the last weight both downwards and upwards, leaving the others
fixed. A sufficiently small perturbation preserves each proper face's first
M strict comparisons and border inequalities, a finite set of conditions.
It also preserves every comparison around the full prefix below h except
the order of u,v. Hence the two full words first differ at their common
first tied rank, at most the displayed binomial bound. Their proper face
words agree through M.

Finally realize these finite strict conditions by prime logarithms. For a
large common scale L, choose shared primes in
[exp(L w_i),2 exp(L w_i)] for i<d, and separate last primes approximating
the two perturbed last weights in the same fashion. Bertrand's postulate
gives them, and the normalized log errors are at most log(2)/L. All strict
conditions persist for sufficiently large L. The finitely many distinct
target weights on each side give disjoint prime intervals when L is large;
all primes can also exceed any prescribed bound. This is the established
finite prime-realization argument, applied simultaneously to the two sides,
not a new prime-distribution theorem. QED.

Let R(n) be the least face dimension sufficient uniformly over all ambient
dimensions and positive independent weights, using n terms of each supplied
face. Theorem 4 and Theorem 5 imply, for sufficiently large n,

    2 floor(log_(27/4) n) <= R(n) <= 1+floor(log_2 n).

Indeed binomial(3k,k)<=(27/4)^k follows by bounding its term in
(1+1/2)^(3k). Choose k=floor(log_(27/4) n), and use Theorem 5 with M=n.
Even all proper faces of dimension 2k then fail to determine the full word.
Thus **R(n)=Theta(log n)**. The leading constant and the exact thresholds
are not determined. The lower bound holds for actual prime bases, not only
abstract weights. Its primes depend on M and can be very large. There is
no claim about a single fixed small prime set, or about entire infinite faces.

Here is an independently checked four-prime witness:

    P_minus = (168323,106937,130069,138373),
    P_plus  = (168323,106937,130069,138389).

Every one of their 14 corresponding nonempty proper coordinate faces has
the same first 1,000 exponent vectors. The full words first differ at rank 10:
the minus word has (0,0,1,1), whereas the plus word has (1,1,0,0). The
arithmetic comparison is

    130069*138373 = 17998037737
      < 168323*106937 = 17999956651
      < 130069*138389 = 18000118841.

The checker verifies primality by trial division and every face word with
both the exact heap and an independent integer-product pointer generator.
The numbers themselves are available to an arithmetic algorithm; the
obstruction concerns exponent-pattern-only reconstruction.

## Prior-art boundary and impact gate

The proof package uses established tools. The following primary comparisons
are essential; absence of an exact phrase in them is not novelty evidence.

| Primary source and inspected location | Established ingredient and remaining comparison |
| --- | --- |
| [Wagner's 2003 thesis](https://ti.inf.ethz.ch/ew/research/phdtheses/wagner03.pdf), Chapter 3, Observations 3.3, 3.13--3.14 and Theorem 3.15 | Divisor boxes, a finite product-bounded ground set, volume estimates and product-bounded difference-matrix encodings are prior art. Our dimension-independent two-point bound and finite face reconstruction need comparison with the full encoding argument and its subsequent uses. |
| [Maclagan, Boolean Term Orders and the Root System B_n](https://arxiv.org/abs/math/9809134), Section 3 and Lemma 4.1 | Complement symmetry, primitive/flippable pairs and changing a coherent order across one comparison hyperplane are established. The lower construction's Boolean core is part of this framework. Its ordinary-monomial cutoff and quantitative face information claim need an exact comparison. |
| [Christian--Conder--Slinko, Flippable Pairs and Subset Comparisons](https://www.math.auckland.ac.nz/~slinko/Research/isipta17.pdf), Section 3 and Theorem 3 | The central complementary pair and determining coherent subset orders by comparisons are established. Their Boolean domain differs from ordinary exponent prefixes, which also include repeated powers. Whether their elicitation results imply this face-dimension statement remains to be checked. |
| [Vietri, 2002](https://www3.risc.jku.at/research/theorema/Groebner-Bases-Bibliography/gbbib_files/publication_455.pdf), Property 2.6 and Proposition 6.2 | Fixed-degree restrictions and finite prime realizations are established. These passages were inspected directly. No new prime-realization principle is claimed. |
| [Vasiliev--Pavlov, 2009](https://doi.org/10.1134/S0361768809020030), Sections 5--6, Propositions 4--7 | The previous full-paper audit is preserved in the original research note. Feasible weight cones, reduced comparisons and strict-prefix enumeration are prior art. The new dimension bound must be compared with those constraint reductions, rather than inferred new from the paper's absence of an asymptotic counting theorem. |
| [Adams-Watters, A163823](https://oeis.org/A163823), formula and comments | The planar comparison threshold and enumeration predate this work. They are attribution and a control, not a new theorem. |

The candidate payoff is a uniform information hierarchy explaining when
lower-dimensional fractal/signature data suffice. It concerns growing
dimension, which the earlier fixed-dimension counting manuscript did not
settle. It supplies a way to constrain reconstruction proposals and identifies
which coordinate interactions can first affect a requested rank.

The current assessment is **continue the mathematical comparison**, with no
contribution manuscript yet. A paper would need a defensible new statement
and a useful consequence beyond this elementary reduction. The next questions
are whether the logarithmic hierarchy is already implicit in published
finite-ordering theory, whether its constants or data size admit a sharp
theorem, and whether the actual occurrence/zero-mask encoding yields a better
reconstruction algorithm. A negative answer must be recorded rather than
turning these supporting lemmas into a publication claim.

## Reproducible exact checks

Run from the repository root:

```bash
python3 -B scripts/run_face_reconstruction_checks.py
```

The report is `results/research/face_reconstruction_checks.json`. It checks
integer-product prefixes with an independent pointer generator; the occurrence
partner and return law; all exponent-pair difference bounds; small primitive
wall boxes by direct scalar scores; the known planar thresholds; reconstruction
collisions over explicit strict integer-weight grids; and all proper faces of
explicit rational-weight lower-bound fixtures in dimensions 2,4,6,8,10.
Full fixture order is checked by an independent weighted pointer merge, and
the first tied rank is counted over a rigorously sufficient degree simplex.
The separate four-prime witness checks all proper faces through length 1,000.

The rational fixtures are not presented as prime fixtures. The uniform prime
statement is proved separately above. Neither finite agreement nor a report
marked passed establishes novelty, an optimal constant or a running-time gain.

## One exponent stream determines the normalized ordering

The first iteration was progress: it supplied the finite face hierarchy,
uniform proofs and exact prime witnesses. It did not settle priority or
provide an access improvement. This second iteration returns more directly
to the individual streams. Surprisingly, the rank of the remaining face is
enough to recover its numerical spacings from the *infinite* stream. The
proof below also supplies a finite error certificate with explicit extra
calibration data. The distinction between those two statements matters.

Normalize the observed coordinate's weight to 1, and write the remaining
weights as theta=(theta_1,...,theta_m), m=d-1>=1. Assume 1,theta_1,...,theta_m
are rationally independent. Let

    0=t_0<t_1<t_2<... be the ordered set {theta.x : x in Z_>=0^m},
    M(T)=#{x : theta.x<=T},
    F(T)=#{(a,x) : a+theta.x<=T}.

From the observed exponent stream E, define A(a,r) to be the 1-based
position of the (r+1)-st occurrence of symbol a. The occurrence-count
identity gives, exactly,

    A(a,r)=F(a+t_r),
    A(a+1,0)-A(a,0)=M(a+1),
    Z(A(a,0))=M(a),                                      (SS1)

where Z counts zeroes through the given position. The second equality also
follows by dividing the counted tuples into those with positive clock
coordinate and those with clock coordinate zero. The placement array and
duality are established sequence machinery; (SS1) specifies their use here.

**Theorem 6 (single-stream rigidity).** In every finite dimension d>=2,
the one infinite stream E determines the multiset of normalized remaining
weights, and hence d and the full ordered exponent sequence up to permutation
of the unobserved coordinate labels. Specifically, for every fixed r>=0,

    t_r = lim_(a->infinity)
          [A(a,r)-A(a,0)]/[A(a+1,0)-A(a,0)].             (SS2)

The infinite binary zero mask of E suffices as well. A common positive
scaling of all weights is necessarily invisible.

Here is a finite statement that proves the limit without assuming a
high-order lattice-counting expansion.

**Theorem 7 (finite height certificate).** Suppose alpha is one of the
remaining generator weights. Let p,b be coprime integers, b>=1, with

    abs(alpha-p/b)<=1/b^2,     0<alpha<=H.

For an integer a>=1, an integer K>=max(1,ceil(t_r)) and

    R_a(r)=[A(a,r)-A(a,0)]/[A(a+1,0)-A(a,0)],

one has the explicit bound

    abs(R_a(r)-t_r)
      <= 8/b + 2^m H b/a
         + (2K+1)[M(a+K)-M(a)]/M(a+1).                 (SS3)

All the M-values in the last term are observable in the zero mask through
(SS1). K can be taken to be one plus the largest observed exponent through
A(0,r), with a minimum of 1. The other inputs are a bound on dimension and
a certified rational approximation to one actual remaining generator.
They are additional calibration information, not quantities inferred for
free from a finite word. For a supplied prime system, an integer H and the
rational approximation can be certified from the corresponding prime pair.

**Proof of Theorem 7.** Use right-closed circular intervals to define
the discrepancy D(a) of the phases {theta.x} over theta.x<=a, normalized
by M(a). If abs(alpha-p/b)<=1/b^2, a block of b points in a rotation by
alpha is within circular distance 1/b of a translated uniform b-grid.
For any circular interval, perturbation can move at most six grid points
across its two endpoints; the grid count differs from its length times b
by at most one. The bound eight is therefore safe. Blocks and a remainder
give discrepancy in a row of L points at most 8L/b+b.

Slice the simplex along the alpha coordinate. Let R(a) count the other
m-1 coordinates at height at most a, with R(a)=1 in dimension zero. Summing
the row bounds gives

    D(a)<=8/b+b R(a)/M(a).

Coordinatewise halving maps the R(a) points into the R(a/2) points, with
at most 2^(m-1) preimages. Also each such latter point can be combined with
any alpha-coordinate from 0 through floor(a/(2 alpha)). Thus

    R(a)<=2^(m-1) R(a/2),
    M(a)>=(floor(a/(2 alpha))+1) R(a/2),
    R(a)/M(a)<=2^m alpha/a<=2^m H/a.                    (SS4)

Put tau=t_r, and write tau=L+phi with 0<=phi<1. For an existing face
point y=theta.x<=a, the number of its clock translates entering (a,a+tau]
is exactly

    floor(tau-y)-floor(-y)
      = L + 1_((0,phi])({y}).

The mean of this function on the circle is tau. Its sum over those M(a)
points differs from tau M(a) by at most M(a)D(a). New face points in
(a,a+tau] contribute at most K+1 clock translates apiece; their number is
at most M(a+K)-M(a). Divide by M(a+1), and replace tau M(a)/M(a+1) by tau.
The additional error is at most K[M(a+1)-M(a)]/M(a+1). Since K>=1, both
boundary terms are bounded by the last term in (SS3). Substitution of
(SS4) proves the stated finite bound. QED.

**Proof of Theorem 6.** Put S=sum theta_i. The unit cubes based at the
points counted by M(T) cover the positive simplex of height T and lie in
the simplex of height T+S. Consequently

    T^m/(m! product theta_i) <= M(T)
      <= (T+S)^m/(m! product theta_i).                  (SS5)

For fixed K, [M(a+K)-M(a)]/M(a+1) therefore tends to zero. The irrational
alpha has coprime continued-fraction approximants with unbounded b satisfying
Theorem 7's hypothesis. First let a tend to infinity with such b fixed;
(SS3) bounds the limit superior of the error by 8/b. Then let b grow. This
proves (SS2), for every r, without an assumed uniform Diophantine exponent.

The set of all recovered t_r is the full additive monoid generated by theta.
Its nonzero elements that cannot be written as a sum of two nonzero elements
are exactly the m generator weights. Indeed a vector of total degree at least
two splits, whereas rational independence prevents a generator from being
a sum of other nonnegative integer combinations. This recovers the multiset
theta and its size, so also the entire lattice ordering up to the unavoidable
permutation. Identifying these atoms refers to the exact infinite set; this
argument is not a finite equality-testing algorithm on estimated real numbers.

Finally the deletion self-similarity reconstructs E from its zero mask:

    E(n)=0 at a zero position,
    E(n)=1+E(n-Z(n)) otherwise.

The first entry is zero, so the referenced index is smaller than n. Thus
every finite prefix of E can be reconstructed successively from the same
length of zero-mask data. This established recursion, followed by (SS2),
proves the last claim. QED.

The same rank-inversion proof applies more generally to a locally finite
height set H with no integer differences between distinct heights, provided
its phases become uniformly distributed and
M(a+K)/M(a)->1 for every fixed K. Before inserting the rotation estimate,
the error is at most its interval discrepancy plus
(2K+1)[M(a+K)-M(a)]/M(a+1). Finite additive generators supply these
hypotheses and the explicit calibration bound in (SS3). This identifies
the reusable mechanism; it is not an additional priority claim.

It is essential to compare with ordinal-uniqueness theorems with their
actual hypotheses. Placement-array duality alone does not imply numerical
recovery. For example, put t_r=r sqrt(2) and

    Phi(x)=floor(x)+{x}^2,     t'_r=Phi(t_r).

Phi is strictly increasing and Phi(x+a)=Phi(x)+a for every integer a.
Thus every comparison between a+t_r and b+t_s is preserved by replacing
t with t', giving the same infinite placement order and observed symbol
stream. But the heights change. The new set is not additively closed:
t'_1=4-2 sqrt(2), t'_2=14-8 sqrt(2), and t'_1<2t'_1<t'_2.
Its transformed phases also fail the uniformity hypothesis. This rules
out deducing Theorem 6 from arbitrary placement-array duality alone. It
does not rule out a published theorem with stronger applicable hypotheses.

## Prime bases: an arithmetic rigidity consequence

**Corollary 8.** For systems containing at least three distinct prime bases,
one infinite exponent stream, or its binary zero mask, determines the
observed prime and the set of all other primes. Their labels may be permuted.
The number of bases need not be supplied.

**Proof.** Theorem 6 leaves only a common positive scale lambda between
the two lists of logarithms. If lambda is irrational, the two numbers
1,lambda are rationally independent. Three logarithms of distinct primes
are also rationally independent. All six exponentials of their pairwise
products would nevertheless be primes from the two systems, hence algebraic.
This contradicts the classical six exponentials theorem. Thus lambda is
rational, say A/B with coprime positive integers. For the observed primes,
q^B=p^A forces q=p and A=B by unique factorization. The scale is one, and
the remaining prime sets coincide. QED.

The transcendence input is established mathematics, not a new result here.
A primary real-variable proof is [Pila, 1993, Theorem 3.1, printed pp.
119--120](https://doi.org/10.1017/S1446788700037022). The logarithmic matrix
form and stronger variations are discussed in [Waldschmidt, 2005](https://hrj.episciences.org/86).
The candidate contribution is the recovery of normalized weights from one
stream; the prime corollary is an application of that recovery and the cited
theorem. It supplies no new prime-gap or transcendence estimate.

For two bases, the proof establishes the normalized ratio. The four
exponentials conjecture would also remove the scale ambiguity there, but
the six-exponentials argument alone does not do so. No equivalence with that
conjecture or unconditional two-prime absolute-uniqueness claim is made.
For one base every prime gives the stream 0,1,2,... . For composite bases,
common powers explicitly preserve the entire stream: (2,3,5) and (4,9,25)
are a control showing why the prime hypothesis matters.

## Finite words have no uniform recovery length

**Theorem 9.** For every finite M>=1 there are two increasing triples of
distinct primes, with the same first two primes, whose observed first-prime
exponent streams agree through M but whose full exponent words first differ
at rank 6. The observed coordinate is already active at rank 2.

**Proof.** Start from weights (1,sqrt(2),2 sqrt(2)). Every exact tie between
distinct exponent vectors has the same first coordinate, because 1 and
sqrt(2) are rationally independent. The only relation direction is (0,2,-1).
Consider the finite box {0,...,max(M,7)}^3. A sufficiently small positive or
negative perturbation of the third weight preserves every strict comparison
between its distinct height classes. Inside a tied class, every first
coordinate is equal, so its projected word is unchanged by any internal
reordering. Every possible first-M point lies in this box, independently
of the positive weights. Thus the two first-coordinate words agree through M.

The first five full vectors on both sides are

    (0,0,0), (1,0,0), (0,1,0), (2,0,0), (1,1,0).

The next two are (0,2,0) and (0,0,1) in opposite orders. The perturbation
can be chosen small enough to preserve this statement. Approximate the
first two weights by shared prime logarithms and the perturbed third
weights by separate prime logarithms using the same finite Bertrand-based
argument as Theorem 5. Preserve every strict inter-class inequality in the
box and the comparisons specifying the first seven vectors. No within-class
comparison needs preservation: its first coordinate is constant. This gives
the stated prime triples, with increasing numerical bases. QED.

One exact witness is

    P_minus=(101,683,466483), P_plus=(101,683,466517),
    466483 < 683^2=466489 < 466517.

Their observed first-prime exponent streams agree through 20,000 entries.
The full word changes at rank 6, from (0,0,1) to (0,2,0). Trial division,
the heap and an independent exact product-pointer generator check this.
The report stores a word hash and the first eight full vectors. It does not
claim that the infinite streams coincide or locate their eventual first
mismatch. The uniform theorem uses a new pair of triples for each M.

There is no conflict with Corollary 8: the infinite data are unique, while
no bound depending only on the requested output rank gives a sufficient
input-prefix length uniformly over unbounded prime bases. Likewise (SS3)
charges its calibration, dimension and boundary terms; it is not a uniform
finite recovery guarantee without numerical information or separation.

## Second-iteration comparison and publication decision

The inspection of [Kimberling--Brown, 2004](https://cs.uwaterloo.ca/journals/JIS/VOL7/Kimberling/kimber67.pdf)
was focused on the rank-array definition, equations (13)--(16), Theorem 4.2
and Corollary 4.6. These are two-generator arrays; their exact floor increments
already recover the one irrational parameter. That case and the placement
array are prior art. [Kimberling, 2022, Section 4 and Example 7](https://cs.uwaterloo.ca/journals/JIS/VOL25/Kimberling/kimber16.pdf)
establishes regular position-array duality and its signature example. Those
inspected statements do not supply the present multidimensional recovery
argument. This is a scoped comparison, not evidence that the general result
has never been proved elsewhere. A missing keyword search is not a novelty
certificate. Classical rotation discrepancy and simplex volume estimates
are explicitly proved here to make their role and hypotheses auditable.

The new candidate has a concrete arithmetic consequence and a finite error
statement, unlike a conclusion based only on visual self-similarity. Its
main unresolved gate is priority and significance in higher-dimensional
signature/ordinal-representation theory. The prior finite-face theorem's
priority audit also remains open. The next comparison should test whether
the general one-stream theorem follows from an existing ordinal uniqueness
or interspersion theorem. If so, credit that theorem and retain the result
as background. No contribution manuscript or completed-goal claim follows
from the present internal proof package alone.

The existing checker now includes zero-mask reconstruction, exact placement
and clock-increment identities, circular-interval endpoint controls, finite
rational-lattice bounds, rigorous prime-log intervals, the 20,000-entry
hidden-comparison witness, and non-vacuous quadratic-irrational certificates
at large clock heights. The latter use an exact special rank-difference
formula checked against direct lattice sums; they do not pretend to generate
the enormous full streams. Finite tests do not verify the transcendence
theorem, the infinite limits or priority. These rest on the written proofs
and the stated classical input. No core kernel or earlier manuscript changes.

## Finite recovery under a quantitative separation promise

The second iteration established infinite-data uniqueness, an arithmetic
consequence and a finite height-error estimate. It left three algorithmic
gaps: obtaining the estimate without calibrating an unknown generator,
identifying the additive atoms from estimates, and bounding how many input
symbols must actually be read. The following statement closes these gaps
under an explicit numerical promise. It does not decide priority or give
a practical random-access improvement.

Continue to normalize the observed weight to 1. The supplied class parameters
are an integer m>=1, positive rationals h,H,c with h<=H and c<=1, and an
integer mu>=1. The unknown remaining weights satisfy theta_i in [h,H] and

    abs(z_0 + sum_i z_i theta_i) >= c / ||z||_infinity^mu
    for every nonzero z in Z^(m+1).                      (ND)

In particular the weights are rationally independent. This is a promise
about the source, not something certified by a finite observed mask. The
input consists of exact mask bits, with 1 meaning that the observed exponent
is zero. Unobserved coordinate names cannot be recovered; use their increasing
weight order as the canonical labeling. The first output vector is all zeroes.
Dimension and weight bounds alone cannot replace (ND): the construction in
Theorem 9 can stay in the fixed class m=2, 1<=theta_i<=3, arbitrarily close
to (sqrt(2),2 sqrt(2)). Its arbitrarily long masks still conceal a rank-6
change. The separating form 2 theta_1-theta_2 approaches zero there.

**Theorem 10 (a uniform finite observation bound).** For every n>=2 and every
class just specified, a computable number L of initial mask bits determines
the first n full exponent vectors with this canonical labeling. No actual
generator value or rational approximation to one is supplied. One explicit
choice, using any integer K>=max(1,ceil(H)), is

    B_atom = max(1,2 ceil(K/h)),
    eta = c/B_atom^mu,     delta = c/n^mu,
    epsilon = min(1/4,h/4,eta/16,delta/(8 m n)),
    J = ceil(H)+1,        c_pair = c/J^mu,
    B_0 = ceil(32/epsilon),
    Q = ceil(B_0^mu/c_pair)+1,
    a = ceil(max(1, 4*2^m H Q/epsilon,
                    4 m(2K+1)(K+mH)/epsilon)),
    L = ceil((a+K+1+mH)^(m+1)/((m+1)! h^m)).            (FR1)

These are sufficient bounds, not optimized estimates. For fixed
m,h,H,c,mu,K, they give

    L = O(n^((m+1)(mu+1)^2)).                            (FR2)

The polynomial is in the numerical rank n, not its bit length. Neither
the exponent nor the constants are uniform as dimension, separation
parameters or prime sizes vary.

**Proof: a uniform height estimate without calibration.** Pick any one of
the actual remaining generators alpha. Dirichlet's elementary approximation
argument at integer budget Q gives a reduced fraction p/b with

    1<=b<=Q,   abs(b alpha-p)<=1/Q,
    abs(alpha-p/b)<=1/(bQ)<=1/b^2.

For completeness, place the Q+1 fractional parts of 0,alpha,...,Q alpha
in Q equal intervals, subtract two in the same interval, and then reduce
the fraction. Reduction preserves these displayed inequalities. Since
abs(p)<=Hb+1<=Jb, (ND) applied to this two-coordinate form implies

    abs(b alpha-p)>=c_pair/b^mu.

The choice Q>B_0^mu/c_pair therefore forces b>B_0. The fraction exists
for every admissible source; the algorithm does not have to find it.
The first two terms of (SS3) are now at most epsilon/4 each.

Write S=sum_i theta_i. Covering by unit cubes gives

    T^m/(m! product_i theta_i) <= M(T)
      <= (T+S)^m/(m! product_i theta_i),    T>=0.

Consequently

    [M(a+K)-M(a)]/M(a+1)
      <= (1+(K+mH)/a)^m-1
      <= 2m(K+mH)/a.

For the last inequality, m(K+mH)/a<=1, and
(1+x)^m<=exp(mx)<=1+2mx in this range. The choice of a in (FR1)
ensures the range and makes the last term of (SS3) at most epsilon/2.
Thus, simultaneously for every face height t_r<=K,

    abs(R_a(r)-t_r)<=epsilon.                           (FR3)

**Proof: all necessary observations are in the charged prefix.** Recover E
from its zero mask by the deletion recursion already proved. Let R=M(K),
obtained as the number of zeroes through A(K,0). This identifies precisely
the face indices 0<=r<R with height at most K. Their required placements
A(a,r), and also A(a+1,0), occur no later than

    A(a+K,0)=F(a+K).

The full (m+1)-dimensional cube bound, with product of weights at least h^m
and sum at most 1+mH, gives F(a+K)<=L. Thus the algorithm can read at
most L bits, stopping at the first occurrence of symbol a+K. Set

    D=A(a+1,0)-A(a,0),
    u_r=A(a,r)-A(a,0),    hat_t_r=u_r/D,    0<=r<R.

All these quantities are integers or exact rational numbers computed from
the input. In particular, an oracle for large placement ranks has not been
silently substituted for reading the mask.

**Proof: recover atoms by separated finite tests.** A positive face point
t_r<=K is decomposable exactly when t_r=t_s+t_t for some 0<s,t<r.
For each such triple, the exponent difference has maximum absolute
coefficient at most 2 ceil(K/h). If its corresponding linear form is
nonzero, (ND) bounds its absolute value below by eta. Since each of the
three estimates has error at most epsilon and 6 epsilon<eta,

    t_r=t_s+t_t  if and only if
    abs(hat_t_r-hat_t_s-hat_t_t)<=3 epsilon.             (FR4)

This decides all decompositions by rational comparisons. The indecomposable
points are precisely the m generators: a vector of degree at least two
splits into two positive vectors, and rational independence prevents a
degree-one vector from having another representation. All generators are
present because theta_i<=H<=K. Moreover their height order is their
canonical coordinate order. Replacing each generator by its estimate
therefore gives m positive weights, each within epsilon of its true value.

**Proof: recover the entire requested prefix.** The first n vectors of any
positive-weight order lie in {0,...,n-1}^(m+1), because a coordinate with
value v has v+1 distinct pure powers no later than that point. Any two
distinct vectors in this box have a true height difference at least delta
by (ND). Replacing only the m remaining weights changes that comparison
by at most mn epsilon<delta/2. All comparisons in the box are preserved.
The first n vectors for the true and estimated weights both lie in the
box, so their prefixes agree. Use the positive integer weights

    (D,u_(first atom),...,u_(last atom))

in a standard heap enumerator to produce this prefix exactly. Ties within
the relevant prefix are excluded by the same separation estimate. This
proves the algorithm and observation bound. Finally epsilon^(-1) is
O(n^(mu+1)); hence Q=O(epsilon^(-mu)),
a=O(epsilon^(-(mu+1))), and the cube bound yields (FR2). QED.

After computing the parameters, a direct implementation reads at most L
bits and stores O(L log L) bits for stream reconstruction and placements.
The number R of tested face heights is bounded by

    R<=ceil((K+mH)^m/(m! h^m)).

Atom detection uses O(R^3) exact comparisons; prefix enumeration uses
O(mn log(mn)) heap operations. Their integers have O(log L+log n) bits,
in addition to the supplied rational-parameter precision. These costs can
be enormous. The bounds do not recognize arbitrary invalid masks, decide
the global promise (ND), recover exact real weights from finite data,
or provide a fast unranking kernel. The case n=1 is the known zero vector.

The class includes nontrivial examples: (1,sqrt(2)) satisfies (ND) with
c=1/3, mu=1 by its quadratic norm. More generally
(1,sqrt(2)/v,sqrt(3)/v), integer v>=1, satisfies it with
c=1/[v(v+4)^3], mu=3: multiply the form by v and take its degree-four
algebraic norm, a nonzero integer; each other conjugate is bounded by
(v+4)||z||_infinity. This also checks scales where several monoid points
occur before the second generator, rather than treating early terms as
automatically indecomposable.

For any fixed prime basis, effective linear-form bounds supply some
admissible c and integer mu after normalization. This follows from the
known Baker--Wustholz bound in [Etessami--Stewart--Yannakakis, Theorem 3.1](https://arxiv.org/pdf/1304.5429):
with primes fixed, its exponent is a constant times log(max(e,||z||)).
Delete zero coefficients, take a uniform constant over the finitely many
subsets, and divide by the observed prime's logarithm. The constants depend
on the primes and can be extremely poor. This is an application of known
logarithm separation, not a new bound on prime logarithms or prime gaps.
If the bases are supplied already, their order can of course be generated
directly; Theorem 10 concerns recovery from an otherwise unlabelled mask
with supplied numerical class bounds.

## Focused ordinal-representation comparison

The following comparisons narrow specific candidate reductions. They do
not amount to a complete search of representation theory.

[Molcho--Wise, Theorem 2.1.3.6](https://arxiv.org/pdf/1807.11364)
states the classical Holder uniqueness theorem for an Archimedean ordered
abelian group. Its proof uses rational comparisons of integer multiples.
That theorem assumes the group operation is given. Here the mask names
face ranks without naming their addition operation. Applying it before
recovering that operation would assume information absent from the input.
The periodic deformation example above demonstrates the issue. This
comparison concerns the cited theorem and proof, not the rest of that paper.

[Wakker, 1991, Theorem 3 and Corollary 7](https://personal.eur.nl/wakker/pdfspubld/91.3algpart1jmp.pdf)
require restricted solvability. The natural placement domain fails it:
for weights (1,sqrt(2)), (1,0)<(0,1)<(2,0), but no integer clock level
at fixed face index 0 equals the middle alternative. Thus neither the
named general theorem nor its discrete corollary applies directly.
Their domains also have different structural assumptions. This hypothesis
failure does not exclude other measurement-theoretic proofs of our results.

[Luce--Narens, 1994, Section 2.2](https://sites.socsci.uci.edu/~lnarens/1994/Luce%26Narens_Book%20Chapter_1994.pdf)
already discusses numerical approximation by standard sequences and the
ordered-group framework. Recovering quantities from ordinal data is
therefore an established objective. No general priority claim, or claim
to solve their broader finite-to-continuum problem, follows from our formula.

[Fishburn--Roberts, 1988](https://doi.org/10.1016/0165-4896(88)90046-7)
is a further finite-uniqueness comparison. Its publisher abstract and author
record were inspected; its full theorems have not been checked here.
Finite exact parameter identification and the present recovery of a finite
pattern under separation promises must be compared carefully. The original
Luce--Tukey article, [DOI 10.1016/0022-2496(64)90015-X](https://doi.org/10.1016/0022-2496(64)90015-X),
was likewise inspected only at abstract level. Neither reading limitation
is evidence of novelty. The original finite-face priority audit also remains
open. No contact with authors or publication action is authorized.

## Third-iteration verification and decision

Reproduce the checks with

```sh
python3 -B scripts/run_face_reconstruction_checks.py
./build.sh
python3 tests/test_smoke.py
```

The report retains all previous controls and adds 192 exact parameter-budget
checks, 1,330 nonzero algebraic-norm controls, 48 adversarial accurate
multi-generator snapshots and 192 independently checked additive witnesses.
Their reconstructed prefixes are compared with independently sorted full
boxes using rigorous root intervals. Four quadratic examples execute the
complete placement-statistic decoder at the actual sufficient clock using
the exact formula

    A(a,r)-A(a,0)=r(a+1)+sum_(j=1)^r floor(j sqrt(2)),
    D=floor((a+1)/sqrt(2))+1.

The previous controls checked this formula against direct counts. The
new examples do not generate or read their enormous full masks. A real
512-entry mask is deliberately rejected as insufficient; the decoder must
not silently substitute snapshot tests for adequate raw input. The general
raw-input observation bound rests on the uniform proof above.
Four additional controls run the entire raw-mask reader successfully at
smaller quadratic clocks with independently certified actual height errors.
They check the stopping rank against an independent exact lattice sum and
the output against a full-box oracle. These exercise mask parsing, placement
extraction and decoding together; their smaller clocks are not executions
of the conservative sufficient-clock formula (FR1).

This is progress on the finite-certification gate: atoms, input length and
data-access costs are now explicit under a natural quantitative promise.
It is not yet evidence of a significant new paper. The large bounds limit
the computational claim; exact primary comparisons still limit the novelty
claim. Continue the original mathematics goal, with neither completion nor
publication readiness asserted. Next compare the finite recovery statement
against existing ordinal-reconstruction results, seek a substantially sharper
bound or a matching obstruction, and finish the face-hierarchy priority
comparison before choosing a manuscript contribution. Preserve all earlier
research and the paused ML work.

## A sharp dimension threshold for finite interleaving

This fourth iteration resolves the first iteration's open leading-constant
question, and also the next logarithmic term, up to a bounded additive error.
The earlier sections retain the historical proofs and their assessments.
The result concerns face words, not the different finite zero-mask input
of Theorem 10, and needs no Diophantine separation constants.

Put gamma=3 sqrt(3)/2. All logarithms in this section are natural. For s>=2 set

    q=floor(s/2),     D_s=binomial(s+q-1,q),
    p=floor((s-1)/2),
    B_s=max(2^(s-1),
            1+ceil(binomial(s+p,p)/s),
            1+ceil(D_s/48)).                            (SH1)

The constant 48 is a convenient sufficient constant, not an optimum.

**Theorem 11 (strict halfspaces of a composition layer).** Let X be uniform
on the nonnegative integer vectors of dimension s and total degree q. If
beta is a nonzero real vector with sum_i beta_i=0 and Y=beta.X, then

    E[Y^4]<=12 E[Y^2]^2,
    Pr(Y<0)>=1/48,    Pr(Y>0)>=1/48.                   (SH2)

**Proof.** Write (q)_j for a falling factorial and s^(j) for a rising
factorial, and let u_j=(q)_j/s^(j). Products with j>q are zero. The elementary
factorial-moment identity for the uniform composition distribution is

    E[product_i (X_i)_(a_i)]
      = (q)_(sum a_i) product_i(a_i!)/s^(sum a_i).

Indeed, the generating function for the numerator sum is
product_i(a_i!) t^(sum a_i)/(1-t)^(s+sum a_i); take its coefficient of t^q
and divide by binomial(s+q-1,q). This gives the identity without invoking
any distributional approximation. Expanding powers, using sum beta_i=0,
and putting A=sum beta_i^2, C=sum beta_i^4, gives

    E[Y^2]=(u_1+u_2) A,
    E[Y^4]=(u_1+7u_2+12u_3+6u_4) C
             +(3u_2+6u_3+3u_4) A^2.                  (SH3)

Here A>0 and C<=A^2. Set r=q/s and t=u_2. For s>=2,

    1/3<=r<=1/2,    0<=t<=r^2,
    u_3<=rt,       u_4<=r^2 t.

Thus E[Y^4]/A^2 <= r+(10+18r+9r^2)t. To verify that this is at most
12(r+t)^2, subtract it from that expression. Its derivative with respect
to t is at most 6r+15r^2-10<0 throughout the specified interval. Its
minimum is therefore at t=r^2, where it equals

    r(2r+6r^2+3r^3-1)>0.

The last polynomial is increasing and positive at r=1/3. This proves the
moment bound. Put sigma^2=E[Y^2]>0. Interpolation, or two applications of
Holder's inequality, gives

    E[abs(Y)] >= sigma^3/sqrt(E[Y^4]) >= sigma/sqrt(12).

Since E[Y]=0, E[Y_-]=E[abs(Y)]/2. Cauchy--Schwarz gives
E[Y_-]<=sigma sqrt(Pr(Y<0)), proving the negative bound. Apply the same
argument to -Y for the positive bound. QED.

The factorial moments and the moment method are standard tools. Their
derivation here makes the particular uniform constant and the strict
inequalities checkable; no separate novelty claim for that machinery is made.

**Lemma 12 (lattice points preceding a full-support comparison).** Let
w_1,...,w_s>0, T>0 and sum_i w_i<=2T. Then

    #{x in Z_>=0^s : w.x<T} >= D_s/48.                 (SH4)

If the w_i are not all equal, center them by beta_i=w_i-(sum_j w_j)/s.
Theorem 11 supplies at least D_s/48 degree-q points strictly below
their mean q sum_i w_i/s, which is at most T. If all weights equal w,
the degree-q layer is strictly below T unless qw=T. In that exceptional
case s=2q; the degree-(q-1) layer is below T and has size

    binomial(s+q-2,q-1) = D_s q/(s+q-1) >= D_s/3.

This also proves (SH4) in the constant case, which must not be inferred
from a zero-variance moment argument.

A useful additional bound follows by taking all vectors of degree at most
p=floor((s-1)/2). In each orbit under cyclic coordinate permutations,
the average height is degree times sum_i w_i/s, strictly below T. Every
orbit therefore supplies a strict point below T. Its size is at most s,
so the number of strict points is at least ceil(binomial(s+p,p)/s).

For a primitive mixed-sign normal z with support size s, the lowest tied
pair u=z_+, v=z_- at height T satisfies sum_(i in support) w_i<=2T.
Together with Lemma 2, these strict counts imply that its first tied rank
at a generic crossing is at least B_s. This charges the rank before either
member of the tied pair; it does not count the later member as the first.

**Theorem 13 (sharp order of the interaction rank).** Let J(s) be the least
first tied rank among positive-weight comparisons with primitive normal of
support s and no second independent relation at the crossing. Uniformly as
s grows,

    J(s)=Theta(gamma^s/sqrt(s)).                        (SH5)

The upper-bound constructions can be perturbed to agree for any prescribed
finite length on every proper coordinate face while differing in the full
word. They can also be realized by distinct prime bases, sharing their
first s-1 labeled bases, using the already established prime approximation.

**Proof.** The lower bound follows from (SH4) and Stirling's formula:

    log D_s=s log(gamma)-(1/2)log(s)+O(1).

For the upper bound, put a=floor(s/2), b=ceil(s/2), and use a positive
coordinates versus b negative coordinates in z, all with absolute
coefficient one. Start with weight b/a in the first group and weight 1
in the second group. Their common height is b. Choose a generic point on
this wall within eta<1/(s+1) of these weights. Its common height is at
most b+a eta, while every vector of degree at least b+1 has height at
least (b+1)(1-eta)>b+a eta. Thus the first tied rank is at most
binomial(s+b,b)=Theta(gamma^s/sqrt(s)). The wall can be chosen so that
its only rational relation is generated by z; every proper face is then
independent. The finite strict inequalities specifying any prescribed
proper-face prefixes survive a sufficiently small perturbation of the
last weight to either side. The only full tie at the relevant height is
u,v, so the full words change there. Apply the simultaneous Bertrand-based
prime realization from Theorem 5 to these finite inequalities. This uses
no new prime-distribution theorem. QED.

**Theorem 14 (face dimension to bounded additive error).** Let R_all(n)
be the least face dimension sufficient uniformly in ambient dimension to
determine an n-term full word from the labeled n-term words of all faces
of at most that dimension. The weights are positive and rationally
independent, as in Theorem 4. Define R_prime(n) in the same way restricted
to distinct prime bases. Then both functions satisfy

    R(n)=[log n+(1/2)log log n]/log(3 sqrt(3)/2)+O(1). (SH6)

One explicit sufficient dimension for n>=2 is

    r_sharp(n)=max{s>=2 : B_s<=n},                     (SH7)

capped by the ambient dimension. The set is finite because B_s>=2^(s-1).
This bound is sufficient; it is not asserted to equal the exact finite R(n).

**Proof.** Repeat the common-face cone path of Theorem 4. A crossing that
changes the first n vectors has lowest tied rank at most n. Lemma 12 and
Lemma 2 force B_s<=n for its support size s. Its face of size s would
therefore be among the supplied data if faces through r_sharp(n) were
given. The crossing changes that face's n-term word, a contradiction.
This proves the sufficient dimension for both weight and prime inputs.

Since D_s=Theta(gamma^s/sqrt(s)), a possible crossing satisfies

    s log gamma-(1/2)log s <= log n+O(1).

The existing bound s=O(log n) now gives
s <= [log n+(1/2)log log n]/log gamma+O(1). For the reverse inequality,
Theorem 5 supplies prime obstructions of dimension 2k whenever
binomial(3k,k)<=n, with agreement through M=n on every proper face.
Stirling's formula gives log binomial(3k,k)=2k log gamma-(1/2)log k+O(1).
Taking k=floor([log n+(1/2)log log n]/(2 log gamma)-C) for a sufficiently
large fixed C proves the required lower bound. Consequently
R_prime(n)<=R_all(n) are squeezed to within a bounded additive error of
the same expression. QED.

This improves the previous Theta(log n) statement to a sharp first and
second logarithmic term. In base-two units the leading coefficient is
approximately 0.726, rather than the earlier sufficient coefficient 1.
It explains how many-coordinate interactions can first appear in a finite
smooth-number word. It does not make requesting all faces inexpensive:
their number can still be sum_(j<=r) binomial(d,j), each with n vectors.
No faster random-access algorithm or improved prime-gap bound follows.

## Fourth-iteration primary comparison and publication gate

[Wagner, Chapter 3, Observations 3.13--3.14 and Theorem 3.15](https://ti.inf.ethz.ch/ew/research/phdtheses/wagner03.pdf)
bound simplex volumes and encode cuts using a full set of difference
directions with dimension-dependent factors. Those arguments remain
credited. They do not directly supply the present constant-fraction count
from one full-support comparison when dimension grows. A complete priority
comparison must still include subsequent uses of the encoding argument.

[Maclagan, Section 3 and Lemma 4.1](https://arxiv.org/pdf/math/9809134)
and [Christian--Conder--Slinko, Theorem 3](https://www.math.auckland.ac.nz/~slinko/Research/isipta17.pdf)
provide the Boolean flip and determining-comparison framework. These
inspected statements use subsets, whose coordinates are 0 or 1. Here
repeated powers contribute enough additional points to change the
interaction-rank scale from 2^s to gamma^s/sqrt(s). Their general methods
are not excluded as a possible route to the same conclusion.

[Alon--Huang--Sudakov, introduction](https://arxiv.org/pdf/1104.1753)
credits averaging and partition methods for nonnegative subset sums.
The uniform objects here are weak compositions with repeated coordinates,
not subsets. Neither averaging nor a moment argument is introduced as new;
the candidate contribution is the sharp ordinary-prefix information law.

[Cheng--Basu, Theorem 2.4 and its definitions](https://arxiv.org/pdf/2603.00286)
concern mixed-integer centerpoints under a linear-radius cube condition.
Those specific hypotheses do not give this thin composition-layer estimate
directly. Our centroid q/s is generally not an integer point. The present
lemma therefore does not prove Oertel's integer-centerpoint conjecture.
Only the stated definitions, theorem and introductory comparison were
inspected; no full audit of that paper is claimed. The older full
Vasiliev--Pavlov reading is preserved in the original research note; its
reduced constraints still need a direct comparison with this dimension law.

The checker adds complete small composition-layer moment and strict-count
controls, exact rational envelope checks, independent cyclic orbits and
odd-support rational obstructions in dimensions 3,5,7,9. The latter check
every proper face through 64 terms and the first full change by heap,
independent pointer merge and a sufficient degree simplex. They complement
the retained even-support witnesses. Finite rank brackets use exact
binomial integers, not a fitted asymptotic curve. The uniform proof, rather
than these tests, establishes (SH5)--(SH6).

The mathematical candidate is now substantially sharper and has a coherent
structural consequence for finite interleaving. This supports continuing
the full goal, not declaring it finished. Novelty and publication significance
remain unestablished. The next gate is a theorem-level priority comparison
of the composition halfspace estimate and the sharp dimension law, followed
by a manuscript decision that states exactly what is new. Preserve the
zero-mask results, earlier manuscripts, core kernels and paused ML research.

## Sparse ordering certificates and the exact interaction parameter

The fourth iteration was progress: it proved a sharp uniform law and saved
exact controls, while leaving priority and a certificate consequence open.
This continuation gives that consequence and compares the published
constraint reductions directly. Numerical prime values are not extra input
in the face-information parameter: primes specify the realizing class of
the exponent patterns. Arithmetic access to supplied primes is a different
problem, for which the core toolkit remains unchanged.

For a realizable ordered prefix W=(a_1,...,a_n), let C(W) be the open cone
of positive weights realizing precisely that strict prefix. Its familiar
finite description consists of positivity and

    w.(a_(j+1)-a_j)>0,                   1<=j<n,
    w.(b-a_n)>0,                        b in min(N^d\W). (SC1)

The minimal excluded border is finite: every excluded point dominates a
border point. This cone description, monomial frontiers and finite-weight
realization are established background, not new constructions here.
The support of an oriented comparison normal counts coordinates with a
nonzero coefficient after common factors have been canceled.

**Lemma 15 (small interaction thresholds).** The parameter of Theorem 13
satisfies J(2)=2, J(3)=5 and J(4)=10.

For support three, a nonbinary normal has V(z)>=12 and first tied rank
at least 6 by Lemma 2. A binary normal splits one coordinate against two.
The box gives three strict points below T. One of the two weights is
strictly smaller than T/2 at a generic wall, so its pure square supplies
another point outside the box. This gives first rank at least 5.

For support four, a nonbinary normal has V(z)>=24. In the binary case
the box gives seven strict points. For a two-against-two split, each pair
has a member of weight below T/2, supplying two additional pure squares.
For a one-against-three split, at least two weights are below T/2 and
the smallest is below T/3; their two squares and one cube supply three
additional points. Equality of two weights, or all three in the last
group, would give a second independent relation, excluded at a generic
crossing. Thus the first rank is at least 10. The support-two bound is
immediate. Upper witnesses are the walls

    (1,1),                 height 1,   first rank 2;
    (100,49,51),           height 100, first rank 5;
    (100,103,101,102),     height 203, first rank 10.

At the displayed height each has only the intended tied pair. Generic
perturbation within the wall preserves the strict comparisons to that
height and removes other relations. A sufficient degree simplex verifies
the counts directly. These small observations refine the finite cutoff;
no separate novelty claim for them is needed.

Define B'_s=max(B_s,J(s)) for s=2,3,4, and B'_s=B_s otherwise, and put

    r_cert(n)=max{s>=2:B'_s<=n},    n>=2;
    r_cert(1)=0.                                         (SC2)

**Theorem 16 (delete dense constraints).** For any nonempty realizable
strict prefix W, retain positivity in (SC1) and delete every other
inequality whose normal has support greater than r_cert(n). The remaining
strict inequalities define exactly C(W).

**Proof.** The cone has a strict interior point by the realizability
assumption. Its closure is therefore the full-dimensional polyhedral cone
defined by the non-strict versions of (SC1). Every facet not on a coordinate
boundary has a relative interior point with positive weights and just one
independent rational relation. Crossing that facet changes the prefix.
The primitive normal z has support s and its lowest tied pair occurs
by rank n. Theorems 13--14 and Lemma 15 give B'_s<=n. Thus every such
facet normal has support at most r_cert(n). Coordinate facets are retained
by positivity. Every facet has a defining inequality among the original
finite list, so every facet remains represented after deletion.

For completeness, polyhedral duality expresses every other valid normal
as a nonnegative combination of the facet normals. This follows by
separating that normal from their cone if no such expression exists.
The combination is nonzero because the original inequality is strict at
the interior point. Strict positivity of all the facet inequalities
therefore implies every removed inequality strictly. This also addresses
points on intersections of redundant supporting hyperplanes. QED.

This is a concrete interval certificate for smooth-number ordering.
Generate W at exact positive rational weights, replay through n+1 to
check strict internal and boundary comparisons, and form its border.
For a rational box w_i in [l_i,u_i] with l_i>0, check only the retained
normals using

    min_(w in box) w.z = sum_(z_i>=0) z_i l_i
                           +sum_(z_i<0) z_i u_i > 0.   (SC3)

If these tests pass, every weight in the box has the same prefix. In
particular a box of rigorous prime-log bounds certifies its exponent
ordering. The realization check is essential: arbitrary impossible words
cannot invoke deletion of their high-support constraints. This procedure
still reads or generates the whole n-term word and its border, with at
most O(nd+d) displayed constraints before filtering. It assumes no free
face or placement oracle and claims no new running-time bound. It is a
research prefix certificate, separate from the CLI's independent rank
auditor and its `rank_certified` output.

For example W=(0,e_1,e_2,e_3) has n=4 and r_cert=2. Its three-coordinate
border comparison is implied by two retained comparisons:

    (1,1,-1) = (2,0,-1)+(-1,1,0).

Thus 2w_1>w_3 and w_2>w_1 imply w_1+w_2>w_3. This exact conic identity
also gives a small independent check of the general facet argument.

**Corollary 17 (one exact parameter for faces and facets).** For n>=2 set

    kappa(n)=max{s>=2:J(s)<=n}.

Then R_all(n)=R_prime(n)=kappa(n). It is also the largest support of a
positive-weight facet of an ordered-prefix cone of length n, over all
ambient dimensions. Any finite linear inequality description of all
such cones sometimes needs that many coordinates in a single normal.
Consequently the sharp law (SH6) describes optimal linear interaction
arity as well as face reconstruction.

For the upper bounds, any first changing wall or positive-weight cone
facet has support s with J(s)<=n, exactly as above. For the lower bounds,
take a generic wall attaining J(s)<=n. The two sides differ by this rank;
all proper-face words agree through any prescribed finite length after a
sufficiently small perturbation. Their simultaneous prime realizations
therefore give R_prime(n)>=s. The crossing also supplies a genuine
prefix-cone facet. A finite inequality description of a full-dimensional
cone must contain a normal parallel to each facet normal; replacing it
by normals of smaller support cannot define that facet. This proves the
claimed equalities and optimality.

There is a further consequence for the state fans of sufficiently generic
n-point ideals. Use the published identification with the corner-cut
polyhedron, for example [Babson--Onn--Thomas, Section 5](https://arxiv.org/pdf/math/0207135).
The genericity condition is that every n-element subset of
V_n={v:product_i(v_i+1)<=n} is a quotient basis. The largest support of a
positive-weight state-fan wall, over all ambient dimensions and such
configurations, is also exactly kappa(n).

The upper bound follows by refining a corner-cut cone by its internal
ordered-prefix cones. To obtain a cut boundary at any prescribed n>=J(s),
start with the generic s-coordinate wall at height T and first rank J(s).
Append n-J(s) new coordinate weights in

    (max(T/2,T-min_i w_i),T).

Each new coordinate contributes precisely its unit vector below T: its
square, its sums with old nonzero points and its sums with other new
coordinates lie above T. Choose the weights generically while keeping
the original relation. The first tied rank becomes exactly n, and its
normal still has support s. The two n-point cut sums differ by v-u, so
their common positive-weight facet has that support. The published
generic-configuration identification transfers this statement to the
state fan. This is a quantitative consequence of the free lattice; the
state-polyhedron correspondence itself is prior art.

## A boundary on extending the law to arbitrary quotient algebras

The genericity in the last consequence cannot be dropped. For a length-n
ideal, define A(n) to be the largest support of a noncoordinate state-fan
facet meeting the positive weight cone, uniformly in ambient dimension.
The following elementary boundary uses the established general framework;
it is not asserted as a separate new universal-basis theorem.

**Proposition 18.** For n>=2,

    2 floor(log_2(n+2))-2 <= A(n) <= 2 floor(log_2 n).  (SC4)

The upper bound is an immediate support consequence of the published
Hilbert-zonotope refinement, [Babson--Onn--Thomas, proof of Proposition 3.2
and Theorem 3.5](https://arxiv.org/pdf/math/0207135). Its directions are
differences of two members of V_n. Each member has at most floor(log_2 n)
nonzero coordinates. A state-fan facet must lie in one of the refining
zonotope's wall hyperplanes, giving the sum of those support bounds.

For the lower bound put k=floor(log_2(n+2))-1, and introduce disjoint
squarefree monomials U=x_1...x_k and V=y_1...y_k. Let

    l=n-(2^(k+1)-2),
    S=divisors(U) union divisors(V) union {t,t^2,...,t^l}.

Let J be the monomial ideal of the complement of S, and I=J+(U-V).
There are |S|=n+1 monomials before the relation. Every variable annihilates
both U and V modulo J, so the ideal generated by their difference there
is one-dimensional. The quotient by I therefore has length n.

For either leading choice of U-V, the monomial generators of J together
with this binomial form a Groebner basis. Indeed every nontrivial
Buchberger pair reduces to a nonconstant monomial multiple of the other
socle monomial, hence lies in J. Thus there are exactly two initial
ideals, J+(U) and J+(V). Their positive-weight cones meet at

    sum_i w_(x_i)=sum_i w_(y_i),

whose normal has support 2k. Both sides contain positive weights; the
padding variable does not enter this wall. This proves the lower bound.
For example a length-6 quotient already has a four-coordinate wall,
whereas a free-lattice wall with that support cannot first occur before
rank 10. The missing monomials explain the difference. In particular
the arbitrary-quotient scale is 2 log_2 n+O(1), rather than (SH6).

## Direct primary comparison and fifth-iteration decision

The supplied [Vasiliev--Pavlov paper](https://doi.org/10.1134/S0361768809020030)
was rechecked in Sections 5--6, Propositions 4--7, equations (11)--(16),
and Figures 7--8, including the rendered formulas. Proposition 5 keeps
comparisons with D_k, a subset of the known full chain. Its retained
normals can still have full support. Proposition 7 uses the chain's front
and earlier multiplication comparisons; that is different input from
separately supplied face words. At a primitive coprime wall the two
monomials have no common factor to cancel. These reductions do not
directly yield a support/rank cutoff or (SH6). They remain credited as
the finite-cone and frontier machinery underlying (SC1), not dismissed
because their statements use different wording. This closes the direct
comparison with those named reductions, not every possible deduction
from finite-ordering theory.

The Hilbert-zonotope proof and generic-configuration identification above
were checked in the primary preprint. The general support upper bound in
(SC4) is credited as a direct consequence. It cannot alone give the
smaller free-lattice law, and the socle example explains why the general
model must retain larger supports. The new candidate is the sharp
free-lattice interaction law and its optimal certificate and generic-fan
consequences, not generic point interpolation or the general algorithm.

Exact controls now include the small J(s) witnesses; complete integer
grids comparing full and sparse cones with both accepted and rejected
weights; explicit conic identities; rejection of unrealizable words and
boxes crossing real walls; and prime-log interval certificates checked
against independent exact product and pointer generators. Padding
witnesses check prescribed-size cut boundaries. Every nontrivial
Buchberger pair in the listed socle-quotient fixtures is checked for both
leading choices. These support the written general arguments; no tests
prove novelty or a running-time improvement. Earlier checks are retained.

This is progress toward the full goal: there is now an ordering certificate
consequence and a precise connection to a broader established mathematical
object, together with a counterexample to an overbroad extension. The
publication case is more coherent, but its novelty and significance are
still under assessment. A focused paper would center on the sharp
interaction law and these consequences. It should credit the composition
moment tools, cone machinery, prime realization and state-fan framework,
and should not imply that a new fast unranking algorithm or a theorem for
all finite quotient algebras has been obtained. The remaining priority
gate concerns the composition halfspace estimate, the sharp free-lattice
law and subsequent uses of the cited framework. The goal remains active.
