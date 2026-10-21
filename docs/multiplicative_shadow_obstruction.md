# A limit of individual rational-shadow certificates

The rational-shadow bound suggests replacing a uniform coordinate height by
the smaller product

\[
 \rho(q)=\prod_{i=1}^3(|q_i|+1),\qquad q\in\mathbb Z^3\setminus\{0\}.
\]

This note rules out a proposed shortcut. There are actual finite quotient
ideals, over every field, with an initial cone requiring
\(\rho(q)=\Theta(n^{4/3})\) before the plane \(q^\perp\) meets its
interior. In particular, a uniform \(O(n)\) certificate for each cone is
false. The examples have exactly **seven** initial ideals. This is a
feasibility obstruction, not a larger-fan construction or an improvement
of the unrestricted \(G_3(n)\) bound. Its precise priority is unestablished.

## What the plane must detect

For an open positive initial cone \(C\), put

\[
 m(C)=\min\{\rho(q):q\text{ primitive integral},\quad
                       C\cap q^\perp\ne\varnothing\}.       \tag{M1}
\]

An intersection with the **interior** is essential: a plane supporting
only a wall does not expose the chosen state vertex in that shadow.
The [existing product-threshold argument](rational_shadow_bound.md) bounds
the cells incident to a wall of small primitive product. That is a
different condition. The present obstruction also rules out obtaining
all cells by replacing the incident-wall condition with (M1) at \(T=O(n)\).

## Realizing three independent comparisons

Let \(N\) be an integer matrix with positive determinant and diagonal,
strictly negative off-diagonal entries, and positive inverse. Suppose
each diagonal entry exceeds the absolute values of the other entries
in its column. Define six monomials with exponent vectors

\[
 a_i=N_{ii}e_i,\qquad
 b_i=\sum_{j\ne i}(-N_{ij})e_j\quad (1\le i\le3).
\]

They form an antichain for divisibility. Let \(D\) be the union of their
divisor sets and \(J_D\) the ideal of all monomials outside \(D\). Set

\[
 I_N=J_D+(x^{a_i}-x^{b_i}:1\le i\le3).                    \tag{M2}
\]

All six monomials are distinct socle elements of the monomial quotient.
Thus the three relations are independent, and their multiples by any
variable vanish. Over any field the quotient has length

\[
 n(N)=\sum_iN_{ii}
        +\sum_i\prod_{j\ne i}(-N_{ij})-2.                 \tag{M3}
\]

Indeed \(D\) consists of the origin, the three axis segments, and
the positive-coordinate rectangles in the three coordinate planes;
these pieces are disjoint. Its cardinality is the expression in (M3)
plus three.

For any monomial order the initial standard set is exactly \(D\) with
the larger monomial in each pair removed. To prove this, the proposed
monomial initial ideal is contained in the actual initial ideal and
already has the correct colength. Equality follows. Hence there are
no additional hidden comparisons.

For positive weight \(w\) the three comparisons are the signs of
\(Nw\). Every sign pattern except \((-,-,-)\) occurs: choose a vector
\(y\) with the required signs, and make one of its positive entries
sufficiently large that \(w=N^{-1}y>0\). Strict positivity of every
entry of \(N^{-1}\) proves this assertion. All-negative signs cannot
occur even for a general global monomial order. The strictly positive
integer row vector \(\alpha={\bf1}^{\mathsf T}\operatorname{adj}(N)\)
satisfies \(\alpha N=\det(N){\bf1}^{\mathsf T}>0\); multiplying the
three supposed all-negative comparisons with multiplicities
\(\alpha_i\) contradicts positivity of monomials.
Consequently (M2) has exactly seven global monomial initial
ideals. Its all-positive comparison cone is

\[
 C_N=N^{-1}\mathbb R_{>0}^3\subset\mathbb R_{>0}^3.         \tag{M4}
\]

For any nonzero row \(q\), the plane meets (M4) if and only if
\(p=qN^{-1}\) has both a strictly positive and a strictly negative
entry. Zeros alone do not suffice. This follows by evaluating
\(qN^{-1}y=p\cdot y\) for \(y>0\).

## An explicit unimodular family

Put

\[
 L=\begin{pmatrix}2&-1&-1\\-1&2&-1\\-1&-1&2\end{pmatrix},\quad
 B=\begin{pmatrix}2&-1&0\\0&1&-1\\-3&2&0\end{pmatrix},\quad
 M(t)=tL+B,\quad N(t)=M(t)^3.                             \tag{M5}
\]

Direct multiplication gives \(\det M(t)=1\) and

\[
 M(t)^{-1}=3t^2{\bf1}{\bf1}^{\mathsf T}
 +t\begin{pmatrix}3&0&3\\4&1&4\\5&2&5\end{pmatrix}
 +\begin{pmatrix}2&0&1\\3&0&2\\3&-1&2\end{pmatrix}.        \tag{M6}
\]

Every entry is positive for \(t\ge1\). Expanding the cube gives

\[
 N(t)=9t^3L+t^2C+tE+F,                                   \tag{M7}
\]

\[
 C=\begin{pmatrix}54&-45&0\\6&15&-21\\-69&48&12\end{pmatrix},\quad
 E=\begin{pmatrix}34&-32&10\\25&-11&-8\\-53&49&-14\end{pmatrix},\quad
 F=\begin{pmatrix}5&-5&3\\9&-6&1\\-6&7&-5\end{pmatrix}.
\]

For every integer \(t\ge1000\), all off-diagonal magnitudes are
between \(8t^3\) and \(10t^3\), and the diagonals are between
\(17t^3\) and \(19t^3\). Thus (M2) applies, with

\[
 192t^6\le n(t)\le301t^6.                                \tag{M8}
\]

**Theorem.** For \(t\ge1000\), the cone \(C_{N(t)}\) of the
length-\(n(t)\) ideal (M2) satisfies

\[
 \frac{t^8}{256000000}\le m(C_{N(t)})\le18000t^8.
                                                               \tag{M9}
\]

In particular \(m(C_{N(t)})=\Theta(n(t)^{4/3})\), uniformly in
the field.

### Lower bound

Let \(q\) be integral and visible in (M4), and put
\(p=qN(t)^{-1}\). Unimodularity makes \(p\) an integer row;
it has both signs. Write \(P=\|p\|_\infty\ge1\) and
\(R=\|q\|_\infty\).

For every row \(u\) with both signs,
\(\|uL\|_\infty\ge\|u\|_\infty\): the difference between
its largest and smallest transformed entries is three times the
original range. The column absolute sums of \(B\) are at most five,
so

\[
 \|uM(t)\|_\infty\ge(t-5)\|u\|_\infty.
\]

Moreover \(uM(t)\) still has both signs. Otherwise multiplication
by the strictly positive inverse (M6) would make \(u\) have one
strict sign. Apply this observation three times:

\[
 R\ge(t-5)^3P\ge t^3P/8.                                \tag{M10}
\]

The identity

\[
 N(t){\bf1}=(9t^2+12t+3,\ 6t+4,\ -9t^2-18t-4)^{\mathsf T}
\]

implies
\(\lvert\sum_iq_i\rvert\le(18t^2+36t+11)P\le20t^2P<R/2\).
Hence at least two absolute coordinates of \(q\) are at least
\(R\) and \(R/4\), respectively. In particular

\[
 \rho(q)\ge R^2\bigl(1+\min_i|q_i|\bigr)/4.              \tag{M11}
\]

Suppose first that \(P\le t/1000\). For each column \(j\), the two
integer coefficients \(pL_j,pC_j\) cannot both vanish. Indeed their
common orthogonal line is generated by, respectively,

\[
 (75,84,66),\qquad(111,93,75),\qquad(30,12,21),
\]

all strictly positive rows. A nonzero row on any such line cannot
have both signs. The column absolute sums of \(C,E,F\) are at most
\(129,112,20\). If \(pL_j\ne0\), integrality and (M7) give

\[
 |q_j|\ge9t^3-129t^2P-112tP-20P\ge8t^3.
\]

If \(pL_j=0\), then \(pC_j\ne0\), and

\[
 |q_j|\ge t^2-112tP-20P\ge t^2/2.
\]

Thus every coordinate has magnitude at least \(t^2/2\), and
(M10)--(M11) give \(\rho(q)\ge t^8/512\).

If \(P>t/1000\), (M10) gives \(R>t^4/8000\), so (M11) gives
\(\rho(q)>t^8/256000000\) without any lower bound on the third
coordinate. This proves the lower half of (M9), including vectors
with a zero coordinate and planes passing through cone-boundary rays.
The latter count only when the plane also meets the interior.

### Upper bound

Take \(p=(1,-1,0)\). Then

\[
 q=pN(t)=\bigl(27t^3+48t^2+9t-4,
              -27t^3-60t^2-21t+1,
               21t^2+18t+2\bigr).                       \tag{M12}
\]

It is primitive because multiplication by a unimodular integer matrix
preserves primitivity. Its plane meets the cone, since \(p\) has both
signs. The three factors in \(\rho(q)\) are at most
\(28t^3,28t^3,22t^2\), giving the upper bound in (M9).

## Published comparison and remaining scope

[Onn--Rozenblit, arXiv:1208.5639v1](https://arxiv.org/pdf/1208.5639),
Theorem 2.3 and its proof, bound a projected integer polytope by
\(O_d((e\|W\|_\infty)^{d(d-1)})\), where \(e\) is edge complexity.
Theorem 4.1 and its proof specialize this to matroid bases using unit
exchanges and greedy optimization. They bound all edge directions by
a coordinate cube and refine by its zonotope. In the planar axis
projection with exponent bound \(H=n\), that route gives \(O(n^2)\),
whereas the repository's hyperbolic, rank-sensitive shadow theorem gives
\(O(n\log(n+1))\). This is a comparison of the inspected proofs, not
an exhaustive current-priority conclusion. No identity with a later
journal version is assumed.

[Deza--Manoussakis--Onn, 20 December 2015](https://optimization-online.org/wp-content/uploads/2015/12/5266.pdf),
Section 3.1, Proposition 3.1 and its proof, replaces the cube-generated
zonotope by primitive bounded-coordinate generators. Proposition 3.2
constructs a matching type of lower comparator using a uniform matroid;
Theorem 3.3 sandwiches the maximum between two Euler-zonotope vertex
counts. These results concern arbitrary matroids and bounded criteria.
They supply established direction-refinement machinery, not a checked
hyperbolic shadow bound or a statement about the original-coordinate
finite-quotient count. The full strongest-current literature audit
remains unfinished.

The construction above is a special case of the socle-identification
class already excluded as a large-fan mechanism in
[socle stability](initial_cone_stability.md). The contribution here is its
uniform quantitative obstruction to an individual product certificate.
Neither a failed search nor these two comparisons establish novelty.

To improve \(G_3\), difficult cones would need a **collective** count.
The present family forces the individual threshold above \(n\) while
contributing only seven cells. It therefore leaves open a small total
exceptional count. Do not enlarge this example family to seek many
initial ideals: the existing socle bound already limits that route.
The next proof task is an aggregate charge using shared quotient
dependencies or state geometry, together with the outstanding exact
priority comparison for the shadow theorem. No such aggregate bound
is asserted here.

## Reproduction

```bash
python3 -B scripts/run_multiplicative_shadow_checks.py
./build.sh
python3 tests/test_smoke.py
```

The report saves exact matrices, endpoints, all seven cone polygons and
strict weights, polynomial matrix identities, the primitive witness,
and bounded independent visibility controls. Small quotient models are
also checked through sparse commuting multiplication on every basis
element. Larger cubic-family fixtures use the proved symbolic
construction, without allocating a quotient of their enormous length.
Finite checks support the uniform proof; they establish neither its
priority nor a new unrestricted count or publication endpoint.
