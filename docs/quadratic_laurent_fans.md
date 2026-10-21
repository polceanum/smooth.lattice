# A sharp quadratic fan bound for a quadratic coordinate and an inverse pair

Let \(I\subset K[x,y,z]\) have positive colength \(n\), over any
field, and contain a monic quadratic \(f(x)\) and \(yz-1\).
Then
\[
                 g(I)\le n^2,\qquad V(I)\le3n^3,\qquad
                              |\mathcal U(I)|\le6n^2.      \tag{L1}
\]
Both the quadratic fan and cubic total-cost powers are sharp,
including for local characteristic-zero quotients with \(x^2=0\).
The proof allows arbitrary dependencies and coupling.

Together with the [zero-product theorem](quadric_cross_fans.md),
this closes the whole class with a quadratic coordinate and a
deficient four-monomial box in the other two coordinates:
\[
 f(x)\in I,\quad \deg f=2,\quad f\text{ monic},\qquad
                     \dim_K\langle1,y,z,yz\rangle\le3.    \tag{L2}
\]
The upper bounds hold uniformly over all fields, and characteristic-zero
examples attain quadratic fan and cubic total wall-cost order.
Arbitrary single-defect surfaces
and arbitrary square-zero quotients remain outside this theorem.
No unrestricted \(G_3\) improvement or novelty is asserted.

## Laurent intervals and a translation-invariant rank test

Write \(Y\) for multiplication by \(y\); it is invertible, with
inverse multiplication by \(z\). All initial standard monomials
avoid \(x^2\) and \(yz\). If the \(x\)-layer is nonempty, a standard
set consequently has the form
\[
 \{y^j:-b\le j\le a\}\ \sqcup\
                       \{xy^j:-d\le j\le c\},            \tag{L3}
\]
where negative powers mean powers of \(z\). Put
\[
 \ell=a+b+1,\quad u=c+d+1,\quad
 \ell+u=n,\quad 1\le u\le\lfloor n/2\rfloor,\quad
 a\ge c,\ b\ge d .
\]
The inequalities follow from divisor closure.
Multiplication of every vector of (L3) by \(Y^b\) is invertible.
Thus independence depends only on
\[
 \delta=b-d,\quad 0\le\delta\le M:=\ell-u,
\]
and is equivalent to the independence of
\[
 1,Y,\ldots,Y^{\ell-1},\qquad
                  xY^\delta,\ldots,xY^{\delta+u-1}.       \tag{L4}
\]
For each fixed \(u\), let \(D_u\subseteq\{0,\ldots,M\}\)
be the set passing this rank test. It can be arbitrary for the
argument. For every \(\delta\in D_u\), all translations
\(0\le d<u\) give vector bases in (L3), with
\[
                 b=\delta+d,\quad
                 a=\ell-1-\delta-d,\quad c=u-1-d.         \tag{L5}
\]
No translation invariance of the positive monomial costs is
being asserted; only vector independence is invariant.

If the \(x\)-layer is empty, there are at most \(n\) possibilities:
the \(n\) consecutive Laurent intervals containing the origin.

## Two parity classes give short monotone optimization paths

Fix \(u\ge1\) and normalize \(w_y+w_z=1\), with \(t=w_y\in(0,1)\).
The \(x\)-cost \(uw_x\) is constant for this part of the proof.
Set \(s=b+d=\delta+2d\). Up to an additive term independent of
\(s,\delta\), the remaining basis cost is
\[
        \frac14\bigl((s+1-nt)^2+(\delta-Mt)^2\bigr).       \tag{L6}
\]
This follows by summing the pure-power exponents in the two hooks.
The feasible points are
\[
 \delta\in D_u,\qquad s=\delta+2d,\quad 0\le d<u.          \tag{L7}
\]

Separate (L7) by the common parity \(p\) of \(s,\delta\).
Each of the two resulting finite sets is closed under
coordinatewise minimum and maximum. Indeed \(\delta\) stays in
the allowed one-dimensional set, \(s\) keeps its parity, and
\[
                          0\le s-\delta\le2u-2
\]
is preserved by both operations.

For completeness, the monotonicity needed here has a direct proof.
Suppose \(P\) is a unique minimizer of (L6) at \(t\), and \(Q\)
at \(t'>t\), in one parity set. Their coordinatewise minimum
and maximum are feasible. For every coordinate in which
\(P_i>Q_i\), replacing the pair \((P,Q)\) by
\((P\wedge Q,P\vee Q)\) changes the sum of their costs by a
negative multiple of
\[
                         (t'-t)a_i(P_i-Q_i),
\]
where \(a_i=n\) for \(s\), and \(a_i=M\) for \(\delta\).
Other coordinates do not change that sum.
Optimality is contradicted if any such decrease occurs.
When \(M=0\), the \(\delta\)-coordinate is identically zero.
Thus both coordinates of the unique minimizers are nondecreasing.
This is the elementary sublattice form of monotone comparative
statics. The established general framework is Topkis,
[“Minimizing a Submodular Function on a Lattice,”
Operations Research 26 (1978), 305–321](https://doi.org/10.1287/opre.26.2.305).
The direct argument above is self-contained; no new general
monotone-optimization theorem is claimed.

Their ranges are \(0\le s\le n-2\), \(0\le\delta\le M\).
Each change in one parity set increases a coordinate by at
least two. There are therefore at most
\[
                     1+(n-2+M)/2=\ell                   \tag{L8}
\]
different unique minimizers in each parity set.

At a generic positive objective, an initial standard set is the
unique minimum monomial vector basis. It must in particular be
the unique optimum among (L3) with its fixed \(u\). It is then
unique within its parity class. Consequently at most \(2\ell\)
standard sets occur for that \(u\). Summing gives
\[
 g(I)\le n+2\sum_{u=1}^{\lfloor n/2\rfloor}(n-u)\le n^2. \tag{L9}
\]
The parity argument counts coherent optima, not all vector bases.
There can be cubically many order-ideal vector bases in this
carrier.

The positive adjacency graph is planar, every primitive cost is
at most \(n\), and every initial ideal has at most six minimal
generators, as in the zero-product case. This proves the other
two bounds in (L1), with all parallel wall pieces retained.

There is also a finite enumeration consequence. From commuting coordinate
matrices and a cyclic vector, compute the \(O(n^2)\) normalized ranks in
(L4). Form the lower envelopes of the affine-in-\(t\) costs in each parity
class; their total input size is \(O(n^3)\), and (L8) leaves only \(O(n^2)\)
unique minimizers. Recover their border normal forms and full positive
cones. Straight elimination and planar clipping give
\(O(n^5)\) field-and-geometric-arithmetic operations, excluding bit
complexity. The same bound holds in the zero-product and linear cases.
This is an elementary enumeration bound, not an optimal-time claim.

## An explicit sharp family and an exact fan count

For every \(n\), let \(q=\lceil n/2\rceil\). In
\(\mathbb Q[v]/v^n\), set
\[
               x=v^q e^{3nv},\qquad y=e^v,\qquad z=e^{-v}.
                                                               \tag{L10}
\]
Exponentials are truncated, \(yz=1\), and \(x^2=0\).
The element \(y-1\) has leading term \(v\), so the image algebra
has length exactly \(n\).

Every size-\(n\) order ideal in the carrier (L3) is a vector basis.
Its lower columns are Hermite evaluation conditions of order zero
at the consecutive integers \(-b,\ldots,a\); its upper columns
have order \(q\) at \(3n-d,\ldots,3n+c\).
Every upper node exceeds every lower node, and \(\ell\ge q\).
Thus the [ordered Hermite lemma](nested_jet_matroids.md) applies.
This uses characteristic zero; it is not inferred from a finite
field specialization.

Counting (L5), the number of all order-ideal vector bases in this family
is exactly
\[
 n+\sum_{u=1}^{\lfloor n/2\rfloor}u(n-2u+1)=\Theta(n^3).
\]
For example it is \(405>20^2\) at length \(20\).
This explains why the rank test alone does not give the quadratic
coherent count.

The positive fan is therefore the ordinary first-\(n\) prefix
fan of the two-layer Laurent carrier. A lower hook of size \(j\)
has \(j\) possible prefixes, with walls at
\(w_y/w_z=(j-a)/a\), \(1\le a<j\).
For a fixed \(u\ge1\), both hooks use the same planar ratio.
Their interval fans have exactly
\[
                 \ell+u-\gcd(\ell,u)=n-\gcd(n,u)
\]
common-refinement intervals: their common walls number
\(\gcd(\ell,u)-1\).
Every such pair is realized by some positive \(w_x\):
choose a threshold in each hook's cutoff gap and let their
positive difference be \(w_x\). This also works when \(\ell=u\)
by choosing different thresholds in the same gap.
The case \(u=0\) has \(n\) hooks. Hence exactly
\[
 g(I_n)=n+\sum_{u=1}^{\lfloor n/2\rfloor}
                                  \bigl(n-\gcd(n,u)\bigr)
                           =\Theta(n^2).                 \tag{L11}
\]

The same family has \(V(I_n)=\Theta(n^3)\).
For a lower bound take coprime positive \(a,b\) in
\([n/3,2n/5]\), put \(k=a+b\), and set the upper count to
\(u=n-k<k-1\), for sufficiently large \(n\).
At \(a w_y=b w_z\), the lower prefix of size \(k\) exchanges
\(y^a,z^b\). Choose \(w_x>0\) so the first \(u\) upper
monomials lie below that threshold and the next one above it.
This is an actual wall of cost \(a+b\ge2n/3\).
There are \(\Theta(n^2)\) such coprime pairs, by Möbius
inversion in this fixed-width proportional box: the main term
is its area divided by \(\zeta(2)\), with \(O(n\log(n+1))\)
error. Their directions are distinct. The cubic upper bound
in (L1) completes the assertion.
This sharpness statement concerns the whole wall family, not
the full-support subfamily alone.

## A single bilinear rank defect with a quadratic remaining coordinate

A nonzero relation in \(\langle1,y,z,yz\rangle\) either is
linear in \(y,z\), or after translating those two variables
has the form \(y'z'-c\).
Translations and nonzero coordinate scalings preserve all
monomial initial ideals for global orders: each monomial
expands into itself plus proper divisors, and the inverse
map has the same property.
If \(c=0\), apply the zero-product theorem. If \(c\ne0\),
rescale to \(y'z'=1\) and apply (L1).
The quadratic in \(x\) is unchanged.

For a linear relation, at weights away from its one price tie,
one of \(y,z\) is excluded. Standard sets lie in a two-layer
planar slab of cardinality \(n\), with at most
\(\lfloor n/2\rfloor+1\) shapes. Two such weight regions
cover all open initial cones. This gives \(O(n)\), and hence
the asserted quadratic bound, in the linear case.
The constant relation alone is impossible at positive colength.
This proves the all-field class statement (L2).

## Verification and limits

Run

~~~bash
python3 -B scripts/run_quadratic_laurent_checks.py
~~~

The report is results/research/quadratic_laurent_checks.json.
It tests invertible multiplication, normalized interval ranks,
translation invariance, both parity domains and their monotone
optimization paths, and complete positive fans from independent
border relations. It includes special nonreduced dependencies,
finite characteristics, nonlocal points and an irreducible
quadratic, and checks the exact sharp count.

The rank translation, integer quadratic-cost calculation,
sublattice uncrossing and elementary hook/Farey count are the
inputs to the upper proof. The ordered Hermite construction and
coprime-pair estimate are established tools. Exact priority of
the class theorem and exact count remains unassessed.
Neither the zero-product branch-cardinality proof nor its
all-order-ideal-bases bound is asserted for the inverse case.
The general single-bilinear-defect surfaces, higher-degree
remaining coordinates and unrestricted extremal endpoint remain
open, with all original publication gates unchanged.
