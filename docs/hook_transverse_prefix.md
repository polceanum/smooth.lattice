# Off-axis prefix ranks control hooks with constant axis units

The maximum axis-prefix gap does not force a large fan. A second rank test,
on the first nonzero transverse coefficient, gives a uniform bound for a
central-only hook over every field. In particular a nonempty open locus
inside the maximal-axis-gap stratum has a quadratic polylogarithmic fan
bound, including all positive weight regions and arbitrary higher transverse
coefficients. A dual functional supplies the reverse estimates without
assuming that the singular axis-prefix matrices are invertible.
The same-layer extension gives a coefficient-independent complete list of
quadratically many weights up to logarithms on the regular locus.

This extends a restricted family bound. The unrestricted sharp-growth
question, an actual polynomial lower separation, and exact priority remain
open. The two outer generator images are zero throughout this note.

## The coefficient conditions and bound

Let \(t\ge2\), \(a\ge2t\), \(A=a-t\), \(L=a-2t\), and
\[
 I=(x^2,y^a,z^a,y^tz^t-xf),\qquad n=4at-2t^2.             \tag{TP1}
\]
In the two arm modules write
\[
 f=y^tU+z^tV,\quad
 R_Y=K[y,z]/(y^A,z^t),\quad R_Z=K[y,z]/(y^t,z^A).
\]
Assume the axis polynomials are constant units. For each nonscalar arm,
let its first nonzero transverse order be \(d_Y,d_Z\in[1,t-1]\):
\[
 U=\alpha+z^{d_Y}F_Y(y)+\sum_{j>d_Y}z^jU_j(y),\qquad
 V=\beta+y^{d_Z}F_Z(z)+\sum_{i>d_Z}y^iV_i(z),             \tag{TP2}
\]
where \(\alpha,\beta,F_Y(0),F_Z(0)\ne0\). All higher coefficients
are arbitrary. No claim is made here for a first transverse polynomial
whose constant coefficient vanishes.

For a unit polynomial \(F(q)\) modulo \(q^m\), call \(r\in[0,m]\)
regular if
\[
 D_{m,r}(F)=\det[F_{r+i-j}]_{0\le i,j<m-r}\ne0,           \tag{TP3}
\]
with \(D_{m,m}=1\) and coefficients outside the range interpreted as zero.
The endpoints \(0,m\) are regular. Let \(h_m(F)\) be the longest run of
nonregular indices. Put
\[
 h_Y=\max\{h_A(F_Y),h_L(F_Y)\},\quad
 h_Z=\max\{h_A(F_Z),h_L(F_Z)\},\quad T=t+h_Y+h_Z,          \tag{TP4}
\]
omitting the length-zero module when \(L=0\).

**Theorem.** If both arms satisfy (TP2), then
\[
 D_x(I)=O(aT(1+\log(a+1))),\qquad
 g(I)=O(naT(1+\log(n+1))^2).                              \tag{TP5}
\]
If either arm is scalar, the same conclusions hold with \(T=t\), without
any condition on the other arm beyond its constant nonzero axis unit.
Retain the earlier box ceiling whenever smaller.

When all four sets of determinants in (TP4) are regular, this gives
\(g(I)=O(n^2\log^2(n+1))\), although the original axis gaps are
\(a-2t-1\). For fixed first orders, regularity is a nonempty open condition
over every infinite field and over sufficiently large finite extensions.
Indeed each determinant is a nonzero polynomial in the coefficients of
\(F\): the specialization \(F=q^r\) makes its matrix the identity.
Their product together with \(F_Y(0)F_Z(0)\) is nonzero in every
characteristic. The higher transverse coefficients are unrestricted.

An explicit common point in characteristic zero or characteristic
\(p>2A\) is
\[
                F_Y(q)=F_Z(q)=(1+q)^{A-1}.               \tag{TP6}
\]
For \(m\le A\), the classical binomial Toeplitz product gives
\[
 D_{m,r}=\prod_{j=0}^{m-r-1}
 \frac{(A-1+j)!\,j!}{(r+j)!\,(A-1-r+j)!}\ne0.
\]
See [Krattenthaler, Section 2.3](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11).
The factorial arguments are at most \(2A-2\). This is a classical
determinant evaluation, not a novelty claim.

## Uncertainty outside the scalar transverse tail

First, for every nonzero \(h\in K[q]/(q^m)\),
\[
             \deg h+\deg(Fh)\ge m-1-h_m(F).              \tag{TP7}
\]
If the degree sum were smaller, start at \(r_0=\deg(Fh)+1\) and
move to the next regular index \(r\le r_0+h_m(F)\). The source
\(h\) fits in degrees \(0,\ldots,m-r-1\), and its product is zero
in rows \(r,\ldots,m-1\). This contradicts (TP3).

Now work modulo \((q^m,v^t)\) with
\(W=\gamma+v^dF(q)+\sum_{j>d}v^jW_j(q)\).
For a nonzero \(H\) of transverse valuation \(r<t-d\), the coefficients
at orders \(r,r+d\) give
\[
 (WH)_r=\gamma H_r,\qquad
             F H_r=(WH)_{r+d}-\gamma H_{r+d}.
\]
Thus \(\deg H_r\le\min(\deg_qH,\deg_qWH)\) and
\(\deg(FH_r)\le\max(\deg_qH,\deg_qWH)\). Equation (TP7) implies
\[
             \deg_qH+\deg_q(WH)\ge m-1-h_m(F).            \tag{TP8}
\]
The corresponding maximum-weight sum is at least \(u(m-1-h_m(F))\)
when the long variable has weight \(u>0\). For \(r\ge t-d\),
\(WH=\gamma H\); no such degree uncertainty is asserted. These
exceptions are excluded by actual exchange ranks below.

## Long exchanges have cheap overlap rows

The [bridge partition](hook_bridge_reduction.md) places all but pure
cross-arm exchanges in an inventory of \(O(at)\) directions. Add the
directions \((1,b,c)\) with \(|b|,|c|<2a\) and \(|bc|\le64aT\), obtaining
\(O(aT(1+\log(a+1)))\) directions. Outside this inventory an exchange,
up to transposition, has labels
\[
 e=(0,Y,Z),\quad 2t\le Y<a,\quad t\le Z<2t,\qquad
 f'=(1,R,S),\quad0\le R<t,\quad t\le S<a,
\]
with \(b=Y-R, c=S-Z>0, bc>64aT\). Write
\[
 w=(\lambda,u,1),\quad\lambda=bu-c>0,\quad
 \tau=Yu+Z=\lambda+Ru+S,\quad d=(Y-t)u+Z-t,\quad Q=Ru+S.
                                                               \tag{TP9}
\]
Both \(b,c>64T\). The pure central \(Z\)-columns span the entire
upper tail \(z^{2t}C\), because \(V\) is a unit. Hence \(f'\) would
already be a loop modulo cheaper columns unless
\[
 \tau\le(2t-1)u+a-1,\qquad
                u<\frac{a}{b-2t}<\frac{2a}{b}.           \tag{TP10}
\]
Consequently \(c>32Tu\) and \(c>64T\). In particular \(Q\) is
strictly larger than both
\[
                 (2t-1)u+t-1,\qquad(t-1)u+2t-1.         \tag{TP11}
\]
All upper rows in the two common \(t\)-by-\(t\) boxes are therefore
strictly cheaper than the exchanged pair. No outer-head cost condition
or invertible axis-prefix basis has been used.

## Primal and dual uncertainty give four estimates

At a wall point avoiding other nonparallel label ties, a literal pair is
nonzero and parallel modulo cheaper columns. Its relation has the form
\[
       fH=\xi y^Rz^S+G,\quad \xi\ne0,                    \tag{TP12}
\]
where \(H\) has leading source \(y^{Y-t}z^{Z-t}\), all its other
sources are cheaper, and \(G\) uses cheaper upper rows. In the two arms,
with maximum planar weight denoted by \(\delta\),
\[
 \delta(H_Y)=d,\quad \delta(UH_Y)<Q-tu,\qquad
 \delta(H_Z)<d,\quad \delta(VH_Z)=Q-t.                    \tag{TP13}
\]
Neither source belongs to its scalar transverse tail. If \(H_Y\) did,
then \(UH_Y=\alpha H_Y\), and (TP13) would imply \(\lambda<t\).
But \(e\) would then be parallel to the strictly cheaper upper identity
\((1,Y,Z-t)\). If \(H_Z\) belonged to its scalar tail, (TP13) would
imply \(\lambda>tu\), and \(f'\) would be parallel to the cheaper
pure central column \((0,t+R,S)\). Both contradict eligibility.
Apply (TP8) at length \(A\). With \(F=(Y+R)u+S+Z\), this gives
\[
 F>(a+t-1-h_Y)u+t,\qquad
                       F>a+tu+t-1-h_Z.                  \tag{TP14}
\]

For the reverse estimates choose a linear functional \(\ell\) that
annihilates every cheaper column and is nonzero on both exchanged columns.
Write its row coefficients as \(\ell_Y(i,j),\ell_Z(i,j)\) in local
arm coordinates, and reverse the two rectangles:
\[
 K_Y=\sum\ell_Y(i,j)y^{A-1-i}z^{t-1-j},\qquad
 K_Z=\sum\ell_Z(i,j)y^{t-1-i}z^{A-1-j}.                   \tag{TP15}
\]
The cheap overlap rows in (TP11) make \(K_Y\bmod y^L=K_Y\) and
\(K_Z\bmod z^L=K_Z\). On pure central source columns the transpose
convolution is exactly multiplication by \(U,V\) followed by these
length-\(L\) truncations. Thus
\[
\begin{array}{ll}
 \delta(K_Y)<(a-1)u+t-1-Q,&
 \delta(UK_Y\bmod y^L)=(A-1)u+t-1-d,\\
 \delta(K_Z)=(t-1)u+a-1-Q,&
 \delta(VK_Z\bmod z^L)<(t-1)u+A-1-d.
\end{array}                                                \tag{TP16}
\]
The equality in the first product comes from \(\ell(e)\ne0\);
the second input equality comes from \(\ell(f')\ne0\).
All four vectors are nonzero, since the products are by units.

Neither reversed input is in its scalar transverse tail. For example,
that exception for \(K_Y\) would mean that \(\ell_Y\) is supported
only on transverse row indices \(j<d_Y\). There \(U\) acts as the
scalar \(\alpha\). Its overlap row coefficients vanish by (TP11).
The \(Y\)-part of \(\ell\) alone would then annihilate every cheaper
column and \(f'\), yet be nonzero on \(e\), contradicting their
parallelism. The \(Z\)-part gives the same contradiction with the roles
reversed. This argument is a separating-functional rank test, not an
assumption that a contraction is an algebra quotient.

Apply (TP8) at length \(L\) to (TP16). It yields
\[
 F<(a+2t-1+h_Y)u+3t-2,\qquad
                    F<a+(3t-2)u+2t-1+h_Z.               \tag{TP17}
\]
This replaces the unavailable cheap full-prefix bases.

## The resulting direction count

If \(a\le64T\), the entire difference box already has \(O(aT)\)
directions. Otherwise eliminate \(F\) from (TP14)--(TP17):
\[
 \frac{a-2t+1-h_Z}{a+t-1+h_Y}<u<
 \frac{a+t-1+h_Z}{a-2t+1-h_Y}.
\]
It follows that \(a|u-1|<4T\), \(a-T<F<a+6T\), and
\(|F-(b+c)|<11T\). In particular
\[
                         |b+c-a|<20T.                    \tag{TP18}
\]
There are \(O(aT)\) integer pairs in this strip. Add the controlled
inventory and apply the [capped-coordinate direction theorem](cross_layer_direction_budget.md)
with \(x^2\in I\), retaining all parallel pieces and subdivisions.
This proves (TP5).

If \(U=\alpha\) is scalar, every pure \(Y\)-row with local long index
at least \(t\) appears only in its upper identity and parallel central
column; bridge sources have long index less than \(t\). Those rows are
direct summands. All remaining labels have \(y<2t\), so their cross-layer
direction inventory is \(O(at)\). The parallel pairs add only the direction
\((1,0,-t)\). The transposed proof handles scalar \(V\).

## Same-layer walls and a common complete weight list

Let \(D_{\rm all}(I)\) count all distinct primitive positive wall
directions, including those with first coordinate zero. Then
\[
                 D_{\rm all}(I)=O(aT(1+\log(a+1))).       \tag{TP19}
\]
To prove the additional assertion, retain same-arm and bridge differences,
an inventory of \(O(at)\) directions, and the bounded box differences
with transverse product at most \(64aT\). A remaining same-layer
exchange has a pure \(Y\)-label \(e=(\epsilon,Y,Z)\) and a pure
\(Z\)-label \(f'=(\epsilon,R,S)\), with
\[
 b=Y-R>0,\quad c=S-Z>0,\quad bc>64aT,\quad u=c/b.
\]
For \(\epsilon=0\), both labels are central sources, with
\(2t\le Y,S<a\) and \(t\le R,Z<2t\). For \(\epsilon=1\), both
are upper identities, with \(2t\le Y,S<a\) and \(0\le R,Z<t\);
the product cutoff ensures these long-coordinate inequalities.
Write \(\tau=w\cdot e=w\cdot f'\), \(Q=\tau-\lambda\),
\(d=\tau-t(u+1)\), and \(F=\tau+Q\).

For two sources, the primal relation has both source-arm maxima equal
to \(d\) and both output-arm maxima strictly below \(Q-tu,Q-t\).
For two upper identities, both source-arm maxima are strictly below
\(d\), and both output maxima equal \(Q-tu,Q-t\). In either case
(TP14) follows. The scalar-tail exceptions would make one exchanged
column parallel to a strictly cheaper column on its own arm, just as
in (TP13).

If \(Q\) exceeds the two overlap costs in (TP11), reflect a separating
functional as in (TP15). For two sources the two input inequalities
are strict and both product inequalities are equalities. For two
upper identities the inputs are equalities and the products are strict.
The same degree sums give (TP17). The scalar-tail separating-functional
contradiction applies separately on each arm. Hence \(a|u-1|<4T\).

For two upper identities, the large-product condition already ensures
this cheap-overlap case: \(Q=Yu+Z=Ru+S\), while
\(c>64T\) and \(c=bu>64Tu\).
For two sources there remains the case
\[
 Q\le\max\{(2t-1)u+t-1,(t-1)u+2t-1\}<2t(u+1).
\]
The original source boxes give
\(\tau<\min(au+2t,a+2tu)\), so
\[
 F<(a+2t)u+4t,\qquad F<a+4tu+2t.
\]
Combine these with (TP14). When \(a>64T\),
\[
 \frac{a-3t-1-h_Z}{a+t}<u<
 \frac{a+t}{a-3t-1-h_Y},\qquad a|u-1|<6T.
\]
In every case \(|b-c|<6T\). Dividing by the common divisor preserves
that bound for the primitive direction. There are \(O(aT)\) such
integer pairs, which proves (TP19).

Fix \(t,a\) and an upper bound \(T_0\) on (TP4). All relevant wall
planes now lie in one coefficient-independent arrangement: the exact
same-arm/bridge inventories, the bounded hyperbolic-product inventories,
the cross-layer strips (TP18), and the same-layer strip
\(|b-c|<6T_0\). After primitive normalization it has
\(O(aT_0(1+\log(a+1)))\) planes. One positive integer weight from
each full-dimensional arrangement chamber is therefore a complete
test list for every ideal satisfying these conditions:
\[
 |\mathcal W_{t,a,T_0}|=O(a^2T_0^2(1+\log(a+1))^2),\qquad
             \|w\|_\infty=O(a^2).                        \tag{TP20}
\]
Indeed each chamber is inside one cone of every such fan. Its extreme
rays are cross products of arrangement normals and coordinate-boundary
normals. Original difference entries have absolute value less than
\(2a\), and their first coordinate is zero or one. Primitive ray
generators therefore have height \(O(a^2)\). The sum of three linearly
independent ray generators lies in the chamber interior and has the
same height order. A full-dimensional cone always has such three rays.
The elementary arrangement bound gives the cardinality.

In the transverse-regular case \(T_0=t\), (TP20) is a common complete
list of \(O(n^2\log^2(n+1))\) weights with \(O(\log a)\) bits per
coordinate. Evaluating exact monomial quotient ranks at that list
enumerates every positive initial ideal in the class. This is a
coefficient-independent compilation consequence; no optimal running
time for constructing the arrangement or evaluating the field ranks
is asserted. The original capped-coordinate bound in (TP5) remains
the sharper count comparison when the transverse gaps are large.

## Verification and limits

```bash
python3 -B scripts/run_hook_transverse_prefix_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checks reconstruct original commuting quotient columns, test all-vector
uncertainty by kernel dimensions including the exceptional scalar tails,
verify the transpose-convolution identity, and test both primal relations
and separating functionals on actual cross-layer and same-layer coherent
walls. Primitive integer interior weights check the stated height mechanism.
Scalar inequalities
are checked separately and are not presented as realizability evidence.
The finite same-layer field witnesses exercise the cheap-overlap branch;
the small-upper-threshold branch has scalar checks and the uniform proof,
without a claimed original-field witness in this check suite.

Toeplitz/Padé normality, the binomial determinant, annihilator separation,
and represented-matroid duality are classical inputs. The normality
background remains in [Jedynak--Gilewicz, Sections 1--2](https://doi.org/10.1155/2013/185648).
The linked direction and shadow proofs retain their sources. Exact priority
of this combined transverse-rank bound is unassessed. First transverse
polynomials with zero constant coefficient, large transverse gaps, arbitrary
outer images, and the unrestricted endpoint remain unresolved.
