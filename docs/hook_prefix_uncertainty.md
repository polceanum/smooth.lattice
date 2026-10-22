# A near-quadratic fan bound for prefix-regular hooks

Let \(t\ge2\), \(a\ge2t\), and
\[
 M=(y^a,y^tz^t,z^a),\quad C=K[y,z]/M,\quad n=4at-2t^2.
\]
Consider every compatible flat square-zero extension with generator
images \(f_0,f_1,f_2\), over any field. Impose only the central
prefix-regular condition (CP3) in the
[prefix-minor note](hook_central_prefix.md). The other images and all
remaining central coefficients are arbitrary.

**Theorem.** Uniformly throughout this locus,
\[
 D_x(I_\phi)=O(n(1+\log(n+1))),\qquad
 g(I_\phi)=O(n^2(1+\log(n+1))^2).                         \tag{PU1}
\]
Here \(D_x\) counts distinct actual cross-layer wall directions and
\(g\) counts all monomial initial ideals for global orders. Actual
parallel pieces, recurrent eligibility and wall subdivisions are retained.
The locus is a nonempty open set in every characteristic, over infinite
fields or sufficiently large finite extensions, by the existing primitive
minor argument. Thus the bound applies to the generic compatible hook
family at growing thickness, including \(a=t^2+1\).

The proof adds a degree uncertainty estimate to the
[central-only wall reduction](hook_central_skeleton.md). It confines
every remaining large-product direction to a diagonal strip of width
\(O(t)\). This improves the previous \(O(n^{3/2}a)\) ceiling on the
stated open locus. It does not bound the prefix-degenerate strata or
unrestricted \(G_3\), and does not establish exact priority.

## Prefix minors forbid simultaneous low degrees

Put \(A=a-t\), \(L=a-2t\), and suppose \(a>2t\). Write
\[
 f_1=y^tU+z^tV,\quad
 R_Y=K[y,z]/(y^A,z^t),\quad R_Z=K[y,z]/(y^t,z^A).
\]
The two units and their row interpretation are those of (CS1a).
For a nonzero element of a rectangular monomial quotient, degrees
below always mean degrees of its unique monomial representative.

**Degree lemma.** For every nonzero \(H\in R_Y\),
\[
                    \deg_y H+\deg_y(UH)\ge L-1.            \tag{PU2}
\]
The transposed estimate holds for \(V\) in \(R_Z\).

Write the two degrees as \(D,E\). If either is at least \(L\),
the assertion is immediate. Otherwise multiplication by \(y^t\)
truncates neither representative. If \(D+E<L-1\), put
\(r=E+t+1\). Then \(t\le r<A\). The nonzero polynomial
\(G=y^tH\) has source support
\[
                  t\le y<a-r,\qquad0\le z<t.
\]
Multiplication by \(U\) sends it into rows with \(y<r\), since
\(\deg_y(UG)=E+t=r-1\). But the prefix minor \(\Delta_Y(r)\)
is exactly the invertible map from this source rectangle to rows
\(r\le y<A,\ 0\le z<t\). Its kernel cannot contain \(G\).
This contradiction proves (PU2), in every characteristic.

For the planar weight \((u,1)\), write \(\delta_u\) for maximum
monomial weight. Nonnegative exponents give
\[
 \delta_u(H)+\delta_u(UH)\ge u(L-1),\qquad
 \delta_u(J)+\delta_u(VJ)\ge L-1                         \tag{PU3}
\]
for nonzero \(H\in R_Y,J\in R_Z\). These are direct consequences
of the invertible prefix minors. They are not positivity, stability or
support-matching assertions.

The minors also have a useful elementary axis interpretation. Ordering
source and target rows by \(z\)-degree makes \(\Delta_Y(r)\)
block triangular with \(t\) identical diagonal blocks from \(U(y,0)\).
Thus it is the \(t\)-th power of that univariate Toeplitz determinant.
The other arm uses \(V(0,z)\). Arbitrary off-axis central coefficients
therefore preserve these minors when the two axis polynomials are fixed.
This is the usual truncated Toeplitz/Padé rank mechanism; the determinant
implication itself is an elementary linear-algebra input. For the classical
denominator systems, determinant tables and normality boundary, see
[Jedynak--Gilewicz, Sections 1--2](https://doi.org/10.1155/2013/185648).
That background is not a priority claim for the hook-family bound below.

## Apply the estimate to an actual exchange relation

Use the central-only extension \(I_\circ\) in (CS1). A direction
outside the controlled inventory \(\mathcal H\) of (CS7) is represented,
up to transposition, by a central lower/upper opposite-arm pair
\[
\begin{split}
 e&=(0,Y,Z),\quad2t\le Y<a,\quad t\le Z<2t,\\
 f&=(1,R,S),\quad0\le R<t,\quad t\le S<a,\\
 q&=(1,-b,c),\quad b=Y-R,\quad c=S-Z,\quad bc>16at.
\end{split}                                                \tag{PU4}
\]
In particular \(b,c>16t\), so \(a>16t\). Normalize a generic
wall point by
\[
 w=(\lambda,u,1),\quad\lambda=bu-c>0,\quad
 \tau=Yu+Z=\lambda+Ru+S<a\min(u,1).                        \tag{PU5}
\]
The strict upper bound is the already proved (CS6). Put
\[
 d=(Y-t)u+Z-t,\quad T=Ru+S,\quad
                      F=(Y+R)u+S+Z.                       \tag{PU6}
\]

The literal pair is nonzero and parallel modulo cheaper columns.
Consequently its rank relation can be written
\[
                       f_1H=\alpha y^Rz^S+G,
                    \qquad\alpha\ne0,                    \tag{PU7}
\]
where \(G\) is supported on upper monomials of planar weight less
than \(T\). The source \(H\) has leading weighted monomial
\(y^{Y-t}z^{Z-t}\); every other source monomial has weight less
than \(d\). This follows by moving the strictly cheaper central
columns to the left and cheaper upper identity columns to the right.
Outer columns are zero in \(I_\circ\).

Let \(H_Y,H_Z\) be the restrictions of this source to the two
rectangles in (CS1a). The leading source lies only on the \(Y\) arm,
so
\[
 \delta_u(H_Y)=d,\quad \delta_u(UH_Y)<T-tu,\qquad
 \delta_u(H_Z)<d,\quad \delta_u(VH_Z)=T-t.                 \tag{PU8}
\]
Both restrictions are nonzero: the first contains the leading source,
and the second has a nonzero image containing \(y^Rz^{S-t}\).
The units make both products nonzero as well. Apply (PU3) to get
\[
             F>(a-1)u+t,\qquad F>a+tu-1.                  \tag{PU9}
\]
All ranks here are ranks of original residue vectors. A relation is
not inferred from a support matching.

## Two cheap prefix bases give the reverse inequalities

We next prove
\[
 \begin{split}
 F&\le au+(2t-1)u+3t-2,\\
 F&\le a+(3t-2)u+2t-1.                                  \tag{PU10}
 \end{split}
\]

For the first, put \(\beta=2t-1-Z\in[0,t-1]\) and
\[
                r=a+t-Y+\lfloor\beta/u\rfloor.
\]
Since \(u>c/b\), \(c>16t\), and \(b\le Y\), we have
\(\beta/u<Y/16<Y-2t\). Thus \(t\le r<A\). Every plain
column in \(Q_Y(r)\) is strictly cheaper than \(e\), because its
maximum-cost gap is
\((\lfloor\beta/u\rfloor+1)u-\beta>0\). For the upper prefix,
\[
 \tau-\max_{P_Y(r)}w\cdot v
 \ge F-au-(2t-1)u-3t+2.
\]
If the first inequality in (PU10) failed, the whole actual prefix
basis would be strictly cheaper and would span \(e\), a contradiction.

For the second, put
\[
 r=\lfloor a+t-1-Z-(Y-2t+1)u\rfloor+1.
\]
The bound \(\tau<a\) gives \(r\ge t\). Also
\((Y-2t+1)u>14t\), using \(Y>16t\) and \(u>c/b\).
Together with \(Z\ge t\), this implies \(r<A\). Every column of
\(Q_Z(r)\) is strictly cheaper than \(f\), by the choice of \(r\).
For its upper prefix,
\[
 \tau-\max_{P_Z(r)}w\cdot v
 \ge F-a-(3t-2)u-2t+1.
\]
Failure of the second inequality in (PU10) would make this entire
actual basis cheaper than \(f\), again a contradiction.

## Remaining directions lie in a narrow diagonal strip

Combine the first lower inequality of (PU9) with the second upper
inequality of (PU10), and then the other two. They give
\[
 \rho^{-1}<u<\rho,\qquad
                   \rho=\frac{a+t-1}{a-3t+1}<\frac{17}{13}.
                                                               \tag{PU11}
\]
Since \(a>16t\), these bounds imply
\(a|u-1|<5t\). Moreover (PU9)--(PU10) give
\(a-1<F<a+6t\). Finally
\[
 F-(b+c)=b(u-1)+2Ru+2Z,
\]
whose absolute value is less than
\(5t+(34/13)t+4t<12t\). Hence every remaining actual direction satisfies
\[
                      |b+c-a|<20t.                        \tag{PU12}
\]
This sufficient constant is not optimized. It includes all recurrent
wall intervals; it is a restriction on their distinct directions.

There are \(O(at)\) positive integer pairs in (PU12), since
\(b,c<a\) and each \(b\) allows at most \(40t+1\) values of \(c\).
The two orientations have the same bound. The controlled inventory
\(\mathcal H\) has \(O(at(1+\log(a+1)))\) directions by (CS9).
Thus \(D_x(I_\circ)=O(n(1+\log(n+1)))\). The direction transfer
(CS10) gives the same bound for every \(I_\phi\) on this locus.
The capped-coordinate direction theorem, which counts all parallel
pieces, now proves (PU1). No assumption about their subdivision was added.

## Collective wall-cost consequences and limits

Let \(W_1(I)=\sum_e\|q_e\|_1\) sum primitive one-norm length over
actual undirected wall pieces. Original hook labels have transverse
coordinates smaller than \(2a\), so every primitive direction has
one-norm \(O(a)\). The adjacency graph is planar and has \(O(g)\)
edges. Therefore
\[
 W_1(I_\phi)=O(a n^2(1+\log(n+1))^2)
                          =O(n^3(1+\log(n+1))^2).         \tag{PU13}
\]
For the divisor threshold
\(\nu(q)=\prod_i(q_i^++1)+\prod_i(q_i^-+1)-2\), the uncontrolled
mixed-sign cross-layer directions have \(\nu(q)\le2\|q\|_1\).
The other sign classes already have the cubic polylogarithmic budget
in the [bounded-coordinate wall theorem](bounded_coordinate_wall_shadows.md).
Consequently the full threshold-weighted wall sum on this locus also
has \(\sum_e\nu(q_e)=O(n^3(1+\log(n+1))^2)\).

For \(a=t^2+1\), (PU1) improves the earlier \(n^{13/6}\) box ceiling
to \(O(n^2\log^2(n+1))\) on a nonempty open locus in every
characteristic. A polynomial separation from this hook family must
therefore lie in prefix-degenerate central-image strata. Their fans
need not be bounded by those of a generic specialization; no such
monotonicity is assumed. Arbitrary degenerate central images and the
original unrestricted sharp endpoint remain open.

## Verification and priority boundary

```bash
python3 -B scripts/run_hook_prefix_uncertainty_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker verifies the degree lemma for all vectors in each tested
degree-bounded subspace by independent full-rank tail tests, with random
off-axis coefficients. It tests the scalar enclosure argument with exact
rational weights, reconstructs original quotient columns, and retains a
unit-only excluded-stratum control where the degree lemma fails. These
finite checks support the uniform proof, not an extrapolated growth law.

Prefix minors, Toeplitz/Padé rank conditions, the state-polytope dictionary,
planarity and the linked shadow and sign-class budgets are classical or
previously established inputs. Exact priority of the combined hook-family
bound is unassessed. This conditional family theorem, and its wall-cost
consequence, do not satisfy the original unrestricted extremal goal or
replace its required strongest-current comparison and final paper.
