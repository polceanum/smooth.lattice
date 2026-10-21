# Overlapping tower layers and coherent arithmetic dependencies

The [price-region note](cyclotomic_price_regions.md) leaves a concrete gap:
monomials of high tag degree can enter before the low-degree tower blocks
are complete. This note gives an exact reduction and a quantitative budget
for that gap. It also distinguishes nonregular column dependencies from
arithmetic changes that actually occur in coherent initial ideals.

These are supporting structural results. No count of the intervening
coprime-period fan, polynomial amplification, new unrestricted bound, or
publication significance is established. The original extremal goal remains
active. In particular, the equal-period example below must not be substituted
for the cyclotomic field of three distinct primes.

## The low-degree basis remains an exact background at every weight

Use the general tower from the earlier note:

\[
 F=K[x]/(f),\quad \deg f=D,\qquad
 A=F[T]/(g),\quad \deg_T g=R,\qquad n=DR,
\]

where both polynomials are monic and \(u,v\in F\) are units. The three
polynomial generators are \(X=x,Y=uT,Z=vT\). The algebra need not be a
field or reduced. Let positive weights be \((\alpha,\beta,\gamma)\), and
call \(b+c\) the tag degree of \(X^aY^bZ^c\).

All standard monomials have \(a<D\), by the monic relation for \(x\).
Restrict the column matroid to

\[
 L=\{(a,b,c):0\le a<D,\ b+c<R\}.
\]

It has rank \(n\). Its minimum basis \(C_\rho\), at a strict weight,
is the union of the independent blocks from the earlier slot formula:
for every \(0\le j<R\), select a minimum basis among
\(x^a(v/u)^c\), \(a<D,c\le j\), at costs \(a+\rho c\), where
\(\rho=(\gamma-\beta)/\alpha\). This statement needs no separation
between the costs of successive blocks. The blocks are a matroid direct sum,
so their internal cost orders suffice even when those orders interleave.

Let \(B\) be the unrestricted minimum monomial basis. Then, exactly,

\[
 B=(C_\rho\setminus Q)\cup H,\qquad
 H=\{a\in B:\operatorname{tag}(a)\ge R\},\quad |Q|=|H|.       \tag{O1}
\]

Indeed a low-ground column rejected by greedy on \(L\) is in the span of
earlier low-ground columns. It remains in the span of earlier columns when
the ground is enlarged. Thus \(B\cap L\subseteq C_\rho\); the equal
cardinalities give (O1). This is ordinary minimum-basis sensitivity, not a new
matroid exchange theorem. It supplies a useful exact reduction: high-degree
columns can delete from the known low basis, but cannot introduce a previously
rejected low column.

There is an immediate whole-family corollary. If \(g(T)=T^R\), every
high-tag column is zero, so \(C_\rho\) describes the **entire** positive
fan at every price. The shared slot breakpoints are
\(\rho=j/h\), \(|j|\le D-1\), \(1\le h\le R-1\).
For \(D,R\ge2\) this gives at most \(2(D-1)(R-1)+2<2n\)
initial ideals globally, rather than only in the old cheap-coordinate
regions. The one-dimensional cases have at most two. Thus a coupled
extension polynomial, with surviving high-tag columns, is essential to
any amplification within this tower framework.

## An excess-degree budget controls the unfinished part

Put \(m=\min(\beta,\gamma)\), \(\Delta=|\beta-\gamma|\).
Let \(\mu_x\) be the sum of the \(X\) exponents in \(B\), and let
\(\mu_e\) sum the exponents of the more expensive tag variable. Define

\[
 J_0=\frac{n(R-1)}2,\qquad
 E=\sum_{a\in B}\operatorname{tag}(a)-J_0.
\]

The degree-at-most-\(j\) columns have rank \(D(j+1)\) for \(j<R\).
They span precisely \(F\oplus FT\oplus\cdots\oplus FT^j\).
Consequently, if the tag degrees of any basis are sorted as
\(d_0\le\cdots\le d_{n-1}\), then \(d_i\ge\lfloor i/D\rfloor\).
In particular \(E\ge0\).

Compare \(B\) with the rectangular basis consisting of \(X^a\) times
the \(j\)-th power of the cheaper tag variable, \(a<D,j<R\).
Its tag moment is \(J_0\), its \(X\)-moment is \(n(D-1)/2\), and
its expensive-tag moment is zero. Optimality gives the uniform inequality

\[
 \boxed{\alpha\mu_x+mE+\Delta\mu_e
           \le\frac{\alpha n(D-1)}2.}                    \tag{O2}
\]

No generic-coefficient, separation-of-layers or characteristic hypothesis is
used. This retains information lost by a total-degree cutoff alone.

If \(M=|H|=Dq+s\), \(0\le s<D\), the \(M\) high-degree columns
occupy the last \(M\) positions in the sorted degree list and each has
\(d_i\ge R\). Therefore

\[
 E\ge\frac{Dq(q+1)}2+s(q+1)\ge\frac{M^2}{2D},\qquad
 M^2\le\frac{Dn(D-1)\alpha}{m}.                          \tag{O3}
\]

The reference rectangle already spans at cost at most
\((R-1)m+(D-1)\alpha\). Every chosen column has cost at most this,
so all selected tag degrees satisfy

\[
 \operatorname{tag}(a)\le R-1+
       \left\lfloor\frac{(D-1)\alpha}{m}\right\rfloor.   \tag{O4}
\]

Equations (O1)--(O4) turn the remaining problem into a controlled extension
of a known basis, with explicit numbers of affected columns and possible
extra layers. They do not count the coherent extensions. For balanced
periods of size \(q\), \(D\) has order \(q^2\), \(R\) order \(q\),
and prices \((1,q,q)\) still allow a positive fraction of high columns
in (O3). Thus the central unresolved scaling has not been assumed away.

## Where characteristic dependence can occur

For the pairwise-coprime circuit algebra

\[
 A_{r,s,t}=K[U,V,W]/(S_r(U),S_s(V),S_t(W)),\qquad
 (X,Y,Z)=(UV,UW,VW),
\]

every coherent basis in any of the three cheap-coordinate regions of the
previous theorem has integer index one in the natural integral algebra.
For example, when \(X\) is cheap, the completed blocks are monomial
bases of the two-circuit algebra \(\mathbb Z[U,V]/(S_r,S_s)\), multiplied
by units and by the independent powers of \(W\). The established
unimodularity theorem for two circuit factors makes every block determinant
\(\pm1\). Multiplication by a unit is an integral automorphism. Their
product is again \(\pm1\). The identical argument applies to the other
axes. Hence the regional standard sets agree in every characteristic.
The new conclusion concerns the specified coherent blocks; the two-circuit
unimodularity theorem itself is prior work.

More generally, fix a finite algebra free over \(\mathbb Z\), its integral
reference basis and a strict positive weight. If its characteristic-zero
standard monomials have determinant \(d\ne0\), their standard set survives
modulo a prime \(p\) exactly when \(p\nmid d\). For that implication, their
whole matrix is invertible over \(\mathbb Z[1/d]\). Each rationally rejected
column has coordinates only in earlier chosen columns, and those coordinates
belong to this localization. Reduction modulo \(p\nmid d\) preserves that
containment and preserves independence of every chosen prefix. If \(p\mid d\),
the proposed full basis loses rank and cannot remain the standard set.

This is the finite free-algebra version of the established lucky-prime and
base-change principle, not a new general modular algorithm. It explains why
a nonregular minor among arbitrary columns does not establish a change in an
initial fan: its determinant must occur in a **coherent** monomial basis.

## A coherent defect exists, but need not amplify a fan

Here is a complete small control, and an infinite extension with a uniform
proof. For equal periods three let

\[
 A_1=\mathbb Z[U,V,W]/(U^2+U+1,V^2+V+1,W^2+W+1),
 \quad (X,Y,Z)=(UV,UW,VW).
\]

The pair-product map is invertible on the period-three labels, so it gives
an actual length-eight quotient in every field. Its characteristic-zero
standard set at \(w=(269002,209961,210033)\) is

\[
 D=\{1,X,X^2,Y,Y^2,Z,Z^2,YZ\}.
\]

Its matrix in the tensor power basis has determinant \(2\). The saved
strict border polynomials certify coherence. In characteristic two it has
rank seven, and the same weight selects a different standard set. Exhausting
all 160 length-eight downsets and clipping their exact border cones gives
nine initials in characteristic zero, six in characteristic two. Six of the
rational cones have index one; three have index two. Thus every odd
characteristic has the same nine-state positive fan, by the determinant
criterion, including the nonreduced characteristic-three algebra.

For every odd integer \(k\ge1\), replace the three defining polynomials
by \(S_3(U^k),S_3(V^k),S_3(W^k)\), and retain the pair-product generators.
The length is \(8k^3\). If \(I_1\) denotes the relation ideal above, then

\[
 I_k=I_1(X^k,Y^k,Z^k).                                   \tag{O5}
\]

To verify this, the natural map from the power-substitution quotient onto
the displayed algebra is surjective: the variables are units of period
\(3k\), and \(U^2=XY/Z\), with two invertible modulo \(3k\).
Both algebras are free over \(\mathbb Z\) of rank \(8k^3\); hence the
surjection is an isomorphism. A Gröbner basis remains a Gröbner basis under
this power substitution, since least common multiples and the ordinary
S-pair reductions commute with it. Its standard set inflates to
\(kD+\{0,\ldots,k-1\}^3\); its positive weight cone is unchanged.

The period-\(k\) grading splits the integral algebra into \(k^3\) rank-eight
pieces. The degree map on remainder triples has determinant \(-2\), hence
is invertible modulo odd \(k\). Multiplication by each remainder monomial
is an integral unit, so a lifted basis has the same index as \(D\) in
each piece. Therefore the three index-two cones have cokernel

\[
 (\mathbb Z/2)^{k^3},\qquad \text{index }2^{k^3},           \tag{O6}
\]

while the other six cones remain unimodular. The full fan still has nine
states outside characteristic two and six in characteristic two, for every
odd \(k\). Large coherent determinant indices therefore coexist with a
constant fan. The substituted reduced-basis coefficients also remain the
same; the index growth is not a coefficient-bit lower bound.

This family has repeated periods. It establishes neither a characteristic
change nor a large fan for three distinct coprime primes. It confirms a real
coherent arithmetic mechanism while preventing its size from being confused
with the number of orderings.

## Prior work and the next proof task

The circuit/simplicial dictionary is established by
[Martin--Reiner, Theorem 1 and its proof](https://arxiv.org/pdf/math/0402206).
[Hidaka--Itoh, Theorem 1.2 and Section 5.2](https://arxiv.org/html/2403.10817v3#S5.SS2)
prove the two-circuit unimodularity input.
[Fu--Ma--Wang, Theorem 4.7 and Section 5.5](https://arxiv.org/html/2609.29227v1#S5.SS5)
prove nonregularity of three nonfree circuit factors by explicit minors.
These are prior results; their abstract bases and minors do not assert
coherence for this monomial presentation. Only the cited passages were
inspected, not the complete subsequent literature.
[Arnold, Section 5](https://www3.risc.jku.at/research/theorema/Groebner-Bases-Bibliography/gbbib_files/publication_355.pdf)
is a primary comparison for the established modular-specialization framework;
the reading scope is recorded with the checks.
Arnold fixes homogeneous ideals for that comparison; the finite free-column
argument above is supplied explicitly rather than assuming his hypotheses
cover the present inhomogeneous presentation.

The next requirement is a uniform bound or a growing coherent exchange
mechanism for the high-column extension (O1), using the actual circuit
relations and budget (O2). Counting all bases of a nonregular matroid,
extrapolating the small fan, or inflating a fixed coherent defect does not
meet the main goal. Current priority of the additional tower application
remains unestablished. No new manuscript is started on this evidence alone.

Run from the repository root:

    python3 -B scripts/run_cyclotomic_overlap_checks.py
    ./build.sh
    python3 tests/test_smoke.py

The checks compare the low-block rule with independently reduced full columns,
verify (O1)--(O4) both inside and outside the earlier region, exhaust all small
fan candidates with exact rational cone coverage, and compute the lifted
integer blocks by ordinary monic polynomial division. Finite checks support
the uniform algebraic arguments and have their scope recorded separately.
