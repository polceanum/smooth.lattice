# Complex basis transforms retain the fan but do not preserve stability

There is a dependency-aware complex transform for every finite quotient.
Its limiting support function recovers the whole positive initial fan.
However, even a length-three complex quotient loses the half-plane and
Lorentzian properties when separate monomial variables are combined into
exponent coordinates. At that same length, determinant-transform zeros
can have any prescribed real part while the initial fan stays fixed.
In characteristic two, an actual length-thirteen quotient has a Fano
minor, preventing any stable polynomial with its full monomial-basis
support, even before combining variables.

These are feasibility guards for an analytic route. The determinant
construction, tropical limit and stability obstructions use classical
results; no new extremal exponent or novelty is asserted. They do not
exclude other complex analytic methods. A bound on simultaneous walls
that uses commuting multiplication remains necessary.

## The transform that retains dependencies

Let \(A=K[x,y,z]/I\) have length \(n>0\). Use the universal ground
\[
 E_n=\{a\in\mathbb N_0^3:(a_x+1)(a_y+1)(a_z+1)\le n\},
\]
and its residue-column matroid \(M_I\), of rank \(n\). For each basis
\(B\), write \(\mu(B)=\sum_{a\in B}a\). Over every field, the
ordinary real polynomial
\[
 b_I(u)=\sum_{B\in\mathcal B(M_I)}\prod_{a\in B}u_a             \tag{C1}
\]
retains exactly which monomial sets are vector bases. Its coefficients
are ordinary real ones, not residues in \(K\). In particular, a
characteristic-two quotient is not silently lifted to characteristic
zero. For any positive real coefficients \(c_B\), define
\[
 F_{I,c}(X,Y,Z)=\sum_B c_B X^{\mu_x(B)}Y^{\mu_y(B)}Z^{\mu_z(B)},
 \qquad
 Z_{I,c}(s;w)=F_{I,c}(e^{-s w_x},e^{-s w_y},e^{-s w_z}).          \tag{C2}
\]
The second expression is entire in \(s\) for each fixed real \(w\).
It requires no analytic continuation. If \(m(w)=\min_B w\cdot\mu(B)\),
then
\[
                 \lim_{t\to+\infty}-{1\over t}\log Z_{I,c}(t;w)
                              =m(w).                           \tag{C3}
\]
Indeed, factor out \(e^{-tm(w)}\) in the finite positive sum; the
remaining factor is bounded above and below by positive constants
independent of \(t\). For unit coefficients, more explicitly,
\[
 m(w)-{\log|\mathcal B(M_I)|\over t}
       \le -{\log Z_{I,1}(t;w)\over t}\le m(w).                 \tag{C4}
\]
The positive domains of strict minimization of \(m\) are the initial
cones by the established state-polytope dictionary. Thus the transform
retains the actual fan, but (C3) alone supplies no complexity bound.
The free smooth-number Euler product does not contain the rank data in
(C1).

When \(K=\mathbb C\), put the residue columns into a full-row-rank
matrix \(V\). The classical Cauchy--Binet construction gives
\[
 q_V(u)=\det(V\operatorname{diag}(u_a)V^*)
       =\sum_B |\det V_B|^2\prod_{a\in B}u_a.                   \tag{C5}
\]
It is stable in the separate variables: if every \(\operatorname{Im}u_a>0\),
then for nonzero \(v\),
\[
 \operatorname{Im}\bigl(v^*V\operatorname{diag}(u_a)V^*v\bigr)
                 =\sum_a\operatorname{Im}u_a\,|(V^*v)_a|^2>0.
\]
Hence the matrix is nonsingular. Changing the vector-space basis of
\(A\) multiplies (C5) by a positive constant, so it does not change
its zeros. Replacing \(u_a\) by \(X^{a_x}Y^{a_y}Z^{a_z}\) gives
(C2) with positive determinant weights. Such a monomial substitution
is not a nonnegative linear substitution of variables.

## An actual length-three substitution obstruction

For any \(c\in\mathbb C^*\), set
\[
                   A_c=\mathbb C[x,y,z]/(x^3,y-cx^2,z).          \tag{C6}
\]
The basis \(1,x,x^2\) proves its length is three. Its only monomial
vector bases are
\(\{1,x,x^2\}\) and \(\{1,x,y\}\). All other ground columns
are zero. Its positive fan therefore has exactly two cones, separated
by \(w_y=2w_x\), independently of \(c\).

In those vector coordinates, (C5) is
\[
       u_1u_x(u_{x^2}+|c|^2u_y).
\]
For \(c=1\), its exponent-coordinate image is
\[
                             F(X,Y,Z)=X^3+XY.                  \tag{C7}
\]
This vanishes at \(X=-1+i\), \(Y=2i\), and any upper-half-plane
\(Z\). Thus even an actual complex quotient does not preserve upper
half-plane stability under the required substitution. Likewise
\(X=1+2i\), \(Y=3-4i\), and any right-half-plane \(Z\) give a
zero with every coordinate in the open right half-plane.

Ordinary total-degree homogenization does not repair the Lorentzian
property: it gives \(H=X^3+XYT\). The Hessian of
\(\partial_XH=3X^2+YT\) is
\[
 \begin{pmatrix}6&0&0\\0&0&1\\0&1&0\end{pmatrix},
\]
with two positive eigenvalues. A Lorentzian cubic cannot have that
quadratic derivative. Equivalently, its two-point support
\(\{(3,0,0),(1,1,1)\}\) fails the symmetric exchange condition.
This refutes automatic inheritance, not every possible alternative
homogenization or an inequality formulated in the separate variables.

For \(w=(1,3,1)\), the determinant-weighted entire transform is
\[
 Z_{A_c}(s;w)=e^{-3s}+|c|^2e^{-4s}
            =e^{-3s}(1+|c|^2e^{-s}).                            \tag{C8}
\]
Its zeros are exactly
\[
                   s=\log|c|^2+(2j+1)\pi i,\qquad j\in\mathbb Z.
                                                                    \tag{C9}
\]
Consequently their real part can be any prescribed real number while
the length, exponent labels and full two-cone fan stay fixed. There
is no coefficient-independent vertical strip containing these zeros
for all quotients of a fixed length. This statement concerns the
determinant-weighted transform. Unit basis coefficients remove this
particular coefficient dependence, but not the obstruction in (C7).

## The all-field obstruction is present in an actual quotient

Work over \(\mathbb F_2\), write \(\mathfrak m=(x,y,z)\), and label
\[
 (u_1,\ldots,u_7)=(x^3,y^3,z^3,x^2y,xy^2,x^2z,xz^2).
\]
Take the quotient by
\[
 \mathfrak m^4+(y^2z,yz^2,xyz,
 u_4+u_1+u_2,
 u_5+u_1+u_3,
 u_6+u_2+u_3,
 u_7+u_1+u_2+u_3).                                           \tag{C10}
\]
All degree-three elements are in the socle, so these seven independent
relations form an ideal modulo \(\mathfrak m^4\). There are ten
independent monomials of degree below three and a three-dimensional
cubic quotient. Thus the length is thirteen. The seven surviving
cubic columns, in the basis \(u_1,u_2,u_3\), are
\[
 (1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1).
\]
They represent the Fano matroid, with 28 bases and seven dependent
triples. Every lower-degree column is a coloop, and all other ground
columns are loops. In particular the full basis support contracts to
the Fano basis support.

Brändén's Theorem 6.6 proves that no stable polynomial has Fano
support, allowing arbitrary nonzero coefficients. If a stable
polynomial with the full basis support of (C10) existed, differentiating
in its ten coloop variables and deleting loop variables would produce
one with Fano support, a contradiction. Thus even a weighted stable
encoding of all bases is unavailable uniformly over fields. This
does not claim that the projected fan alone determines the Fano
matroid or cannot admit another analytic representation.

For the unweighted Fano polynomial, an independent exact check also
gives a negative Rayleigh difference. Set
\((u_3,u_4,u_5,u_6,u_7)=(-3,-3,-3,1,2)\). Its remaining bivariate
polynomial is
\[
                     -3u_1u_2-9u_1+11u_2+9,
\]
and \(\Delta_{12}=(-9)(11)-(-3)(9)=-72\). The checker additionally
perturbs all five fixed variables into the upper half-plane and solves
exactly for an upper-half-plane zero, without relying only on this
necessary stability criterion.

## Reproduction and precise prior-art boundary

```bash
python3 -B scripts/run_complex_basis_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker constructs quotient relations and commuting multiplication,
enumerates original monomial bases, independently checks Cauchy--Binet
by polynomial determinant expansion, and uses rational complex arithmetic
for the substitution and Fano zeros. Numerical checks of (C9) are
separately labelled; its uniform statement follows from (C8).

The relevant primary inputs are:

- [Choe--Oxley--Sokal--Wagner, Theorem 8.1 and Corollary 8.2](https://arxiv.org/pdf/math/0202034):
  the complex determinant construction and its weighted basis support.
- [Brändén, Theorem 5.6 and Theorem 6.6](https://arxiv.org/pdf/math/0605678):
  the Rayleigh characterization and the impossibility of stable Fano support.
- [Brändén--Huh, Theorems 2.25 and 3.10](https://arxiv.org/pdf/1902.03719):
  Lorentzian support conditions and the Lorentzian basis polynomial of
  every matroid, including non-complex-representable ones.
- The state dictionary and its scope are retained in the
  [cross-layer direction reduction](cross_layer_direction_budget.md).
  The Fano realization lies inside the already bounded
  [socle-only class](disjoint_socle_circuits.md); it is not a new
  large-fan candidate.

These sources establish the classical boundary inspected here. They
do not constitute a strongest-current novelty audit for a new analytic
fan-count theorem. Such a theorem has not been proved. The current
unrestricted count, the required aggregate direction/incidence bound,
and the original goal remain unresolved.
