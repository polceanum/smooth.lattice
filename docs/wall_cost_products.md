# A logarithm is necessary in the cubic wall-cost budget

The proposed collective estimate is
\(V(I)=O(n^3\operatorname{polylog}n)\), with repeated wall pieces
counted separately as in the [wall-cost reduction](free_wall_threshold.md).
The logarithmic allowance cannot be removed. There is an explicit
infinite family of radical rational point ideals, of lengths \(N=3p\)
for primes \(p\), such that
\[
             g(H_p)=\Theta(N^2\log N),\qquad
             V(H_p)=\Theta(N^3\log N).                    \tag{1}
\]
The walls producing this lower bound have cost exactly \(N/3\).
Thus they lie in the intermediate range left by the
[compulsory-core estimates](wall_coexistence.md).

This is a quantitative consequence of the existing common-coordinate
fan product and the published generic planar count. It does not give
a polynomial separation from the generic three-variable scale, improve
the unrestricted upper bound, or establish novelty. The polylogarithmic
cubic conjecture remains open.

## A sharp planar cost bound

Let \(J\subset K[x,y]\) have positive colength \(m\), over any field.
Order its positive initial cones by increasing ratio \(w_x/w_y\).
For adjacent standard sets let their profile difference, oriented in
that order, be
\[
                 (-k_e a_e,k_e b_e),\qquad
                 k_e,a_e,b_e\in\mathbb Z_{>0},\quad
                 \gcd(a_e,b_e)=1.
\]
The primitive wall cost is \(a_e+b_e\), the two-variable instance of
the exact divisor threshold. Define
\[
                      W(J)=\sum_e(a_e+b_e).              \tag{2}
\]
Then
\[
                         W(J)\le m(m-1).                \tag{3}
\]

Indeed, the \(x\)-moment decreases and the \(y\)-moment increases at
every change. Comparing the optimality inequalities on the two sides
proves these strict signs; equality in either would force equality
of both profiles. Hence the full weighted jump telescopes:
\[
 \sum_e k_e(a_e+b_e)
   =X_{\mathrm{first}}-X_{\mathrm{last}}
      +Y_{\mathrm{last}}-Y_{\mathrm{first}}
   \le 2\binom m2 .
\]
Each coordinate moment of an \(m\)-element downset lies between zero
and \(\binom m2\). Since \(k_e\ge1\), (3) follows. This also covers
a single cone, where the sum is empty. It counts actual profile edges,
not all possible two-variable directions.

The inequality is sharp at every prime length. More explicitly, over
\(\mathbb Q\), put
\[
 f_m(x)=\prod_{j=1}^m(x-j),\qquad
 T_m=(f_m(x),xy-1).
\]
Every initial standard set is one of the \(m\) hooks
\[
 D_k=\{1,x,\ldots,x^k,y,\ldots,y^{m-1-k}\},
                         \qquad 0\le k<m.               \tag{4}
\]
The relation \(xy-1\) excludes all monomials with both exponents
positive. Each hook is independent: evaluation at the positive nodes
gives consecutive Laurent powers, a Vandermonde matrix after
multiplying rows by a nonzero power of the node.

All hooks are coherent. Their profile is
\((k(k+1)/2,(m-1-k)(m-k)/2)\), and the successive differences are
\((j,-(m-j))\), \(1\le j<m\). The slopes
\((m-j)/j\) are strictly decreasing, so every hook uniquely minimizes
on an open interval of positive prices. There are no other downsets
of size \(m\) avoiding \(xy\), so this is the entire fan. Consequently
\[
 W(T_m)=\sum_{j=1}^{m-1}\frac{m}{\gcd(j,m)}
       =\sum_{\substack{d\mid m\\d>1}}d\varphi(d).
                                                               \tag{5}
\]
For prime \(m=p\), this is \(p(p-1)\), attaining (3).

## The exact wall cost of a common-coordinate fan product

Use the [established product construction](extremal_initial_ideals.md#an-exact-product-construction-reduces-a-lower-bound-to-the-plane).
For planar ideals \(J_1\subset K[x,y]\), \(J_2\subset K[x,z]\) of
lengths \(m_1,m_2\), and a common monic \(f(x)\) of degree \(s\),
put
\[
 H=(f(x),yz,yJ_1,zJ_2),\qquad N=s+m_1+m_2 .
                                                               \tag{6}
\]
Its fan is the Cartesian product of the two interval fans in the
independent ratios \(w_y/w_x,w_z/w_x\). A wall from the first factor
has one distinct piece for every cone of the second factor, and
conversely. These are all adjacencies. Translating each branch's
standard set by its tag changes its profile by a fixed vector and
preserves its edge differences. The two families of primitive
directions have different zero coordinates, so do not coincide.
Thus, exactly,
\[
 V(H)=g(J_2)W(J_1)+g(J_1)W(J_2).                         \tag{7}
\]
Together with (3) and the all-fields planar bound
\(g(J)=O(m\log(m+1))\), this proves
\[
 V(H)=O\!\left(m_1m_2\{m_1\log(m_2+1)
                              +m_2\log(m_1+1)\}\right)
     =O(N^3\log(N+1)).                                  \tag{8}
\]
This controls the entire common-coordinate product class, including
nonreduced and dependent planar branches. No bound for arbitrary
three-variable ideals follows.

## An explicit radical family attaining the logarithmic order

In (6) take \(s=m_1=m_2=m\), \(J_1=T_m\), and
\[
                 J_2=Q_m=(f_m(x),z-x^m).
 \qquad
 H_m=(f_m(x),yz,xy^2-y,z^2-x^mz).                       \tag{9}
\]
The ideal \(H_m\) is the radical ideal of the \(3m\) distinct rational
points
\[
 (j,0,0),\quad (j,j^{-1},0),\quad (j,0,j^m),
                             \qquad 1\le j\le m.        \tag{10}
\]
The displayed generators vanish there. The product construction
gives length \(3m\), equal to the number of distinct points, proving
equality with their vanishing ideal and radicality.

On the planar universal ground for length \(m\), the exponents
\(i+mj\) are distinct: \(0\le i\le m-1\). Every maximal evaluation
minor of \(Q_m\) is nonzero by generalized Vandermonde independence
at positive nodes. Thus it has exactly the generic corner-cut fan,
of size \(C_2(m)=\Theta(m\log(m+1))\). Only finitely many monomials
are needed for this conclusion; no uniformity on the whole
nonnegative lattice is asserted.

Equations (5) and (7) now give, for prime \(p\),
\[
 \begin{split}
 g(H_p)&=pC_2(p),\\
 V(H_p)&=p(p-1)C_2(p)+pW(Q_p),\\
 0&\le W(Q_p)\le p(p-1).
 \end{split}                                             \tag{11}
\]
This proves (1). The \((p-1)C_2(p)\) first-branch wall pieces
all have primitive directions \((j,-(p-j),0)\) and cost \(p=N/3\).
Their compulsory-core residual rank is \(2p+1\). This retains the
actual multiplicities and places the logarithmic obstruction in the
uncontrolled middle band. The family has fewer initials than the
generic order \(\Theta(N^2\log^2 N)\); a large aggregate cost does
not itself imply a fan-count separation.

Any proposed uniform estimate \(V(I)=O(N^3\log^C(N+1))\) therefore
needs \(C\ge1\). This requirement leaves the conditional power
\(12/5\) from the wall-cost reduction unchanged.

## Inputs, verification and limits

The monotone profile chain and its telescoping length are classical
planar convex geometry. The hook calculation is ordinary Laurent
interpolation. The product and its all-order initial formula were
already proved in this repository before this cost calculation.
The matching generic planar lower order is due to
Corteel--Rémond--Schaeffer--Thomas, as restated in
[Elizalde--Galván, Theorem 2.3](https://arxiv.org/html/2312.16353v1).
The standard initial-set/profile dictionary is
[Babson--Onn--Thomas, Theorem 2.4](https://arxiv.org/pdf/math/0207135).
These inputs do not establish priority of this particular
wall-cost observation. No current-novelty claim is made.

Run from the repository root:

    python3 -B scripts/run_wall_cost_product_checks.py
    ./build.sh
    python3 tests/test_smoke.py

The checker independently sweeps free planar cuts, enumerates actual
planar quotient bases, and compares their complete cones. It then
uses evaluation columns of the distinct spatial points (10),
recovers full border relations in each product cone, and verifies
the entire positive fan by exact polygon equality and total area.
Its adjacency costs are extracted independently from the complete
spatial profiles and compared with (7), including repeated
directions. Small rational and explicitly verified finite-field
specializations are checked; the asymptotic rational statement is
proved above, not fitted to those counts.
