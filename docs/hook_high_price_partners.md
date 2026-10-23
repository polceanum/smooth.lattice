# A single hook column has linearly many high-price partners

The clipped-price partner estimate cannot extend unchanged to high prices.
An explicit constant-axis central hook of thickness two has one lower
column with \(\Theta(n)\) distinct genuine upper partners in the other arm.
All four separate pure-list orders stay fixed. Every displayed exchange
has price \(\Theta(n)\), and both transverse displacements have order \(n\).
Thus neither a thickness-only partner cap nor a fixed multiple of the
scalar price range controls these walls.

This is a sharp one-column obstruction, not a polynomial fan separation.
The [passive-band theorem](hook_passive_band_fans.md) already bounds this
entire family by \(O(n^2)\). The unrestricted endpoint, collective counting
and exact current novelty remain open.

## Family and the exchange intervals

Let \(A\ge20\), \(a=A+2\), and work over \(\mathbb Q\), or a field of
characteristic \(p>2A\). Put
\[
 I_A=(x^2,y^{A+2},z^{A+2},y^2z^2-xf_A),\qquad
 f_A=y^2(1+y^{A-2}z)+z^2(1+y(1+z)^{A-1}).                 \tag{HP1}
\]
The two outer images vanish. The compatible hook has length
\(n=8A+8\), with arm operators, in local coordinates,
\[
 W_Y=1+q^{A-2}v,\qquad W_Z=1+v(1+q)^{A-1}
 \quad\text{in }K[q,v]/(q^A,v^2).                        \tag{HP2}
\]
Compatibility and length follow from the
[original hook construction](hook_bridge_reduction.md): each central
image is in its allowed arm, and multiplication by either outer syzygy
vanishes in the plain hook quotient.

Fix the original lower label \(E=(0,A,3)\). For each integer
\[
 \left\lceil3A/5\right\rceil\le K\le\left\lfloor3A/4\right\rfloor,
 \qquad h=A-1-K,
\]
use the nonempty open interval
\[
 J_K=\left((K-1)/(A-3),\ K/(A-3)\right).
\]
At every generic \(u\in J_K\), the pair
\[
 E=(0,A,3),\qquad F_h=(1,0,h+2)                           \tag{HP3}
\]
is nonloop parallel modulo all strictly cheaper original columns at
\[
 w=(\lambda,u,1),\qquad
 \lambda=Au+1-h,\qquad \tau=w\cdot E=Au+3.              \tag{HP4}
\]
Here generic excludes the finite additional comparison equalities;
the rank reduction below actually holds throughout \(J_K\).
Consequently each pair has a nonempty open eligibility segment and
is a genuine initial-fan exchange, with primitive direction
\[
                         (1,-A,h-1).                    \tag{HP5}
\]

Indeed \(1/2<u<1\), \(u>3/5\), and
\(\lambda> A/5+2>2\max(u,1)\). Also \(\lambda=O(A)\).
Both transverse magnitudes in (HP5) are \(\Theta(A)\), uniformly in
this range. Their product is \(\Theta(A^2)\), not \(O(at)\) for
the fixed thickness \(t=2\).

## Exact contraction of the cheaper original columns

Use residual rows \(Y_{q,v},Z_{q,v}\), \(0\le q<A\),
\(v\in\{0,1\}\), after the compulsory plain diagram and upper corner
box are contracted. The local long coordinate is \(y\) on Y and
\(z\) on Z. Write
\[
 G_j=\sum_{q=0}^{A-1}\binom{A-1}{q-j}Z_{q,1},            \tag{HP6}
\]
with out-of-range binomial coefficients zero. Put \(e=Y_{A-2,1}\).
The four central bridges are exactly
\[
\begin{array}{ll}
 H_{00}=Y_{0,0}+Z_{0,0}+e+G_0,&
 H_{10}=Y_{1,0}+Z_{0,1}+Y_{A-1,1},\\
 H_{01}=Y_{0,1}+Z_{1,0}+G_1,&
 H_{11}=Y_{1,1}+Z_{1,1}.
\end{array}                                             \tag{HP7}
\]
All four bridges and all eight upper identities with local \(q<2\)
are strictly cheaper than \(\tau\). For the largest Y-core upper
cost the difference is \(h+1-3u>0\); for the largest Z-core upper
cost it is \(h-1-u>0\). The bridge inequality follows from
\(\tau-(3u+3)=(A-3)u>0\).

For \(q\ge2\), the pure Y lower columns are individual rows, since
\(q^{A-2}v\) kills these inputs. All \(Y_{q,0}\) are cheaper;
among \(Y_{q,1}\), precisely \(q<A-2\) are cheaper. Equation
(HP7) additionally kills \(Y_{A-1,1}\), kills \(G_1\), and
identifies \(e=-G_0\). Thus only this one Y row remains relevant.
Any cheap upper Y column outside the core is redundant: its parallel
lower source is cheaper by \(\lambda-2>0\).

For \(q\ge2\), the Z lower columns are
\[
 L_{q,0}=Z_{q,0}+G_q,\quad c(L_{q,0})=2u+2+q,
 \qquad
 L_{q,1}=Z_{q,1},\quad c(L_{q,1})=3u+2+q.               \tag{HP8}
\]
The cheap tail identities are precisely \(q\le K\), because
\((A-3)u+1\in(K,K+1)\). The cheap upper Z heads have \(q<h\):
their costs are \(\lambda+2+q\). Every corresponding lower head,
including \(q=h\), is cheap, since
\(h<K<(A-2)u+1\). Pairing these columns supplies
\(G_2,\ldots,G_{h-1}\). The head of the tied column \(F_h\)
can be eliminated using its cheaper lower source, leaving \(-G_h\).

Any other cheap lower Z head has an exclusive head row: its upper
identity is not cheap and no other remaining column reaches that row.
It cannot participate in a relation supported only in the tail and
may be eliminated without creating another tail relation. Any cheap
upper Z tail has \(q<h-u\), hence \(q\le h-1<K\), and is already
in the cheap tail-identity span. This accounts for **all** cheaper
columns, not just a necessary projection of them.

The original two images are therefore, up to signs, exactly
\(G_0,G_h\) in
\[
 \langle Z_{q,1}:K<q<A\rangle\big/
                       \langle G_1,\ldots,G_{h-1}\rangle. \tag{HP9}
\]
This is a vector-space contraction; it is not asserted to be an
algebra quotient. Original exponent labels are retained throughout.

## Two classical minors certify the exchange

The tail in (HP9) has dimension \(h\). For \(D=A-1\), the
classical binomial Toeplitz evaluation is
\[
 T_h(r)=\det\left[\binom D{r+i-j}\right]_{0\le i,j<h}
       =\prod_{j=0}^{h-1}
          \frac{(D+j)!\,j!}{(r+j)!\,(D-r+j)!}.           \tag{HP10}
\]
See [Krattenthaler, equation (2.17), Section 2.3](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11).
Transposition converts its displayed binomial indices to the convention
in (HP10). The factorial product is the same evaluation, not a new
determinant identity.

The matrices with columns \(G_1,\ldots,G_h\) and
\(G_0,\ldots,G_{h-1}\), restricted to the tail rows, have determinants
\(T_h(K)\) and \(T_h(K+1)\), respectively. They are nonzero.
Every factorial argument is less than \(2A\), which also proves
nonvanishing in the stated positive characteristics.

It follows that \(G_1,\ldots,G_{h-1}\) have rank \(h-1\), while
both \(G_0\) and \(G_h\) extend them to rank \(h\).
Their images in the one-dimensional quotient (HP9) are therefore
nonzero and parallel. This proves the full original exchange assertion.

## The same dual exceptional label at every wall

All Y reference rows except \(e=Y_{A-2,1}\) lie in the cheaper span.
A separating functional nonzero on the tied pair therefore has Y
coefficients supported exactly on that one row. Reversing the Y rectangle
as in the [filtered-annihilator reduction](hook_filtered_annihilators.md)
gives a nonzero scalar multiple of \(q\). For the same operator
\(N_Y=q^{A-2}v\),
\[
          N_Yq=q^{A-1}v\ne0\quad\bmod(q^A,v^2),
 \qquad N_Yq=0\quad\bmod(q^{A-2},v^2).
\]
The full source map has precisely two independent nonzero columns,
labelled \((0,0),(1,0)\), whereas its length-\(A-2\) truncation
vanishes. Thus the dual exceptional set has size two, and the very
same member \((1,0)\) is witnessed at every displayed exchange.
The linear partner growth cannot be assigned to a growing set of
exceptional labels.

## Sharp local growth and the remaining aggregate problem

There are at least \(3A/20-1\) distinct partners in (HP3).
The entire upper Z arm has only \(2A\) labels, so the maximal
one-column partner order is sharply \(\Theta(A)=\Theta(n)\) at
thickness two. This concerns the maximum possible number of distinct
partners, not an exact count of every partner of this example.

All four separate pure-list orders are fixed on \(1/2<u<1\).
Within a Y list, a possible order change has \(u=1/d\) for a positive
integer \(d\); within a Z list it has an integer ratio. None lies
in that open interval. These are orders within each layer and arm;
combined upper/lower orders still depend on \(\lambda\).

For this label \(Y=A,t=2\), the low-wedge estimate
\(O(at^4/Y+t^5)\) is bounded independently of \(A\). Equations
(HP3)--(HP5) disprove its extension to the entire positive domain.
They also exclude bounding every relevant high-price wall by
\(\lambda\le C t\max(u,1)\), or by \(C t(u+1)\), with an
absolute constant \(C\). The transverse perturbations square to zero,
so this obstruction remains inside the already quadratically bounded
passive-band class.

A proof of a near-linear total direction inventory must therefore
amortize these long partner lists across different original labels or
use additional simultaneous-rank structure. Their existence alone
does not refute a quadratic-polylog full-fan bound. Neither the global
count, its sharpness nor the original extremal objective is resolved.

## Verification and priority limits

~~~bash
python3 -B scripts/run_hook_high_price_partner_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker rebuilds the original quotient and commuting multiplication,
tests the complete cheaper-column ranks, compares independent binomial
determinants with the product, checks adjacent original greedy bases,
and reconstructs complete small original fans. Finite-field examples
use characteristics above the factorial bound. They support the uniform
proof, rather than substituting for it.

The determinant evaluation and greedy rank criterion are classical.
The earlier [low-price binomial obstruction](hook_exception_partners.md)
already has an unbounded partner list at thickness two; the present
statement supplies linear growth at prices of order the length and
shows why the clipped low-wedge estimate cannot cover the complement.
Exact current novelty of this quantitative realization is unestablished.
No new unrestricted exponent or polynomial lower separation is claimed.
