# Two flag birth curves bound all crossings in axis-dependent hooks

There is a collective crossing bound when each central arm depends only
on its long variable. It allows every coefficient degeneration, including
nonunit and zero arms, and does not assume regular Toeplitz prefixes.
The two outer images vanish. Put
\[
 A=a-t\ge t\ge2,\qquad n=4at-2t^2\asymp At.
\]
For arbitrary polynomials \(F,G\) in the respective truncated rings,
consider the actual compatible hook with central image
\[
              f_1=y^tF(y)+z^tG(z),\qquad f_0=f_2=0.
                                                               \tag{AB1}
\]
Equivalently its ideal is
\(I=(x^2,y^a,z^a,y^tz^t-x(y^tF(y)+z^tG(z)))\).
Over every field,
\[
 \boxed{\quad g(I)=O(A^2+At^4),\qquad g(I)=O(n^{17/8}).\quad}     \tag{AB2}
\]
The second bound retains the earlier box bound whenever smaller.
In particular \(A\ge c t^2\), for any fixed \(c>0\), gives
\(g(I)=O_c(n^2)\). This covers the critical aspect \(a=t^2+1\),
with arbitrarily large prefix gaps and inverse degrees. No sharpness or
strongest-current priority is asserted. Arbitrary transverse dependence,
nonzero outer maps and the unrestricted extremal endpoint remain open.

The proof counts intersections of full interaction curves. It does not
infer the count from a list of supporting planes. It also gives an exact
minimum-profile oracle using two one-dimensional elimination problems
and \(t^2\) scalar maxima. The same proof applies to the more general
replicated matrix family described below; commuting multiplication is
not needed for that statement.

## A replicated matrix family

Let \(V_Y,V_Z\) have dimension \(A\), and let \(U_q\), \(U'_q\),
\(0\le q<A\), be bases. Choose arbitrary vectors
\(L_q\in V_Y,L'_q\in V_Z\), \(t\le q<A\), and
\(h_i\in V_Y,k_j\in V_Z\), \(0\le i,j<t\).
Use row space
\[
             (V_Y\otimes K^t)\oplus(K^t\otimes V_Z).
\]
With the standard short bases \(e_i,e_j\), give its columns these
original integer labels:
\[
\begin{array}{c|c|c}
 &\text{vector}&\text{label}\\\hline
 Y\text{ upper}&(U_q\otimes e_j,0)&(1,t+q,j)\\
 Y\text{ lower}&(L_q\otimes e_j,0)&(0,t+q,t+j)\\
 Z\text{ upper}&(0,e_i\otimes U'_q)&(1,i,t+q)\\
 Z\text{ lower}&(0,e_i\otimes L'_q)&(0,t+i,t+q)\\
 \text{head}&(h_i\otimes e_j,e_i\otimes k_j)&(0,t+i,t+j).
\end{array}                                                     \tag{AB3}
\]
The ranges are those above. Loops are retained. The represented rank is
\(2At\). Define its fan by the full-dimensional positive linearity
regions of its minimum basis cost. Both bounds in (AB2), with \(n\)
replaced by \(At\), hold for every such matrix over every field.

For (AB1), take \(V_Y=V_Z=K[q]/q^A\), and, using \(r\) for the index,
\(U_r=U'_r=q^r\), \(L_r=q^rF\), \(L'_r=q^rG\),
\(h_i=q^iF\), \(k_j=q^jG\). The
[original bridge model](hook_bridge_reduction.md) gives exactly (AB3)
after the compulsory diagram and upper corner are removed. Their profile
translation is fixed and does not change the fan. All other outer-source
columns are zero.

## Birth functions in one factor

In one factor give \(U_q\) price \(\theta+q\), and \(L_q\) price
\(q\). Let \(E_\theta(\rho)\) be the span of columns of price at
most \(\rho\), and put \(H_{i-1}=\langle h_0,\ldots,h_{i-1}\rangle\).
If \(h_i\in H_{i-1}\), set \(\alpha_i=-\infty\). Otherwise define
\[
 \alpha_i(\theta)=\min\{\rho:h_i\in E_\theta(\rho)+H_{i-1}\}.
                                                               \tag{AB4}
\]
The minimum exists since the upper columns span. The other factor has
birth functions \(\beta_j\), using \(k_j\).

For each finite birth function:

1. It is continuous, nondecreasing and 1-Lipschitz. Indeed it is a
   minimum, over finitely many spanning subsets modulo \(H_{i-1}\),
   of maxima of their prices. Each maximum has these properties, and
   the minimum preserves them.
2. On every interval between consecutive integers it equals one ground
   price, either \(\theta+q\), \(0\le q<A\), or \(q\),
   \(t\le q<A\). Ordering changes only when an upper and lower
   price meet, at an integer in \([-A,A]\). There are \(O(A)\)
   linear pieces, including the two tails.
3. \(\alpha_i(t)\ge t>i\), because every ground price at \(t\) is
   at least \(t\), and \(h_i\notin H_{i-1}\). Also
   \(\alpha_i(\theta)\to-\infty\) as \(\theta\to-\infty\).
   Therefore \(\alpha_i(\theta)-i\) has exactly one zero
   \(z_i<t\). That zero is an integer: a constant piece has value
   at least \(t\), so the crossing occurs on an upper-price piece.

These statements include arbitrary dependencies among the chosen vectors.
An absent birth remains absent at every parameter. No rank is inferred
from a support matching.

## Two flags give the whole tensor prefix rank

The following elementary flag fact is useful. Let vectors \(C_e\in K^t\)
have prices \(c_e\), and vectors \(D_f\in K^t\) prices \(d_f\).
Give the three tensor families prices
\[
 C_e\otimes e_j:c_e+j,\quad
 e_i\otimes D_f:iu+d_f,\quad
 e_i\otimes e_j:iu+j,\qquad u>0.                              \tag{AB5}
\]
Compare the decreasing price flag of the \(C_e\)'s with the coordinate
suffix flag, and do the same for \(D_f\). Two flags have a common
splitting basis: choose a basis in each successive intersection quotient.
Equivalently the mixed differences of their intersection-rank table give
a permutation after completing partial flags. Let \(c_i^*,d_j^*\)
be the assigned birth prices; use \(-\infty\) for absent directions.

For every real threshold \(\rho\), the rank of all columns in (AB5)
of price at least \(\rho\) is exactly
\[
 \#\{(i,j):\max(c_i^*+j,iu+d_j^*,iu+j)\ge\rho\}.          \tag{AB6}
\]
To prove it, write the first prefix span as
\(\sum_j F_C(\rho-j)\otimes\langle e_j,\ldots,e_{t-1}\rangle\).
The replacement of a single \(e_j\) by its suffix span is valid because
\(F_C(\rho-j)\) increases with \(j\). The second prefix has the
transposed expression, and the third uses the two coordinate flags.
All three are coordinate subspaces in the common tensor splitting.
A coordinate appears if and only if one of the three inequalities in
(AB6) holds. Taking the sum of greedy birth prices proves that the
maximum basis cost is the sum of these \(t^2\) maxima.
This is a statement about every full prefix, not separate marginal ranks.

The simultaneous splitting and permutation description are classical;
see [Eriksson--Linusson, Section 1.2](https://fpsac-archive.github.io/FPSAC96/ARTICLES/18Eriksson.pdf).
The tensor-prefix calculation above is given explicitly to fix its scope.
Three arbitrary flags need not admit one simultaneous splitting basis.

## An exact formula for the original minimum cost

Normalize the positive weight to \((\lambda,u,1)\), and put
\[
 \theta=(\lambda-t)/u,\qquad \phi=\lambda-tu.
\]
Choose a minimum basis of the one-factor pure Y columns at parameter
\(\theta\), and of the Z columns at \(\phi\). Replicate them in
all short coordinates. Their union is a minimum basis of the pure
columns in (AB3); write its total original cost as \(P(\lambda,u)\).
A greedy basis preserves the span of every cheaper prefix. Consequently
adjoining the head columns does not require any omitted pure column.

Express the \(h_i,k_j\) in these bases. The usual kernel matrix for
the retained pure basis plus the \(t^2\) head columns is a dual matroid
of rank \(t^2\). Its columns have exactly the tensor form (AB5).
After subtracting the common head cost \(t(u+1)\), their factor
prices are \(u(q+\theta x)\), respectively \(q+\phi x\),
where \(x=1\) for an upper and \(x=0\) for a lower column.

The assigned flag prices in (AB6) are precisely
\(u\alpha_i(\theta)\) and \(\beta_j(\phi)\).
For example a head combination belongs to the cheaper pure span exactly
when every coefficient on a more expensive selected pure column vanishes.
Annihilating the earlier head coordinates identifies its first threshold
with (AB4). This is also the minimum-row pivot interpretation of the
two-flag rank table, including dependent head vectors.
More explicitly, if \(F_{>\rho}\) is the span of the coefficient
vectors on pure basis columns of price greater than \(\rho\), then
\[
 h_i\in E_\theta(\rho)+H_{i-1}
 \quad\Longleftrightarrow\quad
 F_{>\rho}\cap\langle e_i,\ldots,e_{t-1}\rangle
       \subseteq\langle e_{i+1},\ldots,e_{t-1}\rangle.
\]
This follows by taking the annihilator of the cheaper span and the
earlier head vectors. It proves the identification for the entire prefix.

Basis complementation in the dual now gives the exact original minimum
\[
 \boxed{\quad
 M(\lambda,u)=P(\lambda,u)-
 \sum_{i,j}\max\{u(\alpha_i(\theta)-i),\beta_j(\phi)-j,0\}.
 \quad}                                                        \tag{AB7}
\]
Absent terms are omitted. The sum of head-column costs cancels the common
cost and the \(iu+j\) terms of the dual maximum. This cancellation is
needed for (AB7); no head contribution is silently discarded.
The formula initially follows off the finitely many price ties and
extends to every positive weight by continuity.

Within one open factor-order chamber, each finite birth is one fixed
ground label. The maximum summand in (AB7) therefore has three affine
options. Its maximizing option gives a linear comparison normal, and
its gradient gives the original minimum profile. Ordinary elimination
in two \(A\)-dimensional factors constructs these data in \(O(A^3)\)
field operations and \(O(A^2)\) storage, since \(t\le A\).
Sorting uses \(O(A\log A)\) exact cost comparisons. The \(t^2\)
scalar maxima then evaluate the profile. This is a minimum
profile oracle, not a claim about an optimal basis-output algorithm,
bit complexity or fixed-prime smooth-number unranking speed.

## A projective change makes each interaction monotone

Include all factor price lines
\(\lambda=t+ku\), \(\lambda=tu+k\), \(-A\le k\le A\),
and the line \(\Gamma:\lambda=t(u+1)\). Their arrangement has
\(O(A^2)\) cells. They contain every change of the pure basis or a
birth label, and every zero of a finite birth minus its head index.
The only additional walls allowed by (AB7) are
\[
 u(\alpha_i(\theta)-i)=\beta_j(\phi)-j>0.                 \tag{AB8}
\]

Away from \(\Gamma\), the map \((\lambda,u)\mapsto(\theta,\phi)\)
has inverse
\[
          u={\phi-t\over\theta-t},\qquad \lambda=u\theta+t.
                                                               \tag{AB9}
\]
Thus \(\theta,\phi\) are on the same side of \(t\). Restrict the
inverse to \(\lambda,u>0\). Define
\[
              R_i(\theta)={\alpha_i(\theta)-i\over\theta-t},
 \qquad S_j(\phi)={\beta_j(\phi)-j\over\phi-t}.          \tag{AB10}
\]
On a birth piece \(\alpha_i=\delta\theta+q\),
\[
 R_i(\theta)=\delta+{\kappa\over\theta-t},\qquad
 \kappa=q-i+\delta t>0.                                 \tag{AB11}
\]
For \(\delta=0\), \(q\ge t>i\); for \(\delta=1\),
\(q+t>i\). Hence \(R_i\) is strictly decreasing on either side
of \(t\), by continuity across its breaks. The same holds for \(S_j\).

On the lower side, the positive part in (AB8) is
\(z_i<\theta<t,z'_j<\phi<t\). Both ratios run from zero down
to \(-\infty\). On the upper side both births are positive and
each ratio decreases from \(+\infty\) to zero or one, according
to the last piece. Equality \(R_i(\theta)=S_j(\phi)\) therefore
has at most one increasing graph \(\phi=F_{ij}(\theta)\) on each
side, over the common ratio range. Each graph meets any integer
\(\theta\)-line or \(\phi\)-line at most once. It has \(O(A)\)
pieces. On one piece it is the fractional linear graph
\[
 \phi-t={\kappa_Z(\theta-t)\over
                 \kappa_Y+(\delta_Y-\delta_Z)(\theta-t)}.       \tag{AB12}
\]
Its image under (AB9) is a straight line segment or ray. Clipping to
the positive weight region preserves an interval on each such piece.

## Count all curve intersections, including degeneracies

Fix two interaction graphs on the same side. Follow the first graph,
which is increasing in \(\theta,\phi\). Cut it at every integer
\(\theta\) or \(\phi\) where any of the four relevant birth
functions changes formula. There are \(O(A)\) cuts altogether:
each coordinate is monotone and all finite breaks lie in \([-A,A]\).

On one resulting interval, the first graph has (AB12), while the
second graph is tested by
\(R_k(\theta)=S_l(F_{ij}(\theta))\). Substitution gives
\[
 S_l(F_{ij}(\theta))=
 \delta_l+{\kappa_l(\delta_i-\delta_j)\over\kappa_j}
       +{\kappa_l\kappa_i\over\kappa_j(\theta-t)}.        \tag{AB13}
\]
Both sides are a constant plus a constant divided by \(\theta-t\).
They have at most one isolated intersection on this interval, or agree
throughout it. Coincident arcs are merged; their endpoints are among
the already charged cuts. Intersections at cuts contribute \(O(A)\)
as well. Each pair of full interaction curves consequently has
\(O(A)\) isolated intersections and overlap endpoints, not \(O(A^2)\).
This is the collective saving.

There are at most \(t^2\) interaction curves per side. Their own
breaks and their intersections with factor lines cost \(O(At^2)\),
and all pair contacts cost \(O(At^4)\). On \(\Gamma\), each
finite pair gives at most one point,
\(u(\alpha_i(t)-i)=\beta_j(t)-j\), because both coefficients
are strictly positive. Interaction curves cannot coincide with an open
part of \(\Gamma\). Factor lines contribute only their intersections
or a line already merged with \(\Gamma\).

Euler's formula for this planar arrangement, after merging coincidences
and clipping straight pieces to the positive projective triangle, gives
\(O(A^2+At^4)\) cells. It retains all subdivisions at crossings,
higher-valence vertices, recurrent supporting directions and boundary
pieces. Since the arrangement refines (AB7), this proves the first
bound in (AB2) for (AB3) and hence for the actual hooks.

## Combine with the classical box estimate

Every label in (AB3) has first coordinate zero or one and other
coordinates \(O(A)\). Its rank is \(2At\). The integral basis-profile
polytope lies in a box of volume \(O((At)^3A^2)\). Add the integral
unit cube if needed to obtain full dimension; its normal fan refines
the original. The classical three-dimensional lattice-vertex estimate
therefore gives
\[
                        g=O(A^{5/2}t^{3/2}).              \tag{AB14}
\]
This is the same classical box input used in the
[hook bridge note](hook_bridge_reduction.md), following
[Babson--Onn--Thomas, Corollary 2.5](https://arxiv.org/pdf/math/0207135).
It is not a new lattice-volume theorem.

Write \(N=At\). The \(A^2\) term is at most \(N^2\), and
\[
 \min\{At^4,A^{5/2}t^{3/2}\}
       =\min\{Nt^3,N^{5/2}/t\}\le N^{17/8}.
\]
Split at \(t=N^{3/8}\) to verify the last inequality. Combining
the two bounds proves the second part of (AB2). No sharp exponent is
inferred. For \(A\ge c t^2\), the first bound directly gives
\(g=O_c(A^2t^2)=O_c(n^2)\).

## Verification and scope

~~~bash
python3 -B scripts/run_hook_axis_birth_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker uses original quotient replay, independent entire-prefix
birth tests, two-flag intersection ranks, full tensor-prefix ranks,
minimum profiles and whole original cone patches. Exact rational curve
segments test monotonicity, coincident arcs and all pair contacts against
an independent original-weight line-intersection calculation. Unit,
nonunit, zero and small-characteristic cases are included.

The earlier prefix-regular theorem already covers the generic compatible
hook family. This proof also covers every degeneration in the stated
axis-dependent class, and the larger replicated matrix family.

Greedy prefix preservation, dual complementation, two-flag splitting,
planar arrangements and the lattice-volume bound are classical inputs.
Exact strongest-current priority of the combined crossing bound and
profile oracle remains unestablished. Off-axis terms need not retain
the repeated tensor structure; a basis change is not asserted to remove
them. This result does not settle sharp unrestricted \(G_3\), construct
an infinite actual polynomial excess or complete the original paper goal.
