# Quadratic fans for arbitrary scalar-axis hooks

Let \(K\) be any field, \(A=a-t\ge t\ge2\), and let \(F,G\in K[q]/q^A\)
be arbitrary, including nonunits and zero. For the actual ideal
\[
 I=(x^2,y^a,z^a,
       y^tz^t-x(y^tF(y)+z^tG(z))),\qquad n=4at-2t^2,
\]
the number \(g(I)\) of monomial initial ideals satisfies
\[
                         g(I)=O(A^2t^2)=O(n^2).                \tag{LW1}
\]
The constant is independent of the coefficients, degrees, initial orders
and characteristic. Both outer images are zero and each central arm
depends only on its long variable. No regular-minor hypothesis is imposed.

There is also a coefficient-dependent list of \(O(n^2)\) positive integer
weights, of height \(O(a^2)\), meeting every full-dimensional positive
initial cone. Thus exact basis reconstruction at this many weights
recovers the entire fan. The list is not asserted to be common to all
coefficient choices. Arbitrary transverse dependence, nonzero outer
images and unrestricted \(G_3(n)\) remain unresolved.

The proof uses classical rank-two polynomial lattice reduction. Its
additional step is a short multiplier reduction modulo \(q^t\), followed
by a coarse partition of the moving pivot degree. This retains all
actual rank dependencies, including nongeneric cancellations. The
previous [birth-curve bound](hook_axis_birth_curves.md) applies to a
larger replicated matrix class; the improvement here uses the scalar
polynomial shift structure and is not asserted for that larger class.

## The exact factor problem

Fix one operator \(W\in K[q]/q^A\). In \(V=K[q]/q^A\), the pure sources are
\[
 U_j=q^j\quad(0\le j<A),\qquad
 L_j=q^jW\quad(t\le j<A),
\]
with prices \(\theta+j\) and \(j\), respectively. The heads are
\(h_i=q^iW\), \(0\le i<t\). Write \(E_\theta(\rho)\) for the span of
pure sources priced at most \(\rho\), and define the actual head birth
\[
 \alpha_i(\theta)=\min\{\rho:
       h_i\in E_\theta(\rho)+\langle h_0,\ldots,h_{i-1}\rangle\}.
                                                                  \tag{LW2}
\]
If the head is already in the earlier-head span, set \(\alpha_i=-\infty\).
Each finite birth is a pure-source price. On every open unit interval
it has one fixed label
\[
              \alpha_i(\theta)=\delta_i\theta+q_i,
              \qquad \delta_i\in\{0,1\},\quad 0\le q_i<A.       \tag{LW3}
\]
The entire price order, pure greedy basis and birth assignment can change
only at integers in \([-A,A]\). These are actual elimination statements,
not ranks inferred from a support pattern.

For \(w=(x,y,z)>0\), set
\[
                  \theta={x-tz\over y},\qquad
                  \phi={x-ty\over z}.
\]
The established two-flag tensor-prefix identity gives the original
minimum-basis cost, up to its fixed linear translation, as
\[
 P_Y(w)+P_Z(w)-\sum_{i,j=0}^{t-1}
       \max\{,y(\alpha_i(\theta)-i),
                   z(\beta_j(\phi)-j),0\,\}.                   \tag{LW4}
\]
Here \(P_Y,P_Z\) are the pure replicated minimum-basis costs and
\(\beta\) is the second factor's actual birth sequence. All absent
options are omitted. See [the full prefix-rank proof](hook_axis_birth_curves.md)
and [the original-profile compiler](hook_axis_fan_compilation.md).

## A rank-two lattice and its monotone pivot

Choose the representative \(\deg W<A\), and form
\[
 \mathcal L_W=\{(u,v)\in K[q]^2:u-Wv\in q^AK[q]\}
       =\langle(q^A,0),(W,1)\rangle_{K[q]}.
                                                                  \tag{LW5}
\]
Its determinant is \(q^A\). For nonintegral \(\theta\), give a pair degree
\[
             d_\theta(u,v)=\max(\theta+\deg u,\deg v).
\]
Zero components have degree \(-\infty\). Cancel the leading term of a
row by a monomial multiple of the other whenever their leading
positions agree. Each step is unimodular and strictly lowers the
modified row's weighted degree. Degrees lie in the discrete set
\(\mathbb Z\cup(\theta+\mathbb Z)\), are bounded below for nonzero
polynomial pairs, and never increase; hence reduction terminates.
The final rows \(g_0,g_1\) have different leading positions. Index them
by their leading coordinates, and let their pivot degrees be \(r,s\).
Then
\[
 d_\theta(g_0)=\theta+r,\qquad d_\theta(g_1)=s,
 \qquad r+s=A.                                                   \tag{LW6}
\]
Indeed the product of the two pivot entries has strictly larger degree
than the other determinant term. The determinant remains a nonzero
constant times \(q^A\), so its degree is \(r+s=A\).

Distinct leading positions imply the predictable-degree identity
\[
 d_\theta(f_0g_0+f_1g_1)
     =\max(\theta+r+\deg f_0, A-r+\deg f_1).                    \tag{LW7}
\]
For nonzero \(f_0,f_1\), their leading weighted degrees have different
fractional parts, so cannot tie. The larger one has a leading term
that the other row cannot cancel. This also proves (LW7) if one
multiplier is zero.

In particular, \(r(\theta)\) is the least degree of the first coordinate
of any lattice vector leading in that coordinate: every such vector
has a nonzero \(g_0\) multiplier and pivot degree at least \(r\).
Increasing \(\theta\) preserves the leading position of every vector
already leading in its first coordinate. Thus
\[
          0\le r(\theta)\le A,\qquad r(\theta)
                         \text{ is nonincreasing}.              \tag{LW8}
\]
The leading comparisons of any fixed polynomial pair change only at
integer \(\theta\). Therefore \(r\) is constant on each open unit
interval. No choice of reduced basis changes these minimal pivot degrees.

## Short multipliers locate every active birth

Let \(\gamma_i(\theta)\) be the least weighted degree of a lattice pair
\((u,v)\) such that
\[
                 v\bmod q^t=q^i+\sum_{j<i}c_jq^j.              \tag{LW9}
\]
Then
\[
                     \gamma_i=\max(i,\alpha_i),                 \tag{LW10}
\]
where the right side is \(i\) for an absent birth. To prove this,
first restrict to \(\rho\ge i\). A relation witnessing (LW2) is exactly
a congruence \(u=Wv\pmod {q^A}\) with (LW9),
\(\theta+\deg u\le\rho\) and \(\deg v\le\rho\): the coefficients
of \(v\) below \(t\) specify the head and earlier heads; the higher
coefficients specify the lower pure sources. Both polynomials may be
reduced modulo \(q^A\) without changing the congruence or (LW9), and
without increasing their degrees. Conversely every such pair supplies
the required source relation. Condition (LW9) itself forces
\(d_\theta(u,v)\ge\deg v\ge i\), proving (LW10), including absent
heads.

Write an optimal pair in the basis from (LW6). Replace both
multipliers \(f_j\) by \(f_j\bmod q^t\). Its second coordinate modulo
\(q^t\) is unchanged, so (LW9) persists, and (LW7) says its weighted
degree cannot increase. Thus an optimum uses two multipliers of degree
less than \(t\). Its cost belongs to
\[
 \gamma_i\in
 \{\theta+r+k:0\le k<t\}\ \cup
 \{A-r+k:0\le k<t\}.                                           \tag{LW11}
\]
If \(\alpha_i>i\), (LW10) identifies this cost with the actual birth.
Since \(\theta\) is nonintegral and all labels are integers, its price
type must match the corresponding band in (LW11). Consequently
\[
 \boxed{\quad
 \begin{array}{ll}
 \delta_i=1:&q_i\in r(\theta)+\{0,\ldots,t-1\},\\
 \delta_i=0:&q_i\in A-r(\theta)+\{0,\ldots,t-1\}.
 \end{array}\quad}                                             \tag{LW12}
\]
This is an actual-rank statement for every \(W\), not a generic
Toeplitz formula. It is required only for positive gains in (LW4).

## Coarse intervals have only linearly many possible gains

Inside \([-A,A]\), cut at an integer grid of spacing \(t\), and also
at each integer where \(\lfloor r(\theta)/t\rfloor\) changes.
The grid contributes \(O(A/t)\) intervals. Monotonicity (LW8) gives at
most \(\lfloor A/t\rfloor\) changes of pivot bucket, even if a single
jump skips several buckets. There are therefore \(O(A/t)\) intervals
in the common refinement. Each has length at most \(t\), and some
fixed integer \(R\) satisfies \(R\le r<R+t\) throughout it.

Put \(e=q_i-i\). Formula (LW12) bounds every positive gain label on
that entire coarse interval by
\[
 \begin{array}{ll}
 \delta=1:&R-t+1\le e\le R+2t-2,\\
 \delta=0:&A-R-2t+2\le e\le A-R+t-1.
 \end{array}                                                    \tag{LW13}
\]
Intersect these integer bands with the original range
\(-t+1\le e\le A-1\). Together they contain fewer than \(6t\) labels.

The two tails also have \(O(t)\) labels. If \(\theta<-A\), all upper
prices are negative, so every finite birth is below \(i\); there are
no positive gains. If \(\theta>A\), all lowers precede all uppers.
For nonzero \(W\) put \(\nu=\operatorname{ord}_q W\). The lower
span is \(q^{\nu+t}V\), interpreted as zero when \(\nu+t\ge A\).
Modulo that span every head is supported on upper coordinates
\[
                   \nu\le q<\min(A,\nu+t).
\]
The upper minimum basis has precisely the prefix below
\(\min(A,\nu+t)\). Hence every finite head birth has an upper label
in the displayed interval, giving fewer than \(2t\) possible \(e=q-i\).
For \(W=0\) all births are absent. Pure-source orders are constant on
both tails.

## Count complete local arrangements

Take one coarse interval from each arm. Its inverse image in the
positive weight cone is an intersection of linear inequalities:
\(\theta\in(\ell,h)\) means
\(x-tz-\ell y>0\) and \(hy-x+tz>0\), omitting an infinite endpoint.
The analogous inequalities use \(\phi\). Thus each pair gives one
convex cone, possibly empty. There are \(O((A/t)^2)\) pairs. This
argument requires no injectivity of the map \(w\mapsto(\theta,\phi)\).

On a fixed pair the possible positive gain gradients in (LW4) are
\[
          (\delta,e,-t\delta),\qquad
          (\delta,-t\delta,e),                                 \tag{LW14}
\]
drawn from \(O(t)\) labels on each arm. Include their comparisons
with zero and all cross-arm comparisons. This gives \(O(t^2)\) planes.
Also include the factor-price pencils for every integer in either
finite coarse interval, and its endpoints. These add only \(O(t)\)
planes, because each finite interval has length at most \(t\).
Tails add their single finite endpoint.

These planes refine every change in (LW4). Pure bases and birth labels
are fixed between consecutive factor-price pencils. A head gain can
change sign only on an integer pencil: a lower gain has constant
integer sign, and an upper gain vanishes at \(\theta=-e\).
After the sign and price subdivisions, all remaining changes are
comparisons of the positive gradients (LW14), which were included.
Negative gains cannot win against the zero option.

An arrangement of \(O(t^2)\) lines in the projective weight triangle
has \(O(t^4)\) cells. Counting its full subdivision, including contacts,
intersections and parallel subdivisions, gives
\[
          O((A/t)^2)\,O(t^4)=O(A^2t^2)=O(n^2).                 \tag{LW15}
\]
This is an arrangement refinement of the entire minimum function, not
an identification of one supporting plane with one wall. It therefore
counts all actual full-dimensional initial cones. Since \(A\ge t\),
\(n=4At+2t^2\asymp At\), proving (LW1).

## Integer weights and field work

Every normal used above is integral of height \(O(a)\); the coarse
endpoints are integers in \([-A,A]\). Extreme rays are primitive cross
products of two such normals or coordinate normals, so have height
\(O(a^2)\). In every full-dimensional cell choose three independent
extreme rays and sum their primitive nonnegative representatives.
Their sum is a positive integer interior weight of height \(O(a^2)\).
Taking one per nonempty cell proves the stated complete weight-list
consequence. Exact basis reconstruction at these weights recovers all
standard sets; convex merging recovers their normal cones.

The pivot sequence may be computed independently at the \(2A\)
half-integer parameters. A two-row cancellation lowers weighted degree
by at least \(1/2\). Initially the degree sum is
\(A+\theta+\max(\theta+\deg W,0)\), and finally it is \(A+\theta\).
For \(-A<\theta<A\) this permits \(O(A)\) steps. Intermediate
polynomial degrees are \(O(A)\), since each row's weighted degree
never increases. Each step costs \(O(A)\) field operations, so all
pivots cost \(O(A^3)\) field operations and \(O(A)\) field storage.
Zero \(W\) needs no reductions.

Combining this with the [actual incremental factor compiler](hook_axis_linear_oracle.md)
retains \(O(A^3+A^2t^2)\) total field operations and \(O(A^2)\) field
storage before profile traversal. The established generic profile
oracle still uses \(O(t+\log A)\) rational/integer operations and no
later field elimination. Polygon construction, arrangement generation,
coefficient growth and basis output have additional costs; no optimal
bit-complexity or total scalar-runtime claim is made.

## Verification and classical inputs

```bash
python3 -B scripts/run_hook_axis_lattice_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker independently constructs polynomial lattices and finite
source-prefix ranks. It checks determinant and congruence identities,
least pivots by separate finite linear equations, short-lattice versus
original birth costs, matching price-type bands, monotone buckets and
both tails. It includes every polynomial of degree below five over
\(\mathbb F_2\), every polynomial of degree below three over
\(\mathbb F_3\), and rational and finite-field sparse, dense, nonunit
and zero fixtures. Exact complete subdivisions have area one, every
interior facet belongs to its local inventory, profile polygons merge
convexly, and small integer samples reproduce their profiles. Small
fans are independently compared with the entire original quotient fan;
larger original-quotient comparisons are labelled as spot checks.

The classical primary inputs inspected for this argument are:

- [Mulders--Storjohann, Section 2, Lemmas 2.1--2.6 and Theorems 2.1--2.2](https://cs.uwaterloo.ca/~astorjoh/jscpopov.pdf):
  unimodular leading-position reduction and its termination and cost.
- [Labahn--Neiger--Vu--Zhou, Sections 2.5--2.6, Lemmas 2.1 and 2.4](https://arxiv.org/pdf/2202.09329):
  shifted predictable degree and pivot properties. The nonintegral
  rank-two variant used here has the elementary proof (LW6)--(LW8).

Polynomial lattice reduction and predictable degree are not asserted
new. The inspected statements do not themselves give the hook-fan
count or coefficient-dependent complete weight list. This bounded
comparison does not establish strongest-current novelty of their
combination. The generic [superregular theorem](hook_superregular_axis_fans.md)
additionally gives coefficient-independent profiles and a common weight
list on its hypothesis class; neither stronger uniformity is claimed
here. This closes all scalar-axis zero-outer coefficient strata for
the quadratic upper bound, but supplies no unrestricted endpoint,
infinite polynomial lower excess or completed publication assessment.
