# Compile complete axis-hook fans before traversing their profiles

The [two flag birth formula](hook_axis_birth_curves.md) can be compiled
for all positive weights. The field computation then finishes before
fan traversal starts. Exact rational polygons give the entire actual
profile fan; cells of the auxiliary subdivision must be merged.
This is a coefficient-dependent construction over every field, including
nonunit and zero arms. Both outer images vanish and the central arms
depend only on their long variables.

Put \(A=a-t\ge t\ge2\), \(n=4at-2t^2\). For fixed coefficients,
there is a compilation with
\[
 \begin{array}{ll}
 O(A^4)\text{ field operations},&O(A^2)\text{ field storage},\\
 O(t^2+\log A)\text{ scalar operations per generic profile query},&
 \text{no later field elimination}.
 \end{array}                                                   \tag{AC1}
\]
Scalar operations here include exact rational cost comparisons and integer
profile arithmetic. They are not bit-complexity estimates. A complete
fan traversal uses \(O(K)\) patch-oracle calls, where
\[
                         K=O(A^2+At^4)                       \tag{AC2}
\]
is the number of refinement cells. Thus total field work for the whole
fan stays \(O(A^4)\), rather than repeating elimination in each cell.
Polygon construction and merging require additional exact rational
operations; no optimal bound for those operations or basis-output time
is asserted. The earlier \(O(n^{17/8})\) count bound concerns actual
cones, and is not used as a bound on the finer number \(K\).

This also compiles the larger replicated matrix family (AB3). The axis
hook is its actual commuting-quotient specialization. No arbitrary
replicated matrix is silently asserted to be an algebra quotient.
The construction supplies no new unrestricted exponent, asymptotic lower
excess, common coefficient-independent list or current novelty certificate.

## Compile both one-factor price orders

In one factor upper columns have price \(q+\theta\), and lower columns
have price \(q\), with \(0\le q<A\), respectively \(t\le q<A\).
Their order can change only at integral \(\theta\in[-A,A]\).
For each open interval, including both tails, evaluate at one exact
half-integer. Compute a greedy minimum basis and express the head vectors
in it. The maximum-price coefficient flag and minimum-row pivots then
assign the finite head birth prices to their original labels \((x,q)\),
where \(x\in\{0,1\}\). Dependent head directions remain absent.
These are the two-flag data from the birth theorem, not support matchings.

A record contains the interval, pure-basis profile and assigned labels.
Adjacent records with the same data may be merged. Each one-factor
elimination takes \(O(A^3)\) field operations; there are \(O(A)\)
intervals on each arm. The basis and flag data require \(O(A^2)\) field
storage during compilation. The retained records use \(O(A^2)\) integer
label storage. This proves the field bounds in (AC1).

If a one-factor basis contains \(s_x\) upper labels and has sum of
long inputs \(s_q\), its replicated Y profile is
\[
 \left(ts_x,\ t(At+s_q),\
             t^2(A-s_x)+{At(t-1)\over2}\right).            \tag{AC3}
\]
Transpose its last two entries for the Z profile. This avoids rebuilding
all \(At\) repeated column labels in each query.

## Evaluate one whole open patch

For a positive weight \(w=(x,y,z)\), find the two records containing
\[
                      \theta={x-tz\over y},\qquad
                      \phi={x-ty\over z}.                \tag{AC4}
\]
A binary search uses \(O(\log A)\) rational comparisons. For every
head pair \((i,j)\), make the three integer-vector options
\[
 0,\qquad (\delta_Y,q_Y-i,-t\delta_Y),\qquad
              (\delta_Z,-t\delta_Z,q_Z-j),                \tag{AC5}
\]
omitting an absent birth. Select the option with largest scalar product
with \(w\). Subtract their sum from the two pure profiles in (AC3).
The exact minimum formula (AB7) proves that the result is the original
residual minimum profile. It uses \(t^2\) constant-size maxima and
no operation in the coefficient field.

For the actual ideal, add the fixed translated profile
\[
                \left(t^2,{at(a+t-2)\over2},
                               {at(a+t-2)\over2}\right). \tag{AC6}
\]
It is the sum of the compulsory plain hook diagram and upper corner box.
The checker independently compares this formula with the original fixed
monomials. Profiles in the reported fans include (AC6).

The active record endpoints and the inequalities saying that each chosen
option beats the other options give one convex polygon in the positive
projective triangle. They are inequalities between integer normals of
height \(O(a)\). On this entire open polygon the formula is one linear
function, with the returned profile as gradient. The global exact minimum
formula proves agreement on all of it, rather than only at the query.
Queries on a genuine record boundary or a nontrivial comparison are
rejected; exact interior queries are used in the traversal.
A profile-only query has the scalar complexity in (AC1); constructing
its polygon has an additional cost.

## Traverse the refinement exactly

Start with a generic integer weight whose successive coordinates are
large powers. For each newly encountered patch, clip its inequalities
against the projective triangle. Across every non-coordinate edge, the
[exact adjacent-weight construction](oracle_fan_enumeration.md) gives
an integer weight strictly on the other side and sufficiently close to
preserve every strict bounded-height comparison at the edge midpoint.
Birth values are continuous across a factor boundary, so an active
comparison of the next record can vanish on the edge interior only
when its plane coincides with that edge. The perturbation therefore
enters a full adjacent patch rather than an independent active tie.

All normals used here satisfy the same height restriction as original
label comparisons, so the established universal perturbation denominator
applies. The adjacent patch has a different record/option signature.
The signature identifies one convex open cell: its inequalities and its
profile formula are fixed. It can therefore be safely marked visited.
The adjacency graph of a finite subdivision of a convex triangle is
connected, so exhausting the queue visits every full-dimensional patch.
Repeated queries on already visited signatures are allowed.

The birth-curve proof bounds this refinement by (AC2), including all
factor-line and curve contacts, coincident arcs and the singular price
line. A planar subdivision has \(O(K)\) edges, hence \(O(K)\) calls
in the adjacency walk. All calls after compilation use only the retained
integer labels and exact rational geometry. Assertions verify strict
queries, positive polygon areas and exact total doubled area one.
Truncated or partial traversal does not certify completeness.

## Merge by original profiles and recover full normal cones

A refinement signature is not an initial ideal. Its boundaries can record
an irrelevant change of one factor or cancel in the full minimum formula.
Group patches by their returned original profile. The minimum function is
the support minimum of the original integral profile polytope. Therefore
the full region with one gradient is a convex normal region. Different
actual initial ideals at generic weights have different profiles: strict
monomial costs give a unique greedy residue basis, so two available bases
with the same profile could not be distinct strict optima.

For each group, take the convex hull of all patch vertices. Its doubled
area must equal the sum of the group's patch areas. The hull is the full
original normal polygon. The equality both checks the merging operation
and detects accidental gaps or nonconvex grouping. The merged polygons
again have total area one. This produces the actual normal cones and
count, not the size of the factor arrangement.

For example, one length-56 fixture has 55 refinement patches but only
33 actual cones. Other fixtures include a nonunit characteristic-two
case and a zero arm. The report compares every merged polygon with an
independently enumerated original quotient fan at the stated small sizes.
At larger sizes it retains a complete compiled cover and separately
replays original commuting multiplication and selected original profiles;
those checks are not described as an independent full original fan.

## Actual remainder walls also cross inside a factor cell

The compilation finds a four-cone vertex away from every integral factor
price line and from the singular line. Over \(\mathbb F_{1009}\), take
\(t=4,a=17\), zero outer images and
\[
 F(q)=1+2q^5+3q^6+5q^{11},\qquad
 G(q)=1+7q^7+11q^8+13q^{12}.
\]
At
\[
 w_0=(9,17,18)/44,\qquad
 \theta=-63/17,\quad\phi=-59/18,
\]
the two original remainder directions \((1,-9,8)\) and \((1,9,-9)\)
continue on both sides of the vertex. Both are outside the original
known-direction inventory. Their two coherent exchanged pairs are
\[
 (0,10,5)\leftrightarrow(1,1,13),\qquad
 (0,4,9)\leftrightarrow(1,13,0).
\]
At the midpoint of each half-facet of the first line, the entire strictly
cheaper rank is 97; for the other line it is 80. Either column, or both,
raises that rank by exactly one. The single exchange is independently
verified on the minimum face. Exact adjacent weights reconstruct the
original neighboring bases, and their full normal polygons contain the
displayed half-facet. A neighbor transition can include additional
simultaneous exchanges; its symmetric difference is recorded rather than
assumed to have size two. Exactly four original cones meet at the vertex.

This refutes a noncrossing shortcut even inside one open factor-order
cell. It is a finite certificate, not a growing intersection count or
polynomial lower excess. The collective contact bound remains necessary.

## Reproduction and limits

~~~bash
python3 -B scripts/run_hook_axis_compilation_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The report records the exact operators, fields, complete counts, total
areas, hashes of canonical profile polygons and explicit small fan data.
The checks include a non-Toeplitz replicated matrix with a head already
in the pure lower span, without claiming an algebra realization for it.
Fresh field elimination independently checks compiled profile queries.
These independent verification eliminations are additional to the
field-operation bound for the compilation and traversal themselves.

The data permit further study in the still-open intermediate aspect
\(t\ll A\ll t^2\). A finite sequence of counts does not prove a growth
rate. The critical aspect \(A\ge c t^2\) and generic prefix-regular
hooks are already quadratically bounded as stated in the earlier notes.
Arbitrary transverse terms and outer maps need additional simultaneous
rank data; the compilation is not asserted to handle them.

Greedy/dual profile optimization, two-flag splitting and exact normal-fan
walks are classical inputs. Their earlier general forms and mathematical
citations are preserved in the linked notes. The coefficient-dependent
compilation and exact merging are derived here; exact strongest-current
priority is unestablished. Neither this algorithm nor the finite examples
settle the unrestricted sharp-count goal or supply the original final paper.
