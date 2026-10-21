# smooth.lattice

Fixed-prime smooth-number rank/unrank experiments based on recursive exponent lattices.

A `P`-smooth number for a fixed prime set `P = {p1, ..., pk}` is a number of the form

```text
p1^e1 * p2^e2 * ... * pk^ek,  ei >= 0.
```

This repository focuses on **random access**:

```text
Given P and a 1-based rank N, return the exponent vector of the N-th P-smooth number.
```

It is not intended to replace general `y`-smooth counting/listing algorithms, and it is not a prime-generation or factoring algorithm.

## Current contents

Core kernels:

| File | Purpose |
|---|---|
| `src/smooth_3prime_beatty_ranker.cpp` | Fast special ranker for 3-prime fixed sets using a two-prime floor-sum/Beatty-style inner count. |
| `src/smooth_layer_compressed_general.cpp` | Layer-compressed hybrid MITM kernel, strongest for 5-prime cases in our tests. |
| `src/smooth_sums_only_scalable.cpp` | Sums-only/adaptive MITM fallback for other and higher prime counts. |
| `src/smooth_interval_audit_exps_k6.cpp` | Independent interval-rank auditor for `k <= 6` and primes in `{2,3,5,7,11,13}`. |
| `smooth_toolkit.py` | Python CLI/API wrapper choosing kernels automatically. |

Benchmark comparators:

| File | Purpose |
|---|---|
| `benchmarks/dp_pointer_baseline.cpp` | Standard DP/pointer “super ugly number” baseline. |
| `benchmarks/smooth_xplusy_baseline.cpp` | Practical Cartesian-sum / adaptive `X+Y` baseline. |
| `benchmarks/smooth_xplusy_fj_loh_workbench.cpp` | FJ/LOH-style comparator workbench. |

Documentation and paper material:

```text
docs/       claim matrix, certification protocol, hard result table
paper/      quantitative manuscript and separate research notes
results/    selected certified examples and benchmark summaries
```

The active [extremal initial-ideal goal](docs/extremal_initial_ideals.md)
targets sharp growth, or a polynomial separation from the generic scale with
a general consequence, for finite algebras in three variables. It develops the
exponent-lattice connection behind the original interleavings. Uniform local
bounds and exact coefficient comparisons delimit the candidate families.
Generic torus-grid graphs reduce to ordinary geometric cuts; local rectangles
have genuine coherent rank defects, including in cones using all three
generators. Published constraints on neighboring ideals do not supply that
count: an explicit connected family has exponentially many individually
realizable candidate vertices, while one ideal realizes at most two of them.
A count of simultaneous interacting dependencies remains necessary in
three variables. The [dyadic planar proof](docs/planar_dyadic_bound.md)
now gives `G_2(n) = Theta(n log(n+1))` for all fields and finite planar
ideals, with a general matroid profile bound and a planar universal-basis
bound `O(n^(3/2) log(n+1))`. Its exact priority remains under audit.
Consequently the proved product construction `G_3(3m) >= G_2(m)^2`
cannot give a power separation through planar branches. Further
degree-sensitive and grading-sensitive bounds are retained as sharper
restricted comparisons. A [rational-shadow extension](docs/rational_shadow_bound.md)
now gives the all-field upper bound `G_3(n) = O(n^(23/8) log(n+1)^(3/8))`
and universal-basis cardinality `O(n^(85/24) log(n+1)^(3/8))`.
It bounds rational matroid shadows by their normalized integer height,
then uses short integer relations to cover vertices with bounded facet
normals. This is a uniform upper-bound improvement; sharp growth,
current novelty and the goal's publication assessment remain open.
A [rational-plane transfer](docs/rational_plane_transfer.md) extends this
arithmetic mechanism to every fixed dimension and arbitrary hyperbolic
matroid profiles. It gives `G_4(n)=O(n^(33/7)*log(n+1)^(9/35))`,
improving the older repository power `52/11` by `1/77`. It reproduces
the current three-variable bound; the earlier Dey-based powers remain
stronger in dimensions five and above. No current-best or sharp-growth
claim follows from this extension.
An [exact enumeration consequence](docs/oracle_fan_enumeration.md) now
computes all initial sets, their state fan and reduced bases in
`O(n^(47/8)*log(n+1)^(19/8))` field operations, given one reduced
Gröbner basis or commuting coordinate matrices with a cyclic vector.
It applies in the original coordinates over every field. The oracle
reconstruction and greedy conversion are established methods; with
the classical count alone their composition already gives exponent six.
The shadow theorem supplies the further saving of `1/8`. Coefficient
bit costs, current-best priority and practical speed remain unestablished.
A product-threshold refinement now counts every initial cell incident
to a wall with primitive direction product at most `T` by
`O(n*T^(4/3)*log(n+1))`, with an explicit support-oracle reduction.
For `T=O(n)` this is `O(n^(7/3)*log n)` for that part of the fan.
The remaining cells need a separate bound: the known rational-curve
family already has quadratically many interior cells whose every wall
has quadratic primitive product. No unrestricted improvement follows yet.
An [individual-shadow obstruction](docs/multiplicative_shadow_obstruction.md)
also rules out detecting every cone with primitive plane product `O(n)`.
An explicit all-field family has just seven initial ideals, but one cone's
minimum detecting product is `Theta(n^(4/3))`. The uniform proof and
exact quotient checks distinguish interior detection from touching a wall.
This closes a proposed shortcut; a collective count of difficult cones
and the unrestricted extremal problem remain open.
Fixed internal list orders give exact rank-oracle formulas, but a
polylogarithmic cover by two-list domains is impossible even for that
quadratic family. Three-list ranks require a signed correction with
quadratically many terms in the uniform prefix example. These results
identify a limit of the grouping approach; they do not count the remainder.
A generic union of two scaled three-pole curves also has only
`Theta(n^2)` initial ideals, uniformly in the component sizes.
The six-ray cut and quota proof, supported by independent point-quotient
controls, rules out multiplying these component fan counts.
Special node configurations remain outside that result; it does not
change the unrestricted bound.
For an arbitrary number of scaled curve components, an exact generic
weighted-prefix rank formula now retains only the largest
`floor(sqrt(n-1))` component capacities. It gives a matrix-free rank
oracle using classical floor sums. Balanced samples with at least
`sqrt(n)` components have exactly the generic enclosing-surface fan.
A monotone threshold charge now gives `O(n^2 log(n+1))` initial ideals
for **every** balanced component count, closing that whole generic
construction as a polynomial-amplification route. Unequal-size nested
constraints now have an exact decreasing-phase block description.
Only adjacent size drops greater than three can supply extra rank
constraints. All slowly decreasing profiles, including linear slopes
one through three, therefore obey the same `O(n^2 log(n+1))` fan bound.
An infinite three-component family has a local parallelogram from two
independent exchanges at different costs, verified in exact point ideals.
It exposes the missing simultaneous-cost mechanism in the balanced proof;
four local choices alone give no growing lower bound. Counting those
interactions and the unrestricted count remain open.
A uniform area sandwich now gives the complete leading rank law for
linear size profiles, with an explicit finite error and a curvature
threshold of four. Above that threshold, an infinite family loses a
fixed fraction of the rank and has linearly many active capacity phases.
Sufficiently steep gaps also give exact separate layer selection on
fixed carrier domains. This reduces their fan to a common refinement;
an arithmetic wall-curve count bounds its all-tight part by `O(n^2)`
for every fixed linear profile, across all axis orders and aspects.
A short-block certificate and mixed-intersection bound extend the
quadratic count to fixed pooling widths on compact weight regions.
Logarithmic coordinates give a global `O(n^2 log(n+1))` bound for every
fixed margin above the curvature threshold of a linear profile, including
genuine pooling and all positive aspect ratios. At slope four, the exact
critical harmonic curve reduces to complementary-axis rotation words
and a classical shifted-matroid selection model. An explicit infinite
family has rank defects and canonical block widths of order `sqrt(n)`;
the universal drop-three truncation threshold is therefore sharp.
The exact curve itself has only `O(n^(3/2))` strict witnesses for each
fixed last component size. A final-prefix criterion now closes every
fixed subcritical margin as well: its initials are eventually ordinary
surface prefixes, uniformly across all positive aspect ratios, giving
`O(n^2 log(n+1))` states. At slope four, every exception lies in an
explicit harmonic-ratio strip of width `O(n^(-1/4))`. The critical
family still has a `sqrt(n)` rank deficit at the final surface prefix.
Integrating the rank law now supplies a uniform smooth cost approximation.
Relative-aspect integer-profile caps and the classical lattice-vertex bound
count the whole shrinking strip, giving `O(n^2 log(n+1))` initials for
**every** decreasing profile with `|m_s-4(r-s)| <= D`, for fixed `D`.
This includes bounded irregular perturbations and unbounded block widths,
across all positive aspects. The entire class is therefore excluded as a
polynomial-amplification route. For every higher fixed linear slope, a
quadratic energy comparison also closes the entire subcritical side and
critical bands of width `O(r^(-1/2))` in `K`, where `r` counts components.
On the supercritical side a uniform proximity theorem confines each
layer's change to a constant number of boundary monomials, even when a
pooled block grows without bound. It bounds the total pooling saving by
`O(r)` in normalized cost and gives an exact compressed research algorithm
for the generic component model, with a global exchange certificate and
`O(r)` output generators. This algorithm assumes the generic rank rule;
it neither inspects arbitrary point sets nor extends the prime unranking
CLI. Existing efficient nested-allocation methods are established prior
work. The intermediate shrinking supercritical region, unbounded profile
deviations and unrestricted `G_3` remain open; sparse corrections alone
give no bound on the number of initial ideals.
A uniform truncated-line incidence bound now controls the complete
same-type candidate exchange network through the remaining transition
region, across all aspects. Mixed row-column comparisons are still
essential: a length-54 point ideal has two adjacent initials with
identical layer boundary prefixes and identical same-type boundary
orders. Their separate arithmetic network and its interactions with
the simpler walls still require a uniform count.
An exact rotation/discrepancy reduction now gives arithmetic rigidity
for long pooled blocks: for each fixed saturation threshold, blocks
reaching the `1/(c-K)` scale lie near finitely many critical aspects
in the associated quadratic field, with explicit polynomial bounds
on the resonance parameters. At a nonresonant limiting aspect their
width is `o(1/(c-K))`; a field-norm argument gives a uniform power
improvement at a specified algebraic aspect. This identifies an
arithmetic mechanism behind the remaining degeneracy; it does not
supply the missing mixed-wall fan count.
A weighted resonance cover now supplies a count of actual initials
with sufficiently large pooled blocks, using integer support caps of
one common profile polytope. For each fixed saturation fraction it
gives `O(n^2)` initials throughout `r^(-1/2) <= c-K <= r^(-1/3)`,
where `r` is the component count. An explicit tradeoff also allows
a polynomially decreasing saturation fraction on a smaller window.
Smaller growing blocks and the rest of the supercritical transition
remain uncounted; this is a restricted family result, with no new
unrestricted bound or established novelty claim.
A [common-ray reduction](docs/ray_scaling_bound.md) gives a general
explanation of another limit on amplification: a rank-`r` matroid with
positive prices on `t` common rays has at most `binomial(r+t-1,t-1)`
strict states, and direct sums share that same ceiling at their total
rank. A region has at most `t(t-1)` interior facets, giving a finite
traversal with a proved support-call bound. The proof uses classical
greedy invariance and tropical monomial labels; the independent-component
case is a published tropical hyperplane-arrangement theorem. It is
retained as a structural reduction, with no standalone novelty claim.
An actual length-six quotient has two initials with identical pure-axis
counts, identifying precisely why mixed monomials need additional data.
The proposed relative margin of order `1/n` is now refuted by
[actual all-field quotient families](docs/initial_cone_stability.md).
The sharp worst-cone margin is `Theta(n^-2)` in two variables and
`Theta(n^-3)` in three. The examples have only three and seven
initial ideals, respectively, so poor stability does not imply a
large fan. A disjoint-chain argument also bounds the entire class
of arbitrary linear socle identifications by `O(n^(3/2))` in three
variables, using the classical lattice-vertex theorem. The full
count needs an aggregate mechanism that allows individual thin cones.
A shared profile budget now handles propagating relations in finite
monomial covers as well: coordinate annihilation depths `h_i` give
`O(n^(3/2)*sqrt(h_1*h_2*h_3))` initials. For local quotients the
canonical monomial cover makes this an intrinsic test; any polynomial
separation needs depth product larger than `n` by a polynomial factor.
This necessary condition leaves deep relations and multipoint quotients
uncounted, and gives no new unrestricted exponent or novelty certificate.
For [two-part hypersurface truncations](docs/hypersurface_truncations.md),
an exact homogeneous-initial formula now retains the powers created by
propagating relations. A free-module reduction preserves every original
monomial rank; one explicit cubic has sparse commuting operators using
`O(n)` stored coefficients. A simpler monic model now works over the
integers and in every characteristic. Monic integer border certificates
give its complete positive fans at eight specified cutoffs, through
length 225. A [uniform classification](docs/cubic_fan_classification.md)
now proves the boundary cover for **every cutoff and field**. The full
diagonal Markoff cubic family has exactly `6*(1+M_s)` initials, where
`M_s` counts specified Farey fractions. This gives `3*n/pi^2 + O(sqrt(n)*log(n+1))`
initial ideals. Original monomial price-prefix ranks have an exact
`O(sqrt(n))` arithmetic query, and a standard set needs `O(n)` comparisons.
The proof matches a Laurent three-stream model with a second exact
rank formula at two equal prices, then uses cone convexity. This closes
the entire deep-relation family as a polynomial-amplification route.
The surface, freeness and cone mechanisms are established mathematics;
priority of the classification, unrestricted `G_3` and the publication
assessment remain open.
For [higher interaction order](docs/higher_interaction_fans.md), the family
`(x^2+y^2+z^2+x*y*z^r)+m^(s+1)` now has only `O(n)` initial ideals,
uniformly in `r>=2` and field, with `n=(s+1)^2`.
The integral transport theorem below improves the earlier `O(n^(9/4)/r)`
comparison and excludes the entire family as a polynomial-gain route.
Its boundary fan has an exact Farey count, but an all-field cone witness
refutes the proposed equal-price-plane cover for the full ideal.
An explicit infinite subfamily realizes `Z/m` as torsion of a coherent
price-prefix quotient for every integer `m>=2`; its standard set changes
in precisely the characteristics dividing `m` at that prefix.
The algebra remains integrally free. This arithmetic change does not
establish a larger fan: the same subfamily already has an `O(n^2)` bound.
A coefficient-independent rank law now covers coprime binary forms
`x^2+Q(y,z)+x*H(y,z)` when `y,z` have equal prices, over every field.
Matching its cutoff ranks with the streams proves the entire original
region where `z` is strictly cheapest: exactly the boundary's `2*M_(s,r)+4`
initial cones when `s>=r+2`. An original rank query there uses `O((r+1)*(s+1))` exact
arithmetic operations, and one minimum basis uses `O(n)` comparisons.
The formerly remaining regions are now covered by the uniform integral
transport theorem below. Current priority and the main goal remain open.
A [fixed-equation calculation](docs/fixed_equation_bad_primes.md) now shows
that every prime changes a coherent standard set of some truncation of
`x^2+y^2+z^2+x*y*z^3`. A length `O(p^2)` suffices for characteristic `p`,
with explicit positive weights having no relevant ties. The affected prefix
block has cyclic torsion given by classical ternary-tree coefficients;
the whole integer algebra remains free. For any fixed finite list of primes,
one cutoff makes every sampled characteristic undercount a rational prefix
rank, even if their maximum is used. Explicit residue classes and the Chinese
remainder theorem give the simultaneous witness. This does not prove additional
fan cones, an unrestricted count improvement or novelty.
The [signed-face reduction](docs/signed_face_reduction.md) gives an exact
Laurent-block model for the face obtained by deleting the cheapest variable's
square. Its grading already bounds that face's fan linearly in the quotient
length. A [noncoprime whole-tier rank theorem](docs/middle_interaction_fans.md)
now proves the transfer throughout the region where `z` has the middle price,
over every field. Each ordered middle-price region has at most `12*n+7`
original cones and an exact rank query. Together with the cheapest-`z` result,
only `O(n)` cones meet weights where `z` is not largest; generic
characteristic-dependent changes are confined to the `z`-largest region.
An [equal-degree unit rank theorem](docs/equal_degree_unit_fans.md) now
closes that region for the entire `r=2` family: its full fan has
`Theta(n)` initials, all coherent integral prefix quotients are torsion-free,
and every positive-weight rank query takes `O(s+1)` exact arithmetic/floor
operations with constant auxiliary storage. Its universal denominator is
one. This excludes `r=2` as a growth-amplification route and locates the
first exceptional-prime behavior at `r=3`. The higher-order count is now
bounded linearly as well; the unrestricted extremal goal remains open.
For the fixed degree-five
equation, an integral prepared
normal form now classifies the entire `z^2` cone in every characteristic.
The exceptional-prime witnesses all enter this same elementary cone, so their
rank changes do not demonstrate a fan-count gain.

For the remaining higher interactions, an
[exact binomial reduction](docs/higher_seed_binomial_reduction.md) computes
the two-cheap-equal seed ranks, including their integral arithmetic.
It gives a coefficient-independent formula for coprime binary forms and
explicit upper/lower rank bounds for noncoprime forms.
The [integral transport proof](docs/higher_seed_transport.md) now equates
the complete regular and canonical prefix cokernels, including all
prime-power torsion, at every cutoff and `0<u<Z`. A block recurrence and
an exact global pairing prove the identity uniformly. The rank sandwich
therefore gives the same equal-price law for all nonzero binary forms.
It proves the remaining signed-face transfer and the full-family bound
`g(I_(s,r)) <= 25*n+14`, uniformly in `r>=2` and characteristic.
Thus this entire family cannot supply polynomial amplification over the
generic scale. Exact counterexamples still exclude extending the identity
to arbitrary downward staircases or the opposite price regime.
A [binary ideal reduction](docs/three_point_rank_ideals.md) now identifies
each selected binomial block with a homogeneous component of an explicit
ideal generated by products of three fixed linear forms, preserving all
integral torsion. Its Hilbert-function identity is consequently resolved
by transport too. Existing syzygy-gap and `p`-fractal results provide
context; the proof does not require an unestablished application of them.
A coherent block has rank four over every field despite a support matching
of size five, excluding that shortcut. Exact priority of the new theorem,
the unrestricted count and the original publication goal remain open.

An [affine liaison-addition comparison](docs/liaison_fan_composition.md)
composes two three-variable quotient fans with exact length accounting.
The algebraic construction is classical; the fan is their common refinement,
so branch counts need not multiply. Ordinary box-prefix refinements stay
within `O(N^2*log(N+1)^4)`, while this gluing of two arbitrary balanced
nonlinear rectangles has only `O(N^(3/2))` initials in its true length.
These immediate composition routes therefore do not give the desired
power gain. Neither statement changes the unrestricted bound.

An [exact two-interaction calculation](docs/two_interaction_ranks.md)
shows that replacing two binary forms by one repeated form changes
complete-prefix ranks for every odd interaction exponent and every field.
The rank is determined by a binary-polynomial divisibility condition,
with a necessary distinction at the first cutoff and an explicit
characteristic-two example. The independent-form equation has full-rank
support lattice of constant index. It escapes the previous rank identity,
but no increase in the total number of initial ideals is proved.

A [power-carrier bound](docs/power_carrier_bounds.md) now closes the
first-order two-interaction construction: the square of its quadratic
part belongs to the ideal, forcing every initial basis into a bounded
coordinate slab. Classical lattice geometry gives `O(n^2)` initials
for every fixed-degree tangent equation with all pure coefficients
nonzero and perturbations beginning at a fixed fraction of the cutoff,
regardless of the number of terms or characteristic. Polynomial growth
beyond that scale requires an unbounded number of interaction layers.
An auxiliary `O(n^(8/3))` bound holds after suitable generic linear
coordinate changes; it does not bound the original-coordinate `G_3`.
These are applications of established tools with priority unassessed.

A [five-cone construction](docs/initial_cone_height.md) shows why large
necessary integer weights need not indicate a large fan: it requires
quadratic weight height but only two projection planes. The next count
must control those planes collectively.
A [weighted-facet reduction](docs/weighted_facet_incidence.md) now
isolates that aggregate question and covers all boundary cells with
three planes. A local four-cone family shows that even the large-normal
surface-area tail can be quartic. A three-pole rational curve now gives
actual finite quotients with quartic weighted incidence and only
`binomial(n+1,2)` initial ideals, disproving the proposed cubic budget.
This family also needs at least linearly many projection planes to cover
its fan, ruling out a universal polylogarithmic cover. The count needs a
different aggregate estimate. A nine-cone quotient
also retains the difficult weight and projection exponents in an entirely
interior cell.
A [nested-matroid construction](docs/nested_jet_matroids.md) supplies
explicit non-binomial finite quotients with an exact monomial rank rule,
including an infinite family of embedding dimension three. It replaces
coefficient searches by a classical scheduling oracle on a realizable
exponent lattice. The proof also gives a curvilinear family with exactly
the generic positive fan and an infinite family of non-corner-cut standard
sets. Balanced valuations now reduce exactly to bounded-height lattice
prefixes: their fixed-height count is `Theta_h(n^2 log(n+1))` for `h>=2`,
and a classical geometric bound excludes a power separation at any height.
General valuations need a different argument. A family of length `n=M^2`
has a coherent box of `n^(3/2)` points and explicit layer-rank defects.
The [chain-allocation reduction](docs/nested_chain_allocation.md) now turns
the curvilinear basis conditions into linear cumulative constraints and
gives exact initial-cone certificates. Along a valuation-price line each
quota can cease to be tight only once. The stratum where all quotas are
tight has only `Theta(sqrt(n) log(n+1))` initial sets in that large-box
family. A consecutive-block cut certificate now bounds the entire
equal-valuation family by `O(n^2 log(n+1)^4)`, uniformly in its valuation
parameter. This excludes a polynomial separation from that whole family.
A stronger [cheap-prefix theorem](docs/nested_prefix_bound.md) now gives
the same bound for the **entire** explicit nested construction, including
unequal valuations and embedding dimension three, in terms of its actual
quotient length. It also applies to generic truncated power series with
any fixed leading orders in characteristic zero. One common weight list
samples all these fans. This closes the whole nested route to a power
separation; the unrestricted spatial target and exact priority remain
unresolved. Further lower constructions need dependencies outside that
single additive valuation order.
A [rectangle rank comparison](docs/rectangle_rank_structure.md) now shows
that all mixed multiplication-image ranks stay fixed as nonlinear terms
vary, although the initial fans can differ. Those invariants therefore
miss the selected monomial spans needed for this count. For a square
rectangle of length `n=q^2`, nonlinear terms starting in degree
`D>=q` affect an exact direct summand of rank `(2q-D)^2`.
Classical lattice geometry then gives
`g(I_F)=O(1+q^(3/2)*(2q-D)^3)`, ruling out a power separation from
shallow terminal perturbations. General low-order nonlinear terms,
the unrestricted endpoint and the priority assessment remain open.
The [profile-lattice continuation](docs/profile_lattice_bound.md) now
proves `O_d(n log(n+1))` whenever the intrinsic grading lattice has rank
at most two. In three variables this covers every nonzero integer
grading, uniformly in its weights and characteristic. A classical planar
oracle method then needs `O_d(n log(n+1))` greedy calls.
Two-dimensional faces of unrestricted quotient profile polytopes also
have `O_d(n log(n+1))` vertices; this local bound does not count the faces.
Full-rank degree congruences also sharpen the rectangle volume bound: a single nonlinear
homogeneous correction of degree comparable to the side length cannot
give a power separation. Exact priority remains under assessment; the
unrestricted full-rank case is still the main target.
The [generic multiaffine continuation](docs/generic_multiaffine_bound.md)
now excludes another entire construction class: sufficiently general
points on any irreducible multiaffine surface have
`O(n^2 log(n+1))` initial ideals, with one common weight list for the class.
This is a corollary of classical convex-position geometry and generic
interpolation. Explicit reduced point quotients realize the toric model;
special configurations and local rectangle algebras can violate its rank
rule and remain outside the bound. The unrestricted endpoint stays open.
The [bilinear rectangle continuation](docs/bilinear_rectangle_fans.md)
now classifies one nonuniform class exactly. For square side `q=p^a>=2`
in characteristic `p`, `z=x+y+xy` has exactly `binomial(q+1,2)+1`
initial ideals, only linear in its length `q^2`. Unequal prime-power
sides also have an exact count. Translation and a rectangular Pascal
matrix explain the arithmetic; arbitrary column ranks do not survive
translation. Characteristic zero instead gives a structured tensor
interpolation problem and remains outside this torsion classification.
An [exact outer-cone rule](docs/rectangle_boundary_fan.md) now handles
`z=y+x(1+y)^d` in characteristic zero with one common list of at most
`2n` weights for that region. It also exhibits coherent bilinear walls
with primitive direction product `Theta(n^(3/2))`, ruling out a proposed
linear direction bound. The remaining bilinear fan reduces to fixed
rectangular interpolation on convex integer-node prefixes.
A separate [graph-rank reduction](docs/bilinear_cographic_family.md)
handles side lengths one below powers of the field characteristic.
Its entire monomial residue matroid is the dual of a complete bipartite
graph's forest matroid. Standard sets are complements of maximum spanning
trees with explicit periodic representative costs, giving an
`O(n log n)` comparison algorithm for one weight. Counting those trees
as the three monomial weights vary remains open; arbitrary graph lower
bounds cannot be imported without realizing their costs. Neither result
settles the unrestricted extremal objective.
The [ordered-tree continuation](docs/cographic_order_ideals.md) now
identifies the graph duality with its exact prior sources and proves
additional structure for equal periods. Only two fixed tree vertices
can branch, and a four-cycle certificate reduces `q^2` weighted edges
to at most `19q` candidates. This gives an `O(sqrt(n))` compressed
standard set with `O(sqrt(n) log(n))` arithmetic construction cost;
explicit output still has size `n`. A classical corank/box bound also
excludes the equal-period three-circuit extension from polynomial
separation. The equal-period phase count remains conjectural; distinct
prime periods require a different coherence argument.
The [mixed-period price study](docs/cyclotomic_price_regions.md) proves
that their smallest common coordinate box has exactly `n^2` monomials,
so the global small-box shortcut cannot close that family. A general
finite-algebra tower theorem gives an `O(n)` common weight list in each
cheap-coordinate region, over every field. Comparable periods with bounded
price ratios have only `O(n^2)` states by a classical box argument.
The intermediate price regimes remain open; no polynomial separation or
publication conclusion follows from these regional exclusions.
The [nonnegative-grading theorem](docs/nonnegative_graded_fans.md) now
gives a global `O(n)` bound whenever at least two variable degrees are
positive, over every field and independently of the grading weights.
This order is sharp for the combined class. It also gives `O(n)` local
refinements at every non-axial nonnegative weight of an arbitrary finite
quotient. In a deformed tower, additional states must reach its defining
relation; interleaving independent lower blocks alone is insufficient.
The exact novelty comparison and the unrestricted extremal problem remain open.
The same note extends the linear bound to signed gradings with all three
degrees nonzero and gives an exact prefix-rank formula for the original
mixed-period algebra. Its remaining difficulty is a sum of graph cycle
spaces, which cannot be replaced by connectivity of their union. A wider
price region restricts every basis to one period of candidate slots;
their fixed total rank deficit now gives an `O(n^2)` bound. More generally,
all initial sets staying within one period in any coordinate presentation
number only `O(n^2)`, regardless of prices or period balance. A separating
family must draw its growth from bases reaching all three period thresholds.
This does not count those remaining bases or establish novelty.
The [periodic-cover and arithmetic study](docs/periodic_quotient_arithmetic.md)
separates the quadratically bounded lattice cover from its circuit quotient.
An exact coherent basis for periods `(5,7,11)` has index two and changes
in characteristic two. Monomial substitution can amplify that index
exponentially while leaving the number of initial ideals unchanged;
torsion size alone is therefore insufficient for a counting separation.
Exact traces now verify thirteen quotient states inside a single fixed
cover cone for `(5,7,11)`, including states that cross all three period
thresholds. A coefficient-independent upper-shell rule prunes candidate
omissions. A socle-free shell theorem now bounds their common candidate
set by `O(M^(2/3) k^(1/3) log(M+1)^(2/3))`, where `M` is cover length
and `k` its quotient corank. For comparable large periods this is
`O(M^(8/9) log(M+1)^(2/3))`. The same theorem bounds low-rank basis
monomial multipliers in arbitrary finite quotient algebras; dyadic
monomial quotients show that its logarithm cannot be removed uniformly.
These bounds and finite line traces do not count the full refinement fan.
Local eight-column cube circuits provide another field-independent
pruning rule, including in characteristic two. A uniform boundary-layer
argument now bounds its survivors by three times the sum of the cover
projection areas. An infinite thin-cover family shows why a corank-only
cube bound fails even after all three period thresholds are reached.
Those first examples omit a generator. A period-gap construction now makes
all three generators active over every field, reaches all three thresholds,
and still defeats every cube test. Explicit group relations nevertheless
put its entire price region inside a quadratic fan bound. More generally,
for odd coprime periods ordered `r>t>s`, half-gaps `h=(r-s)/2` and
`e=(r-t)/2`, divided by their gcd to give `h0,e0`,
give regional bounds `O(n^2 sqrt(h*h0))` and `O(n^2 sqrt(e*e0))`.
These classical volume corollaries close additional unbounded-price regions;
they do not count the remaining joint shell/circuit refinements or establish
an unrestricted advance. Two further period identities now cover the whole
positive price domain with four common boxes. They give `O(sqrt(H)*n^(7/3))`
initial ideals for comparable periods, where `H` is the largest normalized
half-gap; bounded normalized gaps include fixed gap ratios and coprime
arithmetic progressions. This is a global restricted bound, above the generic
scale and derived from classical profile geometry. One actual cover needs
a box of size `Omega(M^(5/3))` near an explicit arithmetic ray, although its
quotient has only `O(n)` infinitesimal refinements at that ray. A sharper
global count therefore needs rank information beyond these cover boxes.
An exact score-filtration comparison now separates individual monomial
defects from collective basis excess. A coherent family using all three
generators has zero individual defect on every selected monomial but
quadratic total excess; its whole specified sector still has only
`Theta(n log n)` initial ideals, and an exact opposite-sector rule bounds
the whole explicit family by `O(n^2 log n)`. Prefix ranks and determinant valuations
detect the missing collective dependence. These are supporting identities
and obstructions, with no unrestricted counting or novelty conclusion.
An integral cycle-space reduction now extracts graph-only prefix ranks
and isolates the remaining arithmetic matrix. The coherent characteristic-two
witness reduces to the scalar `2`. A known cyclotomic construction gives
the scalar `3`, but its entire coefficient-deletion family for these
periods cannot have monomial downset representatives in the original
embedding. Arbitrary torsion therefore does not resolve the positive-weight
coherence problem. A uniform count of the residual changes remains open.
An explicit downset basis in the same algebra does have index `3`, but
two of its normal forms impose contradictory monomial comparisons.
Flat substitution gives unbounded indices and linear modular rank loss
with the same failure of coherence. Downset structure alone therefore
does not supply the order constraints needed for the count.
The exact box-residual transfer
is identified with established linkage, while every nonbox group cover
has a residual/initial mismatch measured by its number of maximal monomials.
The additive prime-power class is now classified separately: at length
`n=q^k`, every initial ideal is a prime-power box and the sharp maximum
is `binomial(k+d-1,d-1)`. The proof uses established Frobenius reduction;
its sparse universal bases and complete weight list do not establish an
unrestricted advance or a new algorithmic priority claim.

Current research returns to the original occurrence-count and prime-exponent
interleaving observations. The [face reconstruction study](docs/occurrence_reconstruction.md)
proves that face prefixes involving O(log n) coordinates determine an n-term
full prefix, and that logarithmic dimension is necessary in the worst case,
including for prime bases. Exact checks include two four-prime systems whose
proper faces agree through 1,000 terms while their full words differ at rank 10.
The [manuscript audit](results/research/interleaving_support_audit.json) records
the proof review, bounded primary-literature comparison and publication assessment.
The [goal record](results/research/occurrence_reconstruction_goal.json) preserves
the research milestones. The separate ML track is paused;
its results are preserved without a claimed ML advantage.
The [single-stream continuation](docs/occurrence_reconstruction.md#one-exponent-stream-determines-the-normalized-ordering)
also proves normalized-weight recovery from one infinite exponent stream or
its binary zero mask, an explicit finite height-error bound with calibration,
and absolute prime-basis uniqueness for three or more primes using the
classical six exponentials theorem. Arbitrarily long finite stream prefixes
can nevertheless hide a rank-6 ordering change. Priority and an efficient
recovery algorithm remain unestablished.
The [finite-recovery continuation](docs/occurrence_reconstruction.md#finite-recovery-under-a-quantitative-separation-promise)
supplies an explicit sufficient mask length under stated weight bounds and
a Diophantine separation promise, including generator identification and
all observation costs. Its bound is polynomial in numerical rank for fixed
class parameters but can be enormous. This is a finite-information result;
priority, practical recovery and an unranking advantage remain unresolved.
The [sharp dimension continuation](docs/occurrence_reconstruction.md#a-sharp-dimension-threshold-for-finite-interleaving)
determines the worst-case required face dimension up to a bounded additive
error: `(log n + 0.5 log log n) / log(3 sqrt(3)/2) + O(1)`, including for
prime bases. A uniform composition-layer moment bound explains the sharp
interaction-rank scale. Priority and publication significance remain open.
The [certificate continuation](docs/occurrence_reconstruction.md#sparse-ordering-certificates-and-the-exact-interaction-parameter)
turns the law into sparse linear certificates for realizable prefixes and
identifies the same optimal support parameter in generic point-ideal state
fans. An explicit finite quotient shows why the sharper bound does not
extend to arbitrary ideals. Research interval checks preserve the core
CLI and its independent rank auditor.
An anonymous [review manuscript](paper/interleaving_support/main.tex) and
[PDF](output/pdf/interleaving_support.pdf) now collect the sharp law,
sparse certificates and generic-state-fan consequence with complete proofs.
The [audit](results/research/interleaving_support_audit.json) states the
primary-source boundaries and final internal publication assessment. The
sharp interaction law and its consequences support a focused specialist
paper for independent review. Global priority, acceptance and
submission are not settled; no practical unranking speedup is claimed.

The [mathematical research note](docs/interleaving_research.md) targets general
structure and certification for fixed-prime rank/unrank. It separates known
signature-sequence theory, general occurrence-transform laws, a restricted
cycle classification, and an obstruction to extending finite pairwise patterns.
It also derives a theoretical generic unranking baseline from published log
separation and lattice-counting results. New work must add something beyond
that baseline. Exact research checks and witnesses are saved; novelty and prime
distribution applications remain unestablished.

The [completed quantitative goal](results/research/quantitative_goal.json)
required a new quantitative theorem, a verified literature gap and a useful
mathematical consequence. Its original target was Wagner's published
logarithmic gap in the number of three-dimensional corner cuts, the shapes
underlying smooth exponent prefixes.
The existing proofs and closest primary results have been compared at the
theorem level, including the full Vasiliev--Pavlov finite-ordering paper and
Wagner's stronger face-flag bound. The final internal assessment supports a
focused submission under the documented literature comparison and the user's
explicit final-version assumption.

The quantitative manuscript proves the matching pointwise count
C_d(n)=Theta_d(n^(d-1)(log n)^(d-1)) in every fixed dimension d>=2,
including the original three-dimensional target. It counts primitive normal
directions and uses two integer lattice patches to prevent collisions at a
prescribed size. Exact witnesses test the delicate step for 2,280 ordered pairs of
normals at five boundary heights each. The fresh internal proof audit found
no unresolved logical step. The earlier cumulative construction is preserved
in the note. The count also gives order log n for the mean number of
completions of a coordinate slice, addressing Wagner's suggested route. The
[anonymous mathematical manuscript](paper/corner_cut_growth/main.tex) gives
the complete proof, explicit constants and source comparisons separately
from the earlier structural note. It also
proves the same growth order for complete ordered prefixes and a sharp
worst-case prime-size bound: for fixed dimension and lower prime bound,
O(n)-bit primes realize every admissible n-term pattern, and some patterns
require Omega(n) bits. Known finite-ordering, representation and qualitative
prime-realization results are credited explicitly. This supplies no fixed-prime
unranking speedup or new prime-gap estimate.
The [PDF](output/pdf/corner_cut_growth.pdf) is compiled and visually checked.
The [acceptance audit](results/research/corner_cut_growth_audit.json) records
the mathematical contribution and the reading scope of each comparison.
The novelty position is a documented research assessment; no independent
referee assessment or journal acceptance is claimed.

The [next research plan](docs/interleaving_research.md#higher-impact-research-plan-universal-algebraic-complexity)
aims to explain how algebraic relations change the number of monomial initial
ideals of arbitrary finite polynomial algebras. The
[new goal record](results/research/structural_complexity_goal.json) separates
the full extremal ambition from the first literature and proof-feasibility
gates. Success requires a uniform asymptotic advance and a mathematical
consequence; examples or a roadmap alone do not complete it. This is a new
research objective, not an additional theorem of the existing manuscript.

The existing structural results are preserved as a separate note about the
geometry and containment intervals of pyramidal partitions. Applying Edelman's
region-poset theorem gives a cutting-layer formula for their Mobius function
in every dimension, addressing an observation in an Elizalde--Galvan preprint.
It also gives a subpartition-counting recurrence over removable exposed faces,
with at most linearly many contributing faces at each three-dimensional top.
This is a local structural bound, not a claimed faster overall enumeration method.
An infinite family also refutes the coordinate-slice recognition conjecture
in Galvan's 2023 thesis: triangular sections need not fit one cutting plane.
The [paper draft](paper/pyramidal_mobius.tex), exact checks of 27,100 intervals,
and 1,280 slice certificates are saved. These do not meet the new quantitative
goal. The draft is not submission-ready:
the final published version and
later literature still need comparison. The preceding integer prefix-compression
bound is retained as background; no faster rank/unrank method is claimed.

## Build

Requirements:

```text
C++17 compiler, tested with g++
Python 3.10+
Boost.Multiprecision headers, usually included with a normal g++/libboost-dev installation
```

Build the core kernels:

```bash
./build.sh
```

Run smoke tests:

```bash
python3 tests/test_smoke.py
```

## Basic usage

Unrank a 5-prime fixed smooth number:

```bash
python3 smooth_toolkit.py unrank 2,3,5,7,11 1000000000000
```

Audit a certified case:

```bash
python3 smooth_toolkit.py audit 2,3,5,7,11,13 1000000000000 --timeout 300
```

Compute an exact integer value from an exponent vector:

```bash
python3 smooth_toolkit.py value 2,3,5,7,11,13 55,126,27,54,2,52
```

## Certified examples

The current independent auditor supports prime subsets of:

```text
{2,3,5,7,11,13}
```

with `k <= 6`.

Certified result:

```text
P = (2,3,5,7,11,13)
N = 10^12
exps = [55,126,27,54,2,52]
rank_certified = true
certified_count_le = 1000000000000
```

The auditor proves that exactly `N` `P`-smooth numbers are `<=` the returned candidate.

## Important caveats

1. The fastest kernels use floating-point logarithms for search plus exact integer comparison in the final local band. Treat them as high-performance experimental kernels unless independently audited.
2. The proof-grade audit path currently stops at `k <= 6` for primes `{2,3,5,7,11,13}`.
3. Benchmark results depend on hardware, compiler, and memory pressure. Use the included scripts to regenerate local numbers.
4. This repository targets fixed-prime rank/unrank, not general `Psi(x,y)` smooth-number counting.

## Current best method selection

The wrapper uses this default policy:

```text
k = 3: 3-prime Beatty/floor-sum ranker
k = 5: layer-compressed MITM
other k >= 2: sums-only/adaptive MITM
```

This is an implementation policy, not a theorem of optimality.

## Repo status

Research prototype. The intended next step is to make all headline results independently certifiable and to compare against stronger published Cartesian-sum and fixed-dimensional lattice-counting implementations.

## License

No open-source license has been selected yet. Add a license before publishing if you want others to reuse the code.
