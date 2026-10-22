# Counting a capped-coordinate fan by its cross-layer directions

Let \(I\subset K[x,y,z]\) have positive colength \(n\), over any
field, and contain a monic univariate polynomial of degree \(h\) in
\(x\). Let \(D_x(I)\) be the number of distinct undirected primitive
directions of actual positive walls with nonzero \(x\)-coordinate.
Parallel wall pieces are counted only once in \(D_x\). Then
\[
             g(I)=O\!\left(hn(1+\log(n+1))(D_x(I)+1)\right).       \tag{1}
\]
One may replace \(h\) by \(\min(h,n)\). This bound uses the actual
induced wall directions; the weights of the input coefficients do
not suffice, even in the certified full-parameter generic hook model.

At fixed \(h\), a bound \(D_x(I)=O(n\log^c(n+1))\) would close
that class to a polynomial separation above the generic corner-cut
scale. Conversely, \(g(I)\ge c_0n^{2+\epsilon}\) requires
\[
 D_x(I)+1=\Omega_{c_0,h}\!\left(
                   \frac{n^{1+\epsilon}}{1+\log(n+1)}\right).     \tag{2}
\]
There is a characteristic-zero family with \(h=2\),
\(D_x=\Theta(n)\) and \(g=\Theta(nD_x)=\Theta(n^2)\).
Thus (1) is sharp up to a logarithmic factor along that family.

These are a reduction and a necessary condition, not a new
unconditional exponent for the capped class or unrestricted \(G_3\).
Retain the stronger saved \(O_h(n^{7/3})\) capped-class upper bound
whenever it is smaller. Bounding the actual direction inventory is
still unproved in the regime relevant to a polynomial separation.
No novelty or completion of the original extremal goal is asserted.

## A finite polytope sweep inequality

Let \(P\subset\mathbb R^3\) be any polytope. Write \(g_+(P)\) for
its number of vertices strictly minimizing a positive linear
objective. Let \(s_x(P)\) count the vertices of its projection to
the \((y,z)\)-plane. Let \(E_x^+(P)\) count its positive exposed
edge pieces whose endpoint \(x\)-coordinates differ. Then
\[
                         g_+(P)\le s_x(P)+E_x^+(P).              \tag{3}
\]
Flat polytopes are allowed, and a point has one projected vertex.
No lattice or matroid assumption is needed for (3).

Normalize positive objectives to \((t,u,1)\), with \(t,u>0\).
Every strict vertex region is a convex polygonal region in this
quadrant. Walls from edges with zero \(x\)-displacement have
equations independent of \(t\), hence are vertical in \((u,t)\).

Consider a region having a nonvertical lower boundary edge. Charge
it to one such edge. At a relative interior point of that edge,
exactly two vertex regions meet, and precisely one lies above it.
Consequently no edge receives two charges.

For a region with no nonvertical lower boundary, decreasing \(t\)
at a fixed interior \(u\) remains in the region down to zero. To
see this directly, write its inequalities as
\((v'-v)\cdot(t,u,1)>0\). A lower boundary would have positive
coefficient of \(t\); upper boundaries become easier when \(t\)
decreases, while vertical boundaries are unchanged. Thus such a
region reaches an interval on \(t=0\), rather than just a point.
Its projected vertex minimizes the projected objective there.

In a fiber with the same \((y,z)\)-coordinates, only the endpoint
of smallest \(x\)-coordinate can minimize at \(t>0\). Two of the
remaining regions therefore cannot map to the same projected
vertex. They number at most \(s_x(P)\). This proves (3), including
degeneracies and cells touching only a point of the quadrant boundary.

This is an elementary sweep argument. It makes no general assertion
that a three-dimensional polytope has quadratically many vertices in
its largest shadow.

## Apply the existing layered shadow bound

Use the actual residue matroid on the universal ground
\[
 E=\{(a,b,c)\in\mathbb N_0^3:a<h,
                         (a+1)(b+1)(c+1)\le n\}.
\]
It has rank \(n\) and contains every initial standard set. Its
basis-profile polytope \(P\) retains exactly the relevant strict
positive minima. This is the established state-polytope dictionary.

Every actual edge direction is parallel to an exchange of two ground
labels, so a primitive cross-layer direction satisfies
\(0<|q_x|<h\). The
[layered shadow theorem, (LS8)--(LS9)](layer_shadow_budget.md)
gives, uniformly in the other two entries,
\[
 s_q(P),\ p_I(q)=O(hn(1+\log(n+1))).
\]
The same theorem with projection kernel \((1,0,0)\) bounds
\(s_x(P)\) by this quantity. There are exactly \(D_x(I)\)
cross-layer direction classes, so summing their actual parallel
pieces and applying (3) proves (1).

The full polytope may have additional nonpositive vertices. They
are used only to upper-bound the projected shadow and parallel
fibers, as in the existing shadow proof. The sweep counts the
positive fan, retaining all its parallel pieces.

## The order in the direction budget cannot be lowered

The [sharp zero-product family](quadric_cross_fans.md#a-family-attaining-both-powers)
has, for every \(L\ge1\),
\[
 n=4L+2,\quad h=2,\quad g=(L+1)^2.
\]
Its distinct directions are precisely
\[
 (1,-(2L-2b-1),0),\quad
 (1,0,-(2L-2b-1)),\qquad 0\le b<L.
\]
Hence \(D_x=2L\), and each direction has \(L+1\) parallel
pieces. These are actual characteristic-zero quotients, not merely
admissible lattice normals. The existing Hermite construction proves
the assertions uniformly; the present use is a substitution into (1).

## Input arrow weights miss induced walls

For the full independent-parameter hook at thickness two, let
\(G=((5,0),(2,2),(0,5))\). An original coefficient mapping a head
\(G_i\) to a tail \(v\) contributes comparison normal
\[
                              (-1,G_i-v).                       \tag{4}
\]
The entire generic positive fan was already certified over every
characteristic in the [generic-hook note](generic_hook_fans.md).
It has an actual cross-layer wall of primitive direction
\[
                              (-1,5,-1).                        \tag{5}
\]
No allowed parameter in (4) has this weight. For the first head it
would require the forbidden tail \((0,1)\); for either other head
it would require a negative first coordinate. Alternating paths
through selected columns induce (5).

The accompanying checker recovers both adjacent complete support
cones, a strict wall-interior weight, and their endpoint minors. It
replays the existing thickness-two primitive-minor certificate in
every characteristic. Separately, one common integer specialization
makes both endpoint minors nonzero; original quotient border cones
verify the same actual wall over \(\mathbb Q\). Thus (5) is not an
arbitrary element-pair comparison or a wall of an unavailable sparse
support basis. This example disproves the proposed replacement of
\(D_x\) by the set (4); it does not establish an asymptotic excess
in the number of induced directions.

## Reproduction and priority boundary

```bash
python3 -B scripts/run_cross_layer_direction_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker independently enumerates positive normal cells of finite
profile sets and their projected hulls, checks the injective lower-edge
charge in (3), and uses original commuting quotient matrices for the
sharp products and induced wall. Finite tests support the sweep;
the uniform conclusion uses its proof and the existing shadow theorem.
No complete larger-thickness generic fan is claimed.

The state dictionary comes from
[Babson--Onn--Thomas](https://arxiv.org/pdf/math/0207135).
The rational/layered shadow estimates, convex-token inputs, and
exact-priority limitations are recorded in their linked proofs.
Minimum-basis exchange and the planar sweep are classical geometric
and matroid arguments. Exact priority of the combined direction
gate remains unassessed. Neither a targeted literature search nor
the absence of zero determinants in sampled higher-thickness bases
establishes the uniform generic rank theorem or the unrestricted
extremal endpoint.
