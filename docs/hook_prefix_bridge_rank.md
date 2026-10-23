# Actual bridge rank and the necessity of cheaper prefixes

The [hook bridge reduction](hook_bridge_reduction.md) can use the
actual rank of the bridges, rather than their possible cardinality.
For every compatible flat hook map, put
\[
 M=(y^a,y^tz^t,z^a),\quad a\ge2t,\quad t\ge2,\quad
 n=4at-2t^2,\qquad b=\operatorname{rank}H.
\]
The bridge set \(H=H_0\sqcup H_1\sqcup H_2\) and the original
residue columns are exactly those in that reduction. Over every field,
including exceptional and zero coefficient choices,
\[
                 g(I_\phi)=O\bigl(n^2(1+b)^3\bigr).       \tag{R1}
\]
This is a rank-sensitive corollary of saved greedy compression and
two-colour results. Since \(b\le3t^2\), it recovers the previous
\(O(n^2t^6)\) ceiling. It does not improve the worst-case thickness
power, the unrestricted exponent or establish current novelty.
Retain \(O(n^{7/3})\) and \(O(n^{3/2}a)\) whenever smaller.

The rank of the full span's coupling between the two arms is a
different parameter and cannot replace \(b\). An infinite actual
family below has zero full-span coupling, linear coupling in a
cheaper bridge prefix, and genuine cross-arm walls. Its entire fan
has 89 cones independently of length. Thus this is an obstruction
to a proposed encoding, not a polynomial lower separation.

## Rank-adaptive exact compression

On a convex region of fixed internal orders in \(E_Y,E_Z,H\),
greedily compress **all three** lists. Their retained bases have sizes
\(l,l,b\), where \(l=t(a-t)\). A discarded bridge is dependent on
cheaper original bridges, so the same greedy rejection argument
applies to it as to either arm. This step does not replace a bridge
by its separate arm projections.

The restricted ground now has \(2l+b\) columns and rank \(2l\).
Its dual therefore has rank \(b\). Descending dual compression
retains at most \(b\) original labels from each of the three lists.
The resulting matroid, on at most \(3b\) labels, recovers the full
original optimum on the whole region by complementation. Its labels
retain their binary first coordinates. At \(b=0\) the dual is empty
and no further decision remains.

The [two-colour theorem](two_colour_threshold_braids.md) bounds one
fixed such fan by \(O((1+b)^3)\); maxima follow by negating costs
and the first-coordinate parameter. There are \(O(n^2)\) convex
internal-order chambers by the original bridge inventory. A strict
optimality region cannot split on restriction to a convex chamber.
This proves (R1). The matrix and retained labels remain conditional
on the chamber; no common small sampler is inferred.

In particular an actual family satisfying
\(g(I_n)\ge c n^{2+\epsilon}\) needs
\[
                        1+b=\Omega(n^{\epsilon/3}).       \tag{R2}
\]
Large rank is necessary here, not sufficient.

## The relevant interface is the cheaper bridge prefix

Write the residual coordinate space as \(V_Y\oplus V_Z\).
For a subspace \(S\) spanned by any bridge prefix, define
\[
 \kappa(S)=\dim\pi_Y S+\dim\pi_Z S-\dim S
          =\dim S-\dim(S\cap V_Y)-\dim(S\cap V_Z).       \tag{R3}
\]
This is nonnegative, and it vanishes exactly when
\(S=(S\cap V_Y)\oplus(S\cap V_Z)\).
The images of \(V_Y,V_Z\) in \((V_Y\oplus V_Z)/S\) have
intersection dimension exactly \(\kappa(S)\). This follows by
subtracting their dimensions from the dimension of their sum, or
by taking the quotient of \(S\) by its two pure parts.

Consider a pure cross-arm exchange at cost \(\tau\), away from
bridge ties. Let \(S\) span precisely the original bridges of cost
less than \(\tau\). The other cheaper original columns span a pure
Y subspace and a pure Z subspace. Quotienting by either such subspace
cannot increase the intersection dimension of the two arm images:
for \(Q\subset L\), the new intersection is
\(((L\cap U)+Q)/Q\). Apply this successively to the two pure parts.
Consequently the remaining arm intersection has dimension at most
\(\kappa(S)\). In particular,
\[
 \text{a genuine nonloop-parallel cross-arm exchange implies }
                            \kappa(S)\ge1.                \tag{R4}
\]
This statement uses the actual linear relations and the actual cheaper
prefix. It does not require that a vector-space contraction be an
algebra quotient. It also does not bound the number of exchanges
across different prefixes.

As a direct consequence, if every positive-cost bridge prefix splits,
every nonvertical wall has a same-arm or bridge exchange direction.
The saved known-direction budget then gives
\(g(I)=O(n^2\log(n+1))\). This is a classical split-interface
consequence, not a sharp new theorem for general mixed hooks.

## Zero full coupling can hide linear prefix coupling

Use the exact characteristic-\(1000003\), length-32 common fixture
in [the bridge report](../results/research/hook_bridge_checks.json).
Its actual full fan has 89 cones. The twelve bridge columns span the
entire residual space \(V_Y\oplus V_Z\), with arm dimensions six.
Thus
\[
                    \kappa(\operatorname{span}H)=0.      \tag{R5}
\]
Nevertheless its two wall pieces outside the known-direction
inventory have primitive directions \((1,-7,2)\) and \((1,2,-7)\).
The checker reconstructs their original adjacent bases and verifies
(R4) using full original cheaper-column ranks.

For every integer \(s\ge1\), substitute \(y=Y^s,z=Z^s\), keeping
\(x\). The polynomial ring is free over the original ring, with
residue tags \(Y^rZ^q\), \(0\le r,q<s\). This produces an actual
hook extension with
\[
            t_s=2s,\quad a_s=5s,\quad n_s=32s^2.
\]
The full bridge matrix is a direct sum of \(s^2\) original bridge
matrices, with the same residue blocks in both arms. Hence its full
rank is \(12s^2\), its projected ranks are \(6s^2\) each, and (R5)
holds for every \(s\).

At planar weights \(w_Y=w_Z=1\), take the bridge prefix of cost
less than \(5s-1/2\). Neither outer bridge block has appeared.
The prefix consists exactly of the central columns
\[
                 (0,2s+i,2s+j),\qquad i+j\le s-1.
\]
There are \(s(s+1)/2\) of them. They have distinct residue tags,
and each is a copy of the original central head, whose projections
to both arms are nonzero. The prefix and both its projections
therefore all have rank \(s(s+1)/2\). Consequently
\[
       \kappa(S_s)=\frac{s(s+1)}2\ge\frac{n_s}{64}.       \tag{R6}
\]

Optimization in each residue block is the same original optimization
at transformed weight \((w_x,sw_Y,sw_Z)\), plus a common tag cost.
Thus the entire fan still has exactly 89 cones. Its two displayed
unknown directions become \((1,-7s,2s)\), \((1,2s,-7s)\).
They stay outside the enlarged known inventory: a difference equal
to these multiples forces identical tags in both coordinates.
This also shows why large prefix coupling alone proves no fan growth.

## Verification and priority boundary

~~~bash
python3 -B scripts/run_hook_prefix_bridge_rank_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker compares complete original cone patches with the
rank-adaptive dual, including bridge rank zero and one, rational and
small-characteristic inputs, and all 89 cones of the common fixture.
It verifies the two original unknown walls, their prefix interfaces,
and the free power blocks independently at scales one, two and three.
The uniform assertions follow from the proofs above, not extrapolation.

Greedy optimization and parametric matroids are classical; compare
[Eppstein, Section 1.2](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf).
The saved two-colour proof retains its comparison with
[Hausbrandt--Ruzika, Corollary 3.4](https://arxiv.org/html/2503.08178v2#S3).
Neither source is asserted to prove these particular hook statements.
Exact novelty, a useful aggregate over moving prefixes, the unrestricted
extremal endpoint and the final-paper requirements remain open.
