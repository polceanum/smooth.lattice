# Two monomial arm operators have uniformly quadratic hook fans

Multiple central bridges on one chain do not force a large fan.
With zero outer images, suppose each central arm is an arbitrary
polynomial in one monomial. The whole compatible hook fan then has
\(O(n^2)\) initial ideals, uniformly in the two monomials, polynomial
coefficients and field. Nonunit and zero central operators are allowed.
One common quadratic list of integer weights works for every coefficient
specialization with fixed monomial steps.

The new part of the proof keeps both arms. Their steps generate a
two-dimensional lattice, and the actual row components respect its
cosets. Only finitely many label offsets distinguish the column types.
Counting comparisons in lattice coordinates gives a bound independent
of long propagation. This is a represented direct-sum argument, not an
assumption that the ideal has a variable-induced group grading.

This closes the zero-outer single-monomial-per-arm candidate, including
growing steps smaller than the thickness. Arbitrary outer coupling and
central operators with several independent steps remain outside the
theorem. The unrestricted sharp count, a polynomial lower excess, exact
strongest-current novelty and the original final-paper criteria remain open.

## The uniform statement

Use the [original compatible hook](hook_bridge_reduction.md), with
\[
 t\ge2,\qquad a\ge2t,\qquad A=a-t,\qquad n=4at-2t^2.
\]
Take \(f_0=f_2=0\). For local long and short variables \(q,v\), let
\[
 W_Y=F_Y(q^{\nu_Y}v^{\delta_Y}),\qquad
 W_Z=F_Z(q^{\nu_Z}v^{\delta_Z})
 \quad\text{in }K[q,v]/(q^A,v^t),                          \tag{TO1}
\]
where \(\nu_i\ge0\), \(1\le\delta_i<t\), and \(F_i\in K[T]\)
are arbitrary. Then, with an absolute constant,
\[
                              g(I)=O(n^2).                 \tag{TO2}
\]
For fixed \(a,t,\nu_i,\delta_i\), one complete positive integer
weight list has \(O(n^2)\) members of height \(O(a^2)\), over all
fields and all coefficients in (TO1). The actual fan may depend on the
coefficients and characteristic; the list is common, not the fan.

There are useful stronger parameter bounds. If both steps are steep,
\(\nu_i\ge2\delta_i\), then
\[
 g(I)=O\left(\left\lceil\frac t{\delta_Y}\right\rceil^2
             \left\lceil\frac t{\delta_Z}\right\rceil^2\right).
                                                               \tag{TO3}
\]
If \(\delta_Y=\delta_Z=1\) and \(\nu_Y,\nu_Z\ge1\),
\[
             g(I)=O(t^4),\qquad
             g(I)=O(t^2)\ \text{when }\nu_Y=\nu_Z=1.      \tag{TO4}
\]
The corresponding common lists have the same cardinality orders and
height \(O(a^2)\). These bounds do not assert sharpness. The earlier
[paired binomial theorem](hook_paired_chain_fans.md) remains sharper
for its class: steps at least \(t\), binomial units and zero outer maps
have a sharp worst order \(\Theta(t^2)\).

## A finite-offset comparison lemma

Suppose a represented matroid is a direct sum of row components. Give
its columns integer labels in three variables. Suppose every pair in
one component has label difference
\[
 \kappa_r-\kappa_s+(0,m\mathbf a+n\mathbf b),
 \qquad |m|\le M,\quad |n|\le N,                          \tag{TO5}
\]
where there are \(s_0\) fixed offsets \(\kappa_r\) and fixed
integer vectors \(\mathbf a,\mathbf b\in\mathbb Z^2\).
Then at most \(s_0^2(2M+1)(2N+1)\) comparison planes determine
every component's greedy minimum and hence the entire optimum.
Their positive arrangement has
\(O(s_0^4(1+M)^2(1+N)^2)\) chambers. The statement holds over
every field, including coefficient specializations that split components
further. It needs no claim that an individual comparison is coherent.

Indeed only the order within a direct summand affects its minimum
basis. Arrange every nonzero difference in (TO5); on each open chamber
all these orders are fixed. This is the classical greedy and arrangement
argument. One can use only differences of actual potential original
labels, so their heights remain bounded by the original label box.
The lemma also applies when \(\mathbf a,\mathbf b\) are dependent;
then coincident comparisons can reduce the inventory.

## Split off unmarked chains before using lattice cosets

For each arm, partition its rectangle into chains with step
\((\nu_i,\delta_i)\). For nonprimitive steps, the key
\[
                  (\delta_i q-\nu_i j,\ j\bmod\delta_i)
\]
still identifies exactly one chain. Mark every chain containing an
input \((q_0,j_0)\) with \(q_0<t\). Every central bridge has
such inputs. All pure columns stay within one chain, for every
polynomial \(F_i\); all outer columns are zero. Thus each unmarked
chain is an actual direct summand and can be treated separately.
Its upper and lower input differences are integer multiples of its
step, with multiplier of absolute value less than
\(\lceil t/\delta_i\rceil\).

Every input or output of a marked chain has long coordinate at most
\[
 B_i=t-1+\nu_i\left\lfloor\frac{t-1}{\delta_i}\right\rfloor.
                                                               \tag{TO6}
\]
For an input \((q,j)=(q_0,j_0)+k(\nu_i,\delta_i)\),
positive \(k\) is bounded by the short coordinates; negative \(k\)
gives \(q<t\). The same bound applies to the input of a nonzero
lower column whose output belongs to a marked chain, since both lie
on that chain. This argument does not invert \(W_i\).

Discarding unmarked chains is necessary here. Coarse lattice cosets
alone can contain distant chains all along the long arm, even though
no bridge reaches them. Their internal differences are cheap to count
separately; including unrelated chains in a core comparison would lose
the finite coordinate window.

Map a Y row of local coordinates \((q,j)\) to \((q,j)\), and
a Z row of local coordinates \((q,j)\) to \((j,q)\). Set
\[
 \mathbf a=(\nu_Y,\delta_Y),\qquad
 \mathbf b=(\delta_Z,\nu_Z),\qquad
 \Lambda=\mathbb Z\mathbf a+\mathbb Z\mathbf b.             \tag{TO7}
\]
Partition the marked rows by their cosets modulo \(\Lambda\).
An upper column is supported on its one row. A central column at
input \((i,j)\) has Y support at \((i,j)+k\mathbf a\)
and Z support at \((i,j)+k\mathbf b\). Every such support lies
in the same coset. These are actual row direct sums, even when one
coarse component contains several head bridges or some coefficients vanish.

The original column labels have exactly three offset forms:
\[
 \begin{array}{c|c}
 \text{central lower}&(0,t,t)+(0,i,j)\\
 \text{upper Y}&(1,t,0)+(0,i,j)\\
 \text{upper Z}&(1,0,t)+(0,i,j).
 \end{array}                                               \tag{TO8}
\]
Consequently two labels within one marked coset component differ by
one of nine offsets plus a vector of \(\Lambda\). The offset
classes need not be equal in the quotient lattice. That is why this
proof does not silently invoke the separate
[intrinsic grading theorem](profile_lattice_bound.md).

## Steep steps cancel the long coordinate scale

Suppose \(\nu_Y\ge2\delta_Y\), \(\nu_Z\ge2\delta_Z\), and put
\(D=\nu_Y\nu_Z-\delta_Y\delta_Z>0\). For two inputs in one
marked coset, write their difference as
\((y,z)=m\mathbf a+n\mathbf b\). By (TO6),
\(|y|\le B_Y\), \(|z|\le B_Z\). The exact integer coordinates are
\[
 m=\frac{\nu_Zy-\delta_Zz}{D},\qquad
 n=\frac{-\delta_Yy+\nu_Yz}{D}.                            \tag{TO9}
\]
In particular
\[
 |m|\le\frac{\nu_ZB_Y+\delta_ZB_Z}{D}
       \le\frac{3t}{\delta_Y},\qquad
 |n|\le\frac{3t}{\delta_Z}.                               \tag{TO10}
\]
For the first estimate, use \(D\ge3\nu_Y\nu_Z/4\) and
\(B_i\le t+\nu_i t/\delta_i\). Dividing the numerator by
\(\nu_Y\nu_Z\) bounds it by
\(t/\delta_Y+2t/\nu_Y+\delta_Zt/(\nu_Y\nu_Z)\),
at most \(9t/(4\delta_Y)\). The other estimate is symmetric.
Unmarked-chain comparisons satisfy the same bounds with one coordinate zero.

Apply (TO5) with three offsets. There are
\(O(\lceil t/\delta_Y\rceil\lceil t/\delta_Z\rceil)\)
planes, proving (TO3). The large long steps and the determinant cancel;
neither a large interface rank nor a count of support components is
being used as a lower bound.

## Unit short steps and the parallel case

When \(\delta_Y=\delta_Z=1\) and \(\nu_Y\nu_Z>1\), the same
coordinate formula works without both long steps being at least two.
The bounds \(B_i\le(\nu_i+1)t\) give
\[
 |m|\le t\frac{\nu_Y\nu_Z+2\nu_Z+1}{\nu_Y\nu_Z-1}
       \le7t,\qquad |n|\le7t.                            \tag{TO11}
\]
The largest ratio occurs at a smallest positive product: the required
inequality is \(\nu_Z(3\nu_Y-1)\ge4\), true for positive
integers with product greater than one. The symmetric estimate follows
in the same way. Thus (TO5) gives \(O(t^2)\) planes and \(O(t^4)\)
chambers.

If \(\nu_Y=\nu_Z=1\), both steps are \((1,1)\). Marked components
are diagonal cosets, and their input differences are \(k(1,1)\)
with \(|k|\le2t\). Unmarked comparisons have \(|k|<t\).
Only \(O(t)\) planes occur, proving the \(O(t^2)\) part of (TO4).
This permits arbitrary polynomials in \(qv\) on both arms, not just
binomials or invertible operators.

## Finish the all-step quadratic bound and the common list

If both long-to-short ratios are at least two, (TO3) is \(O(t^4)\),
which is \(O(n^2)\) because \(n\ge6t^2\). Otherwise choose an arm
with \(\nu_i/\delta_i<2\). Its earlier
[one-arm monomial-orbit theorem](hook_monomial_orbits.md) has
\(\beta\le1+\nu_i/\delta_i<3\), or \(\beta=1\) for
\(\nu_i=0\), and already gives \(O(n^2)\) with the other arm
arbitrary. These cases prove (TO2) with an absolute constant.

All partitions and potential original labels above depend only on the
steps, \(a,t\), not on the coefficients or field. A complete arrangement
of their comparison planes therefore gives the common list. Original
labels lie in \([0,1]\times[0,2a)\times[0,2a)\), so normals have
height \(O(a)\). Primitive cross products and sums of three independent
extreme rays give positive integer interior weights of height \(O(a^2)\).
No optimal construction-time or total rank-evaluation complexity is claimed.

## Verification and the open boundary

~~~bash
python3 -B scripts/run_hook_two_orbit_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker reconstructs actual compatible quotients, verifies original
row supports in the fixed coset/chain partition, and checks exact integer
lattice coordinates for potential comparisons. Independent field greedy
reconstructs the whole normal cone from the component minima. Complete
original fans test the common inventory; one shared integer list is reused
across different coefficients, fields and nonunit polynomials. The fixtures
include head columns sharing rows, so the one-bridge-per-component
argument is not being reused outside its hypotheses. The largest fixture
checks the original residue matrix without the full commuting-multiplication
replay; the other fixtures include that replay.

There are also actual coherent guards for the open boundary. Add a
nearby step on both arms, taking
\(W_i=1+q^\nu v+q^{\nu+1}v\), still with zero outer maps.
The original fixtures \((t,a,\nu,p)=(3,9,2,101)\) and
\((4,13,3,3)\) have 76 and 154 cones. Their reports give coherent
exchanges between different old row components, with nonzero first
normal coordinate and supporting planes absent from the entire
two-orbit inventory. This verifies a failure at actual faces, not just
a loss of support sparsity. These finite counts establish no asymptotic
lower bound or polynomial excess.

Two further certificates reach the pure cross-arm remainder in the
[original bridge reduction](hook_bridge_reduction.md), beyond its
known-direction inventory. Write \(M=q^t(1+q)v\), work over
\(\mathbb F_{1009}\), and keep both outer images zero:

| \(t,a\) | \(W_Y,W_Z\) | Exchanged original labels | Positive weight |
|---|---|---|---|
| \(4,25\) | \(1+M,1+M\) | \((1,11,1),(0,7,9)\) | \((35/48,23/192,29/192)\) |
| \(6,49\) | \((1+M)^5,\sum_{j=0}^{5}M^j\) | \((1,13,5),(0,7,17)\) | \((9/16,25/96,17/96)\) |

The primitive planes are \((1,4,-8)\) and \((1,6,-12)\).
In each case the columns are nonzero and parallel modulo the **entire**
strictly cheaper original span. Extending either choice to a basis
therefore gives two minimum bases differing by this pure cross-arm pair,
even if other columns have the same cost. The checker reconstructs the
operators by truncated polynomial multiplication and replays the original
quotient. These cases locate a relevant remaining mechanism; two walls
are not an asymptotic separation.

The lattice, direct-sum, greedy and arrangement ingredients are classical;
the exchange/face dictionary is in
[Feichtner--Sturmfels, Section 2](https://arxiv.org/pdf/math/0411260v1).
Exact strongest-current priority of the combined quantitative class result
is unestablished. Outer images can connect different cosets and unmarked
chains. A second nearby central step can also destroy the index constraint.
The unrestricted sharp \(G_3\), an actual polynomial fan excess and the
original research endpoint remain open.
