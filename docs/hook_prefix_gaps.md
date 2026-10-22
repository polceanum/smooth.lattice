# Prefix gaps give a quantitative bound on degenerate hook fans

Let \(t\ge2\), \(a\ge2t\), and
\[
 M=(y^a,y^tz^t,z^a),\quad C=K[y,z]/M,\quad n=4at-2t^2.
\]
Take every compatible flat square-zero extension, over any field, whose
central image has nonzero coefficients of \(y^t,z^t\). The two outer
images and all other central coefficients are arbitrary. Unlike the
[prefix-regular theorem](hook_prefix_uncertainty.md), singular prefix
minors are allowed.

Put \(A=a-t\), \(L=a-2t\). Use the labelled prefix matrices
\(\Delta_Y(r),\Delta_Z(r)\) of (CP1)--(CP2), for \(t\le r\le A\).
For each arm let \(h_Y,h_Z\) be the longest run of consecutive indices
where its determinant vanishes. The endpoints \(t,A\) are regular:
the first determinant is a power of its leading unit coefficient, and
the last is the empty determinant. When \(L=0\), set both gaps to zero.
Every index has a regular index at or above it, at distance at most the
corresponding gap. Define
\[
                         T=t+h_Y+h_Z.                     \tag{PG1}
\]

**Theorem.** Uniformly in these coefficients and fields,
\[
\begin{split}
 D_x(I)&=O(aT(1+\log(a+1))),\\
 g(I)&=O(naT(1+\log(n+1))^2)\\
     &=O\!\left(n^2\left(1+\frac{h_Y+h_Z}{t}\right)
                       (1+\log(n+1))^2\right).           \tag{PG2}
\end{split}
\]
The previous box ceiling is retained whenever smaller. In particular
\(h_Y+h_Z=O(t\log^c(n+1))\) excludes polynomial excess above the
quadratic scale. Conversely \(g(I)\ge c_0n^{2+\epsilon}\) requires
\(h_Y+h_Z=\Omega_{c_0,\epsilon}(t n^\epsilon/\log^2(n+1))\)
at sufficiently large lengths. Large gaps are a necessary gate, not a
lower construction. Vanished leading coefficients and large-gap strata
remain part of the original unrestricted problem.

## Degree uncertainty loses at most one gap

Write \(f_1=y^tU+z^tV\) in the two rectangular row modules
\(R_Y=K[y,z]/(y^A,z^t)\), \(R_Z=K[y,z]/(y^t,z^A)\).
Both multipliers are units. The argument of (PU2) now gives
\[
 \deg_y H+\deg_y(UH)\ge L-1-h_Y,\qquad H\ne0,            \tag{PG3}
\]
and the transposed bound with \(h_Z\).

Indeed, if the degree sum \(D+E\) were smaller, neither representative
would be truncated by multiplying by \(y^t\). Start at
\(r_0=E+t+1\) and take the next regular \(r\). It satisfies
\(r_0\le r\le r_0+h_Y<A\). The nonzero vector \(y^tH\)
still lies in the prefix source rectangle \(t\le y<a-r\), and its
product lies in rows \(y<r\). It contradicts the invertibility of
that regular prefix minor. If the lower bound in (PG3) is negative,
the assertion is immediate from nonnegative degrees.

The maximum-gap loss is sharp for this degree statement. In characteristic
zero or \(p>2L\), choose axis units \((1+y)^{K_Y},(1+z)^{K_Z}\),
where \(0\le K_Y,K_Z\le L\). The axis Toeplitz determinant has
nonzero indices \(r=t,\ldots,t+K_Y\), together with \(A\), so
\[
                  h_Y=\max(0,L-K_Y-1).                   \tag{PG4}
\]
The same holds on the other arm. For \(K_Y\le L-1\), the vector
\(H=z^{t-1}\) attains the lower bound in (PG3), even with arbitrary
off-axis coefficients: they vanish after multiplication by this vector.
For \(d=r-t\le K_Y\), \(s=L-d\), the classical binomial Toeplitz
product is
\[
 \prod_{j=0}^{s-1}\frac{(K_Y+j)!\,j!}
                         {(d+j)!\,(K_Y-d+j)!},
\]
which stays nonzero in the stated characteristics. For \(r-t>K_Y\)
the nonempty matrix is strictly triangular with zero diagonal. This verifies that the gap parameter has an actual
degenerate realization, rather than being only a formal allowance.

## Move each cheap prefix to the next regular one

If \(a\le64T\), there are only \(O(a^2)=O(aT)\) possible
cross-layer directions in the original label box, proving the first
part of (PG2). Hence assume \(a>64T\).

First an outer pure \(Y\)-source exchanged with an upper \(Z\)-column
has transverse \(Z\)-displacement at most \(8T\). To see this, suppose
\(c>8T\), and normalize \(w=(\lambda,u,1)\). The central \(Z\)-tail
span used in (CP8) either already makes the upper column a loop, or gives
\(u<2\). Positivity gives \(u>c/(2a)>4T/a\). Start at
\(r_0=\lfloor2t/u\rfloor\), and move to the next regular \(Y\)
prefix \(r\le r_0+h_Y\). Here \(t\le r_0<A\). Increasing the
prefix index keeps its plain columns cheaper. Its upper-column cost
increases by at most \(h_Yu\). The gap from the exchanged cost is at least
\[
                     S-5t-2h_Y+1>S-5T+1>0.
\]
Thus the whole regular basis is cheaper, a contradiction. Transpose for
the other outer source.

Enlarge the bridge note's \(O(at)\) inventory by
\[
 \{(1,d,e):|d|,|e|<2a,\quad |de|\le64aT\}.              \tag{PG5}
\]
It has \(O(aT(1+\log(a+1)))\) elements by the signed divisor sum.
All outer-source exchanges lie there by the preceding bound. A direction
outside this inventory must therefore exchange a central lower column
with an opposite-arm upper column. Write its labels as in (PU4), with
\(b=Y-R,c=S-Z\). Then \(bc>64aT\), hence \(b,c>64T\).

The common cost is still below both outer heads:
\[
                  \tau=Yu+Z<a\min(u,1).                  \tag{PG6}
\]
For \(\tau\ge au\), the opposite tail either makes the upper column
a loop, or gives \(u<2\). Start at \(r_0=t+\lfloor2t/u\rfloor\)
and move to the next regular \(Y\) prefix. Positivity gives
\(u>64T/a\), so \(r_0+h_Y<T+a/32<A\). Plain columns are
strictly cheaper, and the upper gap is greater than
\(S-7t-2h_Y+1\ge S-7T+1>0\), a contradiction.

For \(\tau\ge a\), the opposite tail either suffices or gives
\(u<a/(b-2t)<4a/(3b)\). Start at \(r_0=t+\lfloor2tu\rfloor\)
and move to the next regular \(Z\) prefix. Here
\(r_0+h_Z<T+a/31<A\). The plain-column gap is positive; the upper
gap is at least
\[
                    S-3tu-2t-h_Z+1>0,
\]
because \(3tu<4at/b<c/16\) and \(2t+h_Z\le2T\).
This proves (PG6). It uses actual full-rank bases at regular indices,
not ranks at singular ones.

All cheaper columns at such an exchanged cost, and the two exchanged
vectors, are consequently independent of both outer images. The literal
rank test and closure argument in (CS8)--(CS10) therefore transfer every
remaining geometric wall locus and distinct direction to the central-only
extension. This step does not identify parallel subdivisions.

## The four inequalities with gap terms

Use the same actual central-only exchange relation and notation
\(F=(Y+R)u+S+Z\) as in (PU6)--(PU8). Applying (PG3) gives
\[
 F>(a-1-h_Y)u+t,\qquad F>a+tu-1-h_Z.                     \tag{PG7}
\]
In the two cheaper-basis arguments of (PU10), move the chosen index to
the next regular index. The first upper-prefix gap loses at most
\(h_Yu\), and the second at most \(h_Z\). Thus
\[
\begin{split}
 F&\le au+(2t-1+h_Y)u+3t-2,\\
 F&\le a+(3t-2)u+2t-1+h_Z.                               \tag{PG8}
\end{split}
\]
The starting indices are in \([t,A]\) as before; the next regular index
always exists within its stated gap. An empty plain rectangle at \(A\)
causes no exception.

Eliminating \(F\) yields
\[
 \frac{a-3t+1-h_Z}{a+t-1+h_Y}<u<
 \frac{a+t-1+h_Z}{a-3t+1-h_Y}.                            \tag{PG9}
\]
Since \(a>64T\), this is contained in
\[
 \frac{a-3T}{a+T}<u<\frac{a+T}{a-3T}<\frac{65}{61}.
\]
It follows that \(a|u-1|<5T\) and \(a-T-1<F<a+6T\).
As \(F-(b+c)=b(u-1)+2Ru+2Z\), its absolute value is less than
\(12T\). Every remaining direction therefore satisfies
\[
                       |b+c-a|<20T.                       \tag{PG10}
\]
There are \(O(aT)\) pairs in this diagonal strip. Adding (PG5) gives
the direction bound in (PG2), and the capped-coordinate direction
theorem gives its fan bound, retaining every parallel piece.

## A longest-gap family has only five cones

Gap size does not supply a lower separation. Set
\(f_0=f_2=0\), \(f_1=y^t+z^t\), and \(a>2t\). This has
\(h_Y=h_Z=L-1\) when \(L\ge1\), but
\[
                              g(I)=5                     \tag{PG11}
\]
over every field. Its degree uncertainty loss is maximal and sharp.

After the fixed coloops are contracted, each pure central \(Y\)-column
is parallel to one upper \(Y\)-column, and similarly on \(Z\).
There are \(t(a-2t)\) such rank-one pairs on each arm. Each of the
\(t^2\) central bridges gives a disjoint rank-two triangle, with
its two upper identity columns and their sum. The outer columns are loops.
After subtracting the common shift cost, the three costs in every triangle
are \(t(u+1),\lambda+tu,\lambda+t\). The pure pair comparisons
are \(\lambda=t\), \(\lambda=tu\). Their combined actual fan has
exactly the five regions
\[
\begin{array}{c|c}
1&\lambda<\min(t,tu)\\
2&u<1,\quad tu<\lambda<t\\
3&u<1,\quad\lambda>t\\
4&u>1,\quad t<\lambda<tu\\
5&u>1,\quad\lambda>tu.
\end{array}
\]
All are nonempty, and the direct-sum minimum basis is different in
each. This proves (PG11). It rules out treating a large gap, by itself,
as evidence of many initial ideals.

## Verification and scope

```bash
python3 -B scripts/run_hook_prefix_gap_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker uses actual quotient matrices to locate regular indices,
checks the predicted binomial-axis gaps with arbitrary off-axis terms,
and verifies (PG3) for entire degree-bounded subspaces. Exact rational
scalar controls check the shifted-prefix cost bounds and enclosure of
the feasible intervals. Full original-field normal cones independently
verify the five-cone longest-gap controls.

The prefix, degree and direction arguments extend the linked proofs;
normal-index gaps and Toeplitz determinant blocks are classical Padé
structure; see [Jedynak--Gilewicz, Section 1](https://doi.org/10.1155/2013/185648).
The binomial product retains the source in the linked prefix-minor note.
Exact priority of the combined quantitative hook bound is
unassessed. The theorem leaves long gaps and zero-leading strata without
a near-quadratic bound, and it does not settle unrestricted \(G_3\),
produce an actual polynomial separation, or complete the original goal.
