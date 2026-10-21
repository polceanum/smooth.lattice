# Algebraic wall thresholds versus free lattice sorting

The [primitive-wall theorem](primitive_wall_spectrum.md) identifies the
minimum algebraic colength for a primitive mixed direction $q$ as

$$
 \nu(q)=A(q)+B(q)-2,\qquad
 A(q)=\prod_i(q_i^++1),\quad B(q)=\prod_i(q_i^-+1).
$$

There is a sharp-order comparison with the original free exponent lattice.
It concerns when an individual comparison can first appear. It does not
count walls coexisting in one quotient, improve unrestricted $G_3$, or
establish novelty. The socle construction and the usual corner-cut/state
correspondence are established inputs to the two sides of the comparison.

## The first free corner-cut wall

Fix $d\ge2$. Let $\mu_d(q)$ be the smallest integer $r$ for which the
positive fan of ordinary size-$r$ corner cuts has a codimension-one wall
of primitive direction $q$, up to sign. A corner cut is the set of the
first $r$ exponent vectors under a strict positive linear cost. Thus its
columns are all independent, as in a sufficiently generic point ideal.

**Theorem.** For every primitive mixed $q$,

$$
 2^{-d}A(q)B(q)<\mu_d(q)\le(2d)^d A(q)B(q).                \tag{F1}
$$

The constants depend only on $d$. The statement includes zero coordinates.

For the lower bound, take a relative interior weight $w>0$ on such a
wall. Exchanging monomials differ by $kq$, with $k\ne0$. Their common
cost is at least $t=w\cdot q^+=w\cdot q^-$: one exponent contains
$q^+$ as a divisor, the other contains $q^-$. Every integer point in

$$
 0\le a_i\le\lfloor |q_i|/2\rfloor
$$

has cost strictly below $t$. Indeed $w\cdot|q|=2t$, and at least one
coordinate of the primitive $q$ is odd. The box has
$\prod_i(\lfloor |q_i|/2\rfloor+1)\ge2^{-d}A(q)B(q)$ points.
Every one precedes the selection threshold, so the cut has strictly
more points. This argument also allows several collinear monomials to
tie at that threshold.

For the upper bound put $u=q^+$, $v=q^-$, and let $p,s$ count their
nonzero coordinates. Set

$$
 w_i^0=\begin{cases}
  1/(p u_i),&u_i>0,\\
  1/(s v_i),&v_i>0,\\
  4,&q_i=0.
 \end{cases}
$$

Then $w^0\cdot u=w^0\cdot v=1$. In a sufficiently small relative
neighborhood of $w^0$ inside $q^\perp$, every price changes by less
than one quarter of its value. The common threshold lies between
$3/4$ and $5/4$. Every monomial of cost at most that threshold has
zero coordinates off the support of $q$ and
$a_i\le(5/3)d|q_i|\le2d|q_i|$ on its support.

Only finitely many exponents occur in this box. Choose the weight off
every extra tie hyperplane $w\cdot(a-u)=0$ for which $a-u$ is not a
multiple of $q$. Such a choice exists: each of these hyperplanes is a
proper subspace of $q^\perp$. At the common threshold the only remaining
exponents are $u,v$. Primitivity makes every integral vector in the line
$\mathbb Rq$ an integral multiple of $q$; nonnegativity of $u+kq$
allows only $k=0,-1$.

All strictly cheaper monomials together with either endpoint give adjacent
corner cuts. Their size is at most

$$
 \prod_{q_i\ne0}(2d|q_i|+1)\le(2d)^d A(q)B(q).
$$

A small perturbation to either side preserves all other threshold
comparisons, proving (F1).

## A sharp quadratic difference in first-occurrence length

Since $A B\le(A+B)^2/4$, (F1) implies

$$
 \max_{\nu(q)\le n}\mu_d(q)=\Theta_d(n^2).                 \tag{F2}
$$

For the lower bound choose $q=(k,-k-1,0,\ldots,0)$, for which
$\nu(q)=2k+1$ and $A B=(k+1)(k+2)$. This works for every fixed
$d\ge2$. The set here consists of directions occurring at algebraic
length at most $n$; for $d\ge3$ the padding theorem makes it the
exact length-$n$ spectrum as well.

In three variables a full-support family is

$$
 q_k=(k^2,-k,-k-1),\qquad
 \nu(q_k)=2k^2+3k+1,\qquad \mu_3(q_k)=\Theta(k^4).           \tag{F3}
$$

The coordinates have gcd one because $k$ and $k+1$ are coprime. The
minimal algebraic realization over every field has just two initial
ideals. Thus relations can make an individual comparison appear at a
quadratically smaller length than free lattice sorting would permit.
This explains a concrete effect of algebraic dependencies on exponent
orderings. It supplies no polynomial separation in the *number* of
orderings and does not resolve the active extremal objective.

The minimum in $\mu_d$ varies the weight. It does not bound the number
of free monomials preceding a comparison at every weight. In particular,
a coordinate absent from $q$ may have arbitrarily small price. No uniform
free-prefix enumeration algorithm or universal sampling lower bound
follows from (F1)--(F3).

## A collective wall budget sufficient for an improved bound

For an actual colength-$n$ ideal in three variables, let $\mathcal E(I)$
be the undirected adjacencies of its maximal positive initial cones.
Each shared wall has primitive direction $q_e$, counted once per
adjacency, even when different walls have the same direction. Define

$$
 V(I)=\sum_{e\in\mathcal E(I)}\nu(q_e),\qquad
 L=1+\log(n+1).
$$

The direction threshold gives $2\le\nu(q_e)\le n$. Combining it
with the existing rational-shadow theorem yields, for every $T\ge1$,

$$
 g(I)=O\bigl(1+nT^{7/3}L^2+V(I)/T\bigr).                  \tag{F4}
$$

Here is the complete reduction. A wall of direction $q$ is dual to
an exposed edge of the basis-profile polytope parallel to $q$. Under
projection along $q$, that whole edge projects to a shadow vertex.
Distinct such edges project to different vertices: the fiber over one
projected vertex is one interval, hence one face. Consequently the
number of adjacencies of direction $q$ is at most the full shadow
vertex count, bounded in the [rational-shadow note](rational_shadow_bound.md)
by $O(n(1+\|q\|_1^{1/3})L)$. The threshold implies
$\|q\|_1\le\nu(q)$. The exact spectrum has
$O(T^2\log(T+1))$ directions of cost at most $T$. For $T\le n$
their total shadow envelope is therefore $O(nT^{7/3}L^2)$.
When $T>n$, use the same bound at $n$, which is smaller.
There are at most $V(I)/T$ other adjacencies. Finally the positive
cell-adjacency graph is connected, so $g(I)\le1+|\mathcal E(I)|$.
This proves (F4), including degenerate profile polytopes and a single
initial ideal.

Optimizing the threshold gives the unconditional parameterized statement

$$
 g(I)=O\bigl(1+nL^2+n^{3/10}V(I)^{7/10}L^{3/5}\bigr).       \tag{F5}
$$

In particular, a uniform bound $V(I)=O(n^3L^C)$ would imply

$$
 G_3(n)=O\bigl(n^{12/5}L^{(7C+6)/10}\bigr).                \tag{F6}
$$

**The cubic bound on $V$ is unproved.** Equations (F4)--(F5) retain
the unknown collective quantity; (F6) is conditional and is not a new
unrestricted exponent. Exact priority of the reduction is unassessed.

The cubic scale is necessary if such a bound is true. In the actual
[three-pole family](ray_scaling_bound.md#exact-consequence-for-axis-carrier-quotient-algebras),
all allocations of $n-1$ nonconstant pure powers occur. An $x/y$ wall
exchanges $x^A$ with $y^B$, where $A,B\ge1$, $A+B\le n$ and
$C=n-A-B$ powers are selected on the third axis. Its primitive cost is
$(A+B)/\gcd(A,B)$. Therefore

$$
 V(I_n)=3\sum_{A,B\ge1,\ A+B\le n}
                   \frac{A+B}{\gcd(A,B)}=\Theta(n^3).       \tag{F7}
$$

The upper bound is immediate. For the lower bound restrict $A,B$ to
$[n/4,n/3]$ and require coprimality. Elementary Möbius inversion gives
$\Theta(n^2)$ such pairs, each with cost at least $n/2$.
Thus the earlier isolated wall construction cannot test sufficiency
of the cubic budget; a genuine proof must retain all simultaneous
dependencies and their shared quotient length.

## Verification and scope

```bash
python3 -B scripts/run_free_wall_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker constructs rational weights in $q^\perp$, enumerates all
monomials through the common threshold, and verifies that only the two
endpoints tie there. It replays the adjacent first-$r$ sets on both sides
by separate exact sorting, checks their downset property, and compares
with the universal half-box lower bound and divisor-box colength.
Both two-support and full-support quadratic families are included.
Finite witnesses check the construction; the theorem uses the proof above.
Additional controls extract all adjacencies from the saved exact fans,
check their primitive costs and the threshold partition in (F4), and
compare the three-pole sum with the independently generated allocation
adjacencies. They do not verify the unproved cubic budget for arbitrary
ideals or infer (F6) from finite data.

The [primitive-wall note](primitive_wall_spectrum.md) records the
edge-ideal and socle-move comparison with Altmann--Sturmfels. The
[universal state framework of Babson--Onn--Thomas, Sections 2 and 5](https://arxiv.org/pdf/math/0207135)
provides the established generic-point/corner-cut dictionary. Those
primary sections and their relevant arguments were inspected; no claim
that the threshold comparison is new follows. A complete strongest-current
priority audit and the collective $G_3$ bound remain unfinished.
The published [Maclagan--Silversmith planar spine theorem, Theorem 2.6](https://wrap.warwick.ac.uk/id/eprint/181583/7/eScholarship%20UC%20item%2074v009bx.pdf)
realizes certain edges over infinite fields. Definition 3.2 describes
the established significant-arrow framework. Those statements and the
Theorem 2.6 proof were checked; they are local edge and planar results,
and supply no proved spatial bound on $V$ in this argument.
