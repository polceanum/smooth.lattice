# A fixed Fitting core and a branch-excess fan bound

Let \(A=K[x,y,z]/I\) be a length-\(n\) quotient over any field, with
\(yz\in I\). Suppose \(W=yA\cap zA\ne0\), and retain the descending
monic invariant factors and notation of the
[minimum-overlap theorem](minimal_branch_overlap.md):
\[
 W=\bigoplus_{i=1}^r K[x]/(d_i),\quad d_{i+1}\mid d_i,\quad
 h=\sum_i\deg d_i,\quad \Phi=\sum_i i\deg d_i,\quad
 g_j=\prod_{i>j}d_i,\quad D=\deg d_1.
\]
Put
\[
 e_Y=\dim yA-\Phi,\quad e_Z=\dim zA-\Phi,\quad
 s=n-2\Phi=e_Y+e_Z+\dim A/(yA+zA)-h.                     \tag{F1}
\]
All terms on the right are nonnegative. Write \(e=e_Y+e_Z\le s\)
for the total planar branch excess.

**Theorem.** The actual ideal has a canonical containment
\[
 I+W\ \subseteq\ J^*:=(g_0,yz,g_1y,g_2y^2,\ldots,y^r,
                              g_1z,g_2z^2,\ldots,z^r).   \tag{F2}
\]
Here \(I+W\) means the inverse image of the ideal \(W\subset A\).
The quotient by \(J^*\) has length \(2\Phi-h\) and one initial
ideal for every global order. Its standard set \(C\), of that size,
is contained in **every** standard set of \(A\). The kernel of
\(A/W\to K[x,y,z]/J^*\) has exactly dimension \(s\).

There is a monic annihilator of \(x\) in \(A\) of degree at most
\(h+D+s\). With absolute constants, the complete fan satisfies
\[
 \boxed{g(I)=O\bigl((h+s)(\Phi+re+e^2)\bigr)
                    =O(n^2+ne^2).}                       \tag{F3}
\]
In particular the entire class \(e=O(\sqrt n)\) has \(g(I)=O(n^2)\),
even when the shared axis has arbitrary extra length. This includes
\(s=O(\sqrt n)\), since \(e\le s\).
Its quadratic order is sharp in characteristic zero. Every reduced
Gröbner basis has at most \(2r+e_Y+e_Z+4\) polynomials.

This extends the minimum-length result to arbitrary nonminimal branches;
it does not settle all \(yz=0\) quotients or unrestricted \(G_3\).
The Fitting containment below is a classical Hilbert--Burch/Smith-form
corollary. The lattice vertex estimate is also established. Exact priority
and significance of the assembled quantitative theorem remain unassessed.
The full-kernel extension below removes the nonzero-overlap assumption
and uses a smaller branch excess. A separate generic reduced-point
calculation treats an important zero-overlap control class.

## The classical coefficient containment

For a finite planar quotient \(L=K[x,y]/J\), let
\(U=\ker(y:L\to L)\). Its Fitting ideals over \(B=K[x]\) are
\(\operatorname{Fitt}_j^B(U)=(q_j)\), with the convention \(q_j=1\)
beyond its number of nonunit invariant factors. Then
\[
                 J\subseteq\sum_{j\ge0}(q_j)y^j.         \tag{F4}
\]
The right side is an ideal because \(q_{j+1}\mid q_j\).

Here is the determinantal proof, to make the prior-art reduction explicit.
A finite-colength planar ideal has a Hilbert--Burch resolution
\[
 0\longrightarrow R^t\xrightarrow{M}R^{t+1}
       \xrightarrow{\text{signed minors}}R\longrightarrow L
       \longrightarrow0,\qquad R=B[y].
\]
The polynomial, possibly nonminimal form is sufficient. This is the
established planar Hilbert--Burch framework, for example
[Conca--Valla, Section 3, Theorem 3.3 and its proof](https://arxiv.org/html/0708.3576v1).
No claim about generic initial ideals or coordinate changes is needed.

Tensor with \(R/(y)=B\). Its first homology is \(U\), by the
two-term resolution for \(R/(y)\); its second homology is zero.
Thus \(M_0=M(x,0)\) has full column rank. Smith operations over
\(B\), lifted unchanged to \(R\), put it in the form
\[
 M_0=\begin{pmatrix}\operatorname{diag}(a_1,\ldots,a_t)\\0\end{pmatrix},
 \qquad a_1\mid\cdots\mid a_t.
\]
The reduced last map has kernel equal to the first \(t\) coordinate
summands, since its nonzero last entry is a nonzero polynomial in
the domain \(B\). Therefore \(U=\bigoplus B/(a_i)\).

Write \(M=M_0+yN\). A coefficient of \(y^j\) in a maximal minor
uses at most \(j\) entries from \(yN\). Every nonzero remaining
product contains at least \(t-j\) distinct diagonal entries of
\(M_0\), hence is divisible by \(a_1\cdots a_{t-j}=q_j\).
Terms using fewer entries are divisible by the same polynomial.
This proves (F4) for all maximal minors and hence for \(J\).
It uses neither a generic perturbation nor an assertion that arbitrary
matrix perturbations preserve colength.

If \(V\subseteq U\) is a \(B\)-submodule, its primary invariant
exponents are componentwise at most those of \(U\). Consequently
\(\operatorname{Fitt}_j(U)\subseteq\operatorname{Fitt}_j(V)\).
One may therefore replace \(U\) by any prescribed \(V\) killed
by \(y\) in (F4). The earlier planar minimum and its rigidity follow
immediately: the resulting canonical row ideal has colength \(\Phi(V)\),
and equality of finite colengths forces equality of the ideals.
This is a stronger classical reduction of that auxiliary statement.

## The simultaneous quotient for both branches

Apply (F4) with \(V=W\) to the two actual branch modules
\(yA\) and \(zA\), with cyclic vectors \(y,z\). Their annihilator
ideals are \((I:y)\cap K[x,y]\) and \((I:z)\cap K[x,z]\).
They are contained in \((g_0,g_1y,\ldots,y^r)\) and its \(z\)
analogue, respectively.

Modulo \(yz\), write any \(f\in I\) as
\(p_0(x)+\sum_{j\ge1}p_j(x)y^j+\sum_{j\ge1}q_j(x)z^j\).
Multiplication by \(y\) shows that its core-plus-\(y\) part belongs
to the first branch annihilator; the other part is killed by \(y\).
Similarly for \(z\). Thus \(g_0\mid p_0\),
\(g_j\mid p_j,q_j\) for \(j<r\), proving \(I\subseteq J^*\).
If an element of \(W\) has representatives \(yF=zG\), their
difference belongs to \(I\). Coefficient separation in (F2) puts
each representative in \(J^*\), so \(I+W\subseteq J^*\).

Monic row division gives the fixed standard set
\[
 C=\{x^i:i<h\}\ \cup\!\bigcup_{0\le j<r-1}
 \{x^iy^{j+1},x^iz^{j+1}:i<\deg g_{j+1}\}.              \tag{F5}
\]
Its size is \(h+2(\Phi-h)=2\Phi-h\). Each row relation has only
proper-divisor tails in \(x\), so the same initial ideal works for
every global order. Initial-ideal inclusion proves that \(C\) is
forced in every parent standard set. This conclusion does not require
the fan of \(A/W\) itself to be constant.

Let \(V=J^*/I\), an actual ideal of \(A\). Then
\(W\subseteq V\), \(\dim V=h+s\), and \(\dim V/W=s\).
Choose a monic annihilator \(p(x)\) of \(V/W\) of degree at most
\(s\), for example its characteristic polynomial; put \(p=1\)
when \(s=0\). Since \(g_0A\subseteq V\),
\[
                  d_1p g_0\in I.                         \tag{F6}
\]
This is the usual annihilator-of-an-extension argument. In particular
every selected pure-axis label has exponent below \(h+D+s\).

## A bounded relative ground, without moving branch matrices

The canonical planar row downset is
\[
 D_W=\{(i,j):0\le j<r,\ 0\le i<\deg g_j\},\quad |D_W|=\Phi.
\]
At any internal positive planar weight, a branch's greedy standard set
\(D_Y\) contains \(D_W\) by (F4) and has \(e_Y\) additional labels.
Any label rejected within that branch is also rejected in the full
merged scan. Thus selected \(Y\)-labels lie in \(yD_Y\); likewise
for \(Z\). This reasoning concerns ordinary greedy bases, not a
generally false matroid assertion for a fixed branch cardinality.

There is a single finite ground for all these choices. If \(a\notin D_W\)
belongs to \(D_Y\), choose a minimal excluded monomial \(u\le a\)
of \(D_W\). The entire translated divisor box \([u,a]\) is outside
\(D_W\), so \(\prod_i(a_i-u_i+1)\le e_Y\). Hence
\[
 D_Y\setminus D_W\subseteq
 \bigcup_{u\in\operatorname{Min}(\mathbb N^2\setminus D_W)}
                   \left(u+E_{e_Y}^{(2)}\right),           \tag{F7}
\]
where \(E_e^{(2)}=\{a:\prod_i(a_i+1)\le e\}\) and \(E_0\)
is empty. There are \(r+1\) excluded generators. The same holds
for \(Z\). Use (F6) for the axis and the two fixed shells (F7)
for the branches, and discard labels outside the ordinary length-\(n\)
hyperbolic ground. Every parent standard set is contained in this
fixed labelled ground and contains \(C\).

Contract \(C\) in its actual residue matroid. The residual rank is
\(h+s\). A ground of size at most
\[
 2h+D+s+(r+1)\bigl(|E_{e_Y}^{(2)}|+|E_{e_Z}^{(2)}|\bigr) \tag{F8}
\]
suffices, before removing overlaps and loops. It is one fixed matroid
for the whole fan. This is a vector contraction; it need not be an
algebra quotient. Formula (F8) is a candidate-column bound, not a
claim about optimal runtime or coefficient bits.

## Moment budget and fan count

Remove the fixed moment of \(C\) from a parent standard set \(S\).
The residual set has size \(h+s\). Every label's \(x\)-exponent is
at most \(h+D+s-1\): (F6) handles the axis, while a branch's own
standard set contains \(D_W\) and has at most \(e_Y\le s\) additions,
so its maximum \(x\)-exponent is at most \(h+e_Y-1\).
Therefore the residual \(x\)-moment is at most
\[
                   X=(h+s)(h+D+s-1).                     \tag{F9}
\]

The fixed boundary \(yD_W\setminus C\) has total \(y\)-moment
\(\sum_{j<r}(j+1)(\deg g_j-\deg g_{j+1})=\Phi\).
For the additions to \(D_W\), fix a vertical column. Its old height
is at most \(r\); if \(k\) new labels occur there, their shifted
\(y\)-exponents sum to at most \(rk+k(k+1)/2\). Summing and using
\(\sum k=e_Y\) gives
\[
 Y=\Phi+r e_Y+e_Y(e_Y+1)/2,\quad
 Z=\Phi+r e_Z+e_Z(e_Z+1)/2                               \tag{F10}
\]
as residual moment bounds in the other two coordinates.

Take the convex hull of the profiles of actual standard sets after
subtracting the moment of \(C\). Each is a distinct exposed vertex:
at an interior generic weight it uniquely minimizes among all monomial
vector bases, hence among this subset. The hull is integral and lies
in \([0,X]\times[0,Y]\times[0,Z]\). Adding the unit cube preserves
at least its number of vertices and makes it full dimensional.
The established Andrews vertex estimate in dimension three gives
\[
                 g(I)=O\bigl(\sqrt{(1+X)(1+Y)(1+Z)}\bigr).
                                                               \tag{F11}
\]
See [Andrews, 1963](https://doi.org/10.1090/S0002-9947-1963-0143105-7)
and [Dillon, Theorem 1](https://travisadillon.com/documents/volume_vs_vertex_in_lattice_polytopes.pdf).
Since \(D\le h\), (F11) implies the first bound in (F3).
Also \(h+s\le n\), \(\Phi\le n/2\), and \(r(r+1)\le n\);
thus \(\Phi+re+e^2=O(n+e^2)\), proving the second.
A basis has at most \(r+e_Y\) nonempty \(Y\)-rows and
\(r+e_Z\) nonempty \(Z\)-rows. Their frontiers, the axis head and
\(yz\), give the stated reduced-basis cardinality bound.

The previous Pascal family has \(s=0\) and quadratic fan growth in
characteristic zero, proving sharpness of the near-minimum class order.
A polynomially supergeneric \(yz=0\) family with
\(g(I)\ge c n^{2+\epsilon}\) must satisfy
\[
                         e=\Omega_{c,\epsilon}
                               (n^{(1+\epsilon)/2}).       \tag{F12}
\]
This is a necessary condition, not a sufficient construction criterion.
The result need not improve existing bounds when \(e\) is comparable
to \(n\). If both branches attain their planar minima, \(e=0\)
and quadratic growth holds for every axis excess \(s\).

## Verification and remaining work

```bash
python3 -B scripts/run_fitting_branch_core_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker verifies finite planar coefficient containments, genuinely
nonminimal branch gluings, coupling through the shared axis, the actual
canonical quotient map, invariant factors recovered from power ranks,
the fixed shell ground, and whole rational fan cones. These are supporting
checks; (F4)--(F11) supply the uniform proof.

The strongest primary comparisons used here are the classical polynomial
Hilbert--Burch representation of
[Conca--Valla](https://arxiv.org/html/0708.3576v1) and the precise
resolution statement in [Oszer, Theorem 2.1](https://arxiv.org/html/2407.07993v4).
Oszer's introduction and Section 6 explicitly warn that unrestricted
matrix perturbations need not preserve finite colength. No such
perturbation assumption is used here. Smith form, Fitting ideals and
the determinantal expansion in (F4) are established algebraic methods;
the planar auxiliary minimum is now classified as a classical consequence,
rather than an independently new theorem.

The priority of the simultaneous canonical quotient and excess-length
fan theorem still requires a theorem-level comparison with the literature.
The unrestricted exponent, an infinite polynomial separating family,
a consequential unrestricted complexity result and the publication goal
remain open. The next counting obligation is the regime of large planar
excess, with the shared \(x\)-action and actual dependence matrices
retained. The original success criteria are unchanged.

## Full branch kernels: a stronger quotient, including zero overlap

The shared module need not be used in both applications of (F4). Put
\[
 K_Y=\ker(y:yA\to yA),\quad K_Z=\ker(z:zA\to zA),\quad
 k_Y=\dim K_Y,\quad k_Z=\dim K_Z,\quad h=\dim W.
\]
Use descending factors \(d^Y_i,d^Z_i\), their weighted lengths
\(\Phi_Y,\Phi_Z\), generator numbers \(r_Y,r_Z\), and row products
\(g^Y_j=\prod_{i>j}d^Y_i\), \(g^Z_j=\prod_{i>j}d^Z_i\).
For an empty kernel set \(r=k=\Phi=0\) and every row product to one.
Define
\[
 \begin{split}
 q_0&=\operatorname{lcm}(g^Y_0,g^Z_0),\quad c=\deg q_0,\\
 e_Y'&=\dim Y-\Phi_Y,\quad e_Z'=\dim Z-\Phi_Z,\quad e'=e_Y'+e_Z',\\
 d&=\dim A/(Y+Z),\quad S=e'+d-c.                         \tag{K1}
 \end{split}
\]
These parameters apply to every positive-length \(yz=0\) quotient,
whether or not \(W\) is zero. The corresponding canonical ideal is
\[
 J^\dagger=(q_0,yz,\ g^Y_jy^j\ (j\ge1),\
                                g^Z_jz^j\ (j\ge1)).       \tag{K2}
\]
Only finitely many generators are necessary. In particular an empty
kernel contributes the generator \(y\) or \(z\), respectively;
it does **not** contribute the constant generator one. If both kernels
are empty, \(q_0=1\) and the whole canonical quotient is zero.

**Full-kernel theorem.** The actual ideal \(M=K_Y+K_Z\) is contained
in \(J^\dagger/I\), and \(I\subseteq J^\dagger\). This canonical
quotient has one global initial ideal and length
\[
 L=c+\Phi_Y-k_Y+\Phi_Z-k_Z.                               \tag{K3}
\]
Its fixed standard set is forced in every parent standard set. Let
\(m=\dim M=k_Y+k_Z-h\). Then
\[
 S\ge0,\quad R:=n-L=m+S,\quad
           \dim(J^\dagger/I)/M=S.                        \tag{K4}
\]
Let \(b=\operatorname{lcm}(d^Y_1,d^Z_1)\), treating a missing factor
as one, and \(D=\deg b\). There is a monic \(x\)-annihilator of
degree at most \(c+D+S\). The full fan satisfies
\[
 \begin{split}
 X&=R\max(0,c+D+S-1),\\
 Y_0&=\Phi_Y+r_Ye_Y'+e_Y'(e_Y'+1)/2,\\
 Z_0&=\Phi_Z+r_Ze_Z'+e_Z'(e_Z'+1)/2,\\
 g(I)&=O\left(\sqrt{(1+X)(1+Y_0)(1+Z_0)}\right)
       =O\bigl(n^2+n(e')^2\bigr).                        \tag{K5}
 \end{split}
\]
Every reduced basis has at most \(r_Y+r_Z+e'+4\) polynomials.
Thus the quadratic class includes \(e'=O(\sqrt n)\), even with zero
overlap and arbitrary axis length. Where the previous shared-overlap
parameters exist, \(e'\le e\), and its old fixed standard set is
contained in the stronger one. Exact priority remains unassessed.

**Proof.** Both kernels are ideals of \(A\): they are \(x\)-stable,
and \(y,z\) kill them. Their intersection is exactly \(W\).
Apply (F4) to the full kernel of each branch. The core-plus-branch
argument used for (F2) gives \(I\subseteq J^\dagger\).
For \(w=yf\in K_Y\), one has \(y^2f=0\), hence the polynomial
\(yf\) belongs to the branch annihilator \((I:y)\cap K[x,y]\).
Its coefficient containment puts this positive-\(y\) polynomial in
\(J^\dagger\). The same holds for \(K_Z\), proving the inclusion
of \(M\). This also handles nonlocal branches and empty kernels.

The characteristic polynomials of \(x\) on \(K_Y\) and on
\(Y/yY\) agree, by the exact kernel-cokernel sequence for \(y:Y\to Y\).
The latter is cyclic over \(B\), so its monic annihilator is exactly
\(g^Y_0\). The axis surjection \(A/(Y+Z)\to Y/yY\), and its
\(Z\)-analogue, show that \(q_0\) divides the axis annihilator.
Consequently \(d\ge c\). Classical planar minimality gives
\(e_Y',e_Z'\ge0\). Monic row division proves (K3), its universal
initial ideal and its fixed standard set. Length bookkeeping gives (K4).

Set \(V=J^\dagger/I\). The polynomial \(q_0\) sends \(A\) into
\(V\); a degree-at-most-\(S\) annihilator \(p\) of \(V/M\) then
sends it into \(M\), and \(bM=0\). Thus \(bpq_0\in I\).
In fact \(D\le m\), since \(b\) is the minimal annihilator of the
\(B\)-module \(M\). Therefore \(c+D+S\le c+R\le n\).

In each branch replace \(D_W\) in (F7) by the downset from its own
kernel. It has \(\Phi_Y\) or \(\Phi_Z\) labels and at most
\(e_Y'\) or \(e_Z'\) additions. If that downset is empty, its
minimal excluded generator is the origin and its shell is simply the
ordinary length-\(e\) planar hyperbolic cross. These two fixed shells
and the axis bound give one fixed residual matroid of rank \(R\)
for the whole positive fan. The residual ground has at most
\[
 k_Y+k_Z+D+S+(r_Y+1)|E_{e_Y'}^{(2)}|
                      +(r_Z+1)|E_{e_Z'}^{(2)}|             \tag{K6}
\]
labels before discarding overlaps and loops.

The residual moments satisfy (K5) by the same column-addition budget
as (F9)--(F10), now separately in the two branches. This remains true
when the canonical quotient is zero and its compulsory set is empty.
Apply the established three-dimensional lattice vertex theorem.
Finally \(R\le n\), \(c+D+S\le n\), \(\Phi_Y,\Phi_Z\le n\),
and \(r_Y,r_Z\le\sqrt{2n}\), so \(Y_0,Z_0=O(n+(e')^2)\).
This proves the second estimate. Counting branch row frontiers gives
the reduced-basis bound. The inequalities comparing old and new cores
follow from \(W\subseteq K_Y,K_Z\) and the PID invariant-factor
embedding criterion. \(\square\)

For generic reduced points off the shared axis, multiplication by \(y\)
is invertible on \(Y\), and similarly on \(Z\). Both kernels vanish;
then \(e'=n\), so (K5) does not close that control family. Its rank
structure admits a different bound, as follows.

## Generic reduced points on the two planes

Let \(m_Y,m_Z\ge1\), \(n=m_Y+m_Z\), and choose sufficiently general
distinct points, with \(m_Y\) on \(z=0,y\ne0\) and \(m_Z\) on
\(y=0,z\ne0\). Over every characteristic their point ideal satisfies
\[
                         g(I)=O(n^2\log^2(n+1)).           \tag{N1}
\]
Sufficiently general means a specified nonempty Zariski-open condition
on this product of planes. A given finite field need not have rational
points in that open set; extension fields are allowed. Nongeneric point
configurations are not included.

On the finite ground
\(E=\{a:\prod_i(a_i+1)\le n,\ a_ya_z=0\}\), divide labels into
the pure axis, the positive-\(y\) plane and the positive-\(z\) plane.
For a generic configuration a set of at most \(n\) labels is independent
exactly when
\[
            |B\cap Y|\le m_Y,\qquad |B\cap Z|\le m_Z.    \tag{N2}
\]
Mixed \(yz\)-labels are loops. This is a truncation of a partition
matroid, not a uniform rank rule on the whole union of planes.

To prove generic existence, every structurally allowed \(n\)-set has
a permissible assignment of its columns to the \(m_Y\) and \(m_Z\)
rows: assign the branch columns to their own rows and fill the remaining
rows with axis columns. Its evaluation determinant is a nonzero
polynomial in the independent point coordinates. Indeed every nonzero
permutation term is a different monomial: at each labelled row, distinct
surviving columns restrict to distinct monomials on that plane. Thus
there is no cancellation, in any characteristic. The finite product of
these determinants is nonzero. Impose the nonzero branch coordinates
and distinctness as further open conditions. Any smaller structurally
allowed set extends to a permitted \(n\)-set by adding unused axis
labels, proving (N2) for all subsets. Support dimensions give necessity.
This is the usual generic-minor and partition-matroid argument.

At a generic positive objective the constant label is always selected.
Therefore both branch capacities cannot simultaneously be full. There
are only three possibilities:

1. If both capacities are slack, the basis is the \(n\) cheapest
   labels of \(E\). Every exchange is feasible, so any other choice
   would fail optimality.
2. If \(Y\)'s capacity is full, its part is the \(m_Y\) cheapest
   positive-\(y\) labels and its complement is the \(m_Z\) cheapest
   labels on the \((x,z)\)-plane, including the axis. Each part can
   be optimized independently without violating the other capacity.
3. The analogous description holds with \(Y,Z\) reversed.

The ground has size
\[
             M=2\sum_{j=1}^n\lfloor n/j\rfloor-n
                    =O(n\log(n+1)).                       \tag{N3}
\]
It lies on the boundary of the positive orthant. The already proved
[boundary-set lemma](generic_multiaffine_bound.md#a-classical-bound-for-points-on-the-orthant-boundary)
bounds the first case by \(2n(M-n)-M+2=O(n^2\log(n+1))\).
Its input is the classical convex-position facet count in
[Leroux--Rademacher, Section 2.5, Proposition 2.13](https://www.math.ucdavis.edu/~lrademac/algebraicksets.pdf),
with perturbation as in Proposition 2.11. This is not a new spatial
\(k\)-set theorem for matroids.

For a full capacity, deleting the fixed \(y\)-tag leaves an ordinary
planar corner cut of size \(m_Y\); its complement is one of size
\(m_Z\). If \(C_2(m)\) denotes that set count, the two full-capacity
cases together give at most \(2C_2(m_Y)C_2(m_Z)\) bases. The
[established planar bound](planar_dyadic_bound.md)
\(C_2(m)=O(m\log(m+1))\) proves (N1). Internal orders
and extra subdivisions are not multiplied in: a pair of selected sets
determines the parent basis. This upper bound is not asserted sharp.

Thus neither small full-kernel excess nor generic reduced points on this
reducible multiaffine carrier can provide a polynomial separation above
the generic \(n^2\log^2 n\) scale. The remaining \(yz=0\) route must
have large full-kernel branch excess and nongeneric dependent ranks.
This narrows a realizable search mechanism while leaving the original
unrestricted objective intact. The assembled quantitative statements
still require an exact priority and significance assessment.

The checker retains the old 78 witnesses and counters separately. Its
full-kernel extension verifies actual quotient maps and shells, including
empty kernels, unequal branches, partial kernel identifications, zero
overlap and nonlocal base changes. The generic-point controls check
every structurally permitted maximal minor on the recorded ground and
compare complete rational cones with a separate capacity-matroid sweep.
For a concrete scope guard, the four distinct rational points
\((1,1,0),(2,2,0),(3,0,1),(4,0,2)\) have two points on each plane,
but \(\{1,x,y^2,xy\}\) has rank three although (N2) would allow rank
four. The relation \(xy-y^2\) vanishes at all four points. Thus
reducedness and a component count alone do not justify (N2).
