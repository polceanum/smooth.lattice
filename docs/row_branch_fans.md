# Arbitrary coupling of fixed row-polynomial branches has quadratic fan size

For finite quotients with $yz=0$, the whole class with row-polynomial
branch annihilators has a quadratic upper bound, uniformly in the row
count, row lengths, field and compatible coupling. This extends the
fixed-core argument of the [minimum-overlap proof](minimal_branch_overlap.md)
without requiring minimum branch length or small excess. Its mechanism
is the same classical lattice-volume estimate; no new geometric theorem,
current novelty, unrestricted $G_3$ bound or completed goal is asserted.

## Actual branch hypothesis and conclusion

Let $A=K[x,y,z]/I$ have finite length $n>0$, with $yz\in I$.
Use the **actual** cyclic branch modules $yA,zA$, and $B=K[x]$.
Assume their annihilators are

$$
 J_Y=(f_0(x),f_1(x)y,\ldots,f_{s-1}(x)y^{s-1},y^s),
 \quad f_{j+1}\mid f_j,\quad f_s=1,                       \tag{1}
$$

and the analogous ideal $J_Z$ with monic polynomials $g_0,\ldots,g_t=1$.
The nonterminal polynomials have positive degree. Zero branches may
be represented by $J_Y=(1)$ or $J_Z=(1)$ and empty row lists. No
irreducibility, locality, reducedness or genericity condition is imposed.
Monomial staircase ideals are included by taking $f_j=x^{a_j}$
for any decreasing positive row lengths $a_j$; equal lengths are allowed.
The successive factors $f_j/f_{j+1}$ need not divide one another.

Write

$$
 m_Y=\sum_{j<s}\deg f_j,\qquad m_Z=\sum_{j<t}\deg g_j,
 \quad a=\dim B1\le n.
$$

**Theorem.** Uniformly over every field and every compatible coupling,

$$
 g(I)=O\!\left((n+1)\sqrt{(m_Y+1)(m_Z+1)}\right)=O(n^2).
                                                               \tag{2}
$$

The maximum order in this class is sharply quadratic, in characteristic
zero. Its total primitive wall cost and primitive one-norm wall length
are $O(n^3)$, with cubic lower order in the same class. These counts
retain repeated wall pieces. They are not actual edge lengths, whose
integer multiplicities may be larger. The union of monic reduced bases
has cardinality $O(n^{8/3})$ by the usual generator packing bound;
no lower bound for a smallest universal basis is asserted.

The hypothesis concerns exact annihilator ideals. Having these monomial
heads for one order, or just the same Hilbert function, does not suffice.
No theorem for arbitrary planar branch deformations is inferred.

## Fixed quotient by the joint annihilator

Let $F(x)$ be the monic annihilator of $1\in A$; its degree is $a$.
Both $f_0,g_0$ divide $F$. The
[actual envelope construction](univariate_branch_gluing.md) gives

$$
 H=K[x,y,z]/(F,yz,yJ_Y,zJ_Z),\quad
 H=B/(F)\oplus y(B[y]/J_Y)\oplus z(B[z]/J_Z).             \tag{3}
$$

Its map onto $A$ is injective on each summand. The coupling kernel
$T$ is a $B$-submodule of $S=\operatorname{ann}_H(y,z)$.
Conversely any such submodule meeting each summand trivially gives
a quotient with these actual endpoints. Thus this proof permits the
full compatible kernel, including a core component, not just diagonal
identifications or independent socle relations.

The envelope length $N=a+m_Y+m_Z$ is at most $2n$. For clarity, put
$W=yA\cap zA$ and $d=\dim A/(yA+zA)$. Multiplication by $y$ maps
$A/(yA+zA)$ onto $yA/y^2A$. Rank-nullity and $W\subset\ker(y|yA)$
give $\dim W\le d$. Therefore $m_Y+m_Z=n+\dim W-d\le n$,
which proves the length bound. These are the previously established
envelope facts, restated to specify the resource parameter.

Put $\ell=\operatorname{lcm}(f_0,g_0)$, using $f_0=1$ for a zero
branch. In the direct decomposition (3), the joint annihilator is

$$
 S=\ell B/(F)\ \oplus\!
 \bigoplus_{j<s}y^{j+1}f_{j+1}B/(f_j)
 \ \oplus\!
 \bigoplus_{j<t}z^{j+1}g_{j+1}B/(g_j).                  \tag{4}
$$

Indeed $y$ sends the $j$-th $Y$ row to the next by reduction modulo
$f_{j+1}$; its kernel in that row is exactly the displayed multiples
of $f_{j+1}$. The $Z$ argument is identical. On the core, simultaneous
annihilation means divisibility by both $f_0$ and $g_0$. The three
components cannot cancel under multiplication by $y$ or $z$, since
their images occupy different branch rows. This proves (4).

Each polynomial relation in (3) has only proper-divisor tails in $x$.
Ordinary row division gives the same envelope standard set $D_H$ for
every global order: the core $x^i$, $i<a$, and rows
$x^iy^{j+1}$, $i<\deg f_j$, and $x^iz^{j+1}$, $i<\deg g_j$.
Adding (4) replaces these row bounds by $\deg\ell$ on the core,
$\deg f_{j+1}$ on a $Y$ row, and $\deg g_{j+1}$ on a $Z$ row.
The quotient $H/S$ likewise has one fixed standard set $C$.

Since $H\twoheadrightarrow A\twoheadrightarrow H/S$ are actual
algebra quotient maps, initial-ideal inclusion gives

$$
                       C\subseteq D_A\subseteq D_H       \tag{5}
$$

for **every** positive initial standard set of $A$. No span of
compulsory columns is being substituted for an algebra ideal here.
The later contraction by those independent columns is only a vector
matroid contraction.

## The total transverse budget telescopes to branch length

Let $E=D_H\setminus C$. Its $Y$ part is exactly

$$
 E_Y=\{(i,j+1,0):\deg f_{j+1}\le i<\deg f_j,\ 0\le j<s\}.
$$

Summing over **all** labels in this part, rather than one basis, gives

$$
 \sum_{u\in E_Y}u_y
   =\sum_{j<s}(j+1)(\deg f_j-\deg f_{j+1})
   =\sum_{j<s}\deg f_j=m_Y.                              \tag{6}
$$

Likewise $\sum_{E_Z}u_z=m_Z$. All other transverse labels are zero.
This remains true for arbitrary row count and arbitrarily reordered
row drops; a bound on the number of invariant factors is unnecessary.

Contract the independent columns $C$ in the actual residue matroid
on $D_H$, retaining $E$. Every parent's positive minimum basis is
its fixed compulsory set together with a minimum basis of this one
residual matroid. Distinct strictly optimal parent standard sets give
distinct exposed residual profiles. The convex hull $Q$ of **all**
residual basis profiles is an integer polytope satisfying

$$
 Q\subseteq[0,\tbinom N2]\times[0,m_Y]\times[0,m_Z].       \tag{7}
$$

For the first coordinate, every residual basis is a subset of the
downset $D_H$ of size $N$, whose total degree sum is at most
$\binom N2$. Equation (6) supplies the other two coordinates even
for residual bases which are not coherent or downsets.

The lattice polytope $Q+[0,1]^3$ is full dimensional, has at least
as many vertices as $Q$, and has volume at most
$(1+\binom N2)(m_Y+1)(m_Z+1)$. The classical three-dimensional
Andrews bound $f_0(P)=O(\operatorname{vol}(P)^{1/2})$ proves (2),
including lower-dimensional $Q$. See
[Dillon, Theorem 1, first page](https://travisadillon.com/documents/volume_vs_vertex_in_lattice_polytopes.pdf)
for a modern author-hosted statement and proof of the established
volume theorem. This is the same geometric input used in the
minimum-overlap result, not a new lattice-polytope theorem.

The positive state polyhedron has $g(I)$ vertices and $O(g(I))$
edges, by the usual three-dimensional planar graph bound (or a distant
truncation). Each actual undirected wall is dual to one such edge.
The primitive threshold theorem gives
$\|q\|_1\le\nu(q)\le n$ for every wall. Thus (2) implies the
cubic bounds stated after it. The existing $O(n^{2/3})$ minimal-head
bound per basis gives the universal-basis union upper bound.

## Sharpness, prior-art boundary, and the remaining obstruction

The [saved Pascal family](univariate_branch_gluing.md#an-explicit-quadratic-fan-and-its-prime-power-collapse)
has $A=K[t]/t^{2m}$ with
$x=t$, $y=t^m(1+t)^m$, $z=t^m(1+t)^{2m}$. Its two actual planar
branch annihilators are $(x^m,y)$ and $(x^m,z)$, so it satisfies (1).
The existing determinant proof gives $\binom{m+2}{2}$ cones in
characteristic zero and cubic primitive wall cost. This proves sharp
orders for the **whole class** above without making a new construction
claim. It is not an all-characteristic lower claim for that same
coefficient family: the known prime-power collapse remains.

The earlier minimum-overlap proof already telescopes the canonical
row ends to the invariant-factor minimum $\Phi$. Here the exact row
lists may be arbitrary, with $\sum\deg f_j$ much larger than that
minimum. The proof simply uses the actual row lengths in (6) and
the envelope's joint annihilator in place of the minimum quotient.
This extends the closure to all these nonminimal endpoints and all
compatible couplings. Its exact quantitative priority is unassessed;
the envelope, annihilator, matroid, telescoping and lattice-volume
ingredients are established or elementary.

For example, take row drops $1,2,\ldots,r$, all powers of $x$.
The branch length is $\sum_{j=1}^r j^2$, whereas the full-kernel
minimum is $\sum_{j=1}^r j(r+1-j)$. Their difference has linear
order in the branch length as $r\to\infty$. Thus the conclusion is
not restricted to the previous small-excess regime, fixed row count,
fixed coordinate degree, or minimum modules. Arbitrary overlap size
or nonconstant surviving $x$-multiples do not evade this bound while
the actual endpoints retain form (1).

A polynomially supergeneric family with $yz=0$ must therefore leave
this row-polynomial endpoint class in at least one actual branch.
This does not mean that every such departure creates a large fan.
General moving planar branch bases, the compatible truncated minors
for square-zero extensions, and unrestricted $G_3(n)$ remain open.
The stored objective and publication criteria are unchanged.

## A general reduction through an actual quotient

There is a useful form of the argument without the row hypothesis.
For **every** finite $A=K[x,y,z]/I$ with $yz=0$, put

$$
 L=\operatorname{ann}_A(y,z),\qquad A^-=A/L.
$$

For the zero quotient use one empty initial standard set. Then

$$
 g(A)\le C(n+1)\sqrt{(m_Y+1)(m_Z+1)}\,g(A^-).             \tag{8}
$$

More precisely, if $k=\dim L$, the first factor can be replaced by
$C\sqrt{(1+k(n-1))(m_Y+1)(m_Z+1)}$.

In particular **every** such algebra whose actual quotient $A^-$
has a single monomial initial ideal has $g(A)=O(n^2)$, with no
row-polynomial assumption on its branches. This is an actual algebra
quotient, rather than contraction by a possibly nonideal vector span.
Equation (8) is a direct lattice-volume reduction, not a resolved
collective bound when $g(A^-)$ grows.

Here is a proof retaining all labels. Fix an initial standard set $D$
of $A$ and its inherited standard subset $C$ for $A^-$. In the normal
forms supported on $D$, an echelon basis of $L$ has precisely the
pivots $D\setminus C$. If $a$ is one of those leading monomials,
multiplying its representative by $y$ and by $z$ gives elements of
$I$. Multiplicativity of a global order implies

$$
 a+e_y\notin D,\qquad a+e_z\notin D.                    \tag{9}
$$

Because $yz=0$, the downset $D$ is an axis and two planar branches.
Its $y$-positive part is $yD_Y$ for a planar downset $D_Y$;
similarly write $zD_Z$. The same row telescoping as (6), now for
arbitrary $D_Y,D_Z$, gives

$$
 \sum_{a\in D\setminus C}a_y\le |D_Y|\le m_Y,
 \qquad
 \sum_{a\in D\setminus C}a_z\le |D_Z|\le m_Z.             \tag{10}
$$

For instance, if the column of $D_Y$ at $x^i$ has height $h_i$, its
one vertically maximal element contributes $h_i$ after multiplying
by $y$; summing heights counts all of $D_Y$. Equation (9) selects
only a subset of these tops. Selected $y$-positive monomials are
independent vectors in $yA$, proving the last inequality. The same
reasoning applies to $z$.

Now fix one full-dimensional initial cone of $A^-$. The set $C$
is constant there. Subtract its moment from every parent initial
profile whose strict cone intersects this cone. By (10), all such
integer profiles lie in
$[0,k(n-1)]\times[0,m_Y]\times[0,m_Z]$, since there are exactly
$k$ remaining monomials, each with $x$-exponent at most $n-1$.
Each is a vertex
of their convex hull: its original strict minimizing weight still
separates it from this restricted list. Applying Andrews after cube
regularization bounds their number by the first factor in (8).
Every parent open cone meets an open cone of the quotient fan,
since finitely many walls cannot cover an open set. Sum over the
$g(A^-)$ quotient cones to prove (8); multiple counting only
overestimates.

The quotient also has $yA^-\cap zA^-=0$. If $yu-zv\in L$, then
$y^2u=z^2v=0$. Thus $yu$ and $zv$ themselves lie in $L$, proving
the assertion. This is an elementary annihilator argument, not a
classification or a counting theorem for all zero-overlap quotients.

In the row-polynomial envelope above,
$\operatorname{ann}_A(y,z)=S/T$. Indeed $yh,zh\in T$ must both be
zero already in $H$, because they occupy a single branch and the
kernel intersects each branch trivially. Hence $A^-=H/S$, which has
the fixed standard set $C$ computed earlier. Thus (2) is also the
single-initial case of (8). The stronger canonical quotient and the
local transverse budget survive arbitrary branch deformation, but
their quotient fan may move; no small number of global pieces is
implied. The prior same-pair recrossing and arbitrary-chain obstructions
remain intact.

## Exact independent controls

```bash
python3 -B scripts/run_row_branch_checks.py
./build.sh
python3 tests/test_smoke.py
```

The [report](../results/research/row_branch_checks.json) records arbitrary
unsorted drops, off-diagonal $K[x]$-linear shears, nonzero core couplings,
and free monic polynomial base changes, over characteristic zero and
several prime fields. Dense simultaneous-kernel elimination independently
checks (4), not just the displayed generators. The fixed quotient map
is independently checked on every envelope column after passing through
the actual coupled quotient. Complete fans and additional strict weights
check the compulsory sets and profile budgets. An exhaustive small
staircase check includes repeated row widths, and the existing Pascal
formula is compared with entire fans and its characteristic collapses.
Asymmetric branch diagrams with surviving $x$-multiples are also checked.
For the general reduction (8), six nonrow forced-socle quotients compute
the actual joint annihilator by dense simultaneous-kernel elimination,
construct the actual quotient matrices, and enumerate both entire fans.
Strict witnesses from both fans check (9)--(10). Their quotient fans
have multiple initial ideals, guarding against extending the
single-initial conclusion without its hypothesis.
The uniform result is the proof above; these finite computations do
not establish its asymptotic order, novelty or the active goal.
