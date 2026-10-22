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
