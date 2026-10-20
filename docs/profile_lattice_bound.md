# Profile lattices, free gradings and finite degree spacings

This continuation gives a uniform bound for all finite quotients with a
nontrivial integer grading in three variables, and a separate restriction
on nonlinear rectangle embeddings with repeated degree spacings. It does
not settle the unrestricted three-variable target. The lattice dictionary
and final volume estimate use established state geometry; exact priority
of the exchange-lattice refinement remains unassessed.

## The grading lattice is the exact profile-difference lattice

Let $I\subset K[x_1,\ldots,x_d]$ have colength $n$. From one reduced
Gröbner basis, form the lattice $L(I)\subseteq\mathbb Z^d$ generated
by all differences of exponents of nonzero terms in each polynomial.
Monomial polynomials contribute no differences. The familiar intrinsic
grading argument in the [planar note](extremal_initial_ideals.md#the-intrinsic-lattice-and-its-relative-area)
works in every dimension: $I$ is homogeneous for $\mathbb Z^d/L(I)$,
and every variable-induced group grading making $I$ homogeneous kills
these differences. In a reduced basis, a homogeneous component without
the head would otherwise be a nonzero ideal element supported on standard
monomials. Consequently $L(I)$ is independent of the chosen order.

Use the finite monomial ground

$$
E=V_n^d\cup\bigcup_i(e_i+V_n^d)\subseteq V_{2n}^d,
\qquad V_m^d=\{a\ge0:\prod_i(a_i+1)\le m\}.
$$

Let $P$ be the convex hull of all rank-$n$ monomial basis profiles on
this ground. Then

$$
\langle\mu(B)-\mu(C):B,C\text{ bases on }E\rangle_{\mathbb Z}
=L(I).                                                       \tag{1}
$$

For one inclusion, every monomial basis has exactly the dimension of
each graded component in its corresponding group degree. Its sum of
degrees is fixed. Thus all profile differences belong to $L(I)$, even
when the grading group has torsion divisible by the characteristic.

For the other inclusion, fix an initial standard basis $B$. Every
minimal head $u$ of its reduced basis lies in $E$, since $u-e_i\in B$
for any coordinate with $u_i>0$. If the coefficient of a tail monomial
$v\in B$ is nonzero, replacing $v$ by $u$ gives a vector basis. Its
profile difference is $u-v$. These differences generate every support
difference of the reduced basis, proving (1). In particular
$\dim P=\operatorname{rank}L(I)$.

Ranks here are integer-lattice ranks, not matrix ranks after reducing
exponent differences modulo the characteristic.

The enlargement of the ground matters: some minimal heads lie outside
$V_n^d$. The proof uses ordinary vector bases as well as standard bases.
It does not assert that every one-column replacement is a coherent
standard set. Initial standard sets still inject into the exposed
vertices by the established greedy/state dictionary.

## A rational shadow can use the exchange lattice

Here is the additional form of the [rational-shadow theorem](rational_shadow_bound.md).
Suppose a rank-$r$ matroid has weights $a(e)\in\mathbb N_0^d$ with
$\prod_i(a_i(e)+1)\le H$. Suppose every feasible basis exchange has
weight difference in a lattice $L$. Let $A$ be an integer rank-two
projection with $AL$ of rank two, and put

$$
c=\operatorname{covol}(AL),\qquad
\eta=\frac{\sum_{i<j}|\det(A_i,A_j)|}{c},\qquad
\ell_H=1+\lfloor\log_2 H\rfloor.
$$

The projected basis-profile polygon has

$$
O_d\left(r\ell_H+\ell_H^{d-1}
        +(r^2H\eta)^{1/3}\ell_H+1\right)                    \tag{2}
$$

vertices. Element weights need not lie in one coset of $L$. Only
exchange displacements must belong to it.

For completeness, the earlier proof changes in exactly one place.
Keep its anchored dyadic boxes, determinant-coordinate parallelograms,
unique box labels and four quadrant sweeps. Convex token paths cross
only $O_d(\ell_H)$ shared boundary lines at most twice per token, so
label changes cost $O_d(r\ell_H)$. In a box whose maximizing minor
has size $m$ and whose scaled determinant is $D$, summed internal
exchange directions now belong to the lattice $TAL$, of covolume
$mc$. Their coordinate budgets are at most $2rW_s,2rW_t$.
The lattice-triangle direction estimate gives

$$
6+C_d(r^2D/c)^{1/3}
$$

directions in that label. Summing the nonconstant terms uses the same
weighted dyadic pair sum as before, with $c$ in the denominator, and
gives the third term in (2). The constant six costs the number of boxes,
$O_d(\ell_H^{d-1})$.

That additive term must be retained: $D/c$ need not be an integer or
at least one. Absorbing six into a fractional determinant budget would
invalidate this extension. The ground-avoiding boundary lines still
work, because $A$ is integer; their construction does not require the
ground points to belong to $AL$.

## Every two-dimensional profile space has a height-free bound

Suppose all profile differences span a rank-two lattice
$L=\mathbb Zu+\mathbb Zv$ in $\mathbb Z^d$. Choose $A$ with rows
$u,v$, and write $p_{ij}=u_iv_j-u_jv_i$. Cauchy--Binet gives

$$
\operatorname{covol}(AL)=\det(AA^{\mathsf T})
 =\sum_{i<j}p_{ij}^2.
$$

Thus

$$
\eta=\frac{\sum_{i<j}|p_{ij}|}{\sum_{i<j}p_{ij}^2}\le1.    \tag{3}
$$

The final inequality uses integrality. The map $A$ is invertible on
the real span of $L$, so it preserves all vertices of the profile
polygon. Equation (2) consequently bounds the whole polygon without
any dependence on the height or orientation of its affine plane.
One may use any two independent rows spanning a sublattice of $L$;
the ratio is unchanged, because both numerator and denominator acquire
the same sublattice factor.

Apply this to (1), with $r=n$ and $H=2n$. For fixed $d$, the additive
$\ell_H^{d-1}$ is $O_d(n\log(n+1))$. Rank zero gives one profile and
rank one gives at most two. We have proved the uniform statement

$$
\boxed{\operatorname{rank}L(I)\le2
       \quad\Longrightarrow\quad g(I)=O_d(n\log(n+1)).}      \tag{4}
$$

In three variables this covers every ideal homogeneous for any
nonzero integer grading, with arbitrary signs and zero variable weights
allowed. It includes all positive weighted gradings, uniformly in their
possibly large numerical weights. The maximum over the broader class
in (4) has order $n\log(n+1)$: embed the known generic planar ideals
by killing the other variables. No matching lower bound restricted to
positive three-variable gradings is asserted.

This explains a structural boundary: a free grading forces the whole
fan to the planar growth scale. A power-separated spatial construction
must have full-rank intrinsic grading lattice. The inference uses
exchange directions, not an assertion that all element weights lie in
the profile plane.

There is also a uniform oracle consequence. Given one reduced basis,
the lattice and monomial columns can be computed by established methods.
For rank-two $L(I)$, project with the two rows above. Classical planar
support-oracle enumeration finds the whole profile polygon in $O(g_P)$
lexicographic greedy calls, where $g_P=O_d(n\log(n+1))$ counts all its
vertices. Filter the resulting normal cones for a strict positive
objective to recover exactly the global initials. In fixed dimension
this is an exact linear feasibility question. It does not require a
weight list that works for every ideal simultaneously. The support-oracle
algorithm is classical; the bound on the number of calls is the changed
quantitative input. Field-operation and coefficient-bit costs are not
being identified with that oracle count.

The established minimal-generator packing and union argument also gives
$|\mathcal U(I)|=O_d(n^{2-1/d}\log(n+1))$ for this class, in particular
$O(n^{5/3}\log(n+1))$ in three variables. This bounds the number of
polynomials in a universal reduced basis; it supplies no smooth-number
unranking speed claim.

The same argument has a local consequence for unrestricted ideals.
Every face of a matroid basis polytope is itself a matroid basis polytope
(a direct sum of minors). For a minimizing objective, order its cost
levels and take the restriction to each cumulative prefix followed by
contraction of the preceding prefix. Minimum bases take the prescribed
rank increment in every level. Hence every two-dimensional face of the
profile polytope on $E$ has $O_d(n\log(n+1))$ vertices: its face matroid
has rank $n$, retains the same hyperbolic ground, and has a rank-two
profile lattice. In three variables, every facet therefore has this
bound, including facets of ideals with full-rank $L(I)$. An interior
ray of the positive initial fan is incident to at most that many maximal
cones. The matroid-face fact and the local normal-fan dictionary are
classical; the quantitative input is (2)--(3).

The facet bound has sharp order over the full class. Take a generic
planar length-$n$ ideal $J$ and $I=(J,z-1)$. The face $\mu_z=0$ is
the planar profile polygon and has at least the known $\Theta(n\log n)$
initial vertices. Replacing $1$ by $z$ adds an independent vertical
profile direction, making this face a facet of a three-dimensional
polytope. This sharpness example concerns a boundary normal; it does
not prove a matching lower bound for wholly interior rays. A local
degree bound does not count the number of rays and does not improve
the unrestricted $G_3$ exponent by itself.

## Full rank: volume must be measured in the grading lattice

If $d=3$ and $L(I)$ has rank three and index $h$, (1) makes $P$ a
full-dimensional lattice polytope in a coset of $L(I)$. All coordinates
of $P$ lie in $[0,n(2n-1)]$. Changing to a lattice basis and applying
the classical three-dimensional lattice-vertex bound gives

$$
g(I)=O(n^3/\sqrt h).                                       \tag{5}
$$

This is a direct normalized-volume corollary, rather than a new
lattice theorem or an improvement when $h=1$. If the rank is smaller,
one must use the relative lattice and (4); dividing a planar volume
by a full ambient index is not valid.

For rectangles with $r,s\ge2$, the filtration bound is sharper. Let

$$
I_F=(x^r,y^s,z-F),\quad n=rs,\quad
F=\alpha x+\beta y+H,\quad \alpha\beta\ne0,\quad H\ne0,
$$

with $F$ reduced modulo $(x^r,y^s)$ and $H\in(x,y)^2$. Put

$$
h=\gcd\{a+b-1:(a,b)\in\operatorname{supp}H\}.
$$

The displayed generators are already reduced for lex order eliminating
$z$. Their support differences give exactly

$$
L(I_F)=\{v\in\mathbb Z^3:v_1+v_2+v_3\equiv0\pmod h\}.       \tag{6}
$$

Indeed $e_z-e_x,e_z-e_y$ generate the primitive sum-zero plane;
the additional sums are $1-a-b$, whose gcd is $h$. This calculation
works in every characteristic. The reference rectangle basis and its
one-column replacements by $z$ witness the exact index and full rank.

The $x$-adic, $y$-adic and ordinary-degree filtrations bound every
monomial vector-basis profile, respectively, by

$$
\mu_x\le n(r-1)/2,\quad \mu_y\le n(s-1)/2,\quad
\mu_x+\mu_y+\mu_z\le n(r+s-2)/2.
$$

The normalized-volume estimate therefore yields

$$
g(I_F)=O\left(n^{3/2}
 \sqrt{\frac{(r-1)(s-1)(r+s-2)}h}\right).                  \tag{7}
$$

For a square rectangle, $n=q^2$, this is

$$
g(I_F)=O(n^{9/4}/\sqrt h).                                 \tag{8}
$$

In particular, if $H$ is any nonzero homogeneous polynomial of degree
$D$ and $D-1=\Omega(q)$, the entire family has $g(I_F)=O(n^2)$.
The coefficients and the field are unrestricted. More generally this
holds whenever the gcd in (6) is $\Omega(q)$. Thus a single nonlinear
homogeneous correction of order comparable with the side length cannot
give a power separation, even when the previous terminal-band rank
compression leaves almost the entire quotient coupled.

Combining (6) with the [high-order direct sum](rectangle_rank_structure.md)
gives, for first nonlinear degree $q\le D\le2q-2$,

$$
g(I_F)=O\left(\frac{q^{3/2}(2q-D)^3}{\sqrt h}\right).        \tag{9}
$$

To justify the division in this component, use its pure-$x,y$ reference
basis. It contains $x,y$ and every term of $H$. Replacements by $z$
again generate (6), so the component profile polytope is genuinely
three-dimensional with that exact lattice. Its rank is $(2q-D)^2$;
the unchanged middle has at most five cones. This avoids an unjustified
ambient-index argument on a lower-dimensional summand.

A proposed square family with $g(I_F)\ge c n^{2+\epsilon}$ must in
particular have $h=O(n^{1/2-2\epsilon})$. Generic mixed low-degree
terms can have $h=1$; they remain open. This is a restriction on possible
constructions, not a proof of the unrestricted upper endpoint.

## Primary-source comparison

[Babson--Onn--Thomas, Proposition 2.2, Theorem 2.4, Corollary 2.5](https://arxiv.org/pdf/math/0207135)
supplies the universal ground, greedy/state dictionary and classical
volume method. The intrinsic grading argument was already used in the
repository's planar comparison. The new proof step to assess is (2)
with the exchange lattice, including its additive box term, and its
uniform consequence (4). No new Smith-normal-form, Cauchy--Binet or
lattice-vertex theorem is claimed.

The available [Katoh--Tokuyama author preprint](https://ics.uci.edu/~eppstein/pubs/notmine/Katoh-Tokuyama-klevel.ps.gz),
pages 1--2 and 9--11, was inspected as rendered pages. Its Theorem 3.1
and Corollary 3.1 concern levels of concave surfaces; Section 4 explicitly
leaves general two-parameter matroids outside the proof. Its numerically
smaller exponent cannot simply replace the present bound. The
[journal article](https://doi.org/10.1007/s00454-001-0086-z) uses different
numbering; its full text was not retrieved, so version identity is not
assumed. This resolves the applicability of the inspected preprint,
not the strongest-current-bound audit.

[Hausbrandt--Ruzika, Theorem 3.2 and Corollary 3.4](https://arxiv.org/html/2503.08178v1)
were rechecked, including the arrangement-count proof: the unrestricted
$O(m^{2p})$ comparison counts cells for $m$ elements and $p$ parameters.
It does not include a product budget or the exchange-lattice factor.
[Maruri-Aguilar, Theorem 5.1 and Section 6](https://rendiconti.dmi.units.it/volumi/37/05.pdf)
constructs the planar Hilbert-zonotope weight list of order $n^2$.
That universal list is distinct from the number of initials of one ideal.
These are scoped comparisons, not a novelty certificate for (2)--(4).

[Przybylski--Klamroth--Lacour, Section 2.1 and Algorithm 1](https://arxiv.org/pdf/1911.08937)
describes the established dichotomic scheme and its support-facet version.
Those parts were inspected. The planar enumeration above uses the same
edge-normal queries and lexicographic endpoint choice; it is not a new
convex-hull algorithm. The checker recovers each complete polygon with
exactly twice as many support queries as vertices.

## Verification

```bash
python3 -B scripts/run_profile_lattice_checks.py
./build.sh
python3 tests/test_smoke.py
```

The controls compare full basis-profile hulls with four complete projected
sweeps, replay each exchange in the finer lattice, and check the normalized
internal direction budgets. Offset element planes and coarse exchange
lattices are included. Rectangle controls verify the exact index, all
fundamental exchanges and the high-order summand, including characteristics
dividing the index. A separate planar support-oracle enumeration is checked
against every full hull. Gram determinants, Plücker coordinates and gcds
of image minors independently verify normalization under nonprimitive row
choices in dimensions two through five. Computation supports these uniform proofs; it does
not establish their priority or complete the main goal.
