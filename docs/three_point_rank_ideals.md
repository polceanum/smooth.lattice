# Binomial prefix ranks and three-point product ideals

The [higher-seed core identity](higher_seed_binomial_reduction.md), (B8),
remains unproved. This note identifies the precise two-variable ideals
behind its block ranks and compares them with existing syzygy-gap and
$p$-fractal research. It also excludes a support-matching shortcut, even
in characteristic zero. There is no new full-family or unrestricted
fan-count bound here.

## An integral homogeneous-component identity

Write $q=1+t$. Let $L\ge1$ and let a finite set of levels have
nonnegative indices $j$ and intervals $0\le\ell_j\le m_j$.
Consider the integer columns
$$
\{t^j q^i\pmod {t^L}:\ell_j\le i\le m_j\}.
$$
Choose $D\ge L-1$ and $D\ge j+m_j$ for every level. In
$S=\mathbb Z[x,y]$ put
$$
J=\left(x^L,y^{D+1},
 x^j(x+y)^{\ell_j}y^{D-j-m_j}\ :\ j\right). \tag{P1}
$$
Use empty levels by omitting them. An empty collection has $D=L-1$.

**Proposition 1.** The cokernel of the displayed column matrix is
integrally isomorphic to $(S/J)_D$. Thus over every field $K$ its rank is
$$
L-\dim_K(K[x,y]/J_K)_D. \tag{P2}
$$
The isomorphism preserves the free rank and every nonunit Smith
invariant factor. The two presentation matrices may have different
numbers of unit pivots.

**Proof.** Homogenization identifies $\mathbb Z[t]/t^L$ with
$S_D/(x^L)_D$, using the reference $t^h\mapsto x^h y^{D-h}$,
$0\le h<L$. At a fixed $j$, consecutive powers of $1+t$ span exactly
$\mathbb Z[t]_{\le m_j-\ell_j}$ after the factor $t^j q^{\ell_j}$.
The change from $1,q,\ldots,q^{m_j-\ell_j}$ to ordinary powers of $t$
is integral triangular with diagonal one. Hence the homogenized span
is exactly
$$
x^j(x+y)^{\ell_j}y^{D-j-m_j}S_{m_j-\ell_j}.
$$
Summing gives $J_D$ modulo $(x^L)_D$. The generator $y^{D+1}$ has
no degree-$D$ contribution; it makes the binary ideal Artinian without
changing this component. This proves the integral quotient identity
and (P2). $\square$

For the negative-power columns $t^j q^{-i}$, choose
$M=\max_j m_j$ and multiply by the integral unit $q^M$ modulo $t^L$.
Use the positive intervals
$$
\ell'_j=M-m_j,\qquad m'_j=M-\ell_j. \tag{P3}
$$
Proposition 1 applies unchanged. This is multiplication in the finite
row model, not localization at a nilpotent coordinate of the original
algebra.

In particular (P1)--(P3) apply to every regular and canonical block
in (B4)--(B6). Their generators use only $x,y,x+y$, the three fixed
points $0,\infty,-1$ on the projective line. This describes products
of linear forms in an Artinian ideal; it does not identify the ideal
with the ideal of a union of fat points.

## The uniform gap is now a Hilbert-function identity

For a fixed higher-seed instance, form $(J_B,D_B,L_B)$ for each regular
block and $(J'_C,D'_C,L'_C)$ for each canonical block. Both sets of
reference blocks partition the same free core module, so
$$
\sum_B L_B=\sum_C L'_C.
$$
Consequently (B8) is equivalent to
$$
\sum_B\dim_K(K[x,y]/J_B)_{D_B}
=\sum_C\dim_K(K[x,y]/J'_C)_{D'_C}. \tag{P4}
$$
This is an equivalence for the actual coherent index bounds, not a
replacement by a weaker question. The equality is still open.

One exact symmetry is available. Every permutation of $x,y,x+y$, up
to signs, is induced by a matrix in $\mathrm{GL}_2(\mathbb Z)$.
For example $x\mapsto y,y\mapsto x$ swaps $x,y$, while
$x\mapsto-x,y\mapsto x+y$ swaps $y,x+y$ up to a sign.
Their homogeneous-degree maps are integral automorphisms.
Thus permuting the three exponents in **every** generator preserves
all component cokernels, including their torsion, in all characteristics.
A matching of the block components under these symmetries, or an exact
sequence relating their unmatched components, would be a valid proof
mechanism. No such uniform matching or sequence is established here.

## A coherent arithmetic block

In the regular core take $s=11,a=2,r=3$, prices $u=1,Z=7$, cutoff
$T=14$, and signed block $k=2$. It has $L=4$ and precisely the
columns, after harmless signs,
$$
1,\quad q,\quad tq,\quad t^2q^2 \pmod {t^4}. \tag{P5}
$$
At $D=4$, Proposition 1 gives
$$
J=(x^4,y^5,y^3,x(x+y)y^2,x^2(x+y)^2)
 =(x^4,y^3,x^2y^2,2x^3y). \tag{P6}
$$
Indeed subtract $xy^3$ from the fourth displayed generator and
subtract $x^4+x^2y^2$ from the fifth. The reverse inclusions follow
from the same integer identities. Therefore
$$
(\mathbb Z[x,y]/J)_4\cong\mathbb Z/2\mathbb Z. \tag{P7}
$$
The block rank is four away from characteristic two and three in
characteristic two. This is an actual coherent seed block, not an
arbitrarily selected minor.

The same component torsion appears in degree two of the classical
three-power ideal $(x^2,y^2,(x+y)^2)$: its remaining coefficient is
$2xy$. Its whole colength is three outside characteristic two and
four in characteristic two. This identifies an elementary shared
arithmetic mechanism, not an isomorphism of the whole graded algebras.
Over characteristic zero, (P6) has four minimal monomial generators.
Thus the original ideal is not literally a three-generator syzygy-gap
cell. A further reduction is needed before applying a theorem for one
triple to every block.

## A coherent counterexample to support matching

In the canonical core take $s=17,a=2,r=3,u=9/17,Z=1,T=17/3$,
free exponent $b=1$ and residue $k=2$. Here $L=5$ and the selected
columns, after unit signs, are
$$
1,\quad q^{-1},\quad tq^{-1},\quad tq^{-2},\quad t^2q^{-4}
\pmod {t^5}. \tag{P8}
$$
The supports admit a matching of all five columns to different rows.
Nevertheless the first three satisfy $1-q^{-1}-tq^{-1}=0$.
The rank is exactly four over every field: choosing columns one,
three, four and five, the minors on rows $(0,1,2,3)$ and $(0,1,2,4)$
are $2$ and $-7$. Their relatively prime values give rank at least
four in every characteristic, while the displayed relation gives
rank at most four.

Hence even a coherent cutoff and characteristic zero do not permit
replacing these ranks by maximum matchings of nonzero entries. The
obstruction does not refute (B8); it refutes this proposed shortcut.

## Exact prior-art boundary

[Teixeira, arXiv:1008.0583v1](https://arxiv.org/pdf/1008.0583),
Proposition 2.2 and proof, gives the colength of a three-form binary
ideal in terms of its degrees and syzygy gap. Its displayed free
resolution also gives the full Hilbert series. Example 3.8 and the
statement of Han's Theorem 3.9 concern three powers of linear forms.
We have not identified all ideals (P1) with those cells, and have not
audited the cited proof of Han's theorem. Formula (P2) alone supplies
neither a uniform syzygy-gap formula nor (P4).

[Hernández--Teixeira, arXiv:1404.5871v3](https://arxiv.org/pdf/1404.5871),
Section 2.2, Definition 2.7, SG1, Lemma 2.8 and Corollary 2.9,
were compared with (P1). The syzygy-gap identities concern a triple of
forms; the threshold application fixes a two-generator ideal. Neither
has been matched to all the multi-generator ideals above. The threshold
theorems and their proofs were not audited or applied.

For current context,
[Caminata--Zerman, arXiv:2510.19355v1](https://arxiv.org/pdf/2510.19355),
Definition 3.1 and Corollary 4.3, concern weak $p$-fractals of fixed
hypersurfaces and Frobenius-power quotients. Our ordinary cutoffs,
graded components and varying ideals have not been converted to
that setting. Their rationality result is not a rank or fan-count
theorem for this project. Only the stated definitions, statements
and the immediate Corollary 4.3 implication were inspected.

Homogenization, consecutive-power basis changes, linear-form
permutations and Hilbert-function calculations are classical tools.
The connection above is a supporting reduction, not by itself a
paper-level novelty claim or a prime-distribution result.

~~~bash
python3 -B scripts/run_three_point_rank_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker compares independent binary homogeneous ideal matrices
against the original binomial columns, including integer Smith
signatures. It verifies the coherent arithmetic block, the
support-matching counterexample and the three-form comparison.
The main extremal objective and the uniform coherent identity remain
open.
