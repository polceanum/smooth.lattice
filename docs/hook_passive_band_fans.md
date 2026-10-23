# Passive transverse layers give a full hook-fan bound

For constant-axis central hooks, a narrow band of off-axis coefficients
near the transverse cutoff gives a small exact encoding of the **entire**
positive fan. This includes high prices and arbitrary long-coordinate
polynomials, without Toeplitz-prefix regularity or a monomial factorization.

Use
\[
 I=(x^2,y^a,z^a,y^tz^t-xf),\quad a\ge2t,\quad t\ge2,\quad
 n=4at-2t^2,\qquad A=a-t,
\]
with zero outer images. In the two rectangular arm modules write
\[
 W_Y=\alpha+N_Y,\quad W_Z=\beta+N_Z,\qquad
 \alpha\beta\ne0,\qquad N_i\in(v)K[q,v]/(q^A,v^t).
\]
Let \(d\) be the smaller transverse order of \(N_Y,N_Z\), using
order \(t\) for zero, and put
\[
                         s=\max(1,t-d).
\]
More generally any \(1\le s\le t\) with \(N_i\in(v^{t-s})\)
can be used. The first definition ensures \(s\le t-1\).
Uniformly over all fields and coefficients,
\[
 \boxed{\quad g(I)=O(a^2s^8)
                  =O\!\left(n^2s^8/t^2\right).\quad}    \tag{PB1}
\]
This is a full-fan parameter bound. For fixed \(s\), it is quadratic
in \(a\), independently of the thickness \(t\). In particular
\(s=O(t^{1/4}\log^c(n+1))\) excludes polynomial excess above
\(n^2\), with the corresponding logarithmic allowance.

For general \(s\) comparable to \(t\), (PB1) recovers the previous
\(O(n^2t^6)\) order. It gives no new worst-case exponent for arbitrary
mixed hooks or unrestricted \(G_3\). Retain all smaller saved bounds,
including the quadratic-coordinate and box ceilings. Exact priority and
sharpness of the new parameter dependence remain unestablished.

## Active and passive short coordinates

The [original bridge model](hook_bridge_reduction.md) contracts the
compulsory plain diagram and the upper corner box. Its residual rows
are
\[
 V_Y=\langle Y_{q,v}:0\le q<A,0\le v<t\rangle,
 \qquad V_Z=\langle Z_{q,v}:0\le q<A,0\le v<t\rangle.
\]
Upper columns are the corresponding identities. Pure lower Y columns
are \(W_Yq^iv^j\), with \(i\ge t\), and analogously for Z.
The central bridges are
\[
 H_{ij}=W_Y Y_{i,j}+W_Z Z_{j,i},\qquad0\le i,j<t,       \tag{PB2}
\]
with original label \((0,t+i,t+j)\). This notation means multiplication
on the two monomial rectangles, not multiplication of row symbols.

Set
\[
 B=[0,s)\cup[t-s,t),\quad m=|B|=\min(t,2s),
 \qquad P=[0,t)\setminus B,
\]
using integer indices. If \(2s<t\), then \(P=[s,t-s)\).
Every \(N_i\) kills an input whose short coordinate lies in \(P\),
and every nonzero output of \(N_i\) has short coordinate in \(B\).
Thus no off-axis term reaches a passive row. For \(2s\ge t\),
\(P\) is empty and the statements below simply retain the old matrix.

## Three exact eliminations on six scalar price sectors

Partition positive weights by the ordering of
\[
                            w_x/t,\ w_y,\ w_z.           \tag{PB3}
\]
There are six open convex sectors. Normalize \(w=(\lambda,u,1)\).
The only separating planes are \(\lambda=t\), \(\lambda=tu\)
and \(u=1\). All the following choices are fixed on each sector.

First, for a pure-arm passive row \(Y_{q,v}\), \(q\ge t,v\in P\),
only its upper identity and its scalar lower source occupy that row.
They form an isolated parallel pair. Keep the cheaper column and
remove the component. Their cost difference is \(\lambda-t\).
On Z it is \(\lambda-tu\).

Second, for \(i,j\in P\), the two passive rows \(Y_{i,j},Z_{j,i}\)
and their three columns form an isolated rank-two component:
\(H_{ij}=\alpha Y_{i,j}+\beta Z_{j,i}\) and the two identities.
Keep its two cheapest columns. The three pairwise cost differences
are \(\lambda-t,\lambda-tu,t(u-1)\), so this choice is fixed by
(PB3). No other column uses these rows.

Third, suppose \(i\in B,j\in P\). The passive row \(Y_{i,j}\)
occurs in exactly two columns: \(H_{ij}\) and its upper identity.
The other part of \(H_{ij}\) lies entirely in active Z rows.
The cheaper of these two columns is always accepted: no earlier
column has a nonzero entry in its passive row. Contract it, eliminating
that row. The other column becomes a nonzero pure Z column, with its
**original label**. For \(i\in P,j\in B\), do the transposed operation.
These leaf choices again use only \(\lambda-t\) or \(\lambda-tu\).

Each contraction is made by eliminating its passive pivot row.
The only other column using that row is its mate, so no other active
column changes. These are cost-compatible vector-matroid operations;
the contracted span is not asserted to be an algebra ideal.

The total removed rank is
\[
 2(A-t)|P|+2|P|^2+2m|P|=2A(t-m).                        \tag{PB4}
\]
The remaining matroid has exactly
\[
 \begin{split}
 r&=2Am\le4As,\\
 |E|&=4Am-m^2=O(as),\\
 |H_{\rm active}|&=m^2\le4s^2.                         \tag{PB5}
 \end{split}
\]
Its pure arm ranks are \(Am\) each, because all active upper identities
remain. The only mixed columns are \(H_{ij}\) for \(i,j\in B\).
Every original strict minimum is the union of the fixed chosen columns
and a minimum of this remaining matroid. Therefore the full cone patches
inside the sector agree, up to a fixed profile translation. This exact
rank and label compression is a further consequence beyond the count.

The original bridge rank can still be \(t^2\): constant axis units
make its Y projections independent. The smaller decision space comes
from the passive components and cost-compatible leaf elimination,
not from replacing that original rank by \(m^2\).

## The retained labels have only a few short levels

All pure Y labels satisfy
\[
                   x\in\{0,1\},\quad0\le y<a,\quad z\in S,
 \qquad S=B\cup(t+B),                                  \tag{PB6}
\]
and all pure Z labels satisfy the transposed condition. For example,
eliminating a passive Z leaf retains either a bridge label
\((0,t+i,t+j)\) or an upper Z label \((1,i,t+j)\), with
\(i\in P,j\in B\); either has its short Z coordinate in \(t+B\).
No exponent labels are changed by row elimination.

The set \(S\) is a union of at most four integer intervals of length
at most \(s\). Hence \(|S-S|=O(s)\). Differences within either pure
arm give only \(O(as)\) comparison planes. Active bridge labels lie
in four translated \(s\)-by-\(s\) corner boxes, so their internal
differences give \(O(s^2)\) planes, absorbed by \(O(as)\).
This also covers the empty-passive case \(m=t\le2s\).

Include the three scalar price planes. There are therefore
\(O(a^2s^2)\) convex chambers fixing all three retained internal
orders and the scalar sector. This inventory is uniform in all
long-coordinate coefficients.

## Count each remaining decision fan

On one such chamber, compress the two pure lists and the active bridges
in increasing cost. They have ranks \(Am,Am,b\), with \(b\le m^2\).
The restricted ground has corank \(b\). Duality and descending
compression leave at most \(3b\) original labels in a fixed rank-\(b\)
matroid, exactly as in the
[actual bridge-rank reduction](hook_prefix_bridge_rank.md).

Its first labels are still zero or one. The saved
[two-colour theorem](two_colour_threshold_braids.md) gives
\(O((1+b)^3)=O(s^6)\) strict regions, also for dual maxima by
negating the costs and first-coordinate parameter. Convex restriction
cannot split a strict region. Summing over the internal-order chambers
proves (PB1). The compressed matrix still depends on the coefficients;
no coefficient-independent sampler of this small cardinality is inferred.

## A necessary transverse width for polynomial excess

Since \(2at\le n\le4at\), an actual family in this constant-axis
central class with \(g(I_n)\ge c n^{2+\epsilon}\) must satisfy
\[
                   s=\Omega\!\left(t^{1/4}n^{\epsilon/8}\right).
                                                               \tag{PB7}
\]
Thus a fixed number of active transverse layers cannot support such an
excess, even with arbitrarily bad long-coordinate prefix gaps and many
incomparable off-axis monomials. Large \(s\) is necessary, not sufficient.
Together with the [clipped-price localization](hook_clipped_threshold_budget.md),
a thin polynomially separating candidate needs both this growing width
and concentration at high prices in the stated thickness regime.

## Verification and priority limits

~~~bash
python3 -B scripts/run_hook_passive_band_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker reconstructs the original quotient and commuting matrices,
checks passive-row incidences and every elimination, and compares whole
original, active and small-dual cone patches in all six price sectors.
It includes nonmonomial, nonunit long coefficients, incomparable short
layers, zero perturbations, small characteristics and rational examples.
Complete original fans provide independent coverage checks.

The scalar-component, leaf and greedy operations are elementary classical
linear and matroid operations. The quantitative count uses the saved
two-colour theorem; its comparison with
[Hausbrandt--Ruzika, Corollary 3.4](https://arxiv.org/html/2503.08178v2#S3)
remains as stated there. The generic greedy setting is also reviewed in
[Eppstein, Section 1.2](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf).
These sources are not asserted to contain this transverse-band theorem.
Exact current novelty, sharpness of (PB1), arbitrary outer images,
the unrestricted extremal endpoint and the final-paper criteria remain open.
