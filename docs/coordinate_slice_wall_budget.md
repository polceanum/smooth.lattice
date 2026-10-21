# A sharp collective budget for walls in a coordinate plane

Let \(I\subset K[x,y,z]\) have positive colength \(n\), over any field.
For a standard set \(B\), write
\(\mu(B)=\sum_{\alpha\in B}\alpha\), and define

\[
 k_-=\min_B\mu_x(B),\qquad k_+=\max_B\mu_x(B),\qquad D=k_+-k_-.
\]

The extrema range over the actual positive initial standard sets.
Let \(p_I(q)\) count undirected positive wall pieces of primitive
direction \(q\), retaining repeated parallel pieces. Put

\[
 V_x^0(I)=\sum_{\substack{\text{undirected actual walls }e\\
                              (q_e)_x=0}}\nu(q_e).
\]

The [primitive threshold](primitive_wall_spectrum.md) has
\(\nu(0,a,-b)=a+b\) for coprime \(a,b>0\).

**Theorem.** Uniformly over fields and finite quotients,

\[
 \begin{split}
 p_I(0,a,-b)&\le D+1,\\
 V_x^0(I)&\le(D+1)\bigl(n(n-1)-k_--k_+\bigr).
 \end{split}                                                   \tag{S1}
\]

In particular, if \(I\) contains a monic univariate polynomial of
degree \(h\) in \(x\), put

\[
 K_h(n)=\left\lfloor\frac nh\right\rfloor\binom h2
                  +\binom{n\bmod h}{2}.
\]

One may replace \(h\) by \(\min(h,n)\). Then

\[
 V_x^0(I)\le(K_h(n)+1)\bigl(n(n-1)-K_h(n)\bigr)
                         =O(hn^3).                         \tag{S2}
\]

For \(h=2\), this is a cubic collective bound for **all** walls
with zero \(x\)-component, with arbitrary coefficients and
nonreduced dependencies. It includes their entire intermediate-cost
range. Cubic order is sharp even when \(x^2\in I\), and also for
radical characteristic-zero point ideals satisfying \(x(x-1)\in I\).

These statements do not control walls with nonzero \(x\)-component.
For unrestricted \(x\)-moment width \(O(n^2)\), (S1) still permits
quartic cost. No unrestricted \(G_3\) improvement, full cubic
polylogarithmic budget, polynomial lower separation or novelty is
asserted.

## One finite polytope retains every actual wall

Let \(P\) be the convex hull of the exponent-sum profiles of all
actual positive initial standard sets of \(I\). Its positive
support function agrees with the usual residue-matroid profile
polytope: each uniquely minimizing basis is a standard set.
Consequently every actual positive wall is dual to the same
exposed profile edge in \(P\). Distinct undirected wall pieces
give distinct edges, including when their primitive directions
agree. This is the established state-polytope dictionary:
[Babson--Onn--Thomas, Theorem 2.4 and the paragraph following its proof](https://arxiv.org/pdf/math/0207135)
explicitly allow the hull of initial staircase profiles itself.

The argument uses this one polytope of the same quotient. It does
not identify a span of compulsory monomials with an algebra ideal.

Every \(n\)-element downset satisfies

\[
 \mu_x+\mu_y+\mu_z\le T:=\binom n2.                        \tag{S3}
\]

Indeed its \(j\)-th monomial in increasing total degree has degree
at most \(j-1\), since the proper divisor chain of a degree-\(d\)
monomial contains a monomial of every degree \(0,\ldots,d-1\).
Sum these bounds.

For each integer \(k\), take the planar polygon

\[
 Q_k=\operatorname{conv}\{(\mu_y(B),\mu_z(B)):\mu_x(B)=k\}.
                                                               \tag{S4}
\]

Empty levels are ignored. A horizontal profile edge at \(x=k\)
remains an exposed edge of \(Q_k\): its two endpoints are present,
and a better point in this polygon would contradict its original
positive exposure. Its exposing planar costs are positive.

For a fixed direction \((0,a,-b)\), the ratio of those planar
costs is fixed. A convex polygon has at most one exposed edge
at this ratio. Hence at most one actual wall of that direction
occurs at any fixed integer \(k\). There are at most \(D+1\)
levels, proving the first part of (S1).

## Perimeter telescoping in every integer level

As the positive planar price ratio increases, the \(y\)-coordinate
decreases and the \(z\)-coordinate increases along the lower convex
chain of \(Q_k\). Each counted horizontal edge has displacement

\[
                  t_e(-a_e,b_e),\qquad
             t_e\in\mathbb Z_{>0},\quad \gcd(a_e,b_e)=1.
\]

Thus its primitive cost \(a_e+b_e\) is at most its total coordinate
variation \(t_e(a_e+b_e)\). The counted edges may omit some edges
of the chain, which only reduces the sum. By (S3),

\[
 \sum_{\substack{e:(q_e)_x=0\\\mu_x(e)=k}}\nu(q_e)
      \le \operatorname{range}_y(Q_k)+\operatorname{range}_z(Q_k)
      \le2(T-k).                                        \tag{S5}
\]

This also holds for a segment or a point. Summing the arithmetic
progression for \(k_-\le k\le k_+\) proves (S1).
This is a collective bound: it sums all directions and their
parallel multiplicities in a level before summing levels.

A degree-\(h\) univariate relation makes every standard exponent
in \(x\) smaller than \(h\). Packing the \(x\)-parallel downset
columns gives \(0\le k_-\le k_+\le K_h(n)\): a column of height
\(\ell\le h\) contributes \(\binom{\ell}{2}\), maximized by
full columns and at most one partial column. This is the
[existing exact moment cap](interleaving_research.md#sixteenth-feasibility-pass-exact-moment-ranges-and-the-remaining-hard-regime).
Maximizing the right side of (S1) on this interval gives (S2).

If \(x\) is nilpotent, there is also the operator bound

\[
 k_+\le \sum_{j\ge1}\operatorname{rank}(M_x^j)\le K_h(n).
                                                               \tag{S6}
\]

For any monomial vector basis, its columns with \(x\)-exponent
at least \(j\) lie in \(x^jA\). Summing these cardinality bounds
gives the first inequality. Jordan block lengths, or the same
convex column packing, give the second. Substituting a smaller
valid bound for \(k_+\) sharpens (S2).

## A linear residual bound for parallel walls at degree two

Assume \(h=2\), and let \(q=(0,a,-b)\) occur. Its compulsory core
\(C_q\) has \(a+b-1\) monomials, all with \(x\)-exponent zero.
Put \(r=n-\nu(q)+1\), the residual rank from the
[fixed-core theorem](wall_coexistence.md).
Every endpoint standard set has at most \(r\) monomials with
positive \(x\)-exponent. Those exponents are all one. It also
has at most \(\lfloor n/2\rfloor\) such monomials, by downward
closure. Therefore

\[
 p_I(q)\le\min\{r,\lfloor n/2\rfloor\}+1.                 \tag{S7}
\]

To apply the integer-level argument, take the convex hull of the
endpoints of all these \(q\)-walls. Every such edge retains its
positive exposure, and only the stated moment levels can occur.
This is a restriction of profiles, not a contraction of an
unjustified algebra quotient. The earlier exact \(r=1\) rigidity
bound can replace (S7) there.

Linear order in \(r\) is attained by the construction below.
This does not assert \(p_I(q)=O(r)\) for arbitrary ideals: the
earlier general planar embeddings give \(\Theta(r\log(r+1))\).
The bounded coordinate here is specifically the zero coordinate
of \(q\).

## Cubic order and linear parallel multiplicity are attained

Use an actual length-\(n\) quotient whose positive fan is the free
prefix fan of the two-layer slab

\[
             \{0,1\}\times\mathbb N_0^2.                 \tag{S8}
\]

Both the explicit nilpotent characteristic-zero realization and
the generic radical two-plane realization are
[already established](nested_jet_matroids.md#balanced-valuations-give-exactly-bounded-height-lattice-prefixes).
For clarity, the nilpotent version for even \(n=2q\) is
\(\mathbb Q[u]/u^n\) with

\[
 \begin{split}
 x&=u^q\exp((Cq+1)u),\\
 y&=u\exp((C+n)u),\\
 z&=u\exp((C+n^2)u),\qquad C=n(1+n+n^2),
 \end{split}                                                   \tag{S9}
\]

where exponentials are truncated. Here \(x^2=0\), and \(y\)
generates the whole algebra. Every \(n\)-element slab downset is
a basis by the existing ordered Hermite argument: if its layer
sizes are \(a\ge b\), then \(a\ge q\); the degree-ordered
monomials in each layer have valuation no larger than their
assigned indices \(0,\ldots,n-1\).

Set \(n=16m\). Independently choose

\[
 1\le k\le2m,\quad 10m\le c\le12m,\quad c-k\text{ odd},
 \qquad a=n-k,\qquad T_{k,c}=(a+c-1)/2.
\]

Take planar costs \(w_y=1,w_z=c\) and

\[
                  w_x=T_{k,c}-k+\tfrac12>0.              \tag{S10}
\]

There are exactly \(a-1\) lower-layer monomials of cost below
\(T_{k,c}\), and exactly two at that cost:

\[
               y^{T_{k,c}},\qquad y^{T_{k,c}-c}z.
\]

There are exactly \(k\) upper-layer monomials below that cost,
namely \(x,xy,\ldots,xy^{k-1}\); their last and first omitted
costs are \(T_{k,c}-1/2\) and \(T_{k,c}+1/2\).
The numerical ranges ensure \(c<T_{k,c}<2c\) and \(k<c\),
so these statements include every monomial in the slab.

Perturbing \(w_z\) on either side swaps just the displayed lower
pair at the selection threshold. Both resulting \(n\)-sets are
actual bases. Hence this is an actual wall with

\[
                   q_c=(0,c,-1),\quad \nu(q_c)=c+1,
                   \quad\mu_x=k.                       \tag{S11}
\]

Different \((k,c)\) give different wall pieces. At each \(c\)
there are \(m\) permitted values of \(k\). Thus

\[
 \begin{split}
 p_I(q_c)&\ge m,\qquad r=n-c\in[4m,6m],\\
 V_x^0(I)&\ge
       m\sum_{c=10m}^{12m}(c+1)
        =m(2m+1)(11m+1)=\Omega(n^3).                    \tag{S12}
 \end{split}
\]

Together with (S2), this proves sharp cubic order for the specified
whole wall family. All walls used in (S12) have intermediate cost
between fixed positive fractions of \(n\).
The same slab fan has radical realizations with \(n/2\) generic
points on each of \(x=0\) and \(x=1\) over \(\mathbb Q\).
No finite-field asymptotic lower assertion is imported from (S9).

## Verification, inputs and the remaining budget

Run:

~~~bash
python3 -B scripts/run_coordinate_slice_budget_checks.py
~~~

The report is results/research/coordinate_slice_budget_checks.json.
Exact complete fans of two compatible growing-thickness hook
deformations test the collective bound with genuinely coupled
nonreduced ranks. The checker uses the existing compatible-tail
construction in the [hook analysis](interleaving_research.md), and
independently replays commuting multiplication matrices and every
carrier successor. Saved capped-jet fans supply independent additional
models. The lower walls are checked by separate free-prefix sorting;
small instances also have exact integer, fraction-free determinant
checks in the actual algebra (S9). Finite counts are not fitted to
an asymptotic law.

The state-polytope correspondence, planar perimeter telescoping,
the exact moment cap, Jordan rank identities and the slab
realizations are established inputs. The planar weighted budget is
proved in the [product-cost note](wall_cost_products.md).
The additional statement here applies that charge to all integer
coordinate-moment levels in one fan, with a sharp all-rank
bounded-coordinate consequence. Its exact priority is unassessed.

The remaining obligation is not to increase these finite
enumerations. Even with \(x^2=0\), walls with nonzero \(x\)-component
are outside (S1)--(S7). For arbitrary three-variable quotients,
the unrestricted moment width can be quadratic. A proof of the
full middle-cost budget must control those directions or replace
the present width parameter by a stronger collective invariant.
The sharp unrestricted endpoint and all publication criteria remain
unchanged.
