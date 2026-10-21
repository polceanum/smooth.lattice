# Identical planar ranks can hide a growing mixed-wall grid

The remaining bounded-coordinate directions are
\((a,b,-c)\) and \((a,-b,c)\), with \(0<a<h\) and \(b,c>0\).
Separate planar rank flags and the Jordan type of multiplication
by the bounded coordinate do not determine their coherent fan,
even over one fixed field.

Here is an explicit infinite pair of actual quotients with \(x^2=0\).
At the same lengths \(n_k=\Theta(k^3)\), their complete planar
column matroids, planar elimination and colon ideals, and \(x\)-Jordan
types agree. Their planar quotients modulo \(x\) agree as well.
Nevertheless their full positive initial-ideal counts
are \(2\) and \((k+1)^2\). The second fan has no zero-\(x\)-component
walls and no walls where \(x\) is the singleton sign. All its walls
belong to the still unresolved sign patterns.

The construction uses classical monomial quotients and independent
socle identifications. It supplies a failure of a proposed data
reduction, not a supergeneric lower family, a new unrestricted
bound, or a novelty certificate. Both families are in previously
controlled binomial and linear-box classes.

## The same two diagrams, with different identifications

Fix \(k\ge1\), put \(m=2k\), and for \(0\le i<m\) define
\[
 r_i=m-i,\qquad z_i=2m r_i+r_i^2,\qquad
 u_i=y^{2i}z^{z_i},\qquad v_i=y^{2i+1}z^{z_i}=yu_i .
                                                               \tag{P1}
\]
Let \(D_1\) be the union of the planar divisor sets of the \(u_i\),
and \(D_0\) that of the \(v_i\). They are downsets with
\(D_1\subset D_0\). Set
\[
 D=(\{0\}\times D_0)\ \sqcup\ (\{1\}\times D_1),
 \qquad J_D=\langle\hbox{monomials outside }D\rangle .
                                                               \tag{P2}
\]
All \(xu_i,v_i\) are distinct maximal monomials of \(D\).
Indeed \(z_i\) strictly decreases with \(i\); no \(v_i\)
divides any \(u_j\). Thus these \(2m\) classes are nonzero
independent socle vectors of \(K[x,y,z]/J_D\).

For a permutation \(\pi\) of \(\{0,\ldots,m-1\}\), define
\[
 I_\pi=J_D+(xu_i-v_{\pi(i)}:0\le i<m).                    \tag{P3}
\]
The \(m\) relations span an \(m\)-dimensional ideal in the socle.
Consequently
\[
 n_k=|D_0|+|D_1|-m
       =\frac{16m^3+9m^2+11m-3}{3}
       =\frac{128k^3+36k^2+22k-3}{3}.                    \tag{P4}
\]
This is an integer and is positive already at \(k=1\), where it
equals \(61\). Every field is allowed, including characteristic two.
The quotient is cyclic, local and finite; its three multiplication
operators commute and \(x^2=0\).

For every global order, a complete standard set consists of
\[
 D\setminus\{xu_i,v_i:0\le i<m\}
 \quad\hbox{and the cheaper member of every pair }
                       \{xu_i,v_{\pi(i)}\}.              \tag{P5}
\]
To prove this, the monomial ideal \(J_D\) and the larger
member of each pair are contained in the initial ideal.
Their proposed complement is a downset of size \(n_k\).
The colength calculation proves equality of initial ideals.
No independent-choice assumption beyond the actual weight
inequalities is being made.

## What all separate planar rank queries see

Write \(J_0,J_1\subset K[y,z]\) for the monomial ideals
complementary to \(D_0,D_1\), respectively. For every \(\pi\),
\[
 I_\pi\cap K[y,z]=J_0,\qquad
                  (I_\pi:x)\cap K[y,z]=J_1.              \tag{P6}
\]
All lower-layer monomials in \(D_0\) remain independent:
each identified upper tip has its own different lower tip.
The same holds for all upper-layer columns indexed by \(D_1\).
Every other monomial in either layer is zero.
This proves (P6), and in fact the complete rank identities
\[
 \operatorname{rk}\{y^az^b:(a,b)\in S\}=|S\cap D_0|,
 \quad
 \operatorname{rk}\{xy^az^b:(a,b)\in S\}=|S\cap D_1|       \tag{P7}
\]
for every finite set \(S\). Hence all planar rank flags, for
every planar cost order and every truncation rank, agree.
These are equal labelled matroids, not just equal dimensions
or equal planar fan counts.

The third planar endpoint used by the
[square-zero extension presentation](square_zero_extension_fans.md)
is also fixed:
\[
 (I_\pi+(x))/(x)=J_0+(v_0,\ldots,v_{m-1}).                \tag{P6a}
\]
Its diagram is \(D_0\setminus\{v_i\}\).
Thus even all three labelled planar systems leave the
identification in (P3) undetermined. These are nonflat extensions:
the diagram in (P6a) is strictly larger than \(D_1\).

Put \(q_1=|D_1|\). Multiplication by \(x\) has rank \(q_1\),
since its nonzero images are exactly the distinct upper-layer
classes. Its square is zero. Its Jordan type is therefore
\[
                  2^{q_1}1^{\,n_k-2q_1},
 \quad q_1=\frac{8m^3+7m-3}{3},\quad
                       n_k-2q_1=3m^2-m+1.               \tag{P8}
\]
This also is independent of \(\pi\).
The information which varies is the joint relation between the
two planar spaces: their intersection has dimension \(m\),
and (P3) changes its identification.

## Two cones versus a grid containing only the missing directions

For \(\pi=\mathrm{id}\), all comparisons in (P5) have the
same primitive direction \((1,-1,0)\).
Thus the whole positive fan has exactly two cones, one wall,
and total primitive wall cost
\[
                  g(I_{\rm id})=2,\qquad V(I_{\rm id})=2. \tag{P9}
\]

Now let \(\pi\) interchange \(2j\) and \(2j+1\), for
\(0\le j<k\). Put
\[
             L_j=z_{2j}-z_{2j+1}=8k-4j-1 .
\]
The corresponding primitive directions are
\[
             q_j^-=(1,-3,L_j),\qquad
             q_j^+=(1,1,-L_j).                           \tag{P10}
\]
They are all primitive and distinct. Normalize \(w_z=1\).
Their wall lines are
\[
                   w_x=3w_y-L_j,\qquad w_x=L_j-w_y .
                                                               \tag{P11}
\]
Within each family the lines are parallel. Every line of one
family meets every line of the other at
\[
       w_y=(L_i+L_j)/4,\qquad w_x=(3L_i-L_j)/4>0.         \tag{P12}
\]
Indeed \(4k+3\le L_j\le8k-1\), so
\(3\min_jL_j>\max_jL_j\).
There are \(k^2\) distinct crossings and no triple crossings.
Each of the \(2k\) lines meets the positive quadrant and is
split into \(k+1\) wall pieces.
Adding the parallel families successively gives exactly
\[
 g(I_\pi)=(k+1)^2,\qquad
              |\mathcal E(I_\pi)|=2k(k+1).                \tag{P13}
\]
Formula (P5) identifies this arrangement with the complete
actual positive fan; each region has a different standard set.

The two costs at index \(j\) are
\(\nu(q_j^+)=L_j+3\) and \(\nu(q_j^-)=2L_j+4\).
Since \(\sum_jL_j=6k^2+k\),
\[
 V(I_\pi)=(k+1)\sum_j(3L_j+7)
                         =2k(k+1)(9k+5)=\Theta(n_k).      \tag{P14}
\]
Every wall has nonzero \(x\)-component and opposite nonzero
\(y,z\) signs. Therefore both previously controlled budgets
\(V_x^0\) and \(V_x^{\rm opp}\) are exactly zero here.
In particular the remaining walls can cross independently
without any wall in either controlled family.

## Consequence and the next joint-rank requirement

No procedure receiving only the two labelled planar rank oracles
(P7), the planar quotient modulo \(x\), and the \(x\)-Jordan type
can recover the full fan or its
cost, regardless of its number of queries. The two inputs are
identical. Their fan-count ratio is
\(\Theta(n_k^{2/3})\), and their total-cost ratio is
\(\Theta(n_k)\).
In particular a common numerical estimate for these inputs
cannot approximate both fan counts within a symmetric
\(o(n_k^{1/3})\) multiplicative factor, or both costs within
\(o(n_k^{1/2})\). This is an information restriction, not a
running-time lower bound for algorithms given the full ideal.

The [two-list rank identity](rational_shadow_bound.md#two-lists-give-a-partial-matching)
retains the missing information in the joint prefix ranks
\[
 R(i,j)=\operatorname{rk}
       \bigl(\hbox{first }i\hbox{ lower columns},
                       \hbox{first }j\hbox{ upper columns}\bigr).
\]
Neither those ranks nor their relative matching are determined
by the two marginal systems. A collective argument must retain
this coupling, or an equally faithful invariant. Neither a
noncrossing assumption for the remaining wall lines nor charging
all their crossings to zero-component walls is valid.

This family is not a candidate for the original power separation.
Its count is only \(\Theta(n_k^{2/3})\). It is binomial, and its
common monic box has sides \(2,2m,3m^2+1\), hence volume
\(O(n_k)\). The previously proved class bounds already apply.
Larger enumeration in this family will not address the missing
unrestricted budget.

## Verification and prior-art boundary

Run

~~~bash
python3 -B scripts/run_planar_layer_obstruction_checks.py
~~~

The report is results/research/planar_layer_obstruction_checks.json.
It constructs the actual quotient columns and commuting multiplication
maps, checks both planar kernels and the Jordan ranks, independently
recovers every full positive cone from border normal forms, and
extracts adjacency directions and costs from the resulting profiles.
Exact cone coverage and rational crossings distinguish the full
fan calculation from selected weights.

Monomial quotients, socle identifications, binomial residue classes,
the state-profile dictionary and line-arrangement counting are
established tools. The
[primitive-wall note](primitive_wall_spectrum.md) records the
relevant edge-ideal and socle comparison, and the
[two-list note](rational_shadow_bound.md) records the established
interleaving context and its deferred precise comparison with Gusfield.
The [extension-minor presentation](square_zero_extension_fans.md)
correctly retains the map between its fixed planar endpoints; the
present pair quantifies the failure of discarding that map.
The particular rank-indistinguishable pair and quantitative
information consequence are proved here; their exact priority
is unassessed. No claim of a new general gluing operation or
completion of the original extremal goal is made.
