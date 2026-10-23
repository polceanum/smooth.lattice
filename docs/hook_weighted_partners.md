# Joint prefix ranks give weighted partner bounds for mixed hooks

The joint input/output flags give a collective bound for two classes of
annihilator exceptions, with arbitrary mixed arm polynomials. The bound
allows the square-root partner example: a lower label at long coordinate
\(Y\) has \(O(at^4/Y)\) dual partners in the middle-price sector.
Its weighted contribution \(Y\) times the number of partners is
\(O(at^4)\). Primal exceptions in the low-price sector have only
\(O(t^4)\) partners per lower label.

The count is of distinct partners and directions. Actual pieces of a
literal wall may recur. A separate established parallel-shadow estimate
controls those pieces once directions have been counted. The other
price orientations, nonscalar uncertainty failures, arbitrary outer
images and unrestricted \(G_3\) remain open.

## Statement

Use the central-only hook with constant axis units
\[
 I=(x^2,y^a,z^a,y^tz^t-xf),\quad a>2t,\quad t\ge2,\quad
 n=4at-2t^2,\quad A=a-t,\quad L=a-2t.
\]
Its rectangular arm operators are
\[
 W_Y=\alpha+N_Y,\quad W_Z=\beta+N_Z,\qquad
 \alpha,\beta\ne0,\quad N_i\in(v)R_A,
 \quad R_A=K[q,v]/(q^A,v^t).
\]
The two mixed polynomials are arbitrary, over every field. Both outer
generator images are zero.

Fix original labels
\[
 E=(0,Y,Z),\quad 3t\le Y<a,\ t\le Z<2t,\qquad
 F=(1,R,S),\quad0\le R<t,\ 2t\le S<a.
\]
Normalize weights as \(w=(\lambda,u,1)\), \(u>1\).
Count a partner when the pair is nonzero and parallel modulo strictly
cheaper original columns on a genuine open cross-layer wall piece.

Let \(P_{\rm pr}(E)\) count partners having a primal relation with
\(N_YH_Y=0\), in \(0<\lambda<t\).
Let \(P_{\rm du}(E)\) count partners having a separating functional
whose reversed Y polynomial satisfies \(N_YK_Y=0\bmod q^L\), in
\(t<\lambda<tu\), with the cheap-overlap condition of (TP11).
Then
\[
 P_{\rm pr}(E)\le4t^4,\qquad
 P_{\rm du}(E)\le16at^4/Y.                              \tag{WP1}
\]
In particular \(YP_{\rm du}(E)\le16at^4\).
These constants are deliberately loose. No regular-prefix or
monomial-times-unit assumption is made.

The earlier [binomial example](hook_exception_partners.md), with
\(t=2\), \(a=2m^2+m+4\), \(Y=m+2\), \(m\ge4\), has at least
\(m-1=\Theta(a/Y)\) such dual partners. Thus the \(a/Y\) dependence
cannot be reduced uniformly at fixed thickness. Optimal dependence on
\(t\) is not asserted.

## A span threshold with bounded column slopes

In a fixed vector space let \(e\) be a fixed vector, let \(C(u)\)
be an increasing family of column sets, and let \(d_j\) be fixed columns
of costs
\[
                         R_ju+S_j,\quad0\le R_j<t.
\]
Define \(b(u)\) as the least cost threshold whose \(d\)-prefix, together
with \(C(u)\), spans \(e\). If \(e\) is already in \(\langle C(u)\rangle\),
it is a loop and is excluded. Before \(e\) is in the full span, there
is no finite threshold.

On its finite, nonloop portion, for \(u_2>u_1\),
\[
                 b(u_2)\le b(u_1)+(t-1)(u_2-u_1).        \tag{WP2}
\]
A representation at the earlier threshold remains available: \(C\)
only grows and each used \(d_j\)'s cost increases by at most the
displayed amount. Once the full span is reached it remains reached;
once \(e\) is a loop it remains a loop. Between changes of \(C\),
the threshold is continuous, being the minimum of finitely many
maximum column costs.

At a ratio avoiding upper-column ties there is one death column:
the column whose addition first spans \(e\). A nonloop parallel pair
\(e,d_j\) modulo the cheaper prefix must use that column.
It changes only when \(C\) changes or two upper costs cross.
This uses joint prefix spans, not just separate column matroids.

## Fix only the own-arm order

For the fixed Y source, put \(p=(Y-t,Z-t)\) and
\(\widehat p=(A-1-p_q,t-1-p_v)\).
Write \(\rho_u(i,j)=iu+j\) for its local monomial weight.
On \(u>1\), the own-arm monomial order changes only at
\[
 \mathcal R_t=\{b/c:1\le c\le b\le t-1\}.
\]
Thus at most \(1+t(t-1)/2\le t^2\) open sectors suffice, independently
of \(a\). Fix one such sector \(J\). All Y input prefixes through
\(p\) and reflected prefixes through \(\widehat p\) are fixed.
The opposite-arm order is allowed to vary throughout \(J\).

Every central bridge \((0,t+i,t+j)\), \(i,j<t\), is cheaper than
\(E\) when \(Y\ge3t\), \(u>1\). After the projections below, the only
remaining lower constraints are columns with first planar coordinate
\(r\in[t,2t)\). Such a column becomes cheaper than \(E\) when
\[
                         u>(s-Z)/(Y-r).
\]
They therefore form an increasing flag \(C(u)\).
In an interval of length \(\Delta\), at most
\[
                              t(Y\Delta+1)             \tag{WP3}
\]
of these events occur, since there are \(t\) possible \(r\)'s.

The upper Z identities have \(0\le R<t\), \(t\le S<a\).
Their distinct crossing parameters in an interval of length \(\Delta\)
number at most
\[
                 \tfrac12t(t-1)\Delta+t-1.             \tag{WP4}
\]
Their positive differences have denominator at most \(t-1\);
for each denominator \(c\) there are at most \(c\Delta+1\)
possible integer numerators. Simultaneous crossings count once.

## Primal exceptions: quotient only by known cheaper columns

Suppose \(p\) is a primal exceptional label. The
[filtered reduction](hook_filtered_annihilators.md) gives a full kernel
relation \(H^*=p+\text{earlier inputs}\) for \(N_Y\).
Contract the following fixed original column set \(Q_J\):

- pure lower Y sources earlier than \(E\);
- upper Y identities whose local weight is at most the weight of \(p\);
- all central bridges.

Every column in this set is strictly cheaper in \(0<\lambda<t\).
Since \(W_YH^*=\alpha H^*\), the relation shows that the image \(e\)
of \(E\) in this fixed vector quotient lies in the upper Z span.
This assertion also follows directly by writing the global source:
its inputs with long index below \(t\) are bridges, and its other
Y inputs and scalar Y images belong to \(Q_J\).
The vector quotient is not declared to be an algebra.

For any covered actual exchange, its primal relation has the same
support restrictions on Y. It therefore projects to a parallel pair
\(e,F\) modulo \(C(u)\) and the cheaper upper Z prefix.
These projected columns are nonzero: only columns already cheaper
in the original matrix were contracted. Thus every covered partner
is a death column for the fixed system above.

Write
\[
                         \Lambda(u)=Yu+Z-b(u).
\]
By (WP2), its increase is at least
\((Y-t+1)(u_2-u_1)\). Hence all candidate ratios for which
\(0<\Lambda<t\) lie in one interval, possibly with jump endpoints,
of diameter at most
\[
                           \Delta_{\rm pr}\le t/(Y-t+1).
\]
Counting (WP3)--(WP4) within it, plus the first death column, gives
at most \(4t^2\) partners in \(J\). Here \(Y\ge3t\).
Summing the at most \(t^2\) own-arm sectors proves the first part
of (WP1). Original eligibility can remove some candidates; it
does not add a death column.

## Dual exceptions: retain the full restricted functional space

Let \(\mathcal S_J\) be the space of reversed Y polynomials
\[
 K\in R_L,\quad N_YK=0\bmod q^L,\quad
             \operatorname{supp}K\subseteq\{\rho\le\rho(\widehat p)\}.
                                                               \tag{WP5}
\]
Keep exactly these Y functionals and every Z functional. This defines
a fixed linear projection of original columns, by evaluation against
a basis of \(\mathcal S_J\) and against the Z coordinate functionals.

Every cheaper pure lower Y column projects to zero: its transpose
product is \(\alpha K\bmod q^L\), and its reflected weight is above
\(\widehat p\). Every cheaper upper Y identity also projects to zero
when \(\lambda>t\). The reflected matching upper identity is already
more expensive than \(E\), and every allowed functional is supported
at least that far up the original row order.

For a covered original separator, cheap overlap gives long degree below
\(L\). Its scalar truncated product and the cheaper pure Y sources
force every coefficient above \(\widehat p\) to vanish; its coefficient
at \(\widehat p\) is nonzero. It therefore belongs to this retained
space, so \(e,F\)
remain nonzero after projection. Original parallelism projects to
parallelism modulo \(C(u)\) and the cheaper upper Z prefix.
This proves the death-column reduction for dual exceptions. A
projected candidate need not give an actual original exchange; only
the necessary direction is used.

For \(t<\Lambda(u)<tu\), both \(\Lambda\) and
\(\Lambda-tu\) are strictly increasing, by (WP2).
All candidate ratios again lie in one interval. Since
\(b(u)\le(t-1)u+a-1\), any candidate satisfies
\[
 u<a/(Y-2t+1).
\]
For two candidate ratios, (WP2) and the price inequalities give
\[
 \Delta_{\rm du}<
 \frac{t\,u_{\max}}{Y-2t+1}
 \le\frac{ta}{(Y-2t+1)^2}\le9ta/Y^2.                   \tag{WP6}
\]
Equations (WP3)--(WP4), plus the first death column, now give at most
\(16at^2/Y\) partners in \(J\). Summing sectors proves (WP1).
This is a general weighted bound; it does not assume that the four
annihilator basis sets determine the needed projections.

## Distinct directions and recurrent pieces

There are \(t\) choices of \(Z\) for each \(Y\). Summing (WP1)
therefore bounds the distinct directions of the covered primal
and dual classes together by
\[
                          O(at^5(1+\log(a+1))).         \tag{WP7}
\]
The short range \(2t\le Y<3t\) can be added separately with
\(O(at)\) directions by its difference box. Transposition gives
the corresponding assertion when the lower long axis is Z and
its coordinate weight is the larger one.

The established [parallel-shadow bound](cross_layer_direction_budget.md)
is \(O(n\log(n+1))\) actual pieces per direction. Consequently
the total covered edge contribution, including subdivisions and
recurrent pieces, is
\[
                          O(n^2t^4(1+\log(n+1))^2).     \tag{WP8}
\]
This does not count the complete fan. In particular, the other
price/axis orientations and the nonscalar branches were not included.
No connectedness of a literal pair's actual eligibility is inferred
from connectedness of the projected candidate interval.

## Verification and priority

~~~bash
python3 -B scripts/run_hook_weighted_partner_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

Checks use original commuting quotient matrices, nonmonomial primal
relations, an arbitrary three-corner dual operator, and the square-root
family. They construct both fixed projections, verify the original
cheaper constraints, enumerate their whole piecewise-linear death
frontiers, and compare actual coherent pairs and adjacent strict minima.
The threshold Lipschitz inequality and event counts are checked
independently on exact field ranks.

Greedy span thresholds and parametric matroid reasoning are classical;
compare [Eppstein, Section 1.2](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf).
The cost-compatible primal projection, restricted dual space and
displayed counting proof are supplied explicitly. Exact priority of
this mixed-hook application remains unassessed. The unrestricted
extremal endpoint and final-paper criteria remain open.
