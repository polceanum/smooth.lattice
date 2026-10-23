# One annihilator label can have square-root many partners

The filtered-annihilator reduction cannot be completed by assigning a
bounded number of opposite-arm partners to each exceptional label.
Even at thickness two, one fixed dual exceptional label participates
in \(\Omega(\sqrt n)\) actual cross-arm exchanges inside one planar
operator-order interval. The opposite arm has regular unit Toeplitz prefixes.

There is also an information obstruction. Replacing that regular unit
by a constant leaves all four labelled annihilator column matroids,
their complete positive greedy-basis profiles, and their kernel spaces
unchanged. Yet the same exchanged label is a loop throughout the
specified sector in the constant case.

These are local obstructions to a proposed aggregate count. They do
not give a superquadratic fan family, improve unrestricted \(G_3\),
or establish exact novelty.

## The actual quotient and its four matroids

Fix \(m\ge3\), put
\[
 A=2m^2+m+2,\quad a=A+2,\quad t=2,\qquad
 n=8A+8=16m^2+8m+24,
\]
and work in characteristic zero or characteristic \(p>2A\).
For a unit \(G(z)\bmod z^A\), take
\[
 \begin{split}
 f_G&=y^2+z^2+y^{m+2}z+yz^2G(z),\\
 I_G&=(x^2,y^a,z^a,y^2z^2-xf_G).                       \tag{BP1}
 \end{split}
\]
This is the flat central-only hook with zero outer images. Its
local arm units and nonscalar operators are
\[
 U=1+q^mv,\quad V=1+vG(q),\qquad
 N_Y=q^mv,\quad N_Z=vG(q)                              \tag{BP2}
\]
in \(R_A=K[q,v]/(q^A,v^2)\). Write \(L=A-2\).

Use the four labelled column matroids \(B,B_P,B_+,B_-\) from the
[filtered-annihilator reduction](hook_filtered_annihilators.md).
On Z, multiplication by \(G\) is an invertible output row map at
lengths \(A,L\). Thus every one of the four represented matroids is
identical to its counterpart for \(G=1\). More explicitly, in the
full map the labels \((i,0)\), \(0\le i<A\), are independent and
all labels with second coordinate one are loops. Restricting the
input or truncating the output gives exactly the corresponding
restricted independent set. Their kernel spaces agree too.
The Y operator is unchanged.

Consequently all four increasing greedy bases on both arms are
independent of \(G\), and are constant at every positive ratio.
Their primal differences are empty. The Z dual difference is empty,
while the Y dual difference is precisely
\[
                 S_D^Y=\{(A-m-2,0),(A-m-1,0)\}.         \tag{BP3}
\]
In particular its size is two, independently of \(m\).

## Many partners in one interval

Take \(G(q)=(1+q)^{A-1}\), and fix the original lower label
\[
                              E=(0,m+2,3).
\]
For \(j=2,\ldots,m\), set
\[
 \begin{split}
 K&=m^2+j,\qquad J=A-1-K,\\
 u_j&=(K-\tfrac12)/m,\qquad
 \lambda_j=(m+1)u_j+1-J,\\
 F_J&=(1,1,J+2).                                      \tag{BP4}
 \end{split}
\]
Then
\[
 m<u_j<m+1,\qquad
 \lambda_j=2j-\tfrac12+(j-\tfrac12)/m,\qquad
                         2<\lambda_j<2u_j.             \tag{BP5}
\]
The costs of \(E,F_J\) agree at \(w=(\lambda_j,u_j,1)\).
Each pair is nonzero and parallel modulo strictly cheaper original
columns, and hence gives an actual coherent wall at nearby generic
points of its line. Its primitive direction is
\[
                       F_J-E=(1,-(m+1),J-1).           \tag{BP6}
\]
The \(m-1\) partners and directions are distinct. They lie outside
the exact same-arm/bridge inventory of the bridge reduction.
Indeed at thickness two a first-component-one comparison in that
inventory with \(q_y<0<q_z\) has either \(|q_y|\le3\),
\(q_z\le1\), or \(|q_y|\ge a-1\), by its arm and three bridge
boxes. Here \(|q_y|=m+1\) and
\(m^2\le q_z\le m^2+m-2\), which violate all three alternatives.

All these weights lie in
\[
               \mathcal S_m=\{m<u<m+1,\ 2<\lambda<2u\}.
\]
For the two-row rectangles, internal order changes occur only at
positive integers or their reciprocals. Thus \((m,m+1)\) is one
common order interval for the four operator scans. In fact their selected basis
sets remain constant even beyond it.

### Reduction at the pair's actual cheaper threshold

Write \(z_i\), \(0\le i<A\), for the upper Z reference row of label
\((1,1,i+2)\), and \(g_i=q^iG\bmod q^A\) in their span.
At an open interval around (BP4), keep
\[
 K<mu+1<K+1,\quad
 \lambda=(m+1)u+1-J,\quad 2<\lambda<2u.                 \tag{BP7}
\]
At the pair's cost \(\tau=(m+2)u+3\), the relevant cheaper Z space is
exactly
\[
       C_{J,K}=\langle z_0,\ldots,z_{J-1},
                            g_1,\ldots,g_K\rangle.     \tag{BP8}
\]
After eliminating other rows, \(F_J\) is \(z_J\) and \(E\) is
\(-g_0\). This is an original-column reduction, not a replacement
of a contracted vector space by an algebra quotient.

Here are the column details. The bridge \((0,2,2)\) has image
\[
 e_Y+g_0+v_{(1,2,0)}+v_{(1,0,2)},
\]
where \(e_Y\) is the upper Y row \((1,m+2,1)\). The bridge and its
two displayed helpers are cheaper than \(E\). The only other
columns with coordinate \(e_Y\) are \(E\) itself and that upper
identity, whose cost is greater by \(\lambda-2>0\). Eliminating
the bridge therefore sends \(E=e_Y\) to \(-g_0\).

The central Z sources \((0,2,2+i)\) are cheaper precisely when
\(i<mu+1\). For \(1\le i\le K\), their short-coordinate-zero upper
identity helpers are cheaper because \(\lambda<2u\), leaving
\(g_i\). The sole additional Y helper, when \(i=1\), is cheaper
as well. Central Z sources with short coordinate one are parallel
to upper \(z_i\)'s that are already cheaper. Upper Z identities
in this row are cheaper precisely for \(i<J\).

Other Y columns have no Z component, except the bridge
\((0,3,2)\), whose Z component is the already cheaper \(z_0\).
Their Y coordinates do not contain \(e_Y\). Thus they supply no
additional relation on this Z space. Short-coordinate-zero Z
rows can be eliminated separately. These observations prove (BP8).

### The two classical minors

Put \(D=A-1\). The binomial Toeplitz evaluation is
\[
 \det[\binom D{r+i-j}]_{0\le i,j<s}
       =\prod_{i=0}^{s-1}
          \frac{(D+i)!\,i!}{(r+i)!\,(D-r+i)!}.          \tag{BP9}
\]
The product is positive over the rationals in the cases below.
Every factorial argument is less than \(2A\), so it is also
nonzero in characteristic \(p>2A\).
This is the classical evaluation in
[Krattenthaler, Section 2.3, (2.17)](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11).
The same formula with size \(s_0-r\), \(0\le r<s_0\), at
\(s_0=A,L\) proves the claimed regularity of every unit Toeplitz
prefix index; the empty endpoint determinant is one.

Since \(J+K=A-1\), after eliminating \(z_0,\ldots,z_{J-1}\)
there are \(K+1\) rows. Columns \(g_1,\ldots,g_K\) have rank \(K\):
their last \(K\) rows give (BP9) with \((r,s)=(J,K)\).
Adding \(g_0\) gives the nonzero determinant with
\((r,s)=(J,K+1)\). Adding \(z_J\) also gives full rank by the
first minor. Hence (BP8) has codimension one and both \(g_0,z_J\)
are nonzero modulo it. This proves the claimed parallelism and
nonzero quotient rank.

All inequalities used are strict on an open segment of the line
(BP7). Avoiding the finite other label-tie lines yields genuine
adjacent strict minima. A multiple intersection at the specific
rational witness does not need to be assumed generic.

## The same exceptional label at every wall

Take a functional on the Z row space annihilating \(C_{J,K}\)
and nonzero on \(g_0,z_J\). Give \(e_Y\) the opposite value of
\(g_0\), and every other Y row value zero. The bridge evaluations
cancel, and the preceding column reduction shows that this
functional annihilates all cheaper original columns.

The reversed Y polynomial is a nonzero scalar times
\[
                         q^{A-m-1}.                   \tag{BP10}
\]
It is killed by \(N_Y\) at length \(L\), but its full product is
a nonzero multiple of \(q^{A-1}v\). Its leading label is the
same member \((A-m-1,0)\) of (BP3) at every one of the \(m-1\)
walls. The common overlap rows are cheaper: from \(\lambda<2u\)
and \(m\ge3\), their maximum cost is strictly below \(\tau\).
Thus these are precisely dual exceptions covered by (FA8).

This gives \(\Omega(\sqrt n)\) partners for one exceptional label
in one order interval at fixed thickness two. A bound depending
only on thickness for the number of partners per label is false,
even with a regular opposite-arm multiplier.

## The square-root scale is sharp within this binomial class

The same prefix calculation gives a matching local upper bound.
Allow arbitrary integers \(A>\nu+2\), \(\nu\ge3\), replace \(m\)
in the Y multiplier by \(\nu\), and keep \(G=(1+q)^{A-1}\).
Fix \(E_\nu=(0,\nu+2,3)\) and the price sector
\(u>1,\ 2<\lambda<2u\). In one unit interval
\(r<u<r+1\), \(r\) a positive integer, it has at most
\[
                  1+\min(\nu,A/\nu)\le1+\sqrt A       \tag{BP11}
\]
distinct upper Z partners.

To see this, set \(K=\lfloor\nu u+1\rfloor\), away from the
internal source-order ties. For a tied \(F_J\),
\[
 \lambda=(\nu+1)u+1-J.
\]
The original-column reduction above applies unchanged. The cheaper
space has rank \(\min(A,J+K)\). If \(J+K\ge A\), the exchanged
columns are loops. If \(J+K\le A-2\), they are independent modulo
the cheaper space: eliminate \(z_J\) too and use (BP9) on rows
\(J+1,\ldots,J+K+1\), columns \(g_0,\ldots,g_K\).
Thus an eligible pair requires and, under the displayed strict
conditions, is equivalent to
\[
 J+K=A-1,\qquad
             \lambda=(\nu+1)u+K+2-A.                  \tag{BP12}
\]
The price sector gives \(J>0\). If \(K\ge A-1\), the cheaper
space spans every Z row, so \(E_\nu\) is a loop. Every eligible
pair therefore has the nondegenerate ranges above.

There are at most \(\nu\) values of \(K\) in a unit interval.
The inequalities \(2<\lambda<2u\) and
\((K-1)/\nu<u<K/\nu\) additionally imply
\[
 \frac{A\nu}{2\nu+1}<K<
                   \frac{A\nu-\nu-1}{2\nu-1}.
\]
This interval has length less than
\(2A\nu/(4\nu^2-1)\le A/\nu\), so it has at most
\(1+A/\nu\) integer values. Each \(K\) fixes one \(J\), proving
(BP11). At a source-order endpoint, a genuine cross-layer wall
extends to one side unless it is only a multiple intersection;
it supplies no additional full wall with a new \(J\).

The family (BP4) attains \(m-1\) partners with \(\nu=m\) and
\(A=2m^2+m+2\). Hence the worst local order in this binomial
class is \(\Theta(\sqrt n)\). This is not an upper bound for
arbitrary mixed operators or for the complete fan.

## Identical annihilator data can suppress all those exchanges

For the control \(G=1\), all four labelled matroids and kernel
spaces above are unchanged. Nevertheless \(E\) is a loop modulo
cheaper columns throughout \(\mathcal S_m\). The exact relation is
\[
 v_E=v_{(0,2,2)}-v_{(1,2,0)}-v_{(1,0,2)}-v_{(1,1,2)}. \tag{BP13}
\]
All four columns on the right are strictly cheaper when
\(2<\lambda<2u\), \(u>m\), \(m\ge3\).
Thus no exchange of \(E\) with any partner occurs in this sector.

The missing information is the position of the multiplication image
relative to the upper coordinate flags. Output row equivalence
preserves all annihilator-column data but does not preserve those
flags. Joint input/output prefix ranks, such as the two minors in
(BP9), retain this information. The earlier
[real marginal-rank examples](real_marginal_rank_fans.md) already show
a different information loss for four planar quotient rank systems;
this is a specific obstruction for the new annihilator reduction.

## Verification and scope

~~~bash
python3 -B scripts/run_hook_exception_partner_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker reconstructs original commuting quotient columns, verifies
the unit output equivalence and the complete four basis sets, compares
integer determinants with the classical products, checks the entire
remaining Z row space modulo cheaper columns, constructs an original
separator, and checks adjacent strict minima. It also verifies the
control identity and its cheaper helpers. No complete fan count is
inferred from the selected walls.

The binomial determinant and unit row equivalence are classical.
Exact priority of this local obstruction remains unassessed.
A successful collective count must charge the changing joint prefix
ranks and partner interactions rather than just the sizes of the
four basis differences. The unrestricted extremal target, polynomial
lower separation, arbitrary outer images and final-paper criteria
remain open.
