# Fan-preserving padding and an obstruction to absolute-volume reuse charges

Every finite quotient admits explicit univariate thickenings with the
same full positive initial fan and a fixed translation of its state
polyhedron. Formula (F3) gives the initial ideals for every global
order, and (F4) lifts a universal Gröbner basis with \(M\) polynomials
to one with at most \(M+d\), in \(d\) variables. Lengths can grow
arbitrarily while these fans remain fixed.

A same-pair crossing count cannot be bounded by quotient length divided
by the divisor-box volume of either upper monomial, even with any fixed
polylogarithmic allowance. The counterexamples have \(x^2=0\), work over
every field, and retain crossings in the unresolved intermediate
primitive-cost range. At every marked crossing the four incident bases
have a common independent set of size \(n-2\).

This is a proof-route obstruction. The constructions are still
independent socle circuits, a class already known here to have
subgeneric fans. They give no polynomial separation for \(G_3\),
new unrestricted upper bound, or novelty certificate.

## Exact padding of an arbitrary socle relation space

Let \(D\subset\mathbb N_0^d\) be a finite nonempty downset, let
\(E\subseteq\max D\), and let \(W\subset K^E\) be any relation
space of dimension \(\rho\). Over any field set
\[
 I=J_D+(f:f\in W),\qquad n=|D|-\rho .
                                                               \tag{P1}
\]
The support rows may overlap and cancel. They need not be disjoint
circuits. All monomials in \(E\) are in the socle of the monomial
quotient, so \(W\) itself is an algebra ideal there.

For an arbitrary integer \(v\ge0\), form
\[
 D_v=D+[0,v],\quad E_v=v+E,\quad
 I_v=J_{D_v}+(x^vf:f\in W),\quad n_v=|D_v|-\rho ,
                                                               \tag{P2}
\]
where the box and the sum are sets of lattice points.
Equivalently \(D_v\) is the divisor closure of \(v+D\).
Its maximal points are precisely \(v+\max D\); thus the translated
relation space is again a socle ideal of the same rank.

Let \(k=|E|-\rho\). The residue-column matroid of (P1) is the
direct sum of the coloops \(D\setminus E\) and the rank-\(k\)
quotient matroid \(K^E/W\). The same decomposition holds after padding,
with identical quotient columns on translated labels. Every monomial
vector basis therefore corresponds by
\[
 \Psi_v(B)=(D_v\setminus E_v)\ \cup\ (v+(B\cap E)).
                                                               \tag{P3}
\]
These are all downsets, since only maximal points can be omitted.
Translation preserves the order of the \(E\)-labels for every global
monomial order, so (P3) also matches all initial standard sets.

Writing \(\mu(B)=\sum_{\alpha\in B}\alpha\), the exact profile identity is
\[
 \mu(\Psi_v(B))=\mu(B)+\Delta_v,\qquad
 \Delta_v=\mu(D_v)-\mu(D)-\rho v.                       \tag{P4}
\]
Indeed every variable basis uses exactly \(k\) of the \(|E|\) tip
labels. Consequently the complete finite basis-profile polytopes,
and their state polyhedra after adding the positive orthant, differ
only by this fixed translation.

The entire positive fan, its full cones and adjacency, all primitive
wall directions, their integer edge multiples, total primitive cost,
and incidence of compact exposed faces are unchanged. The initial
ideals themselves have the different labels (P3), and the lengths
usually differ. No invariance of reduced-basis degrees, term counts,
or minimal monomial generators is asserted.

If \(C\) is the intersection of all initial standard sets, then
\[
 C_v=(D_v\setminus E_v)\cup(v+(C\cap E)),\qquad
                         n_v-|C_v|=n-|C|.              \tag{P5}
\]
The same formula holds for the intersection over the bases incident
to any one exposed face. These are vector-space cores. They need
not be algebra ideals.

## Exact length without enumerating the padded box

For \(S\subseteq\{1,\ldots,d\}\), let \(\pi_S D\) delete the
coordinates in \(S\). The empty-dimensional projection of a nonempty
set has size one. Then
\[
 |D_v|=\sum_{S\subseteq[d]}
              |\pi_S D|\prod_{i\in S}v_i.              \tag{P6}
\]
Adding a length-\(v_i\) lattice interval along coordinate \(i\)
increases the size of every nonempty fiber by \(v_i\), giving
\(|D+[0,v_i]e_i|=|D|+v_i|\pi_{\{i\}}D|\).
Iterating this identity proves (P6). It is an elementary downset
count, not a new general mixed-volume theorem.

## Monomial padding also preserves propagating relations

The invariance does not require socle support. Let \(D\) be any finite
downset and suppose \(J_D\subseteq I\). This includes arbitrary
origin-supported finite quotients, since their coordinates are nilpotent.
For \(D_v=D+[0,v]\), set
\[
                         I_v=J_{D_v}+x^v I.             \tag{G1}
\]
There is an injection of polynomial modules
\[
 K[x]/J_D\ \xrightarrow{\cdot x^v}\ K[x]/J_{D_v}.
\]
Indeed \(v+\alpha\in D_v\) if and only if \(\alpha\in D\).
Its image has monomial labels \(v+D\), and the other labels
\(H=D_v\setminus(v+D)\) are independent of this image. The full
relation space of \(I_v/J_{D_v}\) is exactly the image of \(I/J_D\),
including all nonconstant multiples.

Thus every monomial vector basis corresponds by
\[
 B\longmapsto H\cup(v+B),\qquad
 n_v=n+|D_v|-|D|.                                      \tag{G2}
\]
Such a new basis is a downset precisely when \(B\) is a downset.
The finite basis-profile polytopes differ by
\[
                  \mu(D_v)-\mu(D)-(|D|-n)v.             \tag{G3}
\]
All the fan, cost and common-core conclusions above follow. This is
an exact module embedding with added coloop columns; it is not a
claim that the two quotient algebras are isomorphic.
Equivalently, for every global monomial order,
\[
                 \operatorname{in}(I_v)=
                     J_{D_v}+x^v\operatorname{in}(I).   \tag{G4}
\]

## Univariate padding applies to every finite quotient

There is an all-field version without an origin-support hypothesis.
Let \(I\subset K[x_1,\ldots,x_d]\) have positive finite length \(n\).
Choose monic univariate annihilators \(f_i(x_i)\in I\), of degrees
\(d_i>0\). Such annihilators always exist by finite-dimensional
linear algebra. Write
\[
 D=\prod_i[0,d_i-1],\quad S=\prod_i d_i,\quad
 F_v=\prod_i f_i^{v_i},\quad
 I^{[v]}=(f_i^{v_i+1}:1\le i\le d)+F_v I
                                                               \tag{F1}
\]
for integers \(v_i\ge0\). Put \(s_i=d_iv_i\),
\(D^{[v]}=\prod_i[0,d_i(v_i+1)-1]\), and
\(H=D^{[v]}\setminus(s+D)\). Then
\[
 \begin{split}
 n^{[v]}&=n+S\left(\prod_i(v_i+1)-1\right),\\
 B^{[v]}&=H\cup(s+B),\\
 \mu(B^{[v]})&=\mu(B)+
        \mu(D^{[v]})-\mu(D)-(S-n)s.                    \tag{F2}
 \end{split}
\]
Here \(B\) runs through the initial standard sets of \(I\), and
\(B^{[v]}\) runs through exactly those of \(I^{[v]}\), with identical
full positive cones. In particular these ideals have the same
number of initials, primitive wall directions, edge multiples,
and total primitive wall cost. Their state polyhedra are translates.
Their global common-core deficits, and the deficits at corresponding
faces, agree.
More explicitly, for every global monomial order,
\[
 \operatorname{in}(I^{[v]})=
 (x_i^{d_i(v_i+1)}:1\le i\le d)
                         +x^s\operatorname{in}(I).     \tag{F3}
\]

To prove this, work first modulo \(J=(f_i^{v_i+1})_i\).
Its polynomial basis consists of
\[
 Q_{b+du}=x^b\prod_i f_i^{u_i},\qquad
           0\le b_i<d_i,\quad0\le u_i\le v_i .
 \]
Successive monic division proves the basis statement in every
characteristic. Multiplication by \(F_v\) embeds \(K[x]/(f_1,\ldots,f_d)\)
as the block \(u=v\); its kernel before passing to this quotient is
exactly \((f_1,\ldots,f_d)\).
The relations added in (F1) lie solely in this block, with quotient
columns identical to the original monomial columns \(x^b\bmod I\).
All other \(Q\)-columns are coloops.

Each \(Q_\alpha\) is \(x^\alpha\) plus proper-divisor terms.
The unitriangular change preserves positive weight prefixes and
global monomial-order prefixes. Hence the polynomial-column greedy
labels are the original monomial greedy labels. In the block \(u=v\)
the common price shift \(w\cdot s\) leaves the original greedy order
unchanged; every halo label is always selected. This proves (F2).
There is also a direct leading-ideal check. The right side of (F3)
is contained in the left side by the leading terms of the supplied
generators in (F1). Its complement is exactly the set in (F2),
of size \(n^{[v]}\). Equality of colengths proves (F3).

Unlike (G1), this argument does not assert equality of the full
monomial vector-basis matroids outside positive objectives.
Arbitrary labelled subsets can change rank under the unitriangular
polynomial basis. The positive state-polyhedron and fan conclusions
follow from the initial-standard-set correspondence.
Splitting, separability and locality are unnecessary.
For example start from \((x-1,y-2,z-3)\), choose
\(f_x=(x-1)^2,f_y=y-2,f_z=z-3\), and \(v=(1,0,0)\).
The result is \(((x-1)^3,y-2,z-3)\).
Its ordinary columns \(x^2,x^3\) are independent, whereas the
polynomial columns \(Q_{(2,0,0)}=(x-1)^2\) and
\(Q_{(3,0,0)}=x(x-1)^2\) coincide in the quotient.
This guard works also in characteristic two.

Formula (F1) gives a constructive realization with the supplied
generators of \(I\), multiplied by the single factored polynomial
\(F_v\), and \(d\) additional polynomials. Formula (F2) represents
the padded fan by the original fan, one fixed profile translation
and an explicit common halo. Expanding all new standard sets, degrees,
term counts and coefficient bits has its own output cost. This
operation preserves fan complexity; it does not amplify it.
If \(G\) is a Gröbner basis of \(I\) for any fixed global order,
then
\[
            \{f_i^{v_i+1}:1\le i\le d\}\ \cup\ F_vG     \tag{F4}
\]
is a Gröbner basis of \(I^{[v]}\) for that order, by (F3).
In particular a universal Gröbner basis with \(M\) polynomials
lifts to one with at most \(M+d\) polynomials.
These lifted bases need not be reduced. This is a factored
generator-count consequence, with the same expansion and bit-cost
limitations just stated.

## Intermediate-cost polynomial reuse with large divisor boxes

Use the [existing exact recrossing family](disjoint_socle_circuits.md),
with \(s\ge1\), \(B_0=4s^2\), and tips
\[
 \begin{split}
 p_i&=(0,2i+\mathbf1_{i\text{ even}},
                   B_0-i^2+\mathbf1_{i\text{ odd}}),\quad 1\le i\le2s,\\
 e&=(1,1,B_0),\qquad f=(1,0,B_0+1).
 \end{split}                                                \tag{P7}
\]
Let \(D\) be their divisor downset and impose the two disjoint
relations \(x^e=\sum_{i\ {\rm even}}x^{p_i}\) and
\(x^f=\sum_{i\ {\rm odd}}x^{p_i}\). The original length is
\[
 n_s=(32s^3+30s^2+10s+12)/3.
\]
Pad by \(v=(0,A,B)\). The projection counts in (P6) are
\[
 |\pi_yD|=8s^2+4,\quad |\pi_zD|=4s+4,\quad
                         |\pi_{\{y,z\}}D|=2,
\]
so the new length is exactly
\[
 n_{s,A,B}=n_s+(8s^2+4)A+(4s+4)B+2AB.                 \tag{P8}
\]

The two translated upper labels have the same threshold functions as
before, in weights \(u=w_y/w_z\), \(t=w_x/w_z\):
\[
 T_e(u)=\max_{i\ {\rm even}}(2iu-i^2),\qquad
 T_f(u)=\max_{i\ {\rm odd}}(2iu-i^2).
                                                               \tag{P9}
\]
There are exactly \(2s-1\) persistent ordinary crossings of this
same upper pair, at
\[
                  (u_i,t_i)=(i+\tfrac12,i^2+i),
                          \quad1\le i<2s.              \tag{P10}
\]
One parity has the unique active index \(i\), the other \(i+1\),
in a neighborhood of this crossing. To verify all other indices
are strictly below the two maxima, use
\[
 t_i-(2ju_i-j^2)=(j-i)(j-i-1)\ge0
                         \qquad(j\in\mathbb Z),         \tag{P11}
\]
with equality only at \(j=i,i+1\).
All four upper-membership patterns occur locally. Their common basis
intersection has size \(n_{s,A,B}-2\). Globally it has size
\(n_{s,A,B}-2s\). The complete fan still has exactly \(4s\) cones.

The primitive directions of the two crossing walls are
\[
 q_i=(1,-2i,i^2),\quad q_{i+1}=(1,-2(i+1),(i+1)^2),
 \qquad \nu(q_i)=2i^2+2i+1.                            \tag{P12}
\]
Both have opposite nonzero \(y,z\) signs. Padding has changed
neither their primitive costs nor their profile edge lengths.

Take first \(A=s,B=0\). Then
\[
 \begin{split}
 N_s&=(56s^3+30s^2+22s+12)/3,\\
 |\operatorname{Div}(v+e)|&=2(s+2)(4s^2+1),\\
 |\operatorname{Div}(v+f)|&=2(s+1)(4s^2+2).
 \end{split}                                                \tag{P13}
\]
Both divisor volumes divided by \(N_s\) tend to \(3/7\).
Nevertheless the same upper pair has
\(\Theta(N_s^{1/3})\) crossings. For a positive fraction of the indices
these walls have cost \(\Theta(N_s^{2/3})\), hence lie strictly in
\[
                   N_s^{3/5}<\nu(q)\le N_s-N_s^{3/7}    \tag{P14}
\]
for all sufficiently large \(s\).
Thus a proposed bound by \(N_s\) divided by either absolute
divisor volume, with any fixed power of \(\log N_s\), is false
even in the intermediate range. Near-full actual common face cores,
persistence and convexity of the two threshold curves do not repair it.

There is also a version whose divisor volumes occupy asymptotically
all of the length while retaining intermediate cost. Take
\[
                 s=k^4,\qquad A=k^5,\qquad B=0.
 \]
Then \(N=\Theta(k^{13})\), the two divisor-volume ratios tend to one,
and the upper pair has \(\Theta(N^{4/13})\) crossings. A positive
fraction have cost \(\Theta(N^{8/13})\), whose exponent is strictly
between \(3/5\) and one. The global common-core deficit is exactly
\(2k^4\), while the deficit at each marked face is two.

## A sharp support budget within the circuit class

For two disjoint socle circuits in two \(x\)-layers, each with one
distinguished \(x\)-exponent-one label and respectively \(m_1,m_2\)
\(x\)-exponent-zero labels, their positive
threshold functions are maxima of at most \(m_1,m_2\) affine
functions. Each affine piece occupies one interval. The combined
breakpoints make at most \(m_1+m_2-1\) open intervals, and on each
interval two distinct affine functions have at most one ordinary
crossing. Therefore persistent ordinary crossings of the fixed upper
pair number at most
\[
                            m_1+m_2-1.                 \tag{P15}
\]
The parity construction has \(m_1=m_2=s\) and attains this bound.
Equivalently, its count is the combined circuit rank minus one.
This is elementary envelope counting. It does not extend to arbitrary
square-zero extensions or replace the joint-minor problem by a socle
model.

## Verification and consequence for the open problem

Run

~~~bash
python3 -B scripts/run_initial_fan_padding_checks.py
~~~

The report is results/research/initial_fan_padding_checks.json.
Small cases use independent quotient columns and commuting multiplication,
exhaustive basis profiles, full border cones, and translated full-fan
comparisons. General overlapping relation spaces and characteristic
cancellations, propagating relations, nonlocal quotients and inseparable
univariate factors are included. The latter compare original monomial
and polynomial-prefix scans, rather than identifying arbitrary-subset
matroids. Large-parameter records verify exact
lengths, divisor volumes, factorized crossing witnesses and primitive
cost inequalities; they do not claim to enumerate those enormous
quotient fans.
The lifted Gröbner bases are also checked by ordinary polynomial
S-pair division in every tested original cone.

The state-polyhedron dictionary is established in
[Babson--Onn--Thomas, Theorem 2.4](https://arxiv.org/pdf/math/0207135).
Socle moves are established in
[Altmann--Sturmfels, Proposition 1](https://www.math.fu-berlin.de/altmann/PAPER/graph.pdf);
that passage was inspected. The padding proof is elementary direct-sum
linear algebra and downset extrusion. Exact priority of this packaging
and the explicit padded families remains unassessed.
The inclusion-and-colength verification of (F3) is an established
leading-ideal method. Related homogeneous biliaison criteria appear in
[Migliore--Nagel, Definition 8.1 and Lemma 8.4](https://academicweb.nd.edu/~jmiglior/MN14-webpg.pdf);
those passages were inspected. Their regularity and homogeneous
codimension hypotheses are not assumed for this finite padding
construction, whose length and leading-ideal identities are proved
directly. This comparison does not establish exact priority.

Absolute divisor volume cannot control weighted reuse. A viable
collective argument must retain relative support geometry or another
constraint that survives removal of compulsory padding. This result
does not refute a budget using primitive costs and the full compatible
rank data. The independent-socle route remains closed as a candidate
for supergeneric growth; the entire fan here is only \(4s\).
All original unrestricted, exact-priority and honest-paper requirements
remain unchanged.
