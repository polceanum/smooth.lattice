# Capped rank filtering for planar fat-point algebras

The generic planar first-jet construction with a square-zero third
coordinate cannot give a polynomial separation above the generic
three-variable scale. More generally, fix a jet order \(r\ge2\). For
generic distinct planar points, adjoin **any** element of the radical as
the third coordinate. The resulting length-\(n\) quotient satisfies

\[
 g(I)=O\!\left(r^2n^2\log(n+1)(1+\log r)\right).             \tag{J1}
\]

In particular the entire fixed-\(r\) class has
\(g(I)=O_r(n^2\log(n+1))\), uniformly in the radical element. This is
a class exclusion, not an improved unrestricted \(G_3\) bound, an
asymptotic lower separation, a sharpness claim for this class, or a
novelty certificate. Special point configurations and arbitrary ideals
containing \(x^r\) are outside the conclusion.

## Filtering inside an actual monic carrier box

Let \(A=K[x_1,\ldots,x_d]/I\) have length \(n\). Suppose \(I\) contains
a monic univariate polynomial of degree \(H_i\) in each \(x_i\), and put

\[
 D=\prod_i\{0,\ldots,H_i-1\}.
\]

These relations have heads \(x_i^{H_i}\) for every positive weight.
Consequently every standard set lies in \(D\). Greedy selection on the
actual residue columns in \(D\) gives the full standard set; this
restriction is justified by the carrier relations, not by replacing
the residue matroid with a uniform one.

Let \(s_D\) be the largest cardinality of an anchored subbox of \(D\)
whose columns do not span \(A\), taking \(s_D=0\) if none exists.
Let \(M_D\) be the supremum, over generic positive weights, of the
number of terms of \(D\) scanned through the last greedy selection.
Then the earlier unrestricted-box argument has the capped version

\[
 1+s_D/d^d\ \le M_D\ \le1+(d+1)^d s_D.                    \tag{J2}
\]

For the upper bound let the last selected monomial have cost \(\tau>0\).
The subbox with coordinate counts

\[
 b_i=\min\{H_i,1+\lfloor\tau/((d+1)w_i)\rfloor\}
\]

is strictly cheaper than \(\tau\), so does not span. Every term scanned
has coordinate count at most
\(\min\{H_i,1+\lfloor\tau/w_i\rfloor\}\le(d+1)b_i\).
This proves the upper bound, with harmless extra one as stated.
For the lower bound take a nonsaturating subbox with caps \(a_i\) and
prices \(w_i=1/(a_i+1)\). Terms of cost below one lie in this box.
Its subbox with caps \(\lfloor a_i/d\rfloor\) has all costs below one
and at least \(\prod_i(a_i+1)/d^d\) terms. A small generic perturbation
preserves separation from every term of \(D\) outside the original
box. Completion must follow these terms. The case \(n=1\) has
\(H_i=1\) if minimal univariate degrees are used, hence \(s_D=0\).
For nonminimal carriers at \(n=1\), greedy still selects \(1\) first
and every anchored box spans, so the same statement holds.

This extension permits noncyclic coordinates: \(D\) is finite even
when unrestricted scanning would run arbitrarily far along an axis.
Its quality still depends on actual ranks.

**Ordered relative-prefix consequence.** In three variables suppose
greedy finishes within \(M\) terms of \(D\), and \(H_x=h\).
If \(|D|\ge M\), take the first \(M\) terms. Let \(\ell_i\) be their
axis counts. The cut contains the subbox with counts
\(\lceil\ell_i/3\rceil\) and is contained in the box with counts
\(\ell_i\). Thus

\[
 M\le\prod_i\ell_i\le27M.
\]

Dyadic enclosing boxes have product \(O(M)\), and there are
\(O(\log(M+1)(1+\log\min(h,M)))\) possible charts. Intersect each
chart with the original caps. Its \(O(M)\) difference vectors give
\(O(M^2)\) central-arrangement chambers; on each chamber the entire
ordered prefix is fixed. Applying actual-column greedy therefore gives

\[
 g(I)=O\!\left(M^2\log(M+1)
                   (1+\log\min(h,M))\right).              \tag{J3}
\]

If \(|D|<M\), order the whole box, which has \(O(M)\) differences;
the same bound follows. Ties may be avoided by a perturbation within
each initial cone. This is the existing
[relative-box ordered-prefix argument](nested_jet_matroids.md),
including its caps. It is a direct classical geometric consequence,
not a new lattice theorem.

The construction also gives one common finite weight list depending
only on the caps and \(M\), with the cardinality in (J3). Choose a
positive rational point in every chart chamber and clear its
denominator. Every quotient obeying the same completion bound is
sampled by this list. At each weight, capped heap enumeration needs
only its first \(M\) monomials, followed by selection using the actual
residue columns. The list is independent of the coefficients; field
operations and coefficient bit costs are separate.

Combining (J2)--(J3), a fixed cap \(h\), with the other two caps at
most \(n\), gives \(g(I)=O_h((s_D+1)^2\log(n+1))\).
Thus \(g(I)\ge n^{2+\varepsilon}\) requires
\(s_D=\Omega_h(n^{1+\varepsilon/2}/\sqrt{\log(n+1)})\).
A noncyclic coordinate by itself does not avoid this necessary
rank-defect test.

## A uniform rectangular interpolation condition

Fix distinct affine points \(P_j=(a_j,b_j)\), \(1\le j\le m\), and
let

\[
 B=K[y,z]/\bigcap_{j=1}^m(y-a_j,z-b_j)^r,\qquad
 n=m\binom{r+1}{2}.                                      \tag{J4}
\]

Use coefficients in the truncated Taylor expansions at the points;
ordinary derivatives divided by factorials are not needed.
For every infinite field \(K\) there is a nonempty Zariski-open set
of point configurations such that, simultaneously for
\(1\le A,C\le rm\),

\[
 A,C\ge r,\quad
 \lfloor A/r\rfloor\lfloor C/r\rfloor\ge m
 \quad\Longrightarrow\quad
 \{y^az^b:0\le a<A,\ 0\le b<C\}\text{ spans }B.           \tag{J5}
\]

Here is an elementary proof sufficient for this coarse bound. Put
\(q=\lfloor A/r\rfloor\), \(s=\lfloor C/r\rfloor\), and specialize
the points to any \(m\) points of a \(q\) by \(s\) Cartesian grid.
The univariate Chinese remainder maps for powers
\((y-a)^r\) and \((z-b)^r\), and their tensor product, identify the
rectangle of degrees \(<rq,<rs\) with the product of the local
rectangular jet rings
\(K[u,v]/(u^r,v^r)\). Project onto the chosen \(m\) factors and then
onto \(K[u,v]/(u,v)^r\). This proves surjectivity at the
specialization, in every characteristic. A maximal minor is therefore
a nonzero polynomial in the point coordinates. There are finitely many
rectangles in (J5); the product of one nonzero minor for each, together
with distinctness conditions, remains nonzero. Its nonvanishing locus
proves the assertion over an infinite field.

The grid witnesses for different rectangles need not be the same
configuration. The finite product of nonzero determinant polynomials
is the justification for simultaneous genericity. Over a finite field,
use (J5) as an explicit rank hypothesis; a nonempty open set need not
have a rational point over that particular field.

## Arbitrary radical third coordinates

Take any \(x\in\operatorname{rad}B\), allowing arbitrary coefficients
and higher Taylor terms separately at the points, and let \(I\) be the
kernel of \(K[X,Y,Z]\to B\), with images \((x,y,z)\).
The map is surjective because \(y,z\) already generate \(B\).
There are actual all-weight carrier relations

\[
 X^r,\qquad
 \prod_j(Y-a_j)^r,\qquad
 \prod_j(Z-b_j)^r\ \in I.
\]

Use \(D=[0,r)\times[0,rm)^2\). Consider a strict free prefix \(S\)
in this box and its \(X=0\) part \(S_0\). Downward closure gives
\(|S_0|\ge|S|/r\). The two-dimensional cut \(S_0\) has axis counts
\(\ell_y,\ell_z\); its half-size rectangle with counts
\(A=\lceil\ell_y/2\rceil\), \(C=\lceil\ell_z/2\rceil\)
lies in \(S_0\), and

\[
 AC\ge |S_0|/4.
\]

If \(|S|>16r^3m\), then \(AC>4r^2m\). Since \(A,C\le rm\),
neither can be below \(r\). Also
\(\lfloor A/r\rfloor\ge A/(2r)\) and similarly for \(C\), so (J5)
shows that this rectangle spans \(B\). Therefore

\[
 M_D\le \min\{|D|,16r^3m+1\}\le32rn+1.                   \tag{J6}
\]

Apply (J3) with \(h=r\). Since \(r\le n\), this gives (J1).
No rank property of \(x\) was used beyond \(x^r=0\).

For first jets, \(r=2\), a slightly smaller constant is available:
\(M_D\le\min\{|D|,72m+1\}\).
Indeed \(|S|>72m\) gives \(AC>9m\).
The caps \(A,C\le2m\) force \(A,C\ge2\), and
\(\lfloor A/2\rfloor\ge A/3\) proves (J5).
The asymptotic count remains \(O(n^2\log(n+1))\).

This closes generic planar fat points with an arbitrary radical
coordinate as a power-separation route, including the square-zero
first-jet probes. It does not close all square-zero-coordinate ideals:
the quotient by \(X\) and the two surviving actions can have ranks
that do not arise from (J4)--(J5).

## Necessary limits and prior art

Rectangular interpolation does not mean that every \(n\)-column set is
independent. At two distinct first-jet points, the six ordinary
quadratic monomials have rank at most five: the square of the line
through the two points vanishes to order two at both.

Generic position in (J5) is also substantive. At nine collinear
first-jet points the \(6\) by \(6\) rectangle meets the numerical
condition in (J5), but has rank \(20<27\). The sixteen independent
multiples of the squared line inside that rectangle already bound
its rank by \(36-16=20\). No assertion for arbitrary configurations
is inferred from generic rank.

Even a large capped defect is insufficient for a lower bound. For
every \(m\ge2\), take \(A=K[t]/(t^{2m})\) and
\((x,y,z)=(t^m,t,t)\). With caps \(D=[0,2)\times[0,2m)^2\),
one has \(s_D=m^2\), but exactly four initial ideals. The cheaper of
\(y,z\) eliminates the other. If its price is \(w\), comparison with
\(w_x=mw\) chooses either its pure chain of length \(2m\), or two
chains of length \(m\), one multiplied by \(x\).

To check \(s_D\), a box with \(x\)-count one and other counts \(B,C\)
has rank \(\min(2m,B+C-1)\); its largest nonsaturating size is
\(m^2\). With \(x\)-count two it spans exactly when \(B+C\ge m+1\);
otherwise its size is at most \(2\lfloor m^2/4\rfloor\).
Taking \(w_y=W,w_z=W+1\), with \(W>(2m)^2\), and making \(x\)
more expensive than all these terms, scans
\((2m)(2m-1)/2+1\) terms before completing. Hence \(M_D=\Theta(m^2)\)
despite four initials. Capped rank defects are necessary tests, not
fan lower bounds.

Bigraded fat-point interpolation is established prior art.
[Carlini--Catalisano--Oneto, arXiv:1711.06193v2](https://arxiv.org/pdf/1711.06193)
formulate this same rectangular Hilbert-function problem in
Definition 1.1; their Theorem 3.10 treats specified low bidegrees and
Theorem 4.6 gives the triple-point function over \(\mathbb C\).
Their introduction states the earlier double-point solution.
Lemma 2.3 and Remark 2.4 give the specialization principle used here.
The definitions, introductory theorem statements, and this short
specialization proof were inspected; the full numerical theorems'
proofs were not audited here. The coarse tensor-grid proof above
supplies (J5) independently in arbitrary characteristic. It does not
claim a new interpolation theorem or solve the intermediate
bigraded Hilbert-function problem for higher multiplicity.

Chinese remaindering, tensor Hermite interpolation, nonzero-minor
genericity and the dyadic prefix count are classical inputs.
Priority of their particular fan-count consequence is unassessed.
The unrestricted extremal goal and the collective middle wall-cost
problem remain open.

## Reproducible exact controls

Run:

~~~bash
python3 -B scripts/run_capped_jet_checks.py
~~~

The report is results/research/capped_jet_checks.json. It records
actual commuting multiplication maps and cyclic monomial replays,
complete finite quotient fans with exact rational cone areas, a
separate exhaustive downset comparison at small length, tensor-grid
rank witnesses in several characteristics, the required rectangular
ranks for the specified point fixtures, capped-prefix greedy replays,
and both negative controls above. These are finite proof controls;
no growth rate is fitted to their counts.
