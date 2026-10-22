# Sharp fan changes after one point is inserted on four parallel lines

A single rational point can increase the number of monomial initial
ideals by quadratic order, already on four fixed parallel lines.
There is an exact square-grid versus diagonal-grid description of the
two whole positive fans. It also gives a cubic increase in primitive
wall length without an arithmetic subsequence.

Let \(s\ge1\), \(D=2s-1\), and define integer point sets over
\(\mathbb Q\) by
\[
 S_s=\{(i,j,k):i,j\in\{0,1\},\ 0\le k\le D,
                                      \ i+j+k\text{ odd}\},
 \qquad P_s=S_s\cup\{(0,0,0)\}.
                                                        \tag{1}
\]
Write \(J_s=I(S_s)\), \(I_s=I(P_s)\). These are radical ideals
of lengths \(4s\) and \(4s+1\), with \(I_s\subset J_s\).
All points lie on the fixed curve
\(x(x-1)=y(y-1)=0\). Then
\[
 \begin{split}
 g(J_s)&=(s+1)^2,& g(I_s)&=s(2s+1),\\
 W(J_s)&=(2s+1)s(s+2),&
 W(I_s)&=\frac83s^3+6s^2-\frac23s-2.                    \tag{2}
 \end{split}
\]
Here \(W\) counts every actual undirected wall piece once, weighted
by the one-norm of its primitive mixed integer direction. Consequently
\[
 \begin{split}
 g(I_s)-g(J_s)&=s^2-s-1=\Theta(s^2),\\
 W(I_s)-W(J_s)&=\frac{2s^3+3s^2-8s-6}{3}=\Theta(s^3).
                                                               \tag{3}
 \end{split}
\]
The orders in (3) concern \(s\to\infty\); at \(s=1\) both
differences are negative. Exact priority of this construction and
the update formulation remains unassessed. The fan counts are below
the generic \(n^2\log^2 n\) scale, so they give no unrestricted
polynomial lower separation.

## A classical polynomial basis makes the rank rule explicit

For \(0\le c\le D\), use normalized binary Krawtchouk polynomials
\[
 R_c(k)=\binom Dc^{-1}\sum_{h=0}^c(-1)^h
                         \binom kh\binom{D-k}{c-h}.
                                                               \tag{4}
\]
They have degree \(c\), nonzero leading coefficient, \(R_c(0)=1\),
and, at the integer nodes \(0,\ldots,D\),
\[
                    R_{D-c}(k)=(-1)^kR_c(k).             \tag{5}
\]
Indeed (5) follows by reversing the coefficients of
\((1+t)^{D-k}(1-t)^k\). Their orthogonality follows by multiplying
two such generating functions and summing with binomial weights:
\[
 \sum_{k=0}^D\binom Dk R_c(k)R_d(k)
                  =\frac{2^D}{\binom Dc}\,\mathbf1_{c=d}.
                                                               \tag{6}
\]
These are classical identities, not new orthogonal-polynomial results.

The coordinate annihilators put every standard set in
\[
                   E=\{0,1\}^2\times\{0,\ldots,D\}.
\]
Replace the monomial with label \(a=(a_1,a_2,a_3)\in E\) by
\[
 Q_a(x,y,z)=(1-2x)^{a_1}(1-2y)^{a_2}R_{a_3}(z).
\]
This replacement has a nonzero leading coefficient and only proper
divisors as other monomials. It is therefore triangular for every
global monomial order and preserves the greedy selected labels.
This is a statement about selected bases, not equality of all subset
ranks in the monomial and transformed columns.

Put \(a^*=(1-a_1,1-a_2,D-a_3)\). On \(S_s\), (5) gives
\[
                              Q_{a^*}=-Q_a.             \tag{7}
\]
There are \(N=4s\) pairs. One column from each pair is independent.
To verify this, take the representatives with \(a_1=0\). On the full
grid, tensor orthogonality (6) makes the \(Q_a\) orthogonal for the
positive weights \(\binom Dk\). Restriction to odd parity replaces
the full inner product by half its difference with the parity-twisted
inner product. The latter pairs \(a\) with \(a^*\), which cannot
both be representatives. Thus the restricted Gram matrix is diagonal
with positive entries. This proves independence without a genericity
assumption on the integer points.

In these reference coordinates the columns of \(J_s\) are the pairs
\(e_r,-e_r\), \(1\le r\le N\). At the added point every \(Q_a\)
equals one, so the columns of \(I_s\) are
\[
                            (e_r,1),\quad(-e_r,1).       \tag{8}
\]
A basis of (8) chooses one member of every pair and both members of
exactly one pair. Every such choice is a basis: the double pair
supplies the last coordinate, and the other columns supply the
remaining independent coordinates. The kernel of the actual quotient
to \(J_s\) is the separator line of \((0,0,0)\). Multiplication by
\(x,y,z\) kills this line. No abstract matroid is substituted for a
point algebra in this argument.

## The two complete fans

Normalize a positive weight as
\[
 U=w_x/w_z>0,\quad V=w_y/w_z>0,\quad
                   P=U+V,\quad Q=U-V,
\]
so \(P>|Q|\). Pair endpoints have common center
\(C=(1/2,1/2,D/2)\). Their differences, with first coordinate one,
are
\[
               (1,1,k),\quad(1,-1,k),
                    \qquad k=-D,-D+2,\ldots,D.          \tag{9}
\]
For \(J_s\), select the cheaper member of every pair. Its fan is
exactly the arrangement
\[
 P=a\quad(a=1,3,\ldots,D),\qquad
 Q=b\quad(b=-D,-D+2,\ldots,D)                            \tag{10}
\]
in \(P>|Q|\). No mixed monomial outside \(E\) can supply another
standard set, by the univariate annihilators. All arrangement chambers
have different selected labels.

For \(I_s\), first select the cheaper member of every pair and then
the cheapest unselected member. Thus double the pair for which
\(|P-a|\) or \(|Q-b|\) is smallest. Equivalently its minimum support
function is that of \(J_s\), plus
\(w\cdot C+\tfrac12\min_r|w\cdot g_r|\), with \(g_r\) from (9).
This proves the following full list of its open cones:

* The unit squares \(i<U<i+1\), \(j<V<j+1\), where
  \(i,j\ge0\) and \(i+j\le D-1\).
* For \(i=0,\ldots,D\), the unbounded region
  \(U>i\), \(V>D-i\), with
  \(2i-D-1<Q<2i-D+1\). Omit the lower \(Q\) inequality at \(i=0\)
  and the upper one at \(i=D\).

For a bounded square, the doubled pair has \(P\)-level
\(i+j+1\) when \(i+j\) is even, and \(Q\)-level \(i-j\)
when it is odd. For the unbounded region it has \(Q\)-level
\(2i-D\). These follow by comparing distances to consecutive odd
integers. All other pair choices are fixed on each listed region.
The regions cover the positive quadrant apart from their walls, with
distinct bases. Therefore they describe whole cones, rather than a
selected sample of objectives.

There are \(D(D+1)/2\) bounded squares and \(D+1\) unbounded cones,
giving \(g(I_s)=s(2s+1)\). In (10), the \(2s\) parallel
\(Q\)-lines give \(2s+1\) initial regions. The line \(P=a\)
is split into \(a\) pieces and adds \(a\) regions. Summing the
positive odd \(a\)'s gives \(g(J_s)=(s+1)^2\).

## Exact primitive wall sums

For \(I_s\), the walls are
\[
 \begin{array}{c|c|c}
 \text{line}&\text{number of pieces}&\text{primitive one-norm}\\\hline
 U=i,\ V=i\quad(1\le i\le D)&D+1-i\text{ on each}&i+1\\
 Q=b\quad(b=-D+1,-D+3,\ldots,D-1)&1&|b|+2.
 \end{array}
\]
The first segments end at \(U+V=D+1\); the second rays start
there. Hence
\[
 W(I_s)=2\sum_{i=1}^D(D+1-i)(i+1)
       +\sum_{\substack{|b|\le D-1\\b\text{ even}}}(|b|+2).
                                                               \tag{11}
\]
For \(J_s\), a positive odd \(P=a\) has \(a\) pieces. A
\(Q=\pm a\) has \((D-a)/2+1\) pieces. Both directions have
primitive one-norm \(a+2\). Thus
\[
 W(J_s)=\sum_{\substack{1\le a\le D\\a\text{ odd}}}
             (D+2)(a+2).                                 \tag{12}
\]
Evaluating (11)--(12) proves (2)--(3). Every primitive direction
contains a coordinate equal to one in magnitude; the jump does not
come from changing gcd factors of long normals. Product examples
in the [ideal-step note](ideal_step_incidence.md) also give cubic
logarithmic jumps and quartic chain variation. The present calculation
instead isolates the exact geometry and the quadratic cone-count change.

## A sharp update order in the two-quadratic class

More generally, suppose a length-\(n\) ideal over any field contains
monic scalar quadratics in \(x\) and \(y\). Every standard set has
\(x,y\)-exponents at most one. Every actual wall direction is parallel
to a difference of two labels of an optimal basis exchange. Its
primitive first two coordinates therefore belong to \(\{-1,0,1\}\),
and its last has magnitude at most \(n\). A pure last-coordinate
direction cannot give a positive wall. Up to sign there are
\(O(n)\) possible directions, of one-norm \(O(n)\).

Their lines in the positive weight simplex have \(O(n^2)\) chambers
and \(O(n^2)\) wall pieces. Consequently, uniformly over this class,
\[
                         g(I)=O(n^2),\qquad W(I)=O(n^3).
\]
Both powers and both one-point update orders are sharp by (1)--(3).
In particular an algorithm maintaining an explicit list of maximal
initial cones may need \(\Omega(n^2)\) additional cone records after
one rational point is inserted. This is an output-size lower bound
for the entire fan, not for one Gröbner basis, compressed representations
or smooth-number random access. A reduced basis in this class has at
most six polynomials, by the four possible \(z\)-tower tips and the
two quadratic heads.

This does not improve unrestricted \(G_3\). In particular a cubic
global wall budget remains compatible with these cubic jumps, and
choosing a special ideal chain is a different question from bounding
every step or every chain.

## Reproduction and priority boundary

~~~bash
python3 -B scripts/run_parity_point_fan_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The report compares actual point evaluation columns with the paired
model, reconstructs both whole fans by exact cone clipping and checks
them against a separate border traversal. It checks Krawtchouk degree,
reflection and weighted Gram identities, the corner separator, ordinary
S-pairs, every primitive wall multiplicity and the exact formulas.
Finite controls support the uniform proof; no fitted growth rate is used.

Normalized Krawtchouk identities and their binomial orthogonality are
classical; see [Harrow--Kolla--Schulman, Section 2, (2.2) and Lemma 2.1](https://www.theoryofcomputing.org/articles/v010a003/v010a003.pdf)
and [DLMF, Tables 18.19.1--2](https://dlmf.nist.gov/18.19).
The stated definitions and identities were inspected. The coefficient
reversal and tensor restriction are proved above. The residue-matroid
and positive initial-fan dictionary is established by
[Babson--Onn--Thomas, Theorem 2.4](https://arxiv.org/pdf/math/0207135).
These comparisons do not establish exact priority of this fixed-coordinate
point insertion example or its update bound. The sharp unrestricted
growth question and an actual polynomial lower separation remain open.
