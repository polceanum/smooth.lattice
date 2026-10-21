# The exact colength threshold for a primitive wall

Let $q\in\mathbb Z^d$ be primitive and have both signs. Write

$$
 u=q^+,\qquad v=q^-,\qquad
 A(q)=\prod_i(u_i+1),\quad B(q)=\prod_i(v_i+1),\quad
 \nu(q)=A(q)+B(q)-2.                                      \tag{1}
$$

A wall means a codimension-one wall between two maximal positive
initial-ideal cones. Its direction is the primitive normal $q$, up to
sign, or equivalently the primitive direction of the dual state edge.
It is not a primitive normal of a state facet.

**Theorem 1.** Over every field, the minimum colength of an ideal having
such a wall of direction $q$ is exactly $\nu(q)$. For $d\ge3$,
every length $n\ge\nu(q)$ is possible. At the minimum length the wall
degeneration is necessarily

$$
 J_D+(x^u-\lambda x^v),\qquad \lambda\in K^*,\qquad
 D=\operatorname{Div}(x^u)\cup\operatorname{Div}(x^v),       \tag{2}
$$

where $J_D$ is the monomial ideal complementary to $D$.

This is an exact local realization and rigidity statement. It does not
bound the number of walls that coexist in one quotient. In particular,
it does not improve the current unrestricted $G_3(n)$ bound or finish
the active extremal goal. Exact priority of the quantitative threshold,
rigidity and counting consequences below remains unassessed.

## Lower bound from the wall grading

At a positive weight in the relative interior of a wall, the initial
ideal is homogeneous for $\mathbb Z^d/\mathbb Zq$. This is the
classical edge-ideal degeneration. Its two monomial initial ideals
have the same Hilbert function for this grading and length $n$.

Here is a justification that does not assume the original ideal is
binomial. The common wall has dimension $d-1$, and its linear span is
$q^\perp$. On its relative interior the weight-initial ideal is
constant. Initial forms, taken from one adjacent marked reduced basis,
can contain only comparisons vanishing on that whole wall. Their
exponent differences are integer multiples of the primitive $q$.
Consequently they generate a $\mathbb Z^d/\mathbb Zq$-homogeneous
ideal. Positive-weight initial degeneration preserves colength.

The two standard sets differ. In some grading fiber their selected
monomial subsets therefore differ. Two different exponents in that
fiber differ by $kq$, $k\ne0$. One is divisible by $x^u$ and the
other by $x^v$, possibly after reversing their roles. Since the
standard sets are downsets, $x^u$ and $x^v$ are nonzero in the wall
quotient. Every monomial in $D$ is nonzero as well.

All monomials of $D$ have different grading degrees except $x^u,x^v$.
Indeed the two divisor boxes meet only at 1. Within either box a
nonzero multiple of the mixed-sign vector $q$ cannot be a difference.
Between the two boxes, the coordinate bounds permit a nonzero multiple
only for the pair $u,v$. Thus $D$ represents

$$
             |D|-1=A(q)+B(q)-2                           \tag{3}
$$

different degrees. Nonzero homogeneous vectors in different degrees
are independent, proving $n\ge\nu(q)$.

## Attainment, rigidity and padding

Both $x^u,x^v$ are maximal in $D$, so every variable kills their
classes in $K[x]/J_D$. Their difference is therefore a socle
relation. The quotient (2) has length $|D|-1=\nu(q)$ over every
field. Its only nontrivial monomial residue class is the pair
$\{x^u,x^v\}$; all other monomials of $D$ are independent, and all
monomials outside $D$ vanish. Hence its two standard sets are

$$
                       D\setminus\{u\},\qquad
                       D\setminus\{v\}.                  \tag{4}
$$

The complete positive wall is $q\cdot w=0$. It meets the positive
orthant because $q$ has both signs. Both sides occur.

For rigidity, every grading fiber represented by a proper divisor in
$D\setminus\{u,v\}$ is a singleton even in the full monomial ground.
Adding or subtracting $q$ violates nonnegativity: the proper divisor
is deficient in at least one required coordinate. The common degree
of $u,v$ contains precisely these two monomials. At equality in (3)
these are all nonzero degree spaces, each of dimension one. Every
outside monomial vanishes, and the unique relation in the last space
has two nonzero coefficients. This proves (2).

For $d\ge3$, choose a coordinate axis which is not the support of
either endpoint if that endpoint is a pure power. At most two axes
are forbidden. Let $l$ be the largest exponent on the chosen axis
already in $D$. Add exactly $n-\nu(q)$ successive points on that
axis beyond $l$. Neither endpoint becomes nonmaximal: a pure power
on this axis cannot dominate a monomial with at least two nonzero
coordinates. The enlarged set is a downset, and the same single socle
relation gives length $n$ and the same wall. Padding changes fixed
monomial columns only. It does not claim that arbitrary padding in
two variables has this property.

## The exact universal direction spectrum

For $d\ge3$, let $\mathcal D_d(n)$ be the set of primitive positive
wall directions, up to sign, occurring among length-$n$ ideals over
a fixed field. The field can be arbitrary. Theorem 1 gives the exact,
field-independent description

$$
 \boxed{\mathcal D_d(n)=\{q\text{ primitive and mixed}:
                         A(q)+B(q)\le n+2\}/\{\pm1\}.}    \tag{5}
$$

In every $d\ge2$ this also describes directions occurring at length
*at most* $n$. This sharpens merely taking all differences of the
universal hyperbolic ground: the two divisor budgets must fit together
in one quotient, after the common constant and the endpoint relation
are removed.

There is an exact common-refinement consequence. In the positive
orthant, the hyperplane arrangement

$$
                     \{q^\perp:q\in\mathcal D_d(n)\}       \tag{6}
$$

is the common refinement of **all** length-$n$ initial-ideal fans,
over any fixed field, for $d\ge3$. Every individual wall belongs to
this arrangement by necessity in (5). Conversely the construction
and padding realize the entire positive part of each hyperplane as
the sole wall of an actual two-initial ideal. Thus every hyperplane
in (6) is required in that common refinement. This is a statement
about common fan geometry, not a lower bound on the number of
weights in an arbitrary sampling list.

## Sharp size in three variables

Put $D_3(n)=|\mathcal D_3(n)|$. Then

$$
             \boxed{D_3(n)\sim
                  \frac{3}{2\zeta(3)}n^2\log n.}           \tag{7}
$$

Directions with two nonzero coordinates contribute exactly
$3\sum_{j=2}^n\varphi(j)=O(n^2)$: on a chosen coordinate pair
their positive magnitudes satisfy $a+b\le n$ and $\gcd(a,b)=1$.

For full support choose the singleton sign in three ways and write
$q=(a,-b,-c)$, with $a,b,c\ge1$. Condition (5) becomes

$$
                      a+bc+b+c\le n.                     \tag{8}
$$

Let $T_h(n)$ count triples satisfying (8) with $h\mid a,b,c$.
Möbius inversion gives the primitive triple count
$F(n)=\sum_{h\ge1}\mu(h)T_h(n)$. With $M=\lfloor n/h\rfloor$,
an exact integer formula is

$$
 T_h(n)=\sum_{B=1}^{\lfloor(M-2)/(h+1)\rfloor}
 \left[C_B(M-B)-(hB+1)\frac{C_B(C_B+1)}2\right],\quad
 C_B=\max\left(0,\left\lfloor\frac{M-B-1}{hB+1}\right\rfloor\right).
                                                               \tag{9}
$$

Indeed the scaled inequality is $A+hBC+B+C\le M$, and the bracket
sums the possible positive $A$ over positive $C\le C_B$.

For completeness the leading asymptotic uses only elementary divisor
sums. Define $W(X)=\sum_{BC\le X}(X-BC)$. The usual harmonic-sum
estimate $\#\{BC\le X\}=X\log X+O(X)$, integrated in $X$, gives

$$
                       W(X)=\tfrac12X^2\log X+O(X^2).
$$

For fixed $h$, the extra $B+C$ in (9) does not change the leading
coefficient. Restrict first to $B,C\ge K$, where
$B+C\le2BC/K$. The counts are squeezed between the weighted divisor
sums with product coefficients $h$ and $h+2/K$, up to
$O_{h,K}(n\log(n+1))$ integer-rounding terms. The discarded
rows $B<K$ or $C<K$ contribute $O_{h,K}(n^2)$. Let $n\to\infty$
and then $K\to\infty$ to obtain

$$
       \frac{T_h(n)}{n^2\log n}\longrightarrow\frac1{2h^3}.
$$

Uniform domination follows by dropping $B+C$:
$T_h(n)\le hW(n/h^2)\le Cn^2\log(n+1)/h^3$ when nonzero.
The convergent sum of $h^{-3}$ permits dominated convergence in
Möbius inversion. The Euler product gives
$\sum_h\mu(h)/h^3=1/\zeta(3)$. Therefore
$F(n)\sim n^2\log n/(2\zeta(3))$, proving (7).

## The fixed-dimensional law

The same spectrum has a sharp count in every fixed dimension:

$$
 |\mathcal D_d(n)|\sim
 \frac{\binom{2d-2}{d-1}}{4(d-2)!\zeta(d)}
                 n^2(\log n)^{d-2},\qquad d\ge3.          \tag{11}
$$

For $d=2$ the formula holds for the spectrum at length at most $n$.
The leading power of $n$ is two independently of the fixed dimension.
This counts directions across the whole class, not cells in any one fan.

Here are the details of the extension. On a fixed support of size
$m=s+t$, choose $s$ positive and $t$ negative coordinates, with both
nonempty. For a fixed common divisor $h$, the two divisor volumes are
products of $s$ and $t$ integers from $h+1,2h+1,\ldots$. If
$P_{k,h}(X)$ counts products of $k$ such factors at most $X$, then

$$
 P_{k,h}(X)\sim
      \frac{X(\log X)^{k-1}}{h^k(k-1)!}.                  \tag{12}
$$

For fixed $h$ this follows by induction from the harmonic sum on that
arithmetic progression. The $k=1$ case is $X/h+O_h(1)$; convolution
with another factor integrates
$(\log(X/u))^{k-2}/u$ and supplies $1/(k-1)$.
The lower bound on each factor only affects lower logarithmic terms.

The number of pairs of such products with sum at most $N=n+2$ is
the Stieltjes convolution of the two counting functions. After
dividing by $N^2(\log N)^{m-2}$, (12) gives the limiting integral
$\int_0^1(1-x)\,dx=1/2$. One can first restrict to
$\epsilon N\le A\le(1-\epsilon)N$. The two omitted ends contribute
$O_{s,t,h}(\epsilon N^2(\log N)^{m-2})$ by the product counting
bounds; then let $\epsilon$ tend to zero. Thus the leading coefficient
for fixed $h$ is
$1/(2h^m(s-1)!(t-1)!)$.

For domination uniformly in $h$, replacing $ha+1$ by $ha$ gives
$P_{k,h}(N)\le Q_k(N/h^k)$, where $Q_k$ counts ordinary positive
integer factors. The standard upper divisor estimate therefore bounds
the pair count by $C_{s,t}N^2(1+\log N)^{m-2}/h^m$.
Möbius inversion and dominated convergence multiply the leading
coefficient by $1/\zeta(m)$. Smaller supports contribute
$O_d(n^2(1+\log n)^{d-3})$ when $d\ge3$.
For full support, choose its positive part in $\binom ds$ ways and
divide by two for the unoriented sign. The leading constant is

$$
 \frac1{4\zeta(d)}
 \sum_{s=1}^{d-1}\frac{\binom ds}{(s-1)!(d-s-1)!}
 =\frac{\binom{2d-2}{d-1}}{4(d-2)!\zeta(d)},
$$

by Vandermonde's identity. This proves (11), using the classical
divisor and Möbius tools explicitly rather than extrapolating finite counts.

## Consequence for the current incidence route

Even knowing the exact possible directions does not allow their
individual shadow envelopes to be summed efficiently. In three
variables,

$$
 \sum_{q\in\mathcal D_3(n)}(1+\|q\|_1^{1/3})
                         =\Theta(n^{7/3}\log(n+1)).       \tag{10}
$$

The upper bound follows from (5), $\|q\|_1\le3n$, and (7).
For the lower bound take $n/4\le a\le n/2$ and
$(b+1)(c+1)\le n/4$, with $\gcd(b,c)=1$. These give
$\Theta(n^2\log n)$ primitive triples and have norm at least $n/4$.
The coprime hyperbola count has order $X\log X$, by the same
Möbius/divisor calculation with the convergent $h^{-2}$ sum.

Consequently separately summing the repository's envelope
$O(n(1+\|q\|_1^{1/3})\log(n+1))$ over the universal spectrum pays a
budget of order $n^{10/3}\log^2(n+1)$. This exceeds the existing
unrestricted bound. It is a limitation of that *universal envelope
sum*, not an actual lower bound for shadow sizes, one ideal's fan,
or an enumeration algorithm. A sharper proof must retain which walls
coexist in a quotient, their incidence or their shared rank costs.
The exact local threshold alone cannot supply the missing collective
square-root-incidence budget.

## Sources and verification

[Altmann--Sturmfels, *The graph of monomial ideals*](https://www.math.fu-berlin.de/altmann/PAPER/graph.pdf),
Propositions 1 and 3, already supplies the socle move and the
codimension-one homogeneous edge-ideal framework. Lemma 13 and
Algorithm 14 identify the edge direction from endpoint profiles.
Their statements and displayed proofs were inspected. The construction
in (2) is an instance of that published socle move. The universal
state/greedy dictionary and hyperplane/zonotope refinement method are
established in [Babson--Onn--Thomas](https://arxiv.org/pdf/math/0207135).
Divisor asymptotics and Möbius inversion are classical inputs. These
comparisons do not establish novelty of these statements or the strongest
current unrestricted bound.

Run from the repository root:

```bash
python3 -B scripts/run_primitive_wall_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker verifies the integer formula against separate primitive
triple enumeration, constructs minimal and padded commuting cyclic
quotients, replays both standard sets by independent field elimination,
tests the threshold on all directed adjacencies in the saved exact
fan artifact, and checks the threshold and equality rigidity on all
pairs of three-variable monomial downsets through length eight having
the same relevant grading Hilbert function. Finite controls support
the uniform proofs. Separate congruence-product convolutions and
bounded exhaustive direction cubes check exact counts in dimensions
two through five, without fitting an asymptotic constant. They do not
enumerate all ideals or certify priority. The report records exact
commands and relevant source hashes. No core kernel, certification
range or prime-random-access claim changes.
