# Exact outer cones for nonlinear rectangles

The exact Pascal model from the [bilinear note](bilinear_rectangle_fans.md)
now gives a uniform description on a whole open weight cone. It also
disproves a proposed shortcut: even the characteristic-zero bilinear
family has primitive coherent wall directions whose coordinate product
is not $O(n)$.

These are supporting structural results. They do not settle the full
rectangle fan, improve the unrestricted $G_3$ exponent, or complete the
active goal. The rank argument is elementary Vandermonde interpolation;
the exact comparison with earlier work is unfinished. No new general
interpolation theorem or publication priority is asserted.

## A class with an exact prefix rule

Let $r,s\ge2$, let $d\ge1$ be an integer, and work in characteristic
zero. Put

$$
A=K[x,y]/(x^r,y^s),\quad F=y+x(1+y)^d,\quad
I=(x^r,y^s,z-F),\quad n=rs.                                \tag{R1}
$$

Always reduce polynomial columns in $A$; its dimension is $rs$ even
when the expression for $F$ has large degree. The same conclusions
hold for $F=\beta y+\alpha x(1+\gamma y)^d$ with all three
coefficients nonzero, by independent coordinate scalings.

Consider the open positive weight cone

$$
\mathcal C_0=\{\min(w_x,w_z)>(s-1)w_y>0\}.
$$

Normalize $w_y=1$ and write

$$
w_z=T,\qquad w_x=T+\rho,\quad\min(T,T+\rho)>s-1.          \tag{R2}
$$

For each $0\le t<r$, form the following finite set of slots:

$$
\begin{split}
S_t={}&\{(0,j):0\le j<s\}\\
 &\cup\{(i,j):1\le i\le t,\ 0\le j<\min(d,s)\}.
\end{split}                                               \tag{R3}
$$

Choose the $s$ cheapest slots by the additive cost $\rho i+j$.
Map a selected slot $(i,j)$ to the monomial $x^iy^jz^{t-i}$.

**Theorem 1.** If $w\in\mathcal C_0$, $\rho>-d$, and the selected
slots have strict cutoffs, the union of these $r$ selections is exactly
the global initial standard set of $I$. If $\rho<-d$, the standard
set is the ordinary $x,y$ rectangle. In both cases it depends on
$\rho$, not $T$.
This is a statement about ranks on the relevant prefixes, not about
all subsets of all monomial columns.

### Why these are the only candidates

For $\rho>-d$, the largest term of $z-y-x(1+y)^d$ is $xy^d$.
Every monomial divisible by that term reduces to cheaper monomials;
monomials divisible by $x^r$ or $y^s$ are zero. When $d\ge s$,
the $y^s$ relation already removes every $y$ exponent at least $s$;
the larger, unreduced graph relation need not be part of a reduced
Gröbner basis. Keeping the candidates (R3) is still valid. Some
redundant candidates can occur after a block has full rank and do
not affect its first $s$ independent columns.

Thus it suffices to consider $y^jz^t$ and $x^iy^jz^{t-i}$ with
$i>0,j<\min(d,s)$. Their $x$ degree in $A$ is at most $t$.
For $t<r$, the coefficient of $x^t$ in the latter column is

$$
y^j(1+y)^{d(t-i)}\pmod{y^s}.                              \tag{R4}
$$

The pure $y^jz^t$ columns ($i=0$) have the same formula. Their
$s$ leading coefficients form a basis of $K[y]/(y^s)$, since
$(1+y)^{dt}$ is a unit.

### Exact prefix ranks within a block

Divide (R4) by the common unit $(1+y)^{dt}$. For a cost prefix of
slots, the columns become

$$
y^j(1+y)^{-di},\qquad 0\le j\le J_i.                      \tag{R5}
$$

For each $i$, the $j$ values are an initial interval. Since
$y=(1+y)-1$, the columns with a fixed $i$ have the same span as

$$
(1+y)^{-di},(1+y)^{-di+1},\ldots,(1+y)^{-di+J_i}.
$$

For $i>0$, one has $J_i<d$, so these intervals of integer exponents
are disjoint. The $i=0$ interval is nonnegative, and all others are
negative. If there are $L\le s$ distinct exponents $\lambda$,
the first $L$ Taylor rows have determinant

$$
\det\left[\binom{\lambda_b}{a}\right]_{
                 0\le a<L,\ 1\le b\le L}
=\frac{\prod_{b<c}(\lambda_c-\lambda_b)}
       {\prod_{a=0}^{L-1}a!}\ne0.                         \tag{R6}
$$

Here binomial coefficients at negative integers have their usual
polynomial meaning. Consequently every slot prefix has exactly
$\min(s,L)$ leading-component rank. A prefix containing more than
$s$ columns contains an independent $s$-column subprefix, so the
same statement covers all prefixes. Arbitrary subsets need not obey
this rank rule.

### The blocks do not interfere

Put $m=\min(T,T+\rho)>s-1$. Every monomial of total $x,z$
degree $t$ has cost at least $mt$. If $\rho\ge0$, the pure
$y^jz^t$ columns supply a full leading-component basis by cost
$mt+s-1$. If $\rho<0$, use $x^ty^j$ instead. These also supply
a full leading-component basis by that cost. Some have $j\ge d$
and are not in (R3), but reducing the graph relation expresses them
in strictly cheaper candidates with total $x,z$ degree at most $t$.
Consequently their span is available by the same deadline in the
candidate scan. Since $mt+s-1<m(t+1)$, induct on $t$. Before the next block
starts, the selected columns span every component of $x$ degree at
most $t$. The lower $x$ terms in (R4) are already in the span of
earlier blocks. Later-arriving columns from an earlier block are
dependent. After $t=r-1$, all $rs$ dimensions are filled, before
any candidate with $t\ge r$ can be selected. If $\rho<-d$, the
three original generators have pairwise coprime leading monomials
$x^r,y^s,z$, and hence form a Gröbner basis. This proves the remaining
assertion of Theorem 1. Equality $\rho=-d$ requires a tie-breaking
order and is not an additional open region.

The proof also works in characteristic $p>d(r-1)+s-1$: the exponents
in (R6) remain distinct modulo $p$ and all required factorials are
invertible. No extension to small characteristic is inferred from
this argument; its failures there are consistent with the torsion
classification in the earlier note.

## A linear common weight list for this entire cone

Every slot has $0\le i<r$ and $0\le j<s$. Every comparison
breakpoint belongs to the finite set

$$
\mathcal Q_{r,s}=\{0\}\cup\{\pm k/h:1\le k<s,\ 1\le h<r\}.\tag{R7}
$$

There are at most $2(r-1)(s-1)+1$ distinct such numbers. The graph-head
transition $-d$ is in this set if $d<s$. If $d\ge s$, all $s$
cheapest slots have $i=t$ whenever $\rho<-(s-1)$, yielding the
ordinary $x,y$ rectangle on both sides of that extra transition.
Thus it creates no missing open state. One $\rho$ from every interval
in the complement of (R7) determines all selections for every $d\ge1$.
The number of distinct initial ideals whose cones meet $\mathcal C_0$
is at most

$$
2(r-1)(s-1)+2\le2rs.                                     \tag{R8}
$$

There is an explicit coefficient-independent common list. Include
the mediant of consecutive finite endpoints in (R7). For the two
unbounded intervals use $\rho=-s+1/2$ and $\rho=s-1/2$.
If $\rho=A/B$, take

$$
(w_x,w_y,w_z)=
\begin{cases}(sB+A,B,sB),&A\ge0,\\
             (sB,B,sB-A),&A<0.
\end{cases}                                               \tag{R9}
$$

All denominators are positive. The inclusion of zero in (R7) ensures
no mediant crosses signs. These weights lie in $\mathcal C_0$ and
have height at most $2rs$: internal mediants have denominator at most
$2(r-1)$ and numerator at most $2(s-1)$ in absolute value; the two
outer samples give height $4s-1\le2rs$.
They cover all the initial ideals in (R8), simultaneously for every
integer $d\ge1$ and every nonzero scaling of the three coefficients.
This is a regional list. It does not sample the entire positive fan.
The counting and weight-height steps are elementary arrangement
consequences of the exact rank rule, not new general enumeration tools.

## The bilinear case has an exact regional count

For $d=1$, the slots are the pure $j=0$ horizontal points and the
vertical $i=0$ points. A threshold tie between $\rho i$ and $j>0$
in a block with $t\ge i$ has exactly $i+j-1$ cheaper slots. It changes
the $s$-point cut precisely when $i+j=s$. Therefore all and only the
positive regional walls are

$$
\rho=\frac{s-i}{i},\qquad 1\le i\le\min(r-1,s-1).          \tag{R10}
$$

The ratios are distinct and each changes at least the block $t=i$.
There are exactly $\min(r,s)$ initial ideals meeting the subcone
$\mathcal C=\mathcal C_0\cap\{w_x>w_z\}$.
When crossing (R10) upwards, all blocks $i\le t<r$ exchange

$$
x^iz^{t-i}\quad\longrightarrow\quad y^{s-i}z^t.            \tag{R11}
$$

The difference of the two exposed basis profiles is
$(r-i)(i,-(s-i),-i)$. For $-1<\rho<0$, each block $t<s$
selects all its $t$ horizontal slots and vertical slots
$j=0,\ldots,s-t-1$. Each block $t\ge s$ instead selects the last
$s$ horizontal slots. This is a single state. It equals the first
positive state precisely when $r\le s$; otherwise $\rho=0$ is
one further wall. The state at $\rho<-1$ is the $x,y$ rectangle,
different from all these states because it contains $xy$. Therefore
the entire outer cone has exactly

$$
g_{\mathcal C_0}(I)=
\begin{cases}r+1,&r\le s,\\s+2,&r>s.\end{cases}             \tag{R10a}
$$

Intermediate optimal vector bases at walls are not asserted to be
global monomial initial standard sets.

## Cubic direction products occur in this fixed bilinear family

Take $r=s=q=2a-1$, $a\ge2$, and $i=a$ in (R10). The wall has
primitive normal

$$
h_a=(a,-(a-1),-a).
$$

It is primitive because $\gcd(a,a-1)=1$. Theorem 1 gives two
adjacent global initial cones on its two sides inside $\mathcal C$;
this is a coherent wall of actual length-$n=q^2$ quotients. Its
direction product is

$$
\prod_{j=1}^3(|h_{a,j}|+1)=a(a+1)^2=\Theta(n^{3/2}).       \tag{R12}
$$

In particular no fixed constant $C$ bounds these primitive products
by $Cn$ for the whole characteristic-zero bilinear class. The
linear-direction-box argument that closed the additive nested family
cannot be imported here. This is not a lower bound above the generic
fan scale: the whole outer cone has only $q+1$ initial ideals. Large
coherent exchange directions and many initial ideals remain different
requirements.

## Scope, sources and verification

### The remaining bilinear fan is an explicit polygon-prefix problem

For $d=1$, the translated model gives a useful reduction valid for the
**whole** positive fan. Write $(\alpha,\beta,\gamma)=(w_x,w_y,w_z)$.
If $\gamma>\alpha+\beta$, the standard set is the $x,y$ rectangle.
Otherwise assume $\gamma<\alpha+\beta$, put
$\delta=\alpha+\beta-\gamma>0$, and use the carrier $XY=Z$.
The column of $X^aY^bZ^c$ depends only on
$(U,V)=(a+c,b+c)$. Among all representatives of this node, the
cheapest has $c=\min(U,V)$, hence at least one of $a,b$ is zero.
Its cost is

$$
f(U,V)=\alpha U+\beta V-\delta\min(U,V)
=\max\{\alpha U+(\gamma-\alpha)V,
        (\gamma-\beta)U+\beta V\}.                       \tag{R13}
$$

The full fan is therefore greedy rectangular polynomial interpolation,
with columns $(U^iV^j)_{i<r,j<s}$, on integer-node prefixes in the
convex polygon $f(U,V)\le L$, $U,V\ge0$. This is the fixed Pascal
model from the earlier note, not generic evaluation at independently
chosen nodes. Exact ranks are still needed.

One finite ground suffices. Put $D=r+s-2$ and retain only

$$
H_{r,s}=\{(U,V)\in\mathbb N_0^2:
 \max(U,V)\le D,\quad -(s-1)\le U-V\le r-1\}.             \tag{R14}
$$

Indeed the original local quotient has maximal ideal power
$\mathfrak m^{D+1}=0$, because $F\in(x,y)$. Translation preserves
initials, so no standard monomial has degree exceeding $D$; for a
canonical representative that degree is $\max(U,V)$. The leading
terms $X^r,Y^s$ impose the two difference bounds in (R14). These
are valid exclusions, not an assertion that every remaining column
is independent. The translated and original individual columns need
not have the same ranks; their global initial standard sets do.

There is one more uniform cone. If
$\max(\alpha,\beta)<\gamma<\alpha+\beta$, the function (R13)
is strictly increasing in each coordinate. A rectangular Lagrange
interpolation formula expresses evaluation at $(U,V)$ using only grid
nodes $(i,j)$ with $i<r,j<s$ and $i\le U,j\le V$: if $U<r$,
only $i=U$ is needed, and otherwise every $i<r$ is smaller; the
same holds for $V$. All distinct stencil nodes are strictly cheaper.
Consequently the selected nodes are exactly
$\{0,\ldots,r-1\}\times\{0,\ldots,s-1\}$, mapped back by
$(U,V)\mapsto((U-V)_+,(V-U)_+,\min(U,V))$.
This uses ordinary tensor Newton/Lagrange interpolation and is not a
new lower-set theorem. For background, the
[Neidinger exposition](https://doi.org/10.1137/17M1124188) explicitly
describes the tensor/lower-set construction as classical; its abstract
was inspected here, not its full proof.

The unresolved counting problem is thus concentrated in the polygons
with $\gamma<\max(\alpha,\beta)$ outside the outer regions already
controlled. One cannot extend the separated-block formula there:
at $r=s=3$ and $w=(13,10,5)$, the actual standard set contains $z^3$
in place of the formula's $y^2$. No full-fan count follows from
convexity of the node prefix or from the number of its lattice points.

The [earlier bilinear rank dictionary](bilinear_rectangle_fans.md)
and its primary-source discussion remain in force. Formula (R6) is
ordinary Vandermonde interpolation in the falling-factorial basis;
no arbitrary-subset interpolation characterization or general tensor-code
rank result is used. The new structural step is separating the $x$-degree
blocks by cost and proving the prefix rule through disjoint exponent
intervals. The exact prior-art comparison remains incomplete, so this
note asserts no new general interpolation principle or verified novelty.

```bash
python3 -B scripts/run_rectangle_boundary_checks.py
```

The checker compares the slot model with direct truncated-polynomial
greedy selection in characteristic zero and a qualifying positive
characteristic. Rational border bases, ordinary S-pairs and exact
coverage after intersection with $\mathcal C_0$ check complete regional
fans. Independent leading-component ranks test every slot prefix.
The cubic-direction fixtures also verify actual shared facets of the
rational quotient cones and their profile differences.

The full characteristic-zero rectangle count and the unrestricted
extremal endpoint remain open. The next uniform argument must handle
interacting degree blocks in the remaining weight regions, or find a
different prescribed dependency family with a proved count. The regional
count and obstruction here do not lower the active goal's success bar.
