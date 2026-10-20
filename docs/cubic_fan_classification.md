# Initial ideals of truncated Markoff cubics

Let $K$ be any field, $\mathfrak m=(x,y,z)$, and

$$
I_s=(x^2+y^2+z^2+xyz)+\mathfrak m^{s+1},\qquad n_s=(s+1)^2.
$$

This note proves a uniform classification for this family. It closes
the boundary-cover question in [the truncation note](hypersurface_truncations.md).
The count is linear in length, despite coordinate relation depths of
order $\sqrt{n_s}$. Thus these deep relations do not give the polynomial
separation sought by the unrestricted $G_3$ goal. That goal, its priority
audit and its paper requirements remain open. The classification below
has a proof; its novelty and publication significance have not been
established by the limited comparisons at the end.

## The statements

For coprime integers $0<p<q$, put $\epsilon(p,q)=1$ when both are odd
and $\epsilon(p,q)=2$ otherwise. Define

$$
\mathcal F_s=\left\{\frac pq:\gcd(p,q)=1,\quad
 \epsilon(p,q)(3q+p)\le2(s+1)\right\},\qquad M_s=|\mathcal F_s|.
                                                               \tag{C1}
$$

**Theorem.** For every field and $s\ge3$,

$$
\#\{\text{monomial initial ideals of }I_s\}=6(M_s+1)
 =\frac3{\pi^2}n_s+O(\sqrt{n_s}\log(n_s+1)).              \tag{C2}
$$

The implicit constant is absolute. The same monomial initial ideals
occur in every characteristic. The assertion also holds for
$(a x^2+b y^2+c z^2+dxyz)+\mathfrak m^{s+1}$ with $abcd\ne0$.
For $s=0,1,2$, the counts are respectively $1,1,3$.

Here is an exact selection and rank rule. Permute coordinates so that
$w_x\ge w_y\ge w_z>0$. If $w_x>w_y+w_z$, the standard set is
the monomials of degree at most $s$ with $x$-exponent below two.
In the strict central region normalize

$$
w=(u+\lambda,u,1),\qquad u>1,\quad0<\lambda<1.           \tag{C3}
$$

Include $1,z,\ldots,z^s$. For each $1\le k\le s$, put $L=s-k+1$
and merge the following three already ordered lists using the
displayed additional prices:

$$
\begin{array}{c|c}
\text{monomials}&\text{additional price}\\\hline
x^a y^{k-a},\quad0\le a\le k&\lambda a\\
y^kz^c,\quad1\le c<L&c\\
x^kz^c,\quad1\le c<L&\lambda k+c.
\end{array}                                                \tag{C4}
$$

Take the cheapest $2L$ terms. Their union is the standard set of
$I_s$. Generic prices avoid ties. A limit from a neighboring generic
region handles order refinements at ties.

For $[v]_h=\max(0,\min(v,h))$, let

$$
\begin{split}
R_0(s,T;u,\lambda)
={}&[\lfloor T\rfloor+1]_{s+1}\\
&+\sum_{k=1}^s\min\left(2L,
 [\lfloor(T-uk)/\lambda\rfloor+1]_{k+1}
 +[\lfloor T-uk\rfloor]_{L-1}
 +[\lfloor T-uk-\lambda k\rfloor]_{L-1}\right).
\end{split}                                                \tag{C5}
$$

**Algorithmic consequence.** This is the exact rank, in $S/I_s$,
of all original monomials whose $w$-price is at most $T$. It uses
$O(s)=O(\sqrt{n_s})$ floor/arithmetic operations and constant auxiliary
storage. The standard set itself is produced in $O(n_s)$ comparisons
by merging (C4). These bounds count exact rational arithmetic
operations; they are not bit-complexity or general Gröbner-basis
running-time claims. No field arithmetic is needed to identify the
initial monomial ideal or answer these price-prefix rank queries.
An outer-cone query counts the elementary standard monomials by two
one-dimensional sums, also in $O(s)$ operations. The original
smooth-number CLI and its certification scope are unchanged.

## Boundary blocks and three streams

Put $J_s=(x^2+y^2+xyz)+\mathfrak m^{s+1}$. Highest weight at
$v=(1,1,0)$ gives $\operatorname{in}_v I_s=J_s$: inclusion is
immediate and both quotients have length $(s+1)^2$ by monic division
in $x^2$. Every monomial basis of $J_s$ lifts to a basis of $I_s$.
For a coherent basis, the usual refinement lemma gives a lifted cone
incident to $v$.

The degree-$k$ summand for the $v$-grading, $k\ge1$, has the exact
column model

$$
x^a y^{k-a}z^c\longmapsto t^az^c
\quad\text{in }K[t,z]/(t^2+zt+1,z^L),\quad L=s-k+1.      \tag{C6}
$$

This is an identity of finite column spaces, not localization of the
original nilpotent quotient. The block has dimension $2L$.
If $0<a<k$ and $c>0$, then

$$
t^az^c=-t^{a+1}z^{c-1}-t^{a-1}z^{c-1}.                 \tag{C7}
$$

Both terms have strictly smaller $\lambda a+c$ price when
$0<\lambda<1$. Iteration reaches the three lists in (C4).

We prove that the cheapest $2L$ terms of those lists are independent
over every field. Work for this proof in the separate ring

$$
K[t,z]/(t^2+zt+1)\simeq K[t,t^{-1}],\qquad z=-t-t^{-1}.
$$

The additional equation $z^L=0$ is $(t^2+1)^L=0$ up to a unit.
If none of the third list is selected, the set is
$1,t,\ldots,t^A,z,\ldots,z^C$ with $A+C=2L-1$ and
$C\le L-1$, so $A>C$. Multiplying any proposed relation by $t^C$
gives a polynomial of degree at most $A+C<2L$. Divisibility by
$(t^2+1)^L$ forces it to vanish identically. The negative Laurent
powers of a nonconstant polynomial in $z$ cannot cancel a polynomial
in $t$, so every coefficient is zero.

If the third list is used, all of $1,t,\ldots,t^k$ occur before it.
The selected set is those powers together with
$z,\ldots,z^C,t^kz,\ldots,t^kz^D$, where
$k+C+D=2L-1$ and $C,D\le L-1$. These inequalities imply
$|C-D|\le k$. A relation has the form

$$
P(t)+Q(z)+t^kH(z)=0\pmod{(t^2+1)^L},\quad
\deg P\le k,\quad Q(0)=H(0)=0.                          \tag{C8}
$$

After multiplication by $t^C$, its Laurent range is contained in
$[0,k+C+D]$, so again its degree is less than $2L$ and (C8) is an
exact Laurent identity. Write $t^a=A_a(z)t+B_a(z)$ using (C6).
For $a\ge1$, $A_a$ has degree $a-1$ and leading coefficient
$(-1)^{a-1}$, in every characteristic. The $t$-coefficient of
$t^kH(z)$ has degree at least $k$ unless $H=0$, whereas the
$t$-coefficient of $P$ has degree at most $k-1$. Thus $H=0$;
negative Laurent powers then give $Q=P=0$.

This proves the three-stream rule for $J_s$, and hence (C5) for it.
One may also use infinite second and third lists: their first $2L$
terms never use $z^L$ or $t^kz^L$. If $k\ge L$, the $L+1$
powers $1,t,\ldots,t^L$ and the $L-1$ smaller $z$-powers precede
$z^L$. If $k<L$, use all $k+1$ powers, the $L-1$ smaller $z$-powers
and $t^kz^d$ for $1\le d\le L-k$. This supplies $2L$ terms of
strictly smaller price. This observation will be useful below.

## The original cubic at two equal prices

The remaining issue is whether the lifted cones cover the ordered
central region. It is enough to prove that the same bases minimize
price for $w=(1+\lambda,1,1)$. Cone convexity then adds any
nonnegative multiple of $v$.

Set $U=y^2+z^2$ and $V=yz$. Over any field, and also over
$\mathbb Z$, $K[y,z]$ is free over $K[U,V]$ with basis

$$
1,\ y,\ z,\ z^2,\qquad\text{degrees }\delta=0,1,1,2.    \tag{C9}
$$

For example the three monic relations
$y^2+z^2-U$, $yz-V$, $z^3-Uz+Vy$ form a Gröbner basis over
$\mathbb Z[U,V]$ for graded lexicographic order in $y,z$.
Their leading monomials are $y^2,yz,z^3$, and their S-pairs reduce
to zero. In particular $zf_1-yf_2=f_3$ and
$z^2f_2-yf_3=Uf_2-Vf_1$; the remaining leading pair is relatively
prime. This proves (C9) without division by two.

Modulo the cubic, $U=-x^2-xV$. Replacing $U$ never decreases
ordinary degree when $x$ has degree one and $V$ degree two.
Consequently the exact truncated module is

$$
\bigoplus_{\delta=0,1,1,2}t_\delta\,
 K[x,V]/\langle x^aV^b:a+2b>s-\delta\rangle.             \tag{C10}
$$

At two equal prices, all original binary monomials of each fixed
$y,z$ degree are present together in a price prefix. The homogeneous
change of basis (C9) therefore preserves that prefix's span. We may
replace them by the columns $t_\delta x^aU^iV^j$, at price
$\delta+(1+\lambda)a+2i+2j$.

After $U$ is eliminated these columns are homogeneous of degree
$q=a+2i+j$ in the two auxiliary variables $x,V$. The degree-$q$
summand of (C10) is $K[\xi]/\xi^L$ as a column space, where

$$
L=\min(q,s-\delta-q)+1,\qquad
x^aU^iV^j\longmapsto(-1)^i\xi^j(1+\xi)^i.              \tag{C11}
$$

Only $0\le q\le s-\delta$ occurs. Again this is a finite column
description; no nilpotent variable has been inverted.

For a price bound $T'=T-\delta$, a fixed $j$ includes precisely
an interval $i_j\le i\le\lfloor(q-j)/2\rfloor$. The lower endpoints
$i_j$ are nondecreasing in $j$, the upper endpoints are nonincreasing,
and the nonempty intervals occur for consecutive $j$ starting at zero.
Put

$$
i_0=\max\left(0,\left\lceil\frac{(1+\lambda)q-T'}{2\lambda}\right\rceil\right).
$$

After factoring $(1+\xi)^{i_0}$, the $j=0$ interval spans all
polynomials up to its maximum degree. At each subsequent $j$ the
maximum degree $j+\lfloor(q-j)/2\rfloor-i_0$ increases by zero or
one. A column attaining a new degree is monic. Induction shows that
the entire prefix span is $(1+\xi)^{i_0}$ times all polynomials up
to the maximum degree. Reducing modulo $\xi^L$, a unit times that
space has rank equal to its dimension capped at $L$.

Equivalently, define the $q+1$ prices

$$
\kappa_{q,i}=\max\{2i,(1+\lambda)q-2\lambda i\},
\qquad 0\le i\le q.                                    \tag{C12}
$$

The rank of this summand is
$\min(L,\#\{i:\kappa_{q,i}\le T'\})$. To check the equivalence
directly, the latter set is the integer interval

$$
\max\left(0,\left\lceil\frac{(1+\lambda)q-T'}{2\lambda}\right\rceil\right)
\ \le i\le\ \min(q,\lfloor T'/2\rfloor).
$$

When the prefix is nonempty, its right endpoint equals the maximum
degree before subtracting $i_0$. Indeed, above the initial degree
$\lfloor q/2\rfloor$, degree $d$ is first reached at $j=2d-q$,
whose cheapest column has price $2d$. If the prefix is empty both ranks are zero.
This proves an independent exact price-prefix rank formula for the
original cubic, uniformly over all fields.

## Equality of the two rank formulas

It remains to prove that (C12), summed over the four channels,
gives (C5) at $u=1$. Take irrational $0<\lambda<1$ to avoid
incidental ties. We compare, for each price, the first cutoff at
which an independent column at that price becomes available.

The list (C12) is the merge of two lists:

$$
\begin{array}{ll}
q+\lambda a,&0\le a\le q,\quad a\equiv q\pmod2,\\
q+b,&1\le b\le q,\quad b\equiv q\pmod2.
\end{array}                                                \tag{C13}
$$

If an item is at position $r$ in this merged list, (C11) retains it
at precisely the cutoffs $s\ge\delta+q+r-1$.

For an item of the first list, with parity $\eta=a\bmod2$, its
position is

$$
r=\lfloor a/2\rfloor+1+
       \left\lfloor\frac{\lambda a+\eta}{2}\right\rfloor.
$$

Write $D=\delta+q$. Its price is $D+\lambda a$ and its first cutoff is

$$
D+\psi(a),\qquad
\psi(a)=\left\lfloor\frac{a+\lfloor\lambda a\rfloor}{2}\right\rfloor.
                                                               \tag{C14}
$$

For $D>a$, the four channel degrees supply exactly two copies of this
item: degrees $0,2$ when $D-a$ is positive and even, or the two
degree-one channels when it is odd. For $D=a$ there is one copy.
These match $x^ay^{D-a}$ and $x^az^{D-a}$, counting their common
pure power only once. For $a=0$, they match the two pure $y,z$
powers, with the constant counted once. Their positions in the
infinite three-stream lists give exactly the same first cutoff
(C14): a spine item has position $a+1+\lfloor\lambda a\rfloor$;
an upper-arm item adds $2(D-a)$ to that position.

The second list of (C13) accounts for mixed $y,z$ monomials. Write
$q=b+2l$ with $b\ge1$, $l\ge0$, and $D=\delta+q+b$. Its price is
$D$. Counting the preceding items in (C13) shows that its first
cutoff is

$$
D+\min\left(l,\left\lfloor\frac{(1-\lambda)b}{2\lambda}\right\rfloor\right).
                                                               \tag{C15}
$$

In the boundary lists, a mixed monomial $y^kz^c$ with $k+c=D$
has price $D$ and first cutoff

$$
D+\left\lfloor\min\left(\frac{(1-\lambda)k}{2},
                       \frac{(1-\lambda)c}{2\lambda}\right)\right\rfloor.
                                                               \tag{C16}
$$

This follows by counting the spine and both arms before its price;
the cases $c<\lambda k$ and $c>\lambda k$ give the two alternatives.
We show that (C15) and (C16) give the same multiset of excess
cutoffs for every $D$. Fix an integer $j\ge1$ and put

$$
B_j=\left\lceil\frac{2\lambda j}{1-\lambda}\right\rceil.
$$

In (C16) an excess at least $j$ requires
$k\ge2j+B_j$ and $c\ge B_j$, because
$2j/(1-\lambda)=2j+2\lambda j/(1-\lambda)$.
The number of such pairs is

$$
\max(0,D-2j-2B_j+1).                                    \tag{C17}
$$

If $D=2m$, (C15) comes from $b+l=m$ and $b+l=m-1$.
The numbers with excess at least $j$ are respectively
$[m-j-B_j+1]_+$ and $[m-j-B_j]_+$; their sum is (C17).
If $D=2m+1$, there are two copies of $b+l=m$, giving
$2[m-j-B_j+1]_+$, again (C17). The total multiplicity in either
case is $D-1$. Equality of all these superlevel counts proves
equality of the multisets, including excess zero.

Thus every price has the same multiplicities at every cutoff in
the two constructions. In particular,

$$
\operatorname{rank}_{S/I_s}\{x^ay^bz^c:
 (1+\lambda)a+b+c\le T\}=R_0(s,T;1,\lambda).             \tag{C18}
$$

Both arguments work over every field. No determinant-sign or
support-matching assertion has been used.

## From equal prices to the whole positive fan

A three-stream standard set $B$ is a basis of $I_s$ by refinement
at $v$. Equation (C18) says that it has the minimum possible price
at $(1+\lambda,1,1)$. Its original reduced border inequalities
therefore hold nonstrictly at this vector and at $v$: a nonzero border
coefficient permits a one-column basis exchange, and minimal total
price gives exactly the corresponding inequality.
Adding $(u-1)v$ gives (C3). A border tail of smaller $v$-degree
has a strict inequality when $u>1$; a tail of equal $v$-degree
has the strict boundary inequality inside its sector. This proves
that $B$ is the original standard set on the entire ordered sector.

The proof used irrational $\lambda$, which are dense in each boundary
sector. The border inequalities are linear and continuous, so their
nonstrict endpoint inequalities hold on the closed interval as well.
In particular this proves the previously conjectured endpoint
inequality (13) in the truncation note for every cutoff and field.
All positive order regions follow by coordinate permutation. The
three outer cones follow from the leading square and equal colength.
Monomial initial ideals have open weight regions, so weight ties
cannot conceal an additional initial ideal. Refining ties, and then
counting the selected basis at the original prices, also extends
the rank query to nongeneric positive weights.

## The exact walls and the count

In a degree-$k$ block only two kinds of cross-list ties can change
the cheapest $2L$ terms. A spine item $a$ and a lower-arm item $c$
tie at $\lambda=c/a$. Exactly $a+c-1$ items precede them, so this is
a wall precisely when $a+c=2L$ and $a\le k$.
A lower-arm item $c$ and an upper-arm item $d$ tie at
$\lambda=(c-d)/k$. Exactly $k+c+d-1$ items precede them, so this
is a wall precisely when $k+c+d=2L$. A spine item cannot tie with
the upper arm before it. All these ties involve two items within
one block; a rational slope shared by several blocks is counted once.

For a reduced fraction $p/q$, the smallest integral representative
$(a,c)=m(q,p)$ with $a+c$ even uses $m=\epsilon(p,q)$.
The first kind of wall occurs at some block of cutoff $s$ exactly when

$$
s\ge\frac{m(3q+p)}2-1.
$$

It persists at every larger cutoff by increasing $k$ and keeping $L$
fixed. The second kind can first occur only one cutoff later for the
same representative, so adds no new slope. This proves that the
positive ordered boundary sectors are exactly the intervals cut out
by (C1), and their number is $c_s=M_s+1$.

For completeness we count their coordinate permutations without
assuming that all give different standard sets. In a central standard
set, the pure $s$th power of the largest-price variable is absent,
whereas both other pure $s$th powers are present. This identifies
which variable is largest. Fix it to be $x$. Formula (C14) shows
that the $x$-bearing part of the basis is invariant under swapping
$y,z$ and decreases as $\lambda$ increases. At every wall one such
monomial is lost, by the two tie descriptions above. Thus different
sectors have different $x$-bearing parts, even after swapping $y,z$.

The mixed $y,z$ part is invariant under that swap only in the last
sector. Indeed, if $D\le s$ and $m=s-D+1$, (C16) excludes precisely
the interval $k\ge2m+B_m$, $c\ge B_m$ within $k+c=D$.
When nonempty, this interval differs from its reflected interval,
as their left endpoints differ by $2m>0$. Hence a swap-invariant
mixed part must contain the full degree-$\le s$ $y,z$ triangle.
It does so as $\lambda\to1$. Monotonicity and fixed basis cardinality
show that this happens exactly in the final sector.

For each choice of largest variable, the first $c_s-1$ sectors
therefore give two sets and the final one gives one. There are
$3(2c_s-1)$ central sets. The three outer sets are different: they
contain $xyz$ when $s\ge3$, whereas central sets exclude it.
The result is $6c_s$, proving the exact formula in (C2).
For $s<2$ the quadratic is beyond the cutoff. At $s=2$ its leading
monomial is the square of the largest-price variable, giving the
three stated initial ideals.

Primitive-pair counting gives the asymptotic. The region
$0<p<q$, $3q+p\le X$ has area $X^2/24$. Coprime odd pairs
have density $2/\pi^2$, and coprime pairs of opposite parity have
density $4/\pi^2$. Möbius inversion, with the elementary boundary
error summed over divisors, gives $O(X\log(X+1))$ error for each
count. Applying the odd count at $X=2(s+1)$ and the opposite-parity
count at $X=s+1$ yields

$$
M_s=\frac{(s+1)^2}{2\pi^2}+O(s\log(s+1)).                \tag{C19}
$$

For example, the odd density is
$\tfrac14\sum_{d\text{ odd}}\mu(d)/d^2=2/\pi^2$;
subtracting it from $1/\zeta(2)$ gives the other density.
This is the classical visible-lattice-point estimate, not a new
number-theoretic input.

Finally, diagonal scaling over a field extension takes every cubic
with $abcd\ne0$ in the statement to the displayed cubic, up to a
nonzero scalar. One can take the common square coefficient
$abc/d^2$, choose the required square roots, and adjust a sign if
necessary. Diagonal scaling and field extension preserve every
monomial-column rank and every monomial initial ideal. This proves
the coefficient extension, including characteristic two.

## Scope, prior work and reproduction

The surface itself and its Vieta symmetries are classical;
[Goldman's introduction and Section 5.2](https://math.umd.edu/~wmg/RemarkableFamily.pdf)
were inspected for that identification. No new Markoff-surface
construction or arithmetic consequence for primes is claimed.

Cone inequalities, refinement and convexity use the established
[Fukuda--Jensen--Thomas, Proposition 2.6, Lemma 2.10,
Corollary 2.11, Lemma 2.13 and Corollary 2.14](https://arxiv.org/pdf/math/0509544),
including the proofs inspected in the preceding note. Monic division,
graded freeness, binary polynomial interpolation and Möbius inversion
are classical ingredients. The new object to compare with prior work
is the complete classification and resulting rank oracle for these
specific nonhomogeneous truncations, not those ingredients.

[Dimitrova--He--Robbiano--Stigler, Proposition 2.7 and proof](https://arxiv.org/pdf/1704.03680)
already gives invariance of standard sets under coordinatewise affine
scaling/translation. That proposition and the paper's definitions
were inspected; the diagonal-scaling step above is an instance of it.
The one-basis characterization and distraction statements in Section 3
do not apply to the increasing counts in (C2). These specific comparisons
are not an exhaustive priority audit. Theorem 4.8 and its proof were
also inspected: they transfer fans between ideals complementary in a
grid, but provide no numerical count without an additional identified
complement and its fan. In particular the strongest
current literature on Markoff jets, fat conics and filtered rank formulas
still needs comparison before any novelty or paper claim.

The supporting checker independently compares exact linear algebra
in the original monomial quotient with both rank formulas, the
three-stream basis and the wall classification. It retains all
previous bounded certificates. Agreement is a check on the uniform
proof above, not a substitute for that proof.

~~~bash
python3 -B scripts/run_cubic_fan_checks.py
python3 -B scripts/run_hypersurface_truncation_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

This is an infinite-family closure and an explicit complexity consequence.
It leaves the unrestricted $G_3$ bound, the absence of a polynomial
lower separation and the original acceptance criteria unchanged.
