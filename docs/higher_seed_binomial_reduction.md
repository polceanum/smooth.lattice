# Higher-interaction seeds as finite binomial matrices

The remaining higher-interaction problem is the largest-$z$ region for
$I_{s,r}=(x^2+y^2+z^2+xyz^r)+\mathfrak m^{s+1}$, $r\ge3$.
The full $r=2$ case is already [closed](equal_degree_unit_fans.md).
Here the two-cheap-equal seed is reduced uniformly to explicit
one-variable binomial matrices, including its integral arithmetic.
The resulting rank identity is still unproved. There is no new fan-count
bound in this note.

## A regular-pair rank theorem

Take integers $r>a\ge1$, nonzero homogeneous binary forms $Q,H$ of
degree $a$, and
$$
A_s(Q,H)=K[x,y,z]/(z^a+Q+Hz^r,\mathfrak m^{s+1}).
$$
Write $\operatorname{rk}_{Q,H}(s,u,Z,T)$ for the rank of all original
monomials of price at most $T$, at prices $(u,u,Z)$ with $0<u<Z$.
The quotient length is
$$
n_{s,a}=\sum_{e=0}^{\min(a-1,s)}\binom{s-e+2}{2}. \tag{B1}
$$

For the core model use formal $U,V$ of ordinary degree $a$, and $z$ of
degree one. Its finite parent ground is
$$
E_{s,a}=\{U^iV^jz^c:a(i+j)+c\le s\}.
$$
Both cores below have the reference rows
$V^b z^d$, $ab+d\le s$:
$$
G_{s,a,r}: U=-z^a-Vz^r,
\qquad C_{s,a,r}: U=-\frac{z^a}{1+z^r}. \tag{B2}
$$
All substitutions are truncated by ordinary degree. In particular the
second denominator is a finite integral unit expansion. Denote their
complete-prefix ranks at prices $(au,au,Z)$ by $g_s(T)$ and $c_s(T)$.

Let $p_a(\delta)$ be the coefficient of $t^\delta$ in
$(1+t+\cdots+t^{a-1})^2$, and use rank zero at negative cutoffs.

**Theorem 1.** Over every field,
$$
\sum_\delta p_a(\delta)c_{s-\delta}(T-\delta u)
\ \le\ \operatorname{rk}_{Q,H}(s,u,Z,T)
\ \le\ \sum_\delta p_a(\delta)g_{s-\delta}(T-\delta u). \tag{B3}
$$
If $Q,H$ are coprime, equality holds on the right. In particular their
complete-tier ranks depend on the coefficients only through
coprimality; they can depend on the characteristic. The left side is
the exact rank for $Q=H=y^a$, with $x$ free.

**Proof.** A finite monic reduction gives the $z<a$ reference in (B1).
For a monomial with $z$-exponent at least $a$, rewrite it using
$z^a=-Q-Hz^r$. The $Q$ branch keeps ordinary degree and decreases
$z$-exponent; the $H$ branch increases ordinary degree by $r$.
The reduction terminates on the finite ground and does not decrease
ordinary degree. Its entries are polynomials in the coefficients.
The lowest ordinary form $z^a+Q$ has degree $a$ and is monic in $z$;
its Hilbert function gives (B1). Thus this spanning reference is a
basis, with no exceptional field. The same reduction works integrally.

For a regular pair the homogeneous free decomposition is
$$
K[x,y]=\bigoplus_\nu t_\nu K[Q,H],\qquad
\sum_\nu t^{\deg t_\nu}=(1+t+\cdots+t^{a-1})^2.
$$
Generation from lifts modulo $(Q,H)$ and equality of Hilbert series
prove it. At equal binary prices this change preserves every entire
binary homogeneous space, hence every complete price tier and ordinary
cutoff. Eliminate $Q=-z^a-Hz^r$ to obtain the direct core sum on the
right of (B3).

The generic coefficient pair is coprime. Polynomial specialization of
the fixed $z<a$ reference matrices cannot increase rank; their length
does not change. This proves the right inequality for any pair.

For the lower inequality, extend the field and change the equally
priced binary coordinates so that both pure-$y$ coefficients are
nonzero. Scaling the other coordinate to zero specializes to
$z^a+A y^a+B y^a z^r$, $AB\ne0$.
Rescale $z=Cz'$, $y=Dy'$ with $C^r=A/B$ and $D^a=C^a/A$ to obtain
the coefficient-one canonical equation. Roots can be taken in an
algebraic field extension even in inseparable characteristic, and
ranks are unchanged by extension. Finally use the integral free
decomposition over $K[x^a,y^a]$ with residues $x^e y^f$, $e,f<a$.
This is precisely the direct sum of canonical cores in (B3).
$\square$

### Exact scalar blocks for the regular core

For a reference row $V^b z^d$ put $k=d+(a-r)b$ and
$$
b_0=\max(0,\lceil-k/(r-a)\rceil),\quad
b_1=\lfloor(s-k)/r\rfloor,\quad L=b_1-b_0+1.
$$
Omit empty blocks. Their reference rows have
$b=b_0,\ldots,b_1$. A parent column has
$$
j\in[b_0,b_1],\quad
0\le i\le\left\lfloor\frac{k+(r-a)j}{a}\right\rfloor,
\quad c=k-ai+(r-a)j.
$$
Its column and price are
$$
(-1)^i t^{j-b_0}(1+t)^i\pmod {t^L},\qquad
kZ-a(Z-u)i+(au+(r-a)Z)j. \tag{B4}
$$
These statements follow by expanding the first substitution in (B2):
the term with index $h$ has row $V^{j+h}z^{ai+c+(r-a)h}$,
coefficient $(-1)^i\binom ih$, and the same $k$.
Reference columns $i=0$ are the unit coordinate columns, proving all
kernel and integral-lattice assertions for the model.

At a cutoff $T$, the allowed $i$ at level $j$ form the explicit interval
$$
\max\left(0,\left\lceil
\frac{kZ+(au+(r-a)Z)j-T}{a(Z-u)}\right\rceil\right)
\ \le i\le\left\lfloor\frac{k+(r-a)j}{a}\right\rfloor. \tag{B5}
$$
The block rank is the rank of these finite binomial columns, not their
number capped at $L$. Arbitrary binomial minors need not be units.

### Exact scalar blocks for the canonical core

Fix the free exponent $b$ of $V$ and the residue $0\le k<r$ of
$ai+c$. Put $L=1+\lfloor(s-ab-k)/r\rfloor$ and omit $L\le0$.
For $0\le j<L$ and $0\le i\le\lfloor(k+rj)/a\rfloor$, the columns
and prices are
$$
(-1)^i t^j(1+t)^{-i}\pmod {t^L},\qquad
au b+kZ-a(Z-u)i+rZj. \tag{B6}
$$
The coefficient at row $j+h$ is
$(-1)^{i+h}\binom{i+h-1}{h}$ for $i>0$; for $i=0$ it is the unit
coordinate column. This is the second substitution in (B2), with
$c=k+rj-ai$.

Both block descriptions are row labels verified by the finite monic
references. Neither inverts $z$ in an Artinian ring.
For exact rational prices and cutoffs, fixed $a$ and $r\le s$,
all blocks have at most $1+s/r$ rows and
there are $O_a(rs)$ of them. Generating columns separately and using
field Gaussian elimination therefore gives an exact seed query with
$$
O_a(s^3(1+s/r)^2)\text{ field operations},\qquad
O_a((1+s/r)^2)\text{ field elements of scratch space}. \tag{B7}
$$
The total number of parent columns is $O_a(s^3)$; selecting them uses
$O_a(s^2)$ additional rational arithmetic/floor operations.
Integer binomial generation and bit complexity are separate. This is
a finite algebraic reduction, not a new speed claim for the public CLI.

## The exact identity still needed

For $r>a$ consider the following statement, over every field:
$$
g_s(T)=c_s(T)\qquad
\text{for all }s,\ 0<u<Z,\ T. \tag{B8}
$$
It holds for the equal-degree unit case $r=a$ by the earlier theorem,
but that proof does not establish (B8).
The new blocks display a possible mechanism: the integral change
$t\mapsto-t/(1+t)$ sends
$t^j(1+t)^i$ to $(-1)^j t^j(1+t)^{-(i+j)}$.
It converts the slope and prices in (B4) into those in (B6), but
retains the extra condition $i+j\ge j$ and different block labels.
Dropping that condition or ignoring the labels is not a proof of
rank equality.

Theorem 1 makes (B8) sufficient for a coefficient-independent rank
law for **all nonzero** binary forms of degree $a$, now allowing
characteristic dependence. Conversely, if the regular and canonical
whole-plane ranks agree for all $s,T$, (B3) gives the convolution of
$g_s-c_s$ with $p_a$. Its constant coefficient is one. Induction on
$s$ inverts the convolution and proves (B8). Thus this core identity
is equivalent to the uniform regular/canonical seed equality; it is
not merely another numerical sufficient test.

For $a=2$ it would prove the remaining original signed-face transfer.
Indeed it would equate the seed ranks for
$Q=x^2+y^2,H=xy$ and $Q=x^2,H=xy$. For $Z>X>Y$, write
$$
(X,Y,Z)=(u,u,Z_0)+t(1,1-r,1),\quad
t=(X-Y)/r,\quad u=((r-1)X+Y)/r,\quad Z_0=Z-t.
$$
Here $Z_0-u=Z-X>0$. The filtered lifting and border inequalities
from the [middle-price proof](middle_interaction_fans.md) would then
transfer the face standard sets in this last region. Together with
the existing regions this would give $g(I_{s,r})=O(n)$ uniformly.
This implication is conditional on (B8). Closing this family would
exclude a route to amplification, not resolve unrestricted $G_3$.

### The original seed's integral decomposition

For the coefficient-one quadratic pair
$U=x^2+y^2,V=xy$, the binary ring is integrally free with basis
$1,x,y,x^2$ over $\mathbb Z[U,V]$.
For example $y^2=U-x^2$, $xy=V$, $x^3=Ux-Vy$ and
$y^3=Uy-Vx$ give generation; independence over $\mathbb Q[U,V]$
and torsion-freeness give independence over $\mathbb Z[U,V]$.
The homogeneous-degree changes are unimodular.
Consequently the actual original seed prefix lattice has the exact
direct decomposition
$$
G_s(T)\ \oplus\ G_{s-1}(T-u)^2\ \oplus\ G_{s-2}(T-2u), \tag{B9}
$$
including reference ambient free modules. The canonical equation has
the analogous decomposition with $C$, using the binary-square residues.
Thus matching core Smith invariants would be stronger than (B8) and
would determine all prime-power torsion as well as field ranks.
Finite matching is not its uniform proof.

## Two exact boundaries on the proposed argument

The seed restriction matters. At $a=2,r=3,s=5$, use original weights
$(9,9,4)$ and price cutoff $20$. For the regular original equation
$x^2+y^2+z^2+xyz^3$ the prefix rank is $15$; for the canonical
$z^2+y^2+y^2z^3$ it is $14$, over every field.
Here $u>Z$, outside the remaining largest-$z$ domain.
This rules out a coefficient-universality statement at arbitrary prices.

Nor does an arbitrary divisor-closed selection suffice even in the
weighted core with $s=15,a=2,r=3$.
Take all $z^c$, $0\le c\le13$, all $Uz^c,Vz^c$, $0\le c\le1$,
and $U^2,UV,V^2$. These 21 columns have regular rank $21$ and
canonical rank $20$, over every field. Unit elimination leaves zero
residual matrices, certifying those ranks in every characteristic.
This selection cannot be a positive price prefix: it contains $U^2$
and $z^{13}$ but excludes $Uz^2$, whose exponent is a strict convex
combination of those two and the constant exponent.
Any universal proof of (B8) must use the coherent cutoff, rather than
an arbitrary downward staircase.

## Prior art and evidence

Homogeneous freeness, finite monic division, binomial expansion,
specialization and Smith form are classical. These exact reductions
and their priority have not been established as new research.

There is a specific classical self-similarity available for a possible
proof of (B8). In characteristic $p$, the identity
$(1+t)^p=1+t^p$ gives, by comparing coefficients,
$$
\binom ih\equiv\prod_\nu\binom{i_\nu}{h_\nu}\pmod p
$$
for the base-$p$ digits of $i,h$. For $h<p^m$, the negative-power
coefficients also satisfy
$\binom{-i}{h}\equiv\binom{(-i\bmod p^m)}{h}\pmod p$,
since $(1+t)^{p^m}\equiv1\pmod{(p,t^{p^m})}$.
Thus both sets of entries have an explicit prime-digit structure.
This is classical Lucas arithmetic, not a new fractal theorem.
The shifts, index bounds and coherent cutoff couple the digits;
entrywise self-similarity alone gives no rank identity.
[Rowland, arXiv:1704.05872v2](https://arxiv.org/pdf/1704.05872),
Section 1, Theorem 1 and its displayed proof outline count individual
binomial-coefficient valuations by digit matrices. They do not compute
the ranks of our selected matrices. A carry-aware rank decomposition
remains a proposed proof route, not a proved consequence.

[French, *Transformations Preserving the Hankel Transform*](https://cs.uwaterloo.ca/journals/JIS/VOL10/French/french13.pdf)
was checked at Example 1.2, Definitions 2.3--2.4 and Lemmas 2.9 and
2.14, including their proofs. It constructs unit triangular changes
whose congruences preserve leading Hankel determinants. No identification
of our selected rectangular binomial matrices with those congruences
has been established. Preserved leading determinants alone would not
prove all the ranks in (B8). This reading is a scoped comparison,
not an application or novelty certificate.

~~~bash
python3 -B scripts/run_higher_seed_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker verifies the exact block maps and rank sandwich against
independent original polynomial columns and Macaulay kernels. It compares
the cores throughout every open equal-price interval at specified
cutoffs and characteristics, checks integral Smith signatures at saved
prefixes, and verifies both all-field counterexamples.
Its finite agreement does not prove (B8), the full original transfer,
a new unrestricted count, or novelty. The goal remains active.
