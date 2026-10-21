# Integral transport of coherent binomial prefix quotients

This resolves identity (B8) in
[the higher-seed reduction](higher_seed_binomial_reduction.md) and hence
the Hilbert-function identity (P4) in
[the three-point model](three_point_rank_ideals.md). The earlier notes
preserve the stages at which those statements were open.

**Theorem 1.** For every integer $r>a\ge1$, cutoff $s\ge0$, real prices
$0<u<Z$, and complete price cutoff $T$, the regular and canonical core
prefix cokernels in (B2) are isomorphic as finitely generated abelian
groups. In particular their free ranks and nonunit Smith invariant
factors agree. The original core matrices have the same shape, so
their full Smith normal forms agree. Auxiliary presentations in the
proof can have different unit pivot counts. Over every field,
$$
g_s(T)=c_s(T). \tag{T1}
$$

The proof pairs the changes between whole blocks. Individual blocks
at the same ordinary degree need not have the same rank. No support
matching, division by binomial coefficients, or characteristic-zero
argument is used.

## One family of integral flags

Divide prices by $Z$, so $Z=1$, and write
$$
d=r-a>0,\qquad A=a(1-u)>0,\qquad B=r-A=d+au>0,
\qquad q=1+t.
$$
All spans below are integer spans in $V_L=\mathbb Z[t]/t^L$.
For $k\ge0$ and real $\tau$ put
$$
R_L(k,\tau)=\left\langle t^j q^i:
0\le j<L,\quad 0\le ai\le k+dj,\quad Bj-Ai\le\tau\right\rangle.
\tag{T2}
$$
For $0\le\rho<a$ let $S_L(\rho,\theta)=R_L(\rho,\theta)$ and
$$
F_L(\rho,\theta;h)=S_L(\rho,\theta)
 +\langle q^{-1},\ldots,q^{-h}\rangle,\qquad
Q_L(\rho,\theta;h)=V_L/F_L(\rho,\theta;h). \tag{T3}
$$
Here $h\ge0$ and the negative powers are finite integral unit
expansions. For $h\ge L$ the quotient is zero: multiplying the first
$L$ negative powers by $q^L$ gives the unit triangular polynomial
basis $1,q,\ldots,q^{L-1}$.

**Lemma 2 (shift).** For every integer $e\ge0$,
$$
q^{-e}R_L(\rho+ae,\theta-Ae)
 =F_L(\rho,\theta;\min(e,H)),\qquad
H=\max(0,\lfloor\theta/A\rfloor). \tag{T4}
$$
Indeed after multiplying, the exponent of $q$ ranges from $-e$ to
$\lfloor(\rho+dj)/a\rfloor$, with price $Bj-Ai\le\theta$.
The nonnegative exponents give $S_L$. Every selected negative exponent
$-v$ has $Av\le\theta$, and its column reduces to the axes using
$$
t^j q^{-v}=t^{j-1}q^{-(v-1)}-t^{j-1}q^{-v}\quad(j,v\ge1).
\tag{T5}
$$
Both terms have smaller price, stay within the allowed negative
exponents, and the $v=1$ first term has exponent zero and an allowed
upper bound. Conversely every indicated axis occurs at $j=0$.
All operations are integral, proving equality of spans.

## Canonical blocks become shifted regular blocks

For $k\ge0$, the canonical block of length $L$ has columns
$$
t^j q^{-i},\qquad 0\le j<L,\quad 0\le ai\le k+rj,
\quad rj-Ai\le\tau.
$$
Put
$$
J=\max(0,\min(L-1,\lfloor\tau/r\rfloor)).
$$
**Lemma 3.** Its integral cokernel is the cokernel of
$$
R_L(k+aJ,\tau-AJ). \tag{T6}
$$

**Proof.** The integral involution $t\mapsto-t/(1+t)$ takes the
canonical columns, up to signs, to $t^j q^v$ with
$-j\le v\le\lfloor(k+dj)/a\rfloor$ and $Bj-Av\le\tau$.
For $v<0$, (T5) reduces every interior column $j>-v$ to cheaper
ones, until $v=0$ or $j=-v$. The latter boundary is $(t/q)^j$,
selected exactly for $0\le j\le J$ when $\tau\ge0$.
If $\tau<0$ no negative column is selected and $J=0$, so (T6)
already follows. Assume $\tau\ge0$ in the remainder of this proof.

Multiply by $q^J$. The boundary columns form an integral unit basis
of all polynomials of degree at most $J$. The nonnegative columns
give the part of (T6) with $q$-exponent at least $J$.
Every missing cheap column $t^j q^i$, $0\le i<J$, belongs to this
span: divide it integrally by the monic polynomial $q^J$.
The remainder has degree less than $J$; the quotient is a polynomial
of degree at most $j+i-J<j$. Its terms $t^h q^J$ are available
because $Bh\le Bj\le\tau-A(J-i)<\tau$.
The upper bound for their original $q$-exponent zero is automatic.
Conversely (T6) contains every $t^h$ for $h\le J$, since
$Bh\le BJ\le\tau-AJ$, and contains the transformed nonnegative
columns. Thus the spans agree after integral row automorphisms.
$\square$

The shifted regular block need not be an original block: its
corresponding free exponent can be negative. The following global
pairing handles that boundary rather than discarding it.

## Transport between adjacent degree blocks

Fix $s,u,T$. For each $0\le m\le s$ put
$$
\rho_m=m\bmod a,\quad E_m=\lfloor m/a\rfloor,\quad
L_m=1+\lfloor(s-m)/r\rfloor,\quad
\theta_m=T-um-(1-u)\rho_m.
$$
Write $Q_m(h)=Q_{L_m}(\rho_m,\theta_m;h)$.
Set $D_m=\lfloor(\rho_m+d)/a\rfloor$.

**Lemma 4 (transport).** If $\theta_m\ge B$ and $m+r\le s$, then
$$
Q_m(h)\cong Q_{m+r}(h+D_m)\quad\hbox{for every }h\ge0. \tag{T7}
$$
If $\theta_m\ge B$ and $m+r>s$, all the $Q_m(h)$ are zero.

**Proof.** At level $j=0$ of $S_{L_m}$ the only column is $1$,
which is selected. Its positive levels are $t$ times
$$
R_{L_m-1}(\rho_m+d,\theta_m-B).
$$
The child parameters satisfy
$$
\rho_{m+r}=\rho_m+d-aD_m,\quad L_{m+r}=L_m-1,\quad
\theta_{m+r}=\theta_m-B+AD_m\ge AD_m.
$$
By (T4) this smaller span is $q^{D_m}F_{L_m-1}
(\rho_{m+r},\theta_{m+r};D_m)$.
Quotient the parent by its unit constant and divide its remaining
rows by $t$. The parent axes transform by
$$
\frac{q^{-h}-1}{t}=-\sum_{v=1}^h q^{-v}.
$$
These consecutive sums span the same axis interval as the individual
powers, by an integral triangular change. Multiplying by $q^{-D_m}$
then joins the child axes $1,\ldots,D_m$ to the new axes
$D_m+1,\ldots,D_m+h$. This proves the integral quotient isomorphism.
When $L_m=1$, the selected constant already spans the parent. $\square$

## A global cancellation, including all torsion

Use the cancellative monoid of finitely generated abelian groups
under direct sum, and its group completion. It is cancellative by
the classification theorem. Thus equality of two direct sums in
this group implies an actual isomorphism. Write
$$
\Delta_m(h)=[Q_m(h-1)]-[Q_m(h)]\quad(h\ge1).
$$
It is zero for $h>L_m$. Transport gives equality of these formal
differences along
$$
(m,k)\longmapsto(m+r,k+d),\qquad k=\rho_m+ah, \tag{T8}
$$
whenever $\theta_m\ge B$. If this step leaves the range
$h\le L_m$, the transported difference is zero.

At fixed $m$ put
$$
M_m=\min\left(E_m+1,\,1+\left\lfloor\frac{r-1-\rho_m}{a}\right\rfloor\right),
\quad J_\ell=\max(0,\min(L_m-1,
\lfloor(\theta_m-A\ell)/r\rfloor)).
$$
The normalized original regular blocks have indices
$$
\ell=0,\ldots,M_m-2,\ E_m, \tag{T9}
$$
and the canonical blocks have indices $\ell=0,\ldots,M_m-1$.
For example, at $s=17,a=2,r=3,u=9/17,T=17/3$, the $m=0$
regular block has length six and spans exactly the polynomials
through degree four. Its cokernel is $\mathbb Z$; the canonical
block has zero cokernel. Thus cancellation cannot be required
separately at each $m$.
To verify (T9), write $m=ab+k$ in a regular block: either $b=0,k=m$,
or $b>0,0\le k<d$. For a canonical block $0\le k<r$.
Each interval of $a$ consecutive integers adds exactly one $k$
of the given residue, so the two lists have the same size $M_m$.
By Lemmas 2 and 3 their cokernels are respectively
$$
Q_m(\min(\ell,H_m)),\qquad
Q_m(\min(\ell+J_\ell,H_m)),\quad
H_m=\max(0,\lfloor\theta_m/A\rfloor).
$$

Consequently the difference of the two whole direct sums is, up to
an irrelevant overall sign, $\sum_{m,h}N_m(h)\Delta_m(h)$, where
$N_m(h)=0$ if $Ah>\theta_m$, and otherwise
$$
N_m(h)=\#\{0\le\ell<M_m:\ell+J_\ell\ge h\}
-\max(0,M_m-1-h)-\mathbf1_{E_m\ge h}. \tag{T10}
$$
Only $1\le h\le L_m$ matters. Put
$$
k=\rho_m+ah,\quad b=(m-k)/a,\quad
\lambda=T-um-(1-u)k=\theta_m-Ah.
$$
The same coefficient is a sum of positive terms indexed by
$$
e\ge1,\quad e<L_m,\quad b+e\ge0,\quad
0\le k-ae<r,\quad Be\le\lambda, \tag{T11}
$$
minus one negative term exactly when
$$
b\ge0,\quad k\ge r,\quad\lambda\ge0. \tag{T12}
$$
For completeness, separate the count in (T10) into $\ell\ge h$
and $\ell<h$. In the latter part put $e=h-\ell$.
The condition on $J_\ell$ becomes $e<L_m$ and
$\theta_m-Ah\ge Be$; its other bounds are exactly (T11).
The former part differs from the regular count by
$-\mathbf1_{M_m\le h\le E_m}$, which is (T12).

Here is the exact pairing. A positive term (T11) maps to
$$
(\widetilde m,\widetilde k)=(m+re,k+de). \tag{T13}
$$
Its new $b$ is $b+e\ge0$, its slack is $\lambda-Be\ge0$,
and $\widetilde k\in[re,r(e+1))$. Thus it satisfies (T12),
and $e=\lfloor\widetilde k/r\rfloor$ recovers it uniquely.
Also $e<L_m$ gives $\widetilde m\le s$.
During the first $e$ transport steps the slack before the last
step is at least $B$, so (T7) applies at every step.
If the new flag index exceeds its row dimension, that positive
term is zero by transport. Otherwise (T13) is a negative term
with exactly the same $\Delta$.

Conversely take any negative term $(m,k)$ and put
$e=\lfloor k/r\rfloor\ge1$. Its unique predecessor is
$(m-re,k-de)$. Since $m\ge k\ge re$, its ordinary index is
nonnegative. Its $k$ is at least $ae$, so its flag index is
positive and no larger than the old one; its row dimension is
$L_m+e$. It satisfies (T11), and the same transport applies.
Every negative term is therefore paired once. This cancels the
whole sum in the group completion and proves Theorem 1. $\square$

This proves actual equality of the complete direct sums of cokernels,
not only equality of Euler characteristics of short exact sequences.
The use of the **split** direct-sum group is essential to retain torsion.

## Consequences for the original family

**Corollary 5.** For arbitrary nonzero degree-$a$ binary forms $Q,H$
over any field, with $r>a$, the complete equal-binary-price ranks of
$z^a+Q+Hz^r$ are the common value of the two bounds (B3).
They depend on $a,r,s,u,Z,T$ and the characteristic, but not on the
coefficients, subject to both forms remaining nonzero.

This follows immediately from the proved rank sandwich and (T1).
It does not assert integral coefficient universality for arbitrary
integer forms that can vanish modulo a prime. For the original
pair $Q=x^2+y^2,H=xy$, the integral decomposition (B9) additionally
transfers the complete core cokernel isomorphism to the original
and canonical seeds, including every prime-power torsion factor.

**Corollary 6.** For every $r\ge2$ and field, all positive weights
with $y$ strictly cheapest transfer from the original ideal
$$
I_{s,r}=(x^2+y^2+z^2+xyz^r)+\mathfrak m^{s+1}
$$
to its signed face $J_{s,r}=(x^2+z^2+xyz^r)+\mathfrak m^{s+1}$.
In particular, for $n=(s+1)^2$,
$$
g(I_{s,r})\le25n+14. \tag{T14}
$$

**Proof.** The middle-$z$ region and $r=2$ are already proved in
the [middle-price](middle_interaction_fans.md) and
[equal-degree](equal_degree_unit_fans.md) notes. For the remaining
region $Z>X>Y$ and $r>2$, put
$$
v=(1,1-r,1),\quad t=(X-Y)/r,\quad
u=((r-1)X+Y)/r,\quad Z_0=Z-t.
$$
Then $w=(u,u,Z_0)+tv$ and $Z_0-u=Z-X>0$.
The signed face identity $\operatorname{in}_v I=J$ is proved in
[the signed-face note](signed_face_reduction.md).
A minimum monomial basis of $J$ at $w$ is also minimum at the
seed $(u,u,Z_0)$, since all its bases have the same graded $v$-sum.
It lifts to a basis of $I$, and every original border difference
$\eta$ has $\eta\cdot v\ge0$.
Corollary 5 equates the complete seed ranks for
$(Q,H)=(x^2+y^2,xy)$ and $(x^2,xy)$. Thus this lifted basis is
minimum for $I$ at the seed too, and every nonzero border coefficient
has $\eta\cdot(u,u,Z_0)\ge0$. Their sum is positive at $w$
when $\eta\cdot v>0$; when it is zero, the coefficient is the
graded face coefficient and genericity supplies strictness.
This is the same filtered-basis argument as in the middle region,
now with its missing seed identity proved. No integral unit basis
in the characteristic-dependent face region is presumed.

The signed-grading bound $g(J)\le12n+7$ handles the whole
cheapest-$y$ region; exchange $x,y$ for cheapest $x$.
The exact cheapest-$z$ count is at most $n$: when $s\ge r+2$,
its wall denominators are at most $s$, so
$2M_{s,r}+4\le s(s-1)+4\le n$; the smaller cutoffs have
at most two sets (one for $s<2$). Every monomial initial has a
generic positive representative, so these regions cover the count.
Their sum is (T14). $\square$

The full family is therefore excluded as a source of polynomial
amplification over the generic corner-cut scale. Exceptional primes
still change some bases and ranks; they do not force a large fan.
The unrestricted $G_3$ problem, the strongest-current comparison and
the original impact/publication requirements remain unresolved.

## Prior art and reproducibility

Unit triangular Pascal identities, fractional-linear changes, finite
abelian-group classification, specialization and filtered initial
theory are established tools. The global pairing (T11)--(T13) has
a self-contained proof here; its exact priority is not certified.
This is not a new prime-distribution or $p$-fractal theorem.

[González--Planas-Vilanova, arXiv:2603.29431v1](https://arxiv.org/pdf/2603.29431),
Notation 2.1 and Section 8, Proposition 8.1 and its displayed proof,
give a reflected-transpose identity for ordinary binomial minors,
with a ratio of binomial coefficients. It is attributed there to
Gessel--Viennot. That ratio need not be a unit at a prime; it is
not substituted for the integral block transport above. The full
other determinant formulas and their application to prime ideals
have not been audited or matched to this result.

~~~bash
python3 -B scripts/run_higher_seed_transport_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker tests the span changes and integral transport separately,
replays the exact positive/negative pairing including zero boundary
terms, compares whole core cokernels through independent monomial
columns, and checks original seed ranks and signed-face transfer.
The uniform claims rest on the proofs, not an extrapolation from
the saved finite fixtures.
