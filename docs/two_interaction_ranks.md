# Two interactions and exact prefix ranks

The single-interaction transport theorem in
[higher_seed_transport.md](higher_seed_transport.md) does **not** extend
by replacing two independent binary forms with one repeated form. The
failure persists uniformly over all fields and arbitrarily large
truncations. This note proves an exact rank formula explaining the
failure. It does not prove a larger number of initial ideals.

## The quotient and its columns

Let $K$ be any field, let $Q,H,L\in K[x,y]$ be nonzero homogeneous
quadratics, and let $r=2k+1\ge3$. Set

\[
 s=2r+1,\qquad
 F=z^2+Q+Hz^r+Lz^{r+1},\qquad
 A=K[x,y,z]/(F,(x,y,z)^{s+1}).                         \tag{C1}
\]

The classes of $x^a y^b z^\epsilon$, with $\epsilon\in\{0,1\}$
and $a+b+\epsilon\le s$, form a basis. In particular,

\[
 n=\dim_K A=(s+1)^2=4(r+1)^2.                          \tag{C2}
\]

Here is a direct justification, including possible bad characteristics.
Use $z^2=-Q-Hz^r-Lz^{r+1}$: a low substitution decreases the
$z$ exponent without changing ordinary degree; a high substitution
increases ordinary degree by at least $r$. Thus reduction terminates
modulo the degree cutoff and spans the asserted set. The lowest
homogeneous part of $F$ is the monic quadratic $z^2+Q$. For any
nonzero polynomial multiplier, the lowest homogeneous part of its
product with $F$ is its lowest part times $z^2+Q$, since the polynomial
ring is a domain. The associated graded quotient therefore has degree
$j$ dimension $2j+1$ for $0\le j\le s$. Summing gives $(s+1)^2$,
exactly the size of the spanning set, proving independence.

Write $q=-Q$. For $j\ge1$, the exact normal forms in this basis are

\[
\begin{aligned}
 [z^{2j}]&=q^j-jHq^{j+k-1}z-jLq^{j+k},\\
 [z^{2j+1}]&=q^jz-jHq^{j+k}-jLq^{j+k}z.                \tag{C3}
\end{aligned}
\]

Multiply these expressions by $x^a y^b$ and discard terms of ordinary
degree greater than $s$. A second high substitution has degree at
least $2+2r=s+1$, so no terms containing two high substitutions can
survive. Expanding $(z^2)^j$ with at most one high substitution proves
(C3). No division by $j$, $2$, or a field element occurs. The cases
$z^0,z^1$ are their own normal forms.

## An exact rank formula

Take the complete monomial prefix for prices $(1,1,2)$ at threshold
$T=r+3$:

\[
 E=\{x^a y^b z^c:a+b+2c\le r+3\}.
\]

All its monomials have ordinary degree at most $s$. Then

\[
 \boxed{\dim_K\langle [E]\rangle
       =2r^2+6r+13-\delta,}                            \tag{C4}
\]

where

\[
\begin{array}{c|c}
 r & \delta=1\text{ precisely when}\\ \hline
 3 & H\mid QL\\
 r\ge5\text{ odd} & H^2\mid QL^2 .
\end{array}                                                     \tag{C5}
\]

Otherwise $\delta=0$. Divisibility is in $K[x,y]$, in the actual
characteristic under consideration.

**Proof.** Let $S_d=K[x,y]_d$, and set $S_d=0$ for $d<0$.
The reference columns already in the prefix form an independent
space

\[
 R=K[x,y]_{\le r+3}\;\oplus\;zK[x,y]_{\le r+1},\quad
 \dim R=\binom{r+5}{2}+\binom{r+3}{2}=r^2+7r+13.       \tag{C6}
\]

Equation (C3) gives

\[
 [z^{2j}]=j q^{j-1}[z^2]+(1-j)q^j,\quad
 [z^{2j+1}]=j q^{j-1}[z^3]+(1-j)q^jz.                \tag{C7}
\]

If the left parent monomial belongs to $E$, every binary monomial
multiplier of the indicated $z^2$ or $z^3$ on the right also belongs
to $E$; the remaining reference terms lie in $R$. Consequently the
prefix image modulo $R$ is generated exactly by

\[
 z^2\bigoplus_{d=1}^{r-1}S_d
 \quad\text{and}\quad
 z^3\bigoplus_{e=0}^{r-3}S_e.                          \tag{C8}
\]

The scalar $z^2$ column has already landed in $R$. The domain in
(C8) has dimension $r^2-r$. Write a vector in it as
$z^2\sum A_d+z^3\sum B_e$. Grouping its normal form outside $R$
by ordinary binary degree and by $z$ exponent, cancelling only
nonzero powers of $Q$, gives exactly

\[
\begin{aligned}
 H A_d-QL B_{d-2}&=0 &&(1\le d\le r-1),\\
 L A_d+H B_d&=0 &&(1\le d\le r-2).                   \tag{C9}
\end{aligned}
\]

An out-of-range $B$ is zero. For example, the $L A_{r-1}$ term
is killed by the ordinary degree cutoff; the $H B_0$ term lies
inside $R$. These boundary terms explain both endpoints of (C9).
The signs in (C9) follow from $q=-Q$ and are valid in characteristic
two as well.

Every odd-indexed $A,B$ is zero, starting with $HA_1=0$.
The even-indexed chain is completely determined in the fraction
field by the scalar $B_0$. If $B_0=0$, the whole kernel is zero.
For $r=3$, only $HA_2=QLB_0$ remains, giving the first line of
(C5). For $r\ge5$, the equations for $A_2,B_2$ require
$H^2B_2=-QL^2B_0$. Thus a nonzero kernel requires $H^2\mid QL^2$.

Conversely, this divisibility also implies $H\mid QL$: compare
the multiplicity of each irreducible factor of $H$ in $Q,L$.
Both

\[
 D=QL/H,\qquad C=QL^2/H^2
\]

are therefore quadratics. The assignments

\[
 B_{2j}=(-C)^j B_0\quad(0\le j\le k-1),\qquad
 A_{2j}=D(-C)^{j-1}B_0\quad(1\le j\le k)             \tag{C10}
\]

solve every equation, including the zero odd boundary at $r-2$.
The kernel has dimension exactly one. Subtracting this dimension
from (C8) and adding (C6) proves (C4). $\square$

This calculation gives the full dependency, not only a determinant
test. When $\delta=1$, put $P=\sum_{j=0}^{k-1}(-C)^j$ for
$r\ge5$, and $P=1$ for $r=3$. The extra relation modulo $R$ is

\[
 z^2(D+z)P\pmod R.                                  \tag{C11}
\]

Its exact normal form lies in $R$; subtracting that normal form
gives a polynomial supported entirely in $E$ that vanishes in $A$.

## What changes, and what remains open

For every odd $r\ge3$, compare the two actual cyclic quotients

\[
\begin{aligned}
 F_{\rm independent}&=z^2+x^2+y^2+xy z^r+x^2z^{r+1},\\
 F_{\rm repeated}&=z^2+x^2+y^2+xy z^r+xy z^{r+1}.
\end{aligned}
\]

The first has $\delta=0$ and the second $\delta=1$ over **every**
field. Their lengths and selected monomial set are equal, while
their complete-prefix ranks differ by one. A rank identity that
replaces both high forms by a repeated form is therefore false
uniformly, rather than just at an exceptional small cutoff.

The first equation also avoids the obvious large-index grading
parameter. Four support differences are

\[
 (2,-2,0),\ (2,0,-2),\ (1,1,r-2),\ (2,0,r-1).
\]

Their four maximal minors are $4r,4(r+1),4,2(r+1)$.
For odd $r$ the support-difference lattice has rank three and index
four. Thus the displayed equation has no nonzero integer-weight
grading, and increasing $r$ does not increase this lattice index.
This is a statement about its support lattice; no identification
with the intrinsic profile lattice is needed or asserted here.

There is also a genuine characteristic-sensitive distinction:
$Q=H=x^2+y^2$, $L=x^2+xy$ has $\delta=1$ for $r=3$ over every
field, but for odd $r\ge5$ it has $\delta=1$ exactly in
characteristic two. In characteristic two, $Q=(x+y)^2$ and
$L=x(x+y)$; in every other characteristic $x+y$ and $x$ are
coprime to $x^2+y^2$ over an algebraic closure. The distinction
between the two divisibility tests in (C5) is necessary.

The first two rank profiles stay different for generic positive
weights sufficiently close to $(1,1,2)$: there is a unit price gap
at the prefix boundary, so a sufficiently small perturbation
preserves this finite prefix. Greedy rank selection then gives
different counts of selected standard monomials in that prefix.
This proves a difference of the realized standard sets. It does
**not** compare their total numbers of initial ideals, or prove
that the dependence produces many simultaneous cones.

The relevant next question is quantitative: can coherent equations
of this kind force more than the generic quadratic-in-length
scale of distinct initial ideals? A single lost relation, even in
an infinite family, cannot answer it. The present exact columns
and coupled recurrence provide a mechanism to analyze; arbitrary
matroid examples cannot replace a realizable quotient or its
length budget. No unrestricted $G_3$ bound changes here.

## Prior art and verification

The monomial-column matroid and its relation to initial staircases
are established in
[Babson--Onn--Thomas, Sections 2--3](https://arxiv.org/pdf/math/0207135).
The truncation and divisibility arguments above are elementary;
their exact priority has not been established. The comparison
with the proved single-interaction identity is internal and exact.
It is not evidence of novelty or publication significance. The
unrestricted objective and the unresolved strongest-current
literature comparison remain in force.

~~~bash
python3 -B scripts/run_two_interaction_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker compares terminating recurrence columns with (C3),
checks the independent Macaulay relation kernel, verifies (C4)
against direct selected columns and the smaller coupled system,
and replays the explicit extra relations. It includes the changed
divisibility criterion and characteristic-two example. Finite
agreement supports the implementation; the uniform assertions
rest on the proof above. No fan count is extrapolated from ranks.
