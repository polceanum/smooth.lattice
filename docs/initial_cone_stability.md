# Sharp relative stability of finite-quotient initial cones

The proposed uniform relative-margin bound of order $1/n$ is false.
Actual local quotients over every field can have a cone with margin
of order $n^{-3}$ in three variables, despite having exactly seven
initial ideals. A planar version has margin of order $n^{-2}$ and
only three initial ideals.

These examples attain the universal worst-case exponents from the
classical integer cofactor argument. They close the individual-cone
packing route proposed in the [common-ray note](ray_scaling_bound.md).
They do not refute the desired sharp count or supply a large-fan
lower construction. The unrestricted $G_3$ goal remains open.

## Statement and interpretation

For an open positive initial cone
$C=\{w>0:z\cdot w>0,\ z\in Z\}$, define

$$
 \mu(C)=\sup_{w>0}\min_{z\in Z}
              \frac{z\cdot w}{|z|\cdot w}. \tag{1}
$$

The vectors are the nonzero exponent differences in a reduced
Gröbner basis, and absolute values are coordinatewise. Pure monomial
generators impose no comparison. An empty $Z$ has margin one.
Let $\Gamma_d(n)$ be the infimum of (1) over positive initial cones
of length-$n$ quotients in $d$ variables, over all fields.

**Theorem.** As $n$ tends to infinity,

$$
 \Gamma_2(n)=\Theta(n^{-2}),\qquad
 \Gamma_3(n)=\Theta(n^{-3}). \tag{2}
$$

The lower bounds hold for every field and quotient. Upper-bound
examples exist over every field at every sufficiently large length.
Their complete positive fans have only three and seven cells,
respectively.

For a fixed witness, write $P=z^+\cdot w$ and $N=z^-\cdot w$.
All independent coordinate perturbations $w_i\mapsto e^{\epsilon_i}w_i$,
$|\epsilon_i|<\rho$, preserve that comparison exactly when
$e^{-\rho}P\ge e^\rho N$. Thus the supremal multiplicative stability
radius of a cone is $\operatorname{arctanh}\mu(C)$. Orthant cones
have infinite radius. Formula (2) is therefore also the sharp
worst-case stability-radius order.

## The classical uniform lower certificate

For a length-$n$ standard set, every standard monomial has total
degree at most $n-1$. A minimal initial generator has degree at
most $n$, because a chain of its proper divisors lies in the
$n$-element standard set. Therefore every comparison vector has
coordinates bounded in absolute value by $n$.

In dimension three, primitive extreme rays of the closed positive
cone have coordinates at most $2n^2$, by two-by-two cofactors of
two active comparison or coordinate inequalities. The cone is
pointed and full dimensional. Sum three independent primitive rays.
Every nonzero defining inequality is strict on their sum, since
it is nonnegative on each ray and cannot vanish on three independent
ones. This gives an integer witness with

$$
 1\le w_i\le6n^2,\qquad z\cdot w\ge1.
$$

Hence $\mu(C)\ge1/(18n^3)$. In dimension two, two primitive rays
have coordinates at most $n$; their sum gives
$\mu(C)\ge1/(4n^2)$. These are classical cofactor bounds applied
to relative slack, not a new integer-cone theorem. The earlier
[weight-height note](initial_cone_height.md) records the positive
weight and cofactor comparisons.

## A socle construction that is an actual quotient

Let $D_0$ be a finite monomial downset with maximal monomials $T$,
and let $J$ be the monomial ideal with standard set $D_0$.
The monomials of $T$ span the socle of $K[x]/J$: every variable
kills them. Consequently any linear subspace of their span is an
ideal of that algebra.

In particular, for disjoint pairs $u_i,v_i\in T$, set

$$
 I=J+(u_i-v_i:1\le i\le m). \tag{3}
$$

Its length is $|D_0|-m$. All lower monomials remain independent;
each pair gives one independent residue class. At positive prices
choose the cheaper representative in each pair. The standard set
is $D_0$ with the more expensive representative of each pair removed.
It remains a downset because only maximal elements are removed.

For any such strict choice, the reduced basis consists of the
oriented pair differences and the other minimal excluded monomials
of that standard set. Its tails are standard. Products of a pair
relation by any variable vanish modulo $J$, so there are no hidden
relations among the lower monomials. This construction gives
commuting multiplication on a cyclic quotient in every characteristic.
It is stronger evidence than an unattached representable matroid.

## The planar family

Fix $k\ge3$. Use the four maximal monomials

$$
 u_1=x^{k+2},\quad v_1=x^2y^{k-1},\qquad
 u_2=y^{k+1},\quad v_2=x^{k+1}y.
$$

Their divisor downset has $5k+1$ elements. The quotient (3)
has length $n=5k-1$. The choice keeping $v_1,v_2$ has exact cone

$$
 k s>(k-1)t,\qquad k t>(k+1)s. \tag{4}
$$

Its slack matrix and inverse are

$$
 M_2=\begin{pmatrix}k&1-k\\-k-1&k\end{pmatrix},\qquad
 M_2^{-1}=\begin{pmatrix}k&k-1\\k+1&k\end{pmatrix}. \tag{5}
$$

Both determinants are one, and the inverse is strictly positive.
Every sign pattern of $M_2w$ except $(-,-)$ is realized at positive
weights: give one positive component a large magnitude and negative
components sufficiently small magnitudes before applying the inverse.
The all-negative pattern cannot give a positive weight. Thus there
are exactly three initial ideals.

Writing $q=s/t$, the interval in (4) is
$((k-1)/k,k/(k+1))$. Equating its two relative slacks gives

$$
 \mu(C_k)=\frac1{(k+\sqrt{k^2-1})^2}
          =\Theta(k^{-2}). \tag{6}
$$

The logarithmic interval width is
$\log(k^2/(k^2-1))=\Theta(k^{-2})$. This already rules out a
uniform $1/n$ margin even in a dimension where the sharp global
count is known in the repository.

## The three-variable family

Fix $k\ge3$. Use six maximal monomials

$$
\begin{array}{lll}
 u_1=x^{k+2},&u_2=y^k,&u_3=z^{4k+3},\\
 v_1=x^2y^{k-1}z,&v_2=x^{k+1}z,&v_3=xyz^2.
\end{array} \tag{7}
$$

They form an antichain. Their divisor downset has $12k+5$ elements:
its layers at $z$ exponents $0,1,2$, and $3,\ldots,4k+3$ have
$4k+1$, $4k-1$, $4$, and $4k+1$ elements in total.
Thus (3) has length

$$
 n=12k+2. \tag{8}
$$

The exact reduced basis for the standard set omitting $u_1,u_2,u_3$ is

$$
\begin{gathered}
 x^{k+2}-x^2y^{k-1}z,\qquad y^k-x^{k+1}z,\qquad
 z^{4k+3}-xyz^2,\\
 yz^3,\quad y^2z^2,\quad xz^3,\quad x^2z^2,\quad x^3y .
\end{gathered} \tag{9}
$$

The full cone is $C_k=\{w:M_kw>0\}$, where

$$
 M_k=\begin{pmatrix}
 k&1-k&-1\\
 -k-1&k&-1\\
 -1&-1&4k+1
 \end{pmatrix},\quad \det M_k=1, \tag{10}
$$

$$
 V_k=M_k^{-1}=
 \begin{pmatrix}
 4k^2+k-1&4k^2-3k&2k-1\\
 4k^2+5k+2&4k^2+k-1&2k+1\\
 2k+1&2k-1&1
 \end{pmatrix}. \tag{11}
$$

All inverse entries are positive, so $C_k$ lies wholly inside the
positive orthant. As in the planar case, all seven sign patterns
except $(-,-,-)$ are realizable at positive weights. They give
different standard sets, proving the exact count seven.

The integer witness

$$
 w=V_k(1,1,1)^T=(8k^2-2,\,8k^2+8k+2,\,4k+1)
$$

has all three slacks equal to one. Its largest absolute comparison
budget is $16k^3+16k^2+4k-1$, yielding the lower bound below.

For the upper bound normalize $w_z=1$ and write $w=(s,t,1)$.
The positive inverse columns show $s>2k-1$, $t>2k+1$.
The positive row vector $\lambda=(2k+1,2k-1,1)$ satisfies

$$
 \lambda M_k=(0,0,1),\qquad
 \lambda |M_k|=(4k^2+2k,\,4k^2-2k,\,8k+1).
$$

If all relative slacks were at least $\mu$, multiplying
$M_kw\ge\mu|M_k|w$ by $\lambda$ would give
$1>\mu(16k^3+4k+1)$. Therefore

$$
 \frac1{16k^3+16k^2+4k-1}
 \le\mu(C_k)\le\frac1{16k^3+4k+1}
 =\Theta(k^{-3}). \tag{12}
$$

The closed triangle in the chart $w_z=1$ has vertices
$V_ke_j/(V_ke_j)_z$. Its area in the $(s,t)$ coordinates is
$1/[2(4k^2-1)]$. On it, $s\ge2k-1$, $t\ge2k+1$ and
$s+t\le4k+1$. The Jacobian of $(s,t)\mapsto(\log s,\log t)$
is $1/(st)$, so its logarithmic area satisfies

$$
 \frac2{(4k^2-1)(4k+1)^2}
 \le \operatorname{area}_{\log}(C_k)
 \le\frac1{2(4k^2-1)^2}
 =\Theta(k^{-4}). \tag{13}
$$

Thus an individual logarithmic area floor of order $1/n^2$ also
fails for actual quotients.

### Every sufficiently large length

The exponents in (2) do not rely on an arithmetic subsequence.
In the planar family, for $k\ge7$ add the $q$ monomials
$x^3y^2,\ldots,x^3y^{q+1}$ to $D_0$, where $0\le q\le4$.
In three variables, for $k\ge12$ add
$x^3y,\ldots,x^3y^q$, where $0\le q\le11$.
These are downset extensions; none has an original $u_i$ or $v_i$
as a divisor. The old paired elements remain in the socle.

The quotient length increases by $q$, while the added residue
vectors are forced independent elements. They translate every basis
profile by the same vector and change no initial cone. This covers
all lengths at least $34$ in dimension two and $146$ in dimension
three. Together with the cofactor lower bound, it proves (2).

## A whole-class bound: socle-only relations have small fans

There is also a useful limit on this realization method.
Let $D_0$ be any finite monomial downset in three variables,
$T$ its maximal elements, and $R$ any linear subspace of their
span. For $I=J+R$, write $Q=D_0\setminus T$ and $q=|Q|$.
Then

$$
 g(I)=O(1+q^{3/2})=O(n^{3/2}),\qquad
 n=\dim_K K[x,y,z]/I>0. \tag{14}
$$

All $Q$ residue vectors are independent and forced in every basis.
All remaining choices form a matroid on $T$. For a fixed coordinate
$i$, the proper predecessor chains
$a-e_i,\ldots,a-a_ie_i$ for $a\in T$ lie in $Q$ and are disjoint:
chains could meet only if two maximal elements agreed in the other
coordinates, making them comparable. Consequently

$$
 \sum_{a\in T}a_i\le q. \tag{15}
$$

After removing the fixed core translation, every basis profile lies
in $[0,q]^3$. Every strict initial state injects into a vertex of its
integer profile polytope. The classical three-dimensional lattice-vertex
bound gives $O(q^{3/2})$ vertices when the polytope is full dimensional.
For a flat hull, its area is $O(q^2)$ and its induced integer lattice
has covolume at least one; the planar bound gives $O(q^{2/3})$.
Segments have at most two vertices. Since $n\ge q$, (14) follows.

The same argument in fixed dimension $d$ gives
$O_d(1+n^{d(d-1)/(d+1)})$. In particular, arbitrary linear socle
identifications cannot produce the requested polynomial separation
above the generic scale. The hard cones above are geometric
conditioning examples in a globally small class.

## Mathematical inputs, verification and the next gap

The antichain/socle construction is established background:
[Agnarsson--Epstein, Section 2 and Proposition 5.2](https://arxiv.org/pdf/1801.02644v2).
The monomial-socle description and the zero-dimensional antichain
statement with its preceding argument were inspected.
The vertex exponents in (14) are classical:
[Andrews](https://doi.org/10.1090/S0002-9947-1963-0143105-7);
[Dillon, Theorem 1 and its relative-lattice formulation](https://arxiv.org/pdf/2606.30856).
Neither input is a new theorem here. The precise priority of the
explicit stability families and the combined sharp stability statement
remains unestablished.

The checker verifies polynomial-coefficient matrix inverse identities,
integer commuting cyclic modules, ordinary S-pairs and normal forms
over the rationals and several prime fields. Complete cone coverage
is checked separately from the module construction; selected Gaussian
bases provide another independent check. Padding and disjoint-chain
budgets are also replayed. Finite agreement supports the displayed
uniform proofs; it does not establish priority.

The original counting target and its acceptance criteria are unchanged.
The relative-margin implication remains a correct conditional
statement, but its universal premise is refuted. Even downsets,
commuting multiplication and all-field realizability do not repair it.
A successful sharp count must bound the aggregate contribution of
thin cones, exploit their restricted profile variation, or use a
different global mechanism. It cannot assume a uniform individual
margin or area. No new contribution manuscript is warranted from
this obstruction alone.

~~~bash
python3 -B scripts/run_socle_stability_checks.py
python3 -B scripts/run_ray_scaling_checks.py
./build.sh
python3 tests/test_smoke.py
~~~
