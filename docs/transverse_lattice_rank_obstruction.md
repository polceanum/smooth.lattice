# Transverse head filters require growing polynomial-module rank

The scalar bounded-shift lattice argument does not extend with constant
module rank to all transverse operators. An explicit compatible hook
has \(\lceil t/2\rceil\) separated positive head-birth costs at one
weight, although its entire fan has an \(O(n)\) upper bound.
For a common cost-preserving polynomial-module representation with
constraints modulo \(q^t\), the necessary rank is \(\Omega(t)\);
the natural rank-\(2t\) representation is optimal in order.

This is a uniform obstruction to one proposed compression mechanism,
not a lower bound for all fan algorithms. The costs are actual one-arm
pure-prefix births; their individual occurrence as coherent cross-arm
walls is not asserted. The family is already covered by the
[two-monomial theorem](hook_two_orbit_fans.md), so it supplies no
polynomial excess or new unrestricted fan-count exponent.

## A compatible family with short transverse steps

Work in characteristic zero, or characteristic \(p>2t\), and put
\[
 t\ge2,\qquad D=t^2,\qquad A=(t-1)D+t,\qquad a=A+t,
 \qquad P(T)=(1+T)^{t-1}.
                                                               \tag{TR1}
\]
Use zero outer images and
\[
 I_t=(x^2,y^a,z^a,
        y^tz^t-x(y^tP(y^Dz)+z^tP(z^Dy))).
                                                               \tag{TR2}
\]
Its length is
\[
             n=4At+2t^2=4t^4-4t^3+6t^2=\Theta(t^4).           \tag{TR3}
\]
Both operators are units with constant axis coefficient one.
The transverse dependence has degrees up to \(t-1\).
On each arm the step is \((D,1)\), with \(D\ge2\).
The saved steep two-monomial theorem gives, uniformly in \(t\),
\[
                              g(I_t)=O(t^4)=O(n).              \tag{TR4}
\]
This uses its full arrangement bound, not a count extrapolated from
the finite tests below.

In the Y arm \(R=K[q,v]/(q^A,v^t)\), pure columns are
\[
 U_{r,k}=q^rv^k,\quad 0\le r<A,\qquad
 L_{r,k}=q^rv^kP(q^Dv),\quad t\le r<A.
\]
Give them normalized prices
\(\theta+r+k\rho\) and \(r+k\rho\), respectively.
The heads are \(h_{i,j}=q^iv^jP(q^Dv)\), \(i,j<t\).
Their actual source-cost order is \(i+j\rho\).
For the head \(h_{0,j}\), let \(\alpha_j\) be the first pure-prefix
price spanning it modulo the earlier heads.

These prices occur at the positive original weight
\[
                 w=(\theta+t\rho,1,\rho).
                                                               \tag{TR5}
\]
The common original affine cost \(t+t\rho\) has been subtracted.

## Exact chain-prefix ranks

Multiplication by \(P(q^Dv)\) preserves \(r-Dk\).
Since \(D>t\), the \(t^2\) head inputs lie on distinct chains.
Consequently earlier heads on other chains cannot help span
\(h_{0,j}\). Its chain has rows \(e_0,\ldots,e_m\), where
\(m=t-1-j\); the long cutoff in (TR1) does not shorten it.
In chain coordinates, with \(c_r=\binom{t-1}{r}\),
\[
 U_\ell=e_\ell,\quad 0\le\ell\le m,\qquad
 L_\ell=\sum_{r=0}^{m-\ell}c_r e_{\ell+r},
                         \quad1\le\ell\le m,\qquad
 h=\sum_{r=0}^m c_r e_r.                                      \tag{TR6}
\]
Write \(V=D+\rho\). After subtracting \(j\rho\), the prices are
\(\theta+\ell V\) and \(\ell V\).

For upper and lower prefixes of sizes \(u,s\), respectively,
\[
\begin{split}
 \operatorname{rank}(U_0,\ldots,U_{u-1},L_1,\ldots,L_s)
                          &=\min(m+1,u+s),\\
 \operatorname{rank}(U_0,\ldots,U_{u-1},L_1,\ldots,L_s,h)
                          &=\min(m+1,u+s+1).
\end{split}                                                     \tag{TR7}
\]
For \(u=0\) the assertions follow from triangularity and \(c_0=1\).
For \(u\ge1\), eliminate the first \(u\) upper rows. The necessary
square minors are consecutive binomial Toeplitz minors, with
differences \(u-1\) for the pure sources and \(u\) after adding
the head. Their nonzero product evaluation is classical:
\[
 \det\!\left[\binom{t-1}{k+i-j}\right]_{0\le i,j<r}
       =\prod_{\ell=0}^{r-1}
           \frac{(t-1+\ell)!\,\ell!}
                {(k+\ell)!\,(t-1-k+\ell)!}.                    \tag{TR8}
\]
Here \(0\le k\le t-1\) and the minors used in (TR7) stay within
the indicated chain. All factorial arguments are at most \(2t-2\).
Thus they are nonzero in characteristic zero or \(p>2t\).
See [Krattenthaler, Section 2.3, equations (2.16)--(2.18)](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11).
No new determinant evaluation is asserted.

The head is therefore born when exactly \(m+1\) pure columns have
appeared. Equivalently its birth, minus \(j\rho\), is the
\((m+1)\)-st smallest member of
\[
 \{\theta+\ell V:0\le\ell\le m\}
                  \ \sqcup\ \{\ell V:1\le\ell\le m\}.           \tag{TR9}
\]
For \(0<\theta<V\), the streams interleave as
\(U_0,L_1,U_1,L_2,U_2,\ldots\). Hence
\[
 \boxed{\quad
 \alpha_j=j\rho+
 \begin{cases}
    \theta+(m/2)(D+\rho),&m=t-1-j\text{ even},\\
    ((m+1)/2)(D+\rho),&m=t-1-j\text{ odd}.
 \end{cases}\quad}                                            \tag{TR10}
\]
Every normalized gain \(\alpha_j-j\rho\) is positive.
The source, pure-prefix and earlier-head definitions are retained
throughout; no rank is inferred from the support alone.

## A sharp-order lower bound for a common module representation

Choose
\[
                 \rho=\frac{12t}{12t+1},\qquad\theta=\frac13 .
                                                               \tag{TR11}
\]
Earlier heads with positive long input lie on other chains and cannot
change these births. For the even values \(m=2r\), (TR10) gives
\[
 \alpha_{t-1-2r}=\theta+rD+(t-1-r)\rho,
       \qquad0\le r\le\left\lfloor\frac{t-1}{2}\right\rfloor .
                                                               \tag{TR12}
\]
These \(\lceil t/2\rceil\) costs are separated by \(D-\rho>t-1\).
Their normalized gains are separated by \(D+\rho\).
In particular no fixed number of intervals of width \(Ct\), for
any fixed \(C\), covers them uniformly as \(t\) grows.
Both arm shifts are positive at this weight: the Z shift is
\((\theta-t/(12t+1))/\rho>1/4\). Hence both arms have positive
head gains; the separation is not obtained by making the other arm
inactive.

Suppose all these exact minimum-cost head problems are represented
by one polynomial module with a common shifted maximum degree,
a rank-\(r\) predictable-degree basis and affine constraints depending
only on the vector modulo \(q^t\). The representation preserves the
displayed costs. If the basis row degrees are \(d_1,\ldots,d_r\),
the [short-multiplier argument](hook_scalar_outer_lattice_fans.md)
at zero shift perturbation gives
\[
                    \alpha_j\in\bigcup_{\ell=1}^r
                                        [d_\ell,d_\ell+t-1].   \tag{TR13}
\]
Indeed reducing each basis multiplier modulo \(q^t\) preserves
the constraint without increasing predictable degree.
Each interval contains at most one cost from (TR12), so
\[
                    r\ge\lceil t/2\rceil=\Omega(n^{1/4}).        \tag{TR14}
\]
The same argument applies to an exact common-degree representation
of the normalized gains. Arbitrary head-specific cost translations,
changing the module with the head, or a larger residue cutoff are
different representation models; (TR14) does not bound them.
Allowing intervals of width \(O(t)\) in a bounded-shift version
still requires \(\Omega(t)\) intervals because \(D=t^2\).

The natural representation has rank \(2t\). Regard multiplication
by \(P(q^Dv)\) as the \(t\)-by-\(t\) polynomial matrix \(W(q)\),
and take
\[
 \mathcal L=\{(u,b)\in K[q]^{2t}:u-W(q)b=0\bmod q^A\}.
                                                               \tag{TR15}
\]
Its basis consists of \((q^Ae_j,0)\) and \((W(q)e_j,e_j)\).
Give component \(u_j\) shift \(\theta+j\rho\), and component
\(b_j\) shift \(j\rho\). The head constraint prescribes \(b\)
modulo \(q^t\). For the heads used above its degree minimum is
exactly \(\alpha_j\); all earlier source costs are below the
positive birth. Thus (TR14) is tight in order for the common
polynomial-module representation of the birth costs.

This separates representation rank from actual fan complexity:
linear fan upper bounds can coexist with a polynomially growing
necessary module rank. A transverse proof needs to amortize the
interactions of many row-cost bands, or use a representation outside
this model. Increasing rank alone is not evidence of many initial
ideals.

## Verification and limits

~~~bash
python3 -B scripts/run_transverse_lattice_rank_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker tests all chain-prefix ranks in (TR7), both sides of the
birth threshold in the complete original rectangular arm, and the
independent full vector-lattice representation (TR15). It checks
separated costs and their row-band cover directly. A chain-order
profile oracle is compared with complete original quotient fans.
The latter checks include commuting multiplication and original
coefficient columns.
For that oracle the head grid points have distinct cosets modulo
\(\mathbb Z(D,1)+\mathbb Z(1,D)\): a difference with both coordinates
of absolute value below \(t\) gives integral lattice coordinates of
absolute value below one, since \(D>t\). Thus every marked component
has one coupled head. Its corank-one birth comparison and the chain
prefix ranks give a field-free profile query; the checker compares
the full resulting polygons with original-field enumeration.

The characteristic restriction matters: at \(t=3\) over
\(\mathbb F_2\), \(P(T)=1+T^2\). The \(m=2\) chain is born at
\(2V\), whereas (TR10) predicts \(\theta+V\). This exact scope
guard is included.

The binomial minors are classical, as are predictable degree and
weak Popov reduction; compare
[Labahn--Neiger--Vu--Zhou, Sections 2.5--2.6](https://arxiv.org/pdf/2202.09329).
The global fan bound (TR4) was already proved for the larger
single-monomial-per-arm class. This note supplies no new extremal
fan exponent, infinite polynomial lower excess or strongest-current
novelty certificate. The unrestricted goal and publication assessment
remain incomplete.
