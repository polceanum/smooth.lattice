# Valuation-sensitive bounds for every constant-axis central hook

The first transverse polynomial need not be a unit. Factoring out its
long-variable valuation extends the transverse-prefix theorem to every
central-only hook with constant axis units. The resulting parameter separates
two effects: gaps in a unit Toeplitz table, and a valuation beyond half the
available module length. Nonzero valuation by itself need not weaken the
quadratic polylogarithmic fan bound.

This is a restricted-family theorem. It does not settle unrestricted
\(G_3(n)\), construct a polynomial lower separation, or establish novelty.
The outer images of \(y^a,z^a\) remain zero.

## A sharp conditional uncertainty bound

Work over any field in
\[
 R=K[q,v]/(q^m,v^t),\qquad
 W=\gamma+v^dF(q)+\sum_{j>d}v^jW_j(q),\quad
 \gamma\ne0,\quad 1\le d<t.                              \tag{TV1}
\]
The coefficients above order \(d\) are arbitrary. Let \(\nu\) be the
valuation of \(F\) modulo \(q^m\), with \(\nu=m\) when it is zero.
If \(\nu<m\), put \(s=m-\nu\) and
\(G=q^{-\nu}F\bmod q^s\). This \(G\) is a unit. Let \(h_s(G)\)
be the longest consecutive gap in the regular Toeplitz indices
\[
 \det[G_{r+i-j}]_{0\le i,j<s-r}\ne0,\qquad0\le r\le s.
\]
Both endpoints are regular. Define
\[
 e_m(F)=
 \begin{cases}
 \max\{h_s(G),\,2\nu-m-1\},&\nu<m,\\
 m-1,&\nu=m,
 \end{cases}
 \qquad \ell_m(F)=m-1-e_m(F).                             \tag{TV2}
\]

**Lemma.** Every \(H\ne0\) with transverse valuation less than \(t-d\)
satisfies
\[
                 \deg_qH+\deg_q(WH)\ge\ell_m(F).          \tag{TV3}
\]
For transverse valuation at least \(t-d\), \(WH=\gamma H\);
these scalar-tail exceptions are unchanged.

To prove the lemma, let \(H_r\ne0\) be its first transverse coefficient.
As in the [unit transverse argument](hook_transverse_prefix.md),
\[
 (WH)_r=\gamma H_r,\qquad
 FH_r=(WH)_{r+d}-\gamma H_{r+d}.                          \tag{TV4}
\]
Write \(D=\deg_qH\), \(E=\deg_q(WH)\). Then
\[
 \deg H_r\le\min(D,E),\qquad
 \deg(FH_r)\le\max(D,E)
\]
when the latter product is nonzero.

If \(FH_r\ne0\), put \(h=H_r\bmod q^s\ne0\). Unit uncertainty
at length \(s\) gives
\[
 \deg h+\deg(Gh\bmod q^s)\ge s-1-h_s(G).
\]
Since \(FH_r=q^\nu(Gh\bmod q^s)\),
\[
                         D+E\ge m-1-h_s(G).              \tag{TV5}
\]
If \(FH_r=0\), the exact kernel of multiplication by \(F\) is
\((q^s)/(q^m)\), so \(\deg H_r\ge s\). Therefore
\[
                              D+E\ge2s.                 \tag{TV6}
\]
The minimum of the two bounds is exactly (TV2). If \(F=0\) in this
module, (TV3) reduces to the valid nonnegative lower bound zero.
No nonzero product is assumed for a nonunit multiplier.

The lower bound is sharp for an explicit infinite family, even with
arbitrary fixed higher coefficients \(W_j\). In characteristic zero or
\(p>2m\), take
\[
             F=q^\nu(1+q)^K,\qquad0\le K<m-\nu.           \tag{TV7}
\]
The classical binomial determinant product gives
\(h_s(G)=s-K-1\) when this is positive, and zero otherwise. Hence
\[
                 \ell_m(F)=\min(\nu+K,\,2(m-\nu)).        \tag{TV8}
\]
The vector \(H=v^{t-d-1}\) attains the first value. When the second
is smaller, \(H=q^{m-\nu}v^{t-d-1}\) is nonzero and attains it:
\(FH_r=0\). Every higher transverse term vanishes on both witnesses.
For \(F=0\), the first witness gives equality zero. This proves
sharpness of the conditional degree lemma, not of the fan-count bound.
The determinant evaluation is classical; see
[Krattenthaler, Section 2.3, (2.16)--(2.18)](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11).

## Apply the lemma to all positive fan directions

Take
\[
 I=(x^2,y^a,z^a,y^tz^t-xf),\quad
 t\ge2,\quad a\ge2t,\quad n=4at-2t^2,\quad
 A=a-t,\quad L=a-2t.                                    \tag{TV9}
\]
Write \(f=y^tU+z^tV\) in the two rectangular arm modules. The axis
units \(U(y,0)=\alpha\), \(V(0,z)=\beta\) are nonzero constants.
For each nonscalar arm use its first nonzero transverse coefficient,
at any order \(1\le d<t\), and call its polynomial \(F_Y,F_Z\).
There is no condition on their constant coefficients or valuations.
Set
\[
 E_Y=\max_{m\in\{A,L\},\,m>0}e_m(F_Y),\qquad
 E_Z=\max_{m\in\{A,L\},\,m>0}e_m(F_Z),\qquad
                         T=t+E_Y+E_Z.                  \tag{TV10}
\]

**Theorem.** If both arms are nonscalar, uniformly over all these
coefficients and all fields,
\[
 \begin{split}
 D_{\rm all}(I)&=O(aT(1+\log(a+1))),\\
 g(I)&=O(naT(1+\log(n+1))^2)\\
     &=O\!\left(n^2\left(1+\frac{E_Y+E_Z}{t}\right)
                         (1+\log(n+1))^2\right).
 \end{split}                                             \tag{TV11}
\]
If either arm is scalar, the previous direct-sum proof gives the same
conclusions with \(T=t\), independently of the other arm.
The earlier box ceiling is retained whenever smaller.

Here is why the entire earlier fan argument still applies. All multipliers
\(U,V\) and their length-\(L\) truncations remain units, because their
axis constants are nonzero. Their first transverse coefficients need not
be units. The cheap-overlap argument uses only those axis constants.
The primal and separating-functional arguments exclude vectors in
transverse orders at least \(t-d\), using only the absence of earlier
transverse coefficients. They never use \(F(0)\ne0\).
Replace the conditional uncertainty bound (TP8) by (TV3), at lengths
\(A,L\). Equations (TP14) and (TP17) follow with \(E_Y,E_Z\) in
place of \(h_Y,h_Z\). Thus the same cross-layer enclosure
\(|b+c-a|<20T\) and same-layer enclosure \(|b-c|<6T\) hold outside
the bounded-product inventory. These bound all actual directions, including
all parallel pieces through the existing capped-coordinate theorem.
The scalar-tail rank exclusions are not enlarged to include the kernel of
\(F\); that kernel was accounted for explicitly in (TV6).

For fixed \(t,a,T_0\), the common-arrangement construction also gives
a coefficient-independent complete list of
\[
 O(a^2T_0^2(1+\log(a+1))^2)
\]
positive integer test weights, of height \(O(a^2)\), for every such ideal
with \(T\le T_0\). No optimal running-time assertion is made.

## Nonunit coefficients without a valuation loss

If the unit factors of both first transverse polynomials are prefix
regular at lengths \(A-\nu,L-\nu\) whenever positive, and
\[
                   \nu_Y,\nu_Z\le(L+1)/2,                \tag{TV12}
\]
then \(E_Y=E_Z=0\). Consequently
\[
                g(I)=O(n^2\log^2(n+1)),                 \tag{TV13}
\]
with a common complete list of the same order. This includes valuations
growing proportionally to \(a\), rather than only bounded valuation.
For each fixed allowed valuation and first transverse order, these
regularity conditions define a nonempty open locus in every characteristic,
after a sufficiently large field extension if necessary. Explicit common
points over characteristic zero or \(p>2A\) are
\[
        F_Y=q^{\nu_Y}(1+q)^{A-\nu_Y-1},\qquad
        F_Z=q^{\nu_Z}(1+q)^{A-\nu_Z-1}.
\]
The same binomial determinant product is nonzero at both required lengths.
The higher transverse coefficients remain arbitrary.

Conversely any polynomial excess within the constant-axis central-only
class requires \(E_Y+E_Z\) to exceed the short-gap scale:
\[
 g(I)\ge c n^{2+\epsilon}\ \Longrightarrow\
 E_Y+E_Z=\Omega_{c,\epsilon}
               (t n^\epsilon/\log^2(n+1)).
\]
By (TV2), this requires a long unit-prefix gap or a sufficiently large
valuation excess beyond half a relevant module length. This is a
necessary condition, not a lower construction.

## Verification and priority boundary

~~~bash
python3 -B scripts/run_hook_transverse_valuation_checks.py
python3 -B scripts/run_hook_transverse_prefix_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The valuation checks reconstruct original commuting quotient columns,
test the conditional statement on full degree-bounded subspaces by exact
kernel ranks, and check both sharp witness mechanisms. The fan controls
use original coherent-wall rank tests and dual functionals with genuinely
nonunit first transverse coefficients. Finite controls do not establish an
unrestricted growth law or replace the uniform proof.

Unit Toeplitz uncertainty, valuation kernels, binomial determinants and
the earlier primal/dual fan reduction are the inputs. Exact priority of
this valuation-sensitive assembly remains unassessed. Large effective
defects, arbitrary outer images, the sharp unrestricted endpoint and the
original final-paper requirements remain open.
