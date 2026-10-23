# Monomial orbits give a quadratic hook class with large interface rank

Small interface rank is not necessary for a quadratic fan. A decomposition
into actual row components removes many pure-arm columns without a
rank perturbation. It gives a full-fan bound for one arm polynomial in
one monomial, allowing arbitrary coefficients on the other arm. When
the monomial's long-to-short slope is bounded, the bound is quadratic
independently of thickness.

The result includes a growing family whose bridge and far-interface
ranks both have order \(t^2\), and whose departure rank has order \(n\).
Thus these large ranks do not, by themselves, identify a large fan.
The unrestricted sharp count, an actual polynomial separation, exact
current novelty and the original final-paper criteria remain open.

## A general component bound for compatible hooks

Use the original flat hook model with
\[
 a\ge2t,\quad t\ge2,\quad A=a-t,\quad n=4at-2t^2,
\]
its residual monomial rows \(T_Y\oplus T_Z\), and groups
\(E_Y,E_Z,H\). In the bipartite support graph on Y rows and nonzero
\(E_Y\) columns, join a column to each row where its coefficient is
nonzero. Mark every connected component met by a Y projection of a
bridge. Let \(L_Y\subseteq E_Y\) be the nonzero columns in these
marked components, and put
\[
 \sigma_Y=|\{e_y-f_y:e\in L_Y,\ f\in E_Z\}|.             \tag{MO1}
\]
Using all potential Z labels only enlarges this quantity. Define
\(\sigma_Z\) by transposition. For every compatible map, every field
and arbitrary outer images,
\[
       g(I_\phi)=O\!\left(n^2
             \left(1+\min(\sigma_Y,\sigma_Z)/t\right)^2\right).
                                                               \tag{MO2}
\]

Indeed the unmarked Y components form actual direct summands of the
represented residue matroid: no bridge reaches them and no Z column
uses Y rows. Their minima are fixed by the internal \(E_Y\) order.
The remaining ground is \(L_Y\cup E_Z\cup H\), with a fixed actual
matrix, independent of the weights.

All comparisons within one arm or involving a bridge have \(O(at)\)
distinct planes, by the original arm and three bridge boxes. A
remaining cross-arm comparison has \(\sigma_Y\) possible Y
differences, \(O(a)\) possible Z differences and three possible first
differences. Thus \(O(a(t+\sigma_Y))\) planes fix the entire core
order and the omitted summands' minima. Their positive arrangement
has \(O(a^2(t+\sigma_Y)^2)\) chambers. Since \(at=\Theta(n)\),
this proves (MO2). Each full optimum is reconstructed by adjoining
the selected bases of the omitted components. No vector contraction
is declared to be an algebra quotient.

The support graph gives a convenient canonical decomposition. Any
coarser row partition containing each pure column in one part also
works, if every part met by a bridge is retained. This observation
permits a coefficient-independent partition in the class below.

## One arm polynomial in one monomial

Assume the Y components of both outer images vanish, and its central
arm operator has the form
\[
 W_Y=F(q^\nu v^\delta)\quad\text{in }K[q,v]/(q^A,v^t),
 \quad \nu\ge0,\quad1\le\delta<t,\quad F\in K[T].        \tag{MO3}
\]
There is no unit or nonvanishing assumption on \(F\). All Z components
of all three generator images, including both outer images, are
arbitrary compatible coefficients. Terms of degree at least
\(h=\lceil t/\delta\rceil\) in \(F\) vanish automatically.

For \(\nu>0\) set
\[
 \beta=\min\left\{\left\lceil t/\delta\right\rceil,
                 1+\nu/\delta,
                 \left\lceil A/\nu\right\rceil,
                 A/t\right\},                         \tag{MO4}
\]
and use \(\beta=1\) when \(\nu=0\). Then, with an absolute constant,
\[
                       \boxed{\ g(I_\phi)=O(n^2\beta^2).\ }\tag{MO5}
\]
In particular a bounded ratio \(\nu/\delta\) gives \(O(n^2)\)
throughout this class, even for growing thickness and arbitrary
opposite-arm coefficients. The same holds for bounded \(A/\nu\).
For a given field and fixed \(t,a,\nu,\delta\), one common complete
weight list works for **every** such coefficient choice:
\[
 |\mathcal W|=O(n^2\beta^2),\qquad
                   \|w\|_\infty=O(a^2).                \tag{MO6}
\]
This is a compilation consequence in addition to the fan bound.

## Keep precisely the monomial chains reached from the core

Partition the Y rectangle into chains with step \((\nu,\delta)\).
For example the pair
\[
                  (\delta q-\nu v,\ v\bmod\delta)        \tag{MO7}
\]
identifies a chain, including nonprimitive steps. Each power of the
monomial preserves this pair. Thus every pure Y column, including an
upper identity, is supported in one chain.

Mark every chain containing an input \((q_0,v_0)\) with \(q_0<t\).
The central bridges have such inputs, so all their Y components lie
in marked chains. The two outer Y components are zero. Unmarked
chains therefore split off exactly, for every polynomial \(F\), even
when some coefficients vanish. This is a coarse marking; the actual
support components can be smaller.

If \((q,v)\) lies in a marked chain, then
\((q,v)=(q_0,v_0)+k(\nu,\delta)\) for an integer \(k\).
The short coordinates imply \(|k|\le h-1\). Negative \(k\)
gives \(q<t\). Hence every marked long coordinate lies in
\[
 S=\bigcup_{k=0}^{h-1}(k\nu+[0,t))\ \cap\ [0,A),         \tag{MO8}
\]
using integer indices. A nonzero pure lower Y column has its input
in the same chain as its output, so its original long label lies
in \(t+S\), just as for upper identities.

The Y differences with a Z label are contained in a translate of
\(S+[-2t,2t]\). The number of such integers is \(O(t\beta)\).
There are four elementary bounds: at most \(h\) intervals of length
\(O(t)\); a containing interval of length \(O(t+h\nu)\);
at most \(\lceil A/\nu\rceil\) nonempty intervals when \(\nu>0\);
and a containing interval of length \(O(A)\). Since
\(h\le2t/\delta\), they give exactly (MO4), up to absolute constants.
For \(\nu=0\), \(S=[0,t)\). Equation (MO2), or its same proof
with the coarser chain marking, proves (MO5).

All these chains and label inventories are independent of coefficients.
Arrange their core comparison planes together with the original pure-Y
internal planes. One integer weight in each positive chamber is the
list in (MO6). Original comparison entries have magnitude less than
\(2a\). Extreme rays are cross products of two plane normals and
have primitive integer height \(O(a^2)\); the sum of three independent
extreme rays lies strictly inside a full chamber and has the same
height order. Exact quotient-rank evaluation at the list recovers every
positive initial ideal. No optimal construction-time or bit-complexity
bound for all rank evaluations is asserted.

## An explicit family with large earlier rank parameters

Take zero outer images, \(A\ge2t-1\), and on both arms put
\[
                       W=(1+qv)^{t-1}.                  \tag{MO9}
\]
Work in characteristic zero or characteristic \(p>2t\). Then
\[
 b=t^2,\qquad
 \chi_Y=\chi_Z=\lfloor t^2/4\rfloor,\qquad
                    \rho=(A-1)(t-1).                   \tag{MO10}
\]
Nevertheless (MO5) gives the full bound \(O(n^2)\), since
\(\nu=\delta=1\) and \(\beta\le2\). More generally the chosen
Y arm may be any polynomial in \(qv\), with the Z arm arbitrary.

Here is an exact check of the interface order in (MO10). The invariant
\(d=q-v\) separates its cross-cut matrix into blocks. Only
\(1\le d<t\) contribute. The block has \(d\) high rows and
\(t-d\) low inputs, with entries
\(\binom{t-1}{r-j}\), \(t-d\le r<t\), \(0\le j<t-d\).
Its first \(m=\min(d,t-d)\) rows and columns have determinant
\[
 \prod_{i=0}^{m-1}
 \frac{(t-1+i)!\,i!}{(t-d+i)!\,(d-1+i)!}\ne0.          \tag{MO11}
\]
This is the classical binomial Toeplitz evaluation in
[Krattenthaler, Section 2.3, (2.17)](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11).
All factorial arguments are below \(2t\). The block rank is therefore
\(\min(d,t-d)\), and summing gives \(\lfloor t^2/4\rfloor\).
The bridge rank follows from the injective unit Y projection.

For departure rank, subtract any transverse polynomial \(C(v)\).
If its constant term differs from one, the difference is a unit.
If its constant is one but its \(v\) coefficient is nonzero, the
difference is \(v\) times a unit. Otherwise the first short-degree
coefficient is \((t-1)q\). Images of \(q^iv^j\),
\(i<A-1,j<t-1\), are independent by their successive lowest short
degrees. Taking \(C=1\) attains the resulting lower bound because
\(W-1=qv\) times a unit. This proves the value of \(\rho\).

At \(a=t^{23}\), the interface theorem's named ceiling has order
\(n^2t^6=n^{9/4}\), whereas this whole family is quadratically
bounded. This comparison is to that coarse ceiling, not to every
coefficient-specific prior result or a strongest-current theorem.
Large rank and wide transverse support are still not sufficient for
an actual polynomial excess. The unrestricted problem remains open.

## A second direction can mix the whole long arm

The single-monomial hypothesis cannot be replaced by a bounded cone of
support slopes in this argument. Take \(A\ge2t+2\), zero outer images,
and the same operator on both arms:
\[
 W_\lambda=(1+qv)^{t-1}+\lambda q^2v,\qquad\lambda\ne0. \tag{MO12}
\]
In characteristic zero or \(p>2t\), this actual compatible family has
\[
 b=t^2,\quad \chi_Y,\chi_Z\ge\lfloor t^2/4\rfloor,
 \quad \rho=(A-1)(t-1),\quad
 \sigma_Y,\sigma_Z\ge A-t-1.                            \tag{MO13}
\]
Thus both interface ranks have order \(t^2\), departure has order
\(n\), and both marked long-difference inventories have order \(A\).
Every nonconstant support monomial \(q^iv^j\) nevertheless satisfies
\(j\le i\le2j\).

To prove the support statement, a pure column with input \((q,0)\),
\(t\le q\le A-3\), has nonzero rows \((q,0)\),
\((q+1,1)\), \((q+2,1)\). The next pure column, with input
\((q+1,0)\), shares \((q+2,1)\). Consequently all rows
\((q,0)\), \(t\le q\le A-2\), belong to one component.
The bridge input \((t-1,0)\) meets that component at \((t+1,1)\).
It is therefore marked. Its upper identities have long labels
\(t+q\); comparison with the opposite upper identity of long label
zero supplies \(A-t-1\) distinct differences. Transposition proves
the other arm's statement.

For the interface lower bound, retain exactly the square submatrices
used in (MO11), for all \(1\le d<t\). The original operator preserves
\(d=q-v\), while the added term increases it by one. Ordering these
row and column blocks by \(d\) therefore makes their combined square
matrix block triangular, with the unchanged nonzero (MO11) minors
on the diagonal. Its determinant is their product, independently of
\(\lambda\). This proves both interface lower bounds for every
nonzero parameter in the stated characteristics, without a genericity
assumption. The constant term of \(W_\lambda\) is one, proving the
bridge rank. Its first short-degree coefficient is
\(q((t-1)+\lambda q)\), a multiple of \(q\) by a unit. The same
lowest-short-degree argument as above proves the departure formula;
\(W_\lambda-1=qv\) times a unit attains it.

This is an obstruction to extending the component argument by slope
width alone, and a concrete family retaining the remaining difficulties.
It is not a lower bound on the number of initial ideals. Counting its
actual coherent exchanges, rather than its connected supports or ranks,
is still necessary for the original goal.

## Verification and priority boundary

~~~bash
python3 -B scripts/run_hook_monomial_orbit_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker reconstructs original quotients and commuting multiplication,
compares support components with independent chain keys, verifies whole
original/core cone patches, and constructs complete small comparison
arrangements. It checks actual fan coverage and original integer weights.
Separate exact binomial minors, interface blocks and departure comparisons
support (MO10). General outer-image cases exercise (MO2) independently
of the monomial hypothesis. The mixed-direction controls separately
check the original quotient, shared support paths, the combined integer
minor and departure comparisons in (MO13).

Support decomposition, greedy selection and hyperplane arrangements are
classical; the binomial determinant is explicitly attributed above.
The quantitative class bound and common-list construction are proved
here, but strongest-current exact novelty and sharpness are unestablished.
The single-monomial support condition is substantive: several independent
monomial steps can connect the chains and are not covered by (MO5).
No unrestricted sharp exponent, actual polynomial lower separation or
completed final-paper claim follows.
