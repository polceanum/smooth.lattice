# A linear endomorphism budget for actual quadratic overlap

Let \(A=K[x,y,z]/I\) have length \(n>0\) over any field, and let
\(x\) satisfy a scalar monic quadratic. Put \(B=K[y,z]\),
\(L=B1\), \(U=Bx\), \(W=L\cap U\), and \(h=\dim_KW\). Then
\[
 \boxed{\dim_K\operatorname{End}_B(W)\le n.}                    \tag{O1}
\]
If \(x^2=0\), there is the stronger sharp bound
\[
 \boxed{\dim_K\operatorname{End}_B(W)\le n-h.}                  \tag{O2}
\]
The ideal structure of the actual overlap is essential. These bounds
are not asserted for an arbitrary contracted rank-\(h\) column space.
They constrain the space of module-compatible identifications, not the
number of induced wall directions of one identification.

There is a fan-count consequence. If a radical planar ideal \(R\) of
colength \(s\ge1\) satisfies \(R^rW=0\), put
\(\ell=\binom{r+1}{2}\). Combining (O1) with the existing binary-section
bound gives
\[
 g(I)=O\!\left((1+\ell\sqrt{sn})(1+n^{4/3})\right).             \tag{O3}
\]
In particular bounded \(r,s\) give \(O_{r,s}(n^{11/6})\), uniformly
in the planar endpoints and their compatible coupling. No genericity,
split-point, characteristic or reducedness assumption on \(A\) is used.
If \(h=0\), the section bound applies directly and the annihilation
hypothesis is unnecessary.

This is a structural reduction using classical finite-length duality,
module theory and the saved section theorem. No current novelty claim,
unrestricted improvement or completion of the original goal is made.
The remaining high-overlap models have nontrivial planar action on \(W\);
the result does not replace their actual moving rank matrices.

## The planar identity behind the bound

For a finite-colength ideal \(J\subset B\) and a finite quotient
\(C=B/I_0\) of length \(p\), one has
\[
 \dim_K\operatorname{Hom}_B(J,C)
     =p+\dim_K\operatorname{Hom}_B(C,B/J).                       \tag{O4}
\]
In particular, if \(I_0\subseteq J\) and \(m=\dim B/J\), then
\[
                  \dim\operatorname{Hom}_B(J,C)=p+m.            \tag{O5}
\]

Here is an elementary all-field derivation of the classical duality
being used. For finite \(B\)-modules \(M,N\), the operators
\(y-Y_M,z-Z_M\) on \(B\otimes_KM\) give a free resolution of
\(M\) with ranks \(\dim M,2\dim M,\dim M\). The first operator
is injective by its monic leading coefficient; its quotient is
\(K[z]\otimes M\), on which the second is injective for the same
reason. The commuting-operator Koszul complex is therefore exact.
Applying \(\operatorname{Hom}_B(-,N)\) gives
\[
 H\longrightarrow H\oplus H\longrightarrow H,
 \quad H=\operatorname{Hom}_K(M,N),
\]
with maps
\[
 f\mapsto(Y_Nf-fY_M,Z_Nf-fZ_M),
 \quad(g,k)\mapsto Z_Ng-gZ_M-Y_Nk+kY_M.
\]
The trace pairing with \(\operatorname{Hom}_K(N,M)\) identifies the
dual of the last cokernel with \(\operatorname{Hom}_B(N,M)\).
Consequently
\[
 \dim\operatorname{Ext}^2_B(M,N)=\dim\operatorname{Hom}_B(N,M),
 \qquad \dim\operatorname{Ext}^0-\dim\operatorname{Ext}^1
                         +\dim\operatorname{Ext}^2=0.           \tag{O6}
\]
No residue-field trace or separability is required: this is the
nondegenerate trace pairing between rectangular \(K\)-matrix spaces.
The long exact sequence from \(0\to J\to B\to B/J\to0\),
together with (O6), proves (O4). If \(I_0\subseteq J\), every map
\(C\to B/J\) is determined freely by the image of one, proving (O5).

Let \(V\) be any ideal of \(C\), of dimension \(h\), and lift it to
\(J/I_0=V\). Inclusion followed by precomposition with \(J\to V\)
gives injections
\[
 \operatorname{End}_B(V)\hookrightarrow
 \operatorname{Hom}_B(V,C)\hookrightarrow\operatorname{Hom}_B(J,C).
\]
Since \(\dim B/J=p-h\), (O5) yields the general planar bound
\[
                      \dim\operatorname{End}_B(V)\le2p-h.       \tag{O7}
\]

## Apply it to the actual two-layer algebra

The [four-quotient reduction](quadratic_intersection_fans.md) proves
that \(W\) is an algebra ideal and that
\(p+q-h=n\), where \(p=\dim L\), \(q=\dim U\). Both \(L,U\)
are cyclic planar modules. Thus \(W\) is an ideal in each corresponding
planar quotient, and (O7) gives
\[
 \dim\operatorname{End}_B(W)\le\min(2p-h,2q-h)\le n.
\]
If \(x^2=0\), multiplication by \(x\) surjects \(L\to U\) and
kills \(W\). Therefore \(p-q\ge h\), whence \(2q\le n\).
Using the \(U\)-bound proves (O2).

For fixed actual overlap embeddings in two fixed planar modules, their
module isomorphisms form an open torsor under
\(\operatorname{Aut}_B(W)\), of dimension
\(\dim\operatorname{End}_B(W)\). The quadratic algebra condition can
cut this space further. This is a dimension statement, after scalar
extension if needed, not a count of finite-field points or fan chambers.
Four marginal rank oracles alone do not specify these embeddings.

Likewise the entire square-zero extension space with fixed
\(J\subseteq K_0\) is linear of dimension at most \(n\). In the
[exact extension classification](square_zero_extension_fans.md), put
\(m=\dim B/J\), \(s=\dim B/K_0\), so \(n=m+s\). Formula (O4)
gives
\[
 \dim\operatorname{Hom}_B(J,B/K_0)
       =s+\dim\operatorname{ann}_{B/J}(K_0)\le m+s=n.           \tag{O8}
\]
Equality holds in the flat case \(J=K_0\). Neither (O8) nor (O1)
charges each induced wall to one parameter: the already verified
[alternating-path example](cross_layer_direction_budget.md)
prevents replacing actual wall directions by the input arrow weights.

## Bounded support and thickness force smaller overlap

More generally, suppose a finite planar quotient \(S\) acts on \(W\).
Write its local factors as \(S_i\), with residue degree \(d_i\) over
\(K\) and composition length \(\ell_i\). Define
\[
                         T(S)=\sum_i d_i\ell_i^2.
\]
Then, for any finite \(S\)-module \(W\) of dimension \(h\),
\[
                    \dim_K\operatorname{End}_S(W)\ge h^2/T(S).
                                                                    \tag{O9}
\]
To prove it, put \(h_i=\dim_KW_i\), and let \(t_i,u_i\) be the
residue-field dimensions of its top and socle. Nakayama gives a
surjection \(S_i^{t_i}\to W_i\). Applying the same argument to the
\(K\)-dual module gives \(S_i^{u_i}\to W_i^*\). Hence
\(h_i\le d_i\ell_i t_i,d_i\ell_i u_i\). All maps from the top
to the socle are endomorphisms, so
\[
 \dim\operatorname{End}_{S_i}(W_i)\ge d_it_iu_i
                         \ge h_i^2/(d_i\ell_i^2).
\]
Summing and applying Cauchy--Schwarz proves (O9). Dualizing over \(K\)
preserves the residue-field dimensions used here, including inseparable
residue fields. For semisimple \(S\), the sharper exact expression is
\(\sum_i d_i m_i^2\), and \(T(S)=\dim_KS\).

Taking \(S=B/R^r\), regularity of the affine plane gives
\(\ell_i=\binom{r+1}{2}\) at every closed support point. Thus
\(T(S)=s\ell^2\). Combining (O9) with (O1) gives
\[
                              h\le\ell\sqrt{sn}.                \tag{O10}
\]
The existing [binary-section bound](profile_section_budget.md) is
\(g(I)=O((h+1)(1+n^{4/3}))\). Substitution proves (O3).
For square-zero \(x\), one may retain the sharper
\(h^2\le s\ell^2(n-h)\).

In particular a family with \(g(I)\ge c n^{2+\varepsilon}\) and
\(R^rW=0\) must satisfy
\[
                 s\binom{r+1}{2}^{\!2}
                       =\Omega_c(n^{1/3+2\varepsilon}).         \tag{O11}
\]
This is a necessary condition, not a construction. The denominator
\(T(S)\) cannot simply be replaced by \(\dim S\) for general
nonsemisimple modules: over \(K[y,z]/(y,z)^2\), take a module with
top \(e_0,e_1,e_2\), socle \(f_0,f_1\), and actions
\(ye_0=f_0,ye_1=f_1,ze_1=f_0,ze_2=f_1\), all others zero.
It has dimension five and endomorphism dimension seven, below \(25/3\).
This last module is a scope control for (O9), not an asserted overlap
realization of a new large fan.

## The unresolved generic hooks saturate the linear budget

For the [generic hook model](generic_hook_fans.md), put
\(a=t^2+1\), \(C=B/(y^a,y^tz^t,z^a)\), and
\(m=2ta-t^2\). On the already established open locus where the
\(y^t\) coefficient of \(\phi(y^a)\) and the \(z^t\) coefficient
of \(\phi(z^a)\) are nonzero, the actual overlap is isomorphic to
\(T=(y^t,z^t)C\). As a \(B\)-module it splits exactly as
\[
 T\simeq B/(y^{a-t},z^t)\ \oplus\ B/(y^t,z^{a-t}).              \tag{O12}
\]
For two rectangular cyclic modules, the Hom dimension is the product
of the two minimum side lengths: the image of one is precisely an
element killed by the two source powers. Consequently
\[
 h=2t(a-t),\qquad
 \dim\operatorname{End}_B(T)=h+2t^2=2ta=n-h,
 \quad n=2m.                                                  \tag{O13}
\]
Thus (O2) is saturated on the growing generic hook locus itself,
not only on the socle examples below. The two diagonal endomorphism
blocks are multiplication by elements of the respective rectangular
algebras. Each off-diagonal Hom block has dimension \(t^2\).
For \(t\ge3\), their direct sum is a square-zero ideal of the
endomorphism algebra. Indeed, writing \(A_0=a-t\), an off-diagonal
map into the first rectangle has image divisible by \(y^{A_0-t}\),
and one into the second has image divisible by \(z^{A_0-t}\).
In either composition the other rectangle kills that factor, since
\(A_0-t\ge t\). Thus the quotient by these two blocks is the product
of the two rectangular multiplication algebras. This is an intrinsic
module statement; it does not remove the actual monomial flags or
the central generator image from the fan problem.
The Loewy order of \(T\) at its single support point is exactly
\(a-1=t^2\), so the bounded-thickness consequence does not close
this candidate family. No generic fan count or survival of all
coherent support minors is inferred. Any further improvement must
use the multiplication and block structure, rather than reducing the
linear endomorphism budget alone.

## Sharpness and the cost of unrestricted matrix mixing

Both linear constants are sharp. For (O1), let
\(C=B/(y,z)^h\), \(h\ge2\), and split its \(h\)-dimensional
socle into nonzero complementary subspaces \(V_0,V_1\). For distinct
nonzero scalars \(\lambda_0,\lambda_1\), set
\[
 A=(C/V_1)\times(C/V_0),\qquad x=(\lambda_0,\lambda_1),
\]
with \(y,z\) induced diagonally. These generators give all of \(A\)
by interpolation in \(x\). The two planar cyclic spaces both have
dimension \(\dim C=h(h+1)/2\). Their intersection is the whole
socle of \(A\), of dimension \(h\). Thus
\(n=h^2=\dim\operatorname{End}_B(W)\). This example requires the
two stated nonzero scalars; characteristic-zero fields suffice for
sharpness. It does not assert arbitrary module gluing is compatible
with that quadratic relation.

For (O2), every field and every invertible \(h\)-by-\(h\) matrix
occurs at the exact length \(n=h^2+h\). Put
\[
 L=B/((y,z)^{h+1},y^h),\qquad U=B/(y,z)^h.
\]
The natural map \(\pi:L\to U\) has \(h\)-dimensional kernel
\(V_0\); the socle \(V_1\) of \(U\) also has dimension \(h\).
Both are killed by \(y,z\). Form the square-zero extension ring
\(H=L\oplus U\), with multiplication
\((l,u)(l',u')=(ll',\pi(l)u'+\pi(l')u)\).
For any isomorphism \(\theta:V_0\to V_1\), its graph
\(\{(v,-\theta v)\}\) is an ideal: \(\pi(V_0)=0\), and
\(\theta\) is \(B\)-linear. Its quotient is generated by \(1,y,z\)
and \(x=(0,1)\), with \(x^2=0\). The images of \(L,U\) embed,
their intersection has dimension \(h\), and
\[
 n=\dim L+\dim U-h=h^2+h,
 \qquad \dim\operatorname{End}_B(W)=h^2=n-h.
\]
Conversely, a full open set of arbitrary matrix identifications on an
\(h\)-space requires \(\dim\operatorname{End}_B(W)=h^2\).
Equation (O2) therefore forces \(n\ge h^2+h\). This is an exact
realization threshold for unrestricted overlap identifications in the
square-zero class. It is not a restriction on realizing a single fixed
matrix, even one with all minors nonzero. The sharp construction is
homogeneous and socle-coupled, inside classes already excluded as
large-fan candidates.

## Verification and prior-art boundary

```bash
python3 -B scripts/run_planar_overlap_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checks compare direct monomial syzygy equations with finite-module
Hom calculations, reconstruct both sharp families from actual quotient
columns, and solve endomorphism equations independently of the dimension
formula. Finite tests supplement the uniform proofs above.

Koszul exactness and its differential conventions are classical; see
[Stacks, 062F](https://stacks.math.columbia.edu/tag/062F) and
[0621](https://stacks.math.columbia.edu/tag/0621). The local length
calculation uses the associated-graded statement in
[Stacks, Lemma 10.106.1](https://stacks.math.columbia.edu/tag/00NN).
The trace argument spells out finite-length planar duality rather than
claiming a new duality theorem. The section estimate and its geometric
inputs retain the priority limitations in their linked proof. Exact
priority of the combined overlap threshold and support-thickness
consequence has not been established. The unrestricted extremal target
and actual direction-inventory bound remain open.
