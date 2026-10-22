# Axis-constant hooks retain arbitrary smaller couplings in a weight sector

Long prefix gaps need not have the five-cone fan of the zero-off-axis
example. There is an exact sector reduction: every arbitrary central-only
hook occurs inside the next-thickness axis-constant hook. Within the sector,
a cost-compatible column reduction gives translated child labels and fixed
halo coloops. An explicit
infinite family below has identical axis polynomials and prefix determinants
to a five-cone family, but a linearly growing lower bound for its number
of initial ideals.
This does not give a superquadratic lower separation.

## Exact reduction in the low-price sector

Let \(t\ge1\), \(a\ge2t+1\), and
\[
 I_f=(x^2,y^a,z^a,y^tz^t-xf),\qquad n=4at-2t^2,            \tag{AS1}
\]
where \(f\) is any compatible central image in the planar hook quotient.
No leading-unit or prefix-regular assumption is made on \(f\).
Put \(s=t+1\), choose \(\alpha,\beta\in K^*\), and form
\[
 P_f=(x^2,y^{a+1},z^{a+1},
       y^sz^s-x(\alpha y^s+\beta z^s+yzf)).                \tag{AS2}
\]
This is a compatible flat hook of thickness \(s\), cutoff \(a+1\),
and length \(n+4a+2\). Its two axis units are constant. All its
off-axis coefficients are precisely those of \(yzf\); this is a
bijection between the allowed central images in (AS1) and the off-axis
images in (AS2) with those fixed axis coefficients.

Normalize a positive weight to \(w=(\lambda,u,1)\), and let
\[
                  \mathcal S_s=\{0<\lambda<s\min(u,1)\}.  \tag{AS3}
\]
Define the fixed halo
\[
 H_a=\{(\epsilon,i,j):\epsilon\in\{0,1\},\ 0\le i,j\le a,
                                      \quad ij=0\},
 \qquad |H_a|=4a+2.
\]

**Theorem.** At every generic weight in \(\mathcal S_s\), the initial
standard sets correspond by
\[
                  B_{P_f}=H_a\ \sqcup\ ((0,1,1)+B_{I_f}). \tag{AS4}
\]
Their full normal-cone intersections with this sector agree. Walls,
vertices and recurrent intervals inside it are retained. In particular
the number of distinct child cones meeting the sector is a lower bound
for \(g(P_f)\). This is a sector statement, not equality of the two
complete positive fans.

### Removing the axis terms preserves every cost-prefix space

First replace the central image of \(P_f\) by \(yzf\), obtaining
\(\widetilde P_f\). Use the original complete hook residue-column
matrices, with the same reference basis and labels for both ideals.
A central-source column of label \((0,s+i,s+j)\) changes by
\[
       \alpha v_{(1,s+i,j)}+\beta v_{(1,i,s+j)}.            \tag{AS5}
\]
Every nonzero helper on the right is an unchanged upper identity column.
Its cost minus the source cost is respectively \(\lambda-s\) and
\(\lambda-su\), both strictly negative in (AS3). A helper which
vanishes outside the upper diagram contributes zero.

Induct over the cost levels. Every changed column differs by a combination
of strictly earlier unchanged columns. Thus the two matrices have identical
prefix spans and identical quotient matroids on each tied cost level.
All minimum-weight monomial bases have the same labels. This argument
also covers generic refinements at a tie; continuity gives the closed
normal-cone intersections with the open sector. It does not discard an
expensive helper or change a coefficient based on support matching.

### The resulting ideal is an actual monomial padding

Put \(J=(x^2,y^{a+1},z^{a+1})\). There is an exact ideal identity
\[
 \widetilde P_f=J+yzI_f,\qquad
 (\widetilde P_f:yz)=I_f,\qquad
                  J:yz=(x^2,y^a,z^a)\subset I_f.           \tag{AS6}
\]
All terms outside the \(yz\)-cone belong to the monomial part alone.
The [monomial padding identity](canonical_fan_deflation.md#the-fixed-halo-and-exact-initial-identities)
therefore gives, for every global order,
\[
 \operatorname{in}\widetilde P_f=J+yz\operatorname{in}I_f,
 \quad B_{\widetilde P_f}=H_a\sqcup((0,1,1)+B_{I_f}).       \tag{AS7}
\]
Equivalently the dropped-parent matrix is a direct sum of the halo coloops
and the original child matrix with translated labels. The profiles satisfy
\[
 \mu(B_{P_f})=\mu(H_a)+n(0,1,1)+\mu(B_{I_f})
                         \quad\hbox{in }\mathcal S_s.     \tag{AS8}
\]
This proves the theorem with an actual colon algebra, rather than a
quotient by an unjustified vector-space span.

The axis polynomials of (AS2) are constant units, regardless of \(f\).
The prefix determinants are nonzero only at the two endpoint indices
\(s,a-t\), with the empty-endpoint convention. Their maximum gaps are
\(\max(0,a-2t-2)\). Thus a maximal-gap unit stratum can encode arbitrary
lower-thickness central coupling, including child images with zero leading
coefficients. It cannot be classified from its axis determinants alone.

## Identical axis data with an unbounded fan-count ratio

Take an integer \(a\ge5\) in characteristic zero, or characteristic
\(p>2(a-1)\). Use the thickness-one child
\[
 f_a=y(1+y)^{a-2}+z+z^{a-1},\qquad
 I_a=(x^2,y^a,z^a,yz-xf_a),
\]
and the thickness-two parent
\[
 P_a=(x^2,y^{a+1},z^{a+1},
             y^2z^2-x(y^2+z^2+yzf_a)),\qquad N=8a.       \tag{AS9}
\]
Its axis polynomials and all prefix determinants equal those of
\(P_a^0\), obtained by replacing its central image by \(y^2+z^2\).
The latter has exactly five cones by (PG11). Nevertheless
\[
                         g(P_a)\ge\lfloor(a-1)/2\rfloor, \tag{AS10}
\]
so their ratio is unbounded, already at fixed thickness two.
Consequently any estimate using only the length, axis polynomials and
prefix determinants, and approximating every such fan count within a
multiplicative factor \(R(N)\), must satisfy
\[
                        R(N)^2\ge (N-16)/80               \tag{AS10a}
\]
along this sequence. Indeed both ideals supply identical input data,
whereas their counts have ratio at least \((a-2)/10\). Thus these
data cannot give a uniform constant-factor estimate, even over
\(\mathbb Q\). Equal full prefix ranks are not asserted.

Here is a direct coherent proof. Put \(A=a-1\), and contract the
literal fixed columns of the child. Its upper \(Y\)-columns have local
indices \(0\le j<A\). Its pure central \(Y\)-columns have source
indices \(1\le i<A\), with binomial convolution multiplier
\(U=(1+y)^{A-1}\). After a common affine cost is subtracted, their
costs are
\[
                   j,\qquad i+\delta,
                   \quad\delta=(1-\lambda)/u.             \tag{AS11}
\]
If \(0<\lambda<u<1/2\), a nonempty price prefix consists of upper indices
\(0,\ldots,r-1\) and source indices \(1,\ldots,k\), with \(r\ge k+1\)
whenever a source occurs. Empty prefixes have rank zero.
Its actual rank is \(\min(A,r+k)\). Eliminate the upper identities;
the needed square minor on rows \(r,\ldots,r+q-1\), sources
\(1,\ldots,q\), is
\[
              \det[\binom{A-1}{r-1+i-j}]_{0\le i,j<q},
                    \qquad q=\min(k,A-r).
\]
The [classical binomial Toeplitz product](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11)
is nonzero in the stated
characteristics; its factorial arguments are at most \(2A-2\).
Hence the \(Y\)-minimum basis is exactly the first \(A\) columns in
the merged streams (AS11).

On the \(Z\) arm, all pure source columns are parallel to the matching
upper identities and are more expensive, since \(\lambda<u\).
The one bridge has a nonzero coefficient on the last upper \(Z\)-row,
of cost \(\lambda+A\). This exceeds every \(Y\)-column cost and the
bridge cost when \(u<1/2\). The corank-one greedy/dual argument therefore
selects the bridge, drops precisely that last \(Z\)-column, and keeps
the entire \(Y\)-minimum basis. Thus its changes survive in the whole
child algebra.

For each integer
\[
        2\le h\le A-2,\quad h\equiv A\pmod2,\qquad
        u_h=\frac2{2h+1},\quad\lambda_h=\frac1{2h+1},      \tag{AS12}
\]
one has \(\delta=h\). Exactly \(A-1\) pure \(Y\)-columns lie below
the tied source and upper pair
\[
 i=(A-h)/2,\quad j=(A+h)/2,
 \qquad (0,i+1,1),\quad(1,j+1,0).
\]
Each is nonzero modulo those cheaper columns, and the two are parallel
there. This is an actual coherent wall of primitive direction
\((1,h,-1)\). Perturbing \(\lambda\) to either side exchanges the pair.
The source count decreases by one at each successive eligible \(h\),
so the \(\lfloor A/2\rfloor-1\) walls yield at least
\(\lfloor A/2\rfloor\) distinct states. All these weights lie in
\(\mathcal S_2\), and (AS4) preserves their states and directions in
the parent. This proves (AS10).

The proved lower bound is linear in \(N\); it does not establish a
superquadratic count.
It establishes a limitation of axis-only descriptions inside the long-gap
stratum. It supplies no growth law for arbitrary off-axis coefficients or
an upper bound for the full parent fan outside (AS3).

## Verification and priority boundary

```bash
python3 -B scripts/run_hook_axis_sector_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker reconstructs original commuting quotient columns, verifies
every triangular column identity and the padding row map, and compares
entire finite fan patches within the sector. It separately verifies all
price-prefix ranks used in the lower family and its adjacent original-field
minimum bases. Equal axis determinants are distinguished from equal full
prefix ranks, which are not asserted.

Cost-compatible column elimination, monomial padding, greedy duality and
the binomial determinant product are classical inputs. Their sources and
limits remain in the linked proofs. Exact priority of the assembled sector
recursion and explicit axis-data obstruction is unassessed. The full
unrestricted growth problem, a polynomial lower separation above the
generic corner-cut scale, and the original final-paper criteria remain open.
