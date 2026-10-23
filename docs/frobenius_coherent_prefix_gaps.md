# Frobenius gaps in complete positive cost prefixes

In every fixed positive characteristic, an infinite family of actual
square-zero hooks has complete positive cost prefixes which are independent
in the support matching model but have growing actual nullity. This happens
throughout a coefficient stratum with arbitrary dense first transverse
polynomials, rather than only at a special numerical coefficient choice.
The earlier-head columns and the other arm are retained in the prefix.

Let the characteristic be a prime \(p\), and put
\[
 Q=p^e\ge2,\quad d=\lfloor Q/2\rfloor,\quad t=Q+1,
 \quad A=t+2d,\quad a=A+t.
\]
In \(K[q]/q^A\), take arbitrary \(H,G\) with
\(H_0H_1\ne0\). Consider
\[
 I=(x^2,y^a,z^a,
 y^tz^t-xf),\qquad
 f=y^t(1+zH(y))+z^t(1+yG(z)).                         \tag{F1}
\]
Terms outside the hook vanish, and the usual generator syzygies hold.
Thus the quotient has reference basis \(D\sqcup xD\), where
\[
 D=\{(i,j):0\le i,j<a,\ i<t\text{ or }j<t\},\qquad
 n=4At+2t^2=\Theta(Q^2).
\]
Both outer images are zero. The degrees of \(H,G\) can grow with \(A\),
including \(\lfloor A/2\rfloor\); their other coefficients are unrestricted.

**Theorem.** There is a nonempty open set of positive weights and a
complete cheaper-column prefix such that, after compulsory coloops are
removed, the prefix has
\[
 m=t(A+1)\text{ columns},\qquad
 \operatorname{rank}_{\rm support}=m,
 \qquad \operatorname{rank}_K=m-d.                    \tag{F2}
\]
Every generic weight in this set has a unique support-minimum basis whose
actual rank is at most \(n-d\). Consequently the deficit is
\(\Omega(\sqrt n)\) along an infinite sequence in each fixed positive
characteristic. This is a rank deficit, not a fan-count lower bound.

The statement is uniform in the coefficients under \(H_0H_1\ne0\).
In particular, a rule replacing complete coherent ranks by support matching
fails generically on this first-transverse-degree-one stratum in positive
characteristic. No assertion is made about the full independent-parameter
hook family, where higher transverse and nonzero outer coefficients are
also allowed.

## The complete original prefix

Put \(r=t+d-1\), and start with
\[
 w_0=((Q+t)A,1,A),\qquad
 \tau_0=(Q+t)A+t+r+1/3.                              \tag{F3}
\]
Retain every original monomial column of cost less than \(\tau_0\).
The strict gap to every monomial cost is at least \(1/3\), so membership
is unchanged in a neighborhood of \(w_0\).

The noncompulsory columns in this prefix are exactly the following:
\[
\begin{array}{ll}
 U_i=(1,t+i,0),&0\le i\le r,\\
 L_{i,j}=(0,t+i,t+j),&t\le i<A,\ 0\le j<Q,\\
 L_{i,Q}=(0,t+i,t+Q),&t\le i\le r,\\
 B_{i,j}=(0,t+i,t+j),&0\le i,j<t.
\end{array}                                                   \tag{F4}
\]
The last line includes every earlier central head, with none discarded.
There are no upper Z-arm columns or pure lower Z-arm columns in the prefix.
Zero outer columns remain loops and have no effect on either rank.

Here are direct price checks. A central source has cost
\(t+tA+i+jA\), so its source cutoff is
\(QA+r+1/3\). All Y sources at levels \(j<Q\) occur, since
\((Q-1)A+A-1=QA-1\); at \(j=Q\) the cutoff is \(i\le r\).
The maximum head source cost is \(Q+QA<QA+r+1/3\).
The cheapest Z pure source, with long source index \(t\), exceeds
the cutoff by \(A-r-1/3=d+2/3\).
Upper Y prices relative to their fixed translation have cutoff
\(i+jA<r+1/3\), giving just the first line, since \(r<A\).
The cheapest upper Z column exceeds \(\tau_0\), because
\(tA>t+r+1/3\).

The number of Y columns in the first three lines is
\[
 (r+1)+2dQ+d=t(2d+1).
\]
Adding the \(t^2\) heads gives \(m=t(A+1)\).

## Exact rank, with all heads present

On the Y rectangle \(R=K[q,v]/(q^A,v^t)\), write
\(W=1+vH(q)\). The selected lower columns are
\(Wq^iv^j\) with the indices in (F4). They are independent, because
\(W\) is a unit. Their number is \(2dQ+d\).

For \(t\le i\le r\), the finite inverse identity gives
\[
 q^i=W\sum_{j=0}^{Q}(-1)^j q^iH(q)^jv^j.             \tag{F5}
\]
At every level \(j<Q\), all sources on the right have long index at
least \(t\), and every such source is selected. Frobenius gives
\[
 H(q)^Q=\sum_k H_k^Qq^{kQ},\qquad
 q^iH(q)^Q=H_0^Qq^i\pmod {q^A},                      \tag{F6}
\]
since \(i\ge t\) and \(t+Q\ge A\). Hence the last level in (F5)
also uses only selected sources. All \(d\) upper columns with
\(t\le i\le r\) lie in the selected lower span.

The remaining \(t\) upper columns, with \(i<t\), are independent
modulo that span: every selected lower vector is divisible by \(q^t\).
Thus the Y-prefix nullity is exactly \(d\).

Only heads have a Z component in this prefix. Their Z projections are
the \(t^2\) corner monomials multiplied by the unit \(1+vG(q)\)
in the transposed rectangle. Those projections are independent. Therefore
the heads increase actual rank by exactly \(t^2\). They cannot repair
the Y-prefix dependencies. Compulsory prefix columns are coloops in the
full original residue matrix and contribute equally to both ranks.
This proves the actual-rank equality in (F2), using the complete original
prefix rather than an isolated arm minor.

## A coefficient-independent perfect matching

Match each selected upper column to its own Y row. For every
\(t\le i\le r\), route the displaced lower column along
\[
 (i+\min(j,d),j),\qquad 0\le j\le Q.                 \tag{F7}
\]
For \(j<d\), the next edge uses the nonzero coefficient \(H_1\);
for \(d\le j<Q\), it uses \(H_0\). Distinct paths never share
an original row. All intermediate source columns are selected, and their
end rows have long indices \(t+d,\ldots,t+2d-1\), beyond the selected
last-level sources. Match every lower column off these paths to its own
identity row. This matches all Y prefix columns to distinct rows.

Match the heads to their distinct constant-term Z rows, and compulsory
prefix columns to their own reference rows. This is a perfect support
matching of the entire prefix, establishing the other equality in (F2).

At any strict monomial cost order, a support-greedy basis contains every
column of an independent prefix. The \(d\) actual dependencies persist
after the other basis columns are added. Its actual rank is consequently
at most \(n-d\).

An explicit strict positive integer weight is available if desired.
Let \(M=4a+1\), \(K_0=300aM^2\), and take
\[
                 w_*=K_0w_0+(1,M,M^2).                       \tag{F8}
\]
All nonzero column labels have x degree at most one and y,z degrees
less than \(a\). The perturbation is smaller than the scaled prefix
gap, so the cutoff \(K_0\tau_0\) selects exactly (F4) and its compulsory
columns. Distinct labels of equal \(w_0\)-cost have distinct balanced-base
\(M\) perturbations; unequal integer \(w_0\)-costs cannot be cancelled
by that perturbation. Thus the support minimum is unique.

## Characteristic-zero comparison

Eliminating the identity pivots in the Y prefix leaves a \(d\)-square
matrix, up to a nonzero sign,
\[
                [\,[q^{d+i-j}]H(q)^Q\,]_{0\le i,j<d}.        \tag{F9}
\]
One may also obtain (F9) by transporting a free first-level source
through all \(Q\) lower blocks. In positive characteristic this matrix
is zero by (F6).

In characteristic zero, specialize \(H=(1+q)^D\), \(D\ge1\).
Only coefficients below \(A\) are relevant. The matrix in (F9) becomes
\([\binom{DQ}{d+i-j}]\). The classical binomial Toeplitz product is
\[
 \prod_{j=0}^{d-1}
 \frac{(DQ+j)!\,j!}{(d+j)!\,(DQ-d+j)!}>0.                     \tag{F10}
\]
Therefore this same prefix is independent on a nonempty coefficient-open
set in characteristic zero. This comparison does not prove that every
positive prefix is generically independent there.

## Verification and scope

Run from the repository root:

```bash
python3 -B scripts/run_frobenius_prefix_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker constructs original quotient columns independently, verifies
the complete strict prefix and the explicit perfect matching, computes
its actual rank, checks all dependencies in (F5), and inspects the unique
whole support-minimum basis. Prime fields and a proper characteristic-two
extension are included. Small cases are compared with the existing
independent original-quotient constructor. Characteristic-zero controls
check (F9)--(F10) by exact integer arithmetic.

Frobenius is a classical ring endomorphism; see
[Stacks, Section 33.36](https://stacks.math.columbia.edu/tag/0CC6).
The determinant in (F10) is classical; see
[Krattenthaler, Section 2.3, (2.16)--(2.18)](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11).
The earlier [coherent convolution obstruction](coherent_convolution_rank_gap.md)
already disproves a stratum-independent matching rule using nonzero outer
images. Here the outer images vanish, and the obstruction holds throughout
the stated first-transverse coefficient stratum in each positive
characteristic. Exact priority of this combined construction is unassessed.

No lower bound on the number of actual cones follows from the rank deficit.
For this parameter sequence, \(a=\Theta(\sqrt n)\), so the existing
[box ceiling](hook_bridge_reduction.md) \(O(n^{3/2}a)\) already gives
\(g(I)=O(n^2)\). The family therefore cannot supply a polynomial excess
over the generic fan scale. Its role is to test the proposed rank rule,
including the exact cheaper-prefix and coefficient-generic qualifications.
The unrestricted sharp \(G_3\) problem, a polynomial fan-count separation,
the consequential-theorem requirement, current novelty and the final paper
remain open. The characteristic-zero coherent matching question and actual
winning exchange comparisons require further work.
