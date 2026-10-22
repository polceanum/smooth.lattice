# A linear rank gap at a positively coherent support basis

Perfect support matching does not decide positively coherent basis
availability for compatible square-zero multiplication, even generically
on a linear coefficient stratum. There is an explicit five-parameter
example over every field. Finite free power substitution gives lengths
\(N=240s^2\), for every \(s\ge1\), at which a unique positive support
minimum has actual column rank
\[
                         N-N/240.                         \tag{1}
\]
The entire actual fan is preserved, up to a diagonal change of weights,
as \(s\) grows. Thus this rank gap supplies no growing fan count.

This refutes a stratum-independent coherent matching rule. It does
**not** refute the full independent-parameter generic hook model:
the example sets most of that model's parameters to zero. The complete
generic model remains certified only at thicknesses two and three.
There is no new unrestricted bound for \(G_3\), polynomial lower
separation, or novelty assertion here.
Indeed at the displayed weight the full-parameter support model has
a different basis, cheaper by 4438, so this sparse cone is not being
passed off as a cone of that larger relaxation.

## The actual quotient and its support relaxation

Let
\[
 M=(y^{17},y^4z^4,z^{17}),\qquad C=K[y,z]/M,
 \qquad D=\{(i,j):0\le i,j<17,\ i<4\text{ or }j<4\}.
\]
Then \(|D|=120\). For independent coefficients \(A,B,D_0,E,H\), set
\[
 \begin{split}
 f_0&=A z^{15}+B yz^{13},\\
 f_1&=D_0 yz^{12}+E y^2z^{13}+H y^{14}z,\\
 f_2&=0.
 \end{split}                                                \tag{2}
\]
These are allowed images of the three generators of \(M\) under a
\(K[y,z]\)-linear map \(\phi:M\to C\). Both adjacent syzygies vanish:
\(z^4f_0=y^{13}f_1=0\) and \(z^{13}f_1=y^4f_2=0\) in \(C\).
The quotient
\[
 A_0=K[x,y,z]/(x^2,xM,y^{17}-xf_0,y^4z^4-xf_1,z^{17})          \tag{3}
\]
has the reference basis \(D\sqcup xD\), of length 240. This follows
also by defining multiplication on \(C\oplus xC\): the syzygies
make the reductions independent of the chosen divisor of a monomial.

The plain \(D\)-columns are compulsory coloops. After contracting
them the actual matrix is \([1_D\mid\Phi]\), with the original
monomial labels. Replace its nonzero positions by support matching
to obtain a transversal matroid of rank 120. All assertions about
this support pattern use \(ABD_0EH\ne0\); this torus is nonempty
over every field, including the two-element field.

## One complete positive support cone

Remove from the upper identity columns the following 18 rows:
\[
 R=\{(0,15),(0,16)\}\ \cup\
       \{1,2,3\}\mathbin\times\{12,13,14,15,16\}
       \ \cup\ \{(16,3)\}.                                \tag{4}
\]
Insert the following 18 plain columns:
\[
 \begin{split}
 Q_0&=\{(17+i,j):0\le i\le3,\ 0\le j\le1\},\\
 Q_1&=\{(4,4+j):0\le j\le3\}\cup
       \{(4+i,4+j):1\le i\le2,\ 0\le j\le2\}.
 \end{split}                                                \tag{5}
\]
Together with plain \(D\), these give a set \(\mathcal B\) of 240
monomials. It is a downset and is a basis of the support relaxation.
It is its unique minimum throughout the open cone
\[
 X>0,\qquad 2Z>5Y,\qquad13Y>5Z,\qquad5Z>X+11Y.              \tag{6}
\]
These inequalities imply \(Y,Z>0\). The integer weight
\((X,Y,Z)=(39,43,109)\) is strictly inside.

Here is a finite exact certificate for the claim, rather than a
sampling assumption. In a basis with a perfect matching, the rows
reachable from a new column by alternating paths give exactly its
possible single-column exchanges. Enumerating these paths on the
149 nonzero residual columns gives 46 distinct primitive comparison
normals. All 46 are strictly positive at \((39,43,109)\), with
minimum value 3. Their intersection with the normalized positive
triangle has, in cyclic order, vertices
\[
 (0,2/7,5/7),\quad(0,5/18,13/18),\quad
 (5/14,5/28,13/28),\quad(3/10,1/5,1/2).
\]
Its doubled area is \(11/2520\); its three non-coordinate facets
are exactly those in (6). The report supplies the matching, all
46 normals and an independent minimum-cost assignment dual.
The matroid minimum-basis exchange criterion proves uniqueness on
the whole open cone. No actual determinant is inferred from matching.

## The determinant cancellation and exact rank

Eliminate the selected identity columns. Order the remaining rows
lexicographically as in (4), and the remaining columns lexicographically
as in (5). The resulting matrix \(L\) has order 18. Its determinant
has exactly two nonzero permutation terms. They have opposite signs
and the same product
\[
                         A^7B D_0^8 E H.
\]
Thus \(\det L=0\) over the integral parameter ring and in every
characteristic.

One can see the cancellation without expanding a large determinant.
Repeatedly expand at a row or column with only one nonzero. Fourteen
such expansions leave, with the remaining rows and columns in order,
\[
 \begin{pmatrix}
 D_0&0&A&0\\
 0&D_0&B&0\\
 E&0&0&A\\
 0&E&0&B
 \end{pmatrix}.                                           \tag{7}
\]
The kernel vector \((A,B,-D_0,-E)\) explains the identity. The
eliminated factors multiply to \(A^6D_0^7H\), up to sign.
Deleting zero-based row 2 and column 3 of the original \(L\) gives
a 17-minor with the single signed term
\[
                           A^8D_0^6E^2H.
\]
Hence \(\operatorname{rank}L=17\) whenever \(AD_0EH\ne0\), and
\(\mathcal B\) has actual column rank 239 on the coefficient torus.
This is failure at an exposed support minimum, not merely at an
arbitrary noncoherent minor. It is also not a characteristic-specific
failure of one numerical specialization.

## A uniform linear deficit with the same full fan

For any integer \(s\ge1\), extend (3) by the substitution
\(y=Y_0^s,z=Z_0^s\), retaining \(x\). The target polynomial ring
is free of rank \(s^2\) over the source, with basis
\(Y_0^rZ_0^q\), \(0\le r,q<s\). The quotient has length
\(N=240s^2\). Every monomial column lies in precisely one residue
block, and each block is a copy of the original quotient matrix.
The same is true of the support relaxation.

The expanded set
\[
 \mathcal B_s=\{x^iY_0^{sj+r}Z_0^{sk+q}:
                 x^iy^jz^k\in\mathcal B,\ 0\le r,q<s\}
\]
is therefore a support basis of actual rank \(239s^2\).
At weight \((39s,43,109)\), costs within a residue block equal
\(s\) times the original costs plus the same tag cost per column.
The minimum basis in each block is unique, proving (1) at a full
positive cone for every \(s\).

More generally, at every positive weight \(w\) the optimization
splits into identical copies of the original problem with weight
\((w_x,sw_y,sw_z)\). Thus the entire actual positive fan, and
separately the support fan, is transported by this diagonal map.
In particular neither fan count grows with \(s\). Free power
substitution amplifies the rank discrepancy and length together;
it supplies no lower bound above the generic corner-cut count.

## Structural downsets do not certify availability

There is a useful general support-level fact behind these checks.
Suppose multiplication of column labels induces a partial injection
on rows, with each new neighborhood exactly the shifted old
neighborhood after deleted rows are removed. This holds for \(y,z\)
on the contracted matrices above. Such multiplication preserves
transversal-matroid closure.

Indeed a circuit \(T\) has exactly \(|T|-1\) neighboring rows by
minimal Hall dependence. For any \(a\in T\), a matching of
\(T\setminus\{a\}\) covers all those rows. Restricting that matching
to the surviving rows after multiplication still covers every new
neighbor. The multiplied \(a\) is consequently in the closure of
the multiplied remaining columns. Multiplication by a positive-cost
variable therefore preserves rejection in a greedy scan. Its
accepted monomials, with the compulsory lower diagram, form a downset.

This elementary Hall argument explains why downsets persist even
when determinant availability fails as in (7). It is not a
representation theorem for a support matroid by a polynomial quotient.

## Inputs and reproducibility

Minimum-basis exchange, alternating paths, determinant expansion,
Hall dependence and finite free base change are established tools.
[Reissig, condition (2) and Theorem 3.2](https://arxiv.org/pdf/1305.4277)
concern whole lower triangular block Toeplitz truncations with
independent coefficient-block entries. The additional repetitions
and selected submatrices here are outside that conclusion.
[Gluesing-Luerssen--Rosenthal--Smarandache, Lemma 3.8 and Theorem 3.9](https://math.nd.edu/assets/209984/01580796.pdf)
give symbolic superregularity and finite-field existence for scalar
lower triangular Toeplitz matrices. They do not supply the bivariate
selected-minor rule refuted here. Exact priority of this combined
coherent obstruction is unassessed; these comparisons do not establish
novelty or settle the strongest-current unrestricted bound.

Run from the repository root:

```bash
python3 -B scripts/run_coherent_convolution_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker reconstructs the support cone and an independent assignment
dual, expands both minors exactly by two methods, verifies original
commuting multiplication and finite-field ranks, and independently
builds power-substitution columns. Large lengths use the proved
residue-block decomposition, rather than enumerating large full fans.
The original unrestricted extremal objective and all its priority and
proof requirements remain open.
