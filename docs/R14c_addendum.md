# R14c — Doubt the Machine: arithmetic addendum

**24 November 2026.** This is an AI-assisted, correction-preserving follow-up to the November 23 manuscript. The original paper, data, preregistrations, and predictions are retained as historical evidence. See the [Doubt the Machine protocol](https://github.com/Dojo-1/doubt-the-machine).

## 1. R2 via a finite modular-form congruence

Assume the coefficient identities of Proposition 2.1 and the usual weight-one theta/eta modularity. Put

$$
B(z)=\frac{\theta_1(z)-\theta_2(z)}2=\eta(z)\eta(23z),
\quad V=B^2,\quad
W(z)=\theta_1(2z)B(2z),\quad E(z)=\theta_1(2z)B(z).
$$

Here \(\theta_1=\theta_{1,1,6}\), \(\theta_2=\theta_{2,1,3}\). These are weight-one forms of character \(\chi_{-23}\). Consequently \(V,W,E\) are weight-two forms of trivial character on \(\Gamma_0(46)\). Their coefficients are integral. The involution \((x,y)\mapsto(-x,-y)\) gives \(\theta_1\equiv1\pmod2\) as a formal series, and \(B^2\equiv B(2z)\pmod2\). Thus \(G=(V-W)/2\) is a holomorphic integral-coefficient weight-two form on \(\Gamma_0(46)\). For odd \(n\), \(G_n=v_n/2\) and \(E_n\equiv B_n\pmod2\).

For a series \(F=\sum_n F_nq^n\), let \(P_S F=\sum_{n\bmod92\in S}F_nq^n\). Define

$$
S=\{n\bmod92:n\text{ odd},\ (n/23)=1\},\qquad
S_3=\{n\in S:n\equiv3\pmod4\}.
$$

The series \(D=P_S G-P_{S_3}E\) has integral coefficients.

**Conservative level justification.** For an integer \(m\), a residue projection is a linear combination of the translates \(F(z+j/m)\), with coefficients the discrete Fourier transform of the residue mask. Let \(F\) have trivial character at level \(L\), and let \(\gamma=\left(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right)\in\Gamma_0(Lm^2)\). With \(T_j=\left(\begin{smallmatrix}1&j/m\\0&1\end{smallmatrix}\right)\) choose \(j'\equiv jd^2\pmod m\). Then \(T_j\gamma T_{j'}^{-1}\) belongs to \(\Gamma_0(L)\): integrality follows from \(m^2\mid c\) and \(ad\equiv1\bmod m\), and the lower-left entry is still \(c\). Hence the translates are permuted under the weight-two slash action. The masks \(S,S_3\) are invariant under multiplication by a unit square modulo 92, as are their Fourier transforms. The projected forms are therefore modular with trivial character at level \(L m^2\); rational translations also preserve holomorphy at cusps. This argument avoids assumptions about the minimal level of the twists.

Take \(L=46\), \(m=92\). The level \(389344\) has index \(389344(1+1/2)(1+1/23)=609408\). Sturm's bound at weight two is \(101568\). The exact [conservative certificate](../review/R14c/conservative_check.py) tests all coefficients modulo 2 through that bound, including the zero constant term. It uses the pentagonal-number theorem and carry-free integer packing for the square of \(B\bmod4\); there is no FFT or floating-point rounding. The [raw log](../review/R14c/conservative.log) reports zero nonzero coefficients. Hence Sturm's theorem gives \(D\equiv0\bmod2\).

For an odd prime \(p\ne23\) represented by the nonprincipal class \(2x^2+xy+3y^2\), \(p\in S\) and \(B_p=-1\) (the two nonprincipal reduced classes are inverses, each has two units; equivalently use the weight-one class-character trace). Thus

$$
v_p/2\equiv {\bf1}_{p\equiv3\ (4)}\pmod2,\qquad
\boxed{v_p\equiv p-1\pmod4}.
$$

The congruence holds with either conjugate choice of eigenform since \(v\mapsto-v\) leaves residues 0 and 2 fixed. In particular, a nonprincipal collision prime must satisfy \(p\equiv1\pmod4\).

The smaller bound 2208 and the theta-series implementation are additional cross-checks; the argument above needs only the conservative level. The full 2-adic image, exact \(F_2\), and novelty of this congruence remain separate questions.

## 2. Exact local factor at 3, conditional on the image model

Write \(\ell=3\) and \(G_0=\langle i\rangle\cdot 2.A_5\) of order 240 in \(GL_2(\mathbf F_9)\) with rational determinant. The explicit enumeration in [verify.py](../review/R14c/verify.py) gives 102 rational-trace elements, of which two are scalar (\(\pm I\)). The following is a count for the **full preimage model**, not an identification of the actual Galois image.

For a nonscalar rational-trace residue class, the differentials of the nonrational trace and determinant coordinates are independent over \(\mathbf F_3\). Indeed the nonrational coordinate is a nonzero linear functional on \(\mathbf F_9\), and a linear dependence between these two functionals on the matrix space would force the adjugate matrix to be scalar. Thus each subsequent rational-trace lift, conditional on rational determinant, has probability \(1/3\); these 100 classes contribute \(3(100/240)\) to the scaled density.

For a scalar residue class impose rational trace first. Since 2 is invertible, write \(g=tI+3X\) with \(t\) rational and \(X\) traceless. The determinant is \(t^2-9(a^2+bc)\), where \(a,b,c\in\mathcal O_3\). In the basis \(1,\phi\), the nonrational coordinate of \(a^2+bc\) is
$$
(2a_0+a_1)a_1+b_0c_1+b_1(c_0+c_1).
$$
By an invertible linear change of variables this is the dot product of two independent uniform vectors of length three over \(\mathbb Z/3^m\). Put
$$
H_m=3^m\Pr(U\cdot V\equiv0\bmod3^m)
=\sum_{j=0}^{m-1}\frac{(3^3-1)3^j}{3^{3(j+1)}}+3^{-2m}.
$$
This follows by conditioning on the minimum 3-adic valuation \(j\) of the three entries of \(U\), treating \(U=0\bmod3^m\) separately. In a scalar residue class the determinant-coordinate condition alone is smooth, so its probability in unrestricted lifts is \(3^{-(k-1)}\), the same as that of rational trace alone. Given rational trace, the determinant condition is the dot-product equation modulo \(3^{k-2}\). Consequently the scalar classes contribute \(9(2/240)H_{k-2}\) to \(F_3^{(k)}\) for \(k\ge2\).

$$
F_3^{(k)}=\frac54+\frac3{40}H_{k-2},\qquad
\lim_{m\to\infty}H_m=\frac{13}{12},\qquad
\boxed{F_3=\frac{213}{160}}.
$$

At \(k=1,2\) this agrees with the manuscript's exact \(51/40,53/40\). The \(4/3\) extrapolation in the original open problems is retired **for this model**. The full-preimage assumption, and any joint exceptional-factor entanglement, remain to be established independently. The conditional factor at 3 does not certify the global asymptotic constant.

## 3. Evidence and chronology

| Claim | R14c status |
|---|---|
| Mod-4 R2 law | Mathematical deduction from coefficient identities, modularity, Sturm's theorem and the exact finite certificate; external specialist verification pending |
| \(F_3=213/160\) | Exact full-preimage model limit, not an unconditional Galois image result |
| \(F_2=2.554\) | Empirical input, not an exact limit |
| M1 / M2 | Historical finite-range rejection / exploratory held-out agreement, unmodified |
| Infinitude or asymptotic of collision primes | Open |

Both preregistrations and the large dataset first appear in the same public commit [c25c030](https://github.com/Dojo-1/J0-23-collision-primes/commit/c25c030f8528ee15b2f110ba4029ca8c1cad9b9f) on November 23. The author-reported freeze is not independently timestamped by that Git history. A matching hash proves content integrity, not temporal precedence. The new factor produces a retrospective corrected M2 prediction of 914.4 rather than 914.9 on the historical held-out range; the 909 observation and the original M2 prediction are retained.

## 4. Reproduction and limits

```bash
python3 review/R14c/replay.py
```

See [exact coefficient and group results](../review/R14c/results.json), [point counts](../review/R14c/pointcounts.json), and [forecast arithmetic](../review/R14c/predictions.json). The commit-time cross-check independently recomputed the conservative projection (zero violations, 24288 selected odd split indices), all 11 recorded curve point counts over the base and quadratic extension fields, and the prime sieve and six forecast weight sums (3,562,113 good odd primes, agreement to \(6\cdot10^{-14}\)). It did not regenerate the entire historical coefficient dataset or establish its chronology. No CI execution or Lean certification is claimed. The next proof obligation is the actual 3-adic image, independently of fitting the factor; the next experimental test should freeze an independently timestamped forecast on a genuinely new range.

