# R14b — Hecke eigenvalue collisions for \(J_0(23)\): congruences and density heuristics

**Date:** 2026-11-23. **Extends:** R14 open item *"investigate collision-prime distributions with an explicitly declared analytic hypothesis and adequate sample size"*.
**Code:** [`experiments/collision_primes_r14b.py`](../experiments/collision_primes_r14b.py) · **Tests:** [`tests/test_collision_primes_r14b.py`](../tests/test_collision_primes_r14b.py) (10) · **Ledger:** [`R14b_Collision_Prime_Audit.json`](../results/R14b_Collision_Prime_Audit.json) · **Pre-registrations:** [M1](../results/R14b_PREREG_M1.md), [M2](../results/R14b_PREREG_M2.md).
> **R14c update (24 November 2026):** This dossier remains a historical R14b record. The [R14c addendum](R14c_addendum.md) provides a Sturm proof of R2 and corrects the conditional 3-adic factor to 213/160. The complete 2-adic and 3-adic images remain open in the repository's claim ledger.

## 0. Object

\(f, f^\sigma\in S_2(\Gamma_0(23))\), Hecke field \(K=\mathbb Q(\sqrt5)\), \(\mathcal O=\mathbb Z[\phi]\), \(\phi^2=\phi+1\). Write \(a_p(f)=u_p+v_p\phi\). A **collision prime** is \(p\ne2,23\) with \(a_p(f)=a_p(f^\sigma)\) ⇔ \(v_p=0\) ⇔ R14's \(\Delta_p=5v_p^2=0\) ⇔ \(T_p\in\mathbb Z\) inside \(\mathbb T=\mathcal O\).

**Coefficient engine (new, exact).** With \(\theta_1=\sum q^{x^2+xy+6y^2}\), \(\theta_2=\sum q^{2x^2+xy+3y^2}\) and \(B=(\theta_1-\theta_2)/2=\eta(z)\eta(23z)\):
\[
\boxed{\;v_p=[q^p]\,\eta(z)^2\eta(23z)^2,\qquad u_p=[q^p]\,\theta_2\,\eta(z)\eta(23z).\;}
\]
The combination was solved on \(M_2(\Gamma_0(23))=\langle\theta_1^2,\theta_1\theta_2,\theta_2^2\rangle\), then checked for exact equality against PARI `mfcoefs` for every \(n\le10^4\). It was also cross-checked against **independent** PARI `hyperellcharpoly` point counts at 468 primes, with 0 mismatches. FFT convolution to \(6\cdot10^7\) gives 3,562,113 primes, with rounding residual \(<5\cdot10^{-8}\).

## 1. Theorem A — collision primes are rigid mod 22 [P, conditional on cited inputs]

**Inputs (imported, not re-proved):**
- (I1) Mazur, *Eisenstein ideal*: for prime \(N\), \(\mathbb T/I\cong\mathbb Z/n\) with \(n=\mathrm{num}((N-1)/12)\). For \(N=23\), \(n=11\), and \(T_\ell-1-\ell\in I\) for every \(\ell\ne23\).
- (I2) \(\operatorname{End}(J_0(23))=\mathcal O\) (Hashimoto 2000; LMFDB; R14). Since \(T_2\) has minimal polynomial \(x^2+x-1\), we get \(\mathcal O=\mathbb Z[T_2]\subseteq\mathbb T\subseteq\mathcal O\), so \(\mathbb T=\mathcal O\).
- (I3) Eichler–Shimura: Frobenius on \(T_\ell J\), a free \(\mathcal O\otimes\mathbb Z_\ell\)-module of rank 2, has trace \(T_p\) and determinant \(p\).

**Theorem A.** If \(p\ne2,23\) is a collision prime and \(a=a_p(f)\in\mathbb Z\), then
\[
a\equiv p+1\pmod{11},\qquad a\equiv\begin{cases}1&p=2x^2+xy+3y^2\\0&\text{otherwise}\end{cases}\pmod 2 .
\]
So \(a\bmod 22\) is determined by \(p\).

*Proof, mod 11.* \(T_p=a\in\mathbb Z\), and \(a-1-p\in I\cap\mathbb Z=11\mathbb Z\) by (I1–I2).

*Proof, mod 2.* Complete the square in the model \(y^2+(x^3+x+1)y=-x^5\):
\[
(2y+x^3+x+1)^2=(x^3+x+1)^2-4x^5=(x^3-2x^2-3x-1)(x^3-2x^2+x-1).
\]
Both cubics are irreducible with discriminant \(-23\). The cubic field of discriminant \(-23\) is unique, so \(\mathbb Q(J[2])\) is the Hilbert class field \(H\) of \(\mathbb Q(\sqrt{-23})\), and \(\mathrm{Gal}(H/\mathbb Q)\cong S_3\). The Galois action on \(J[2]\cong\mathbb F_4^2\) is \(\mathbb F_4\)-linear with determinant the mod-2 cyclotomic character, which is 1. Its image therefore lies in \(SL_2(\mathbb F_4)\cong A_5\) and has order dividing 6. Every such element has order 1, 2 or 3, and so has trace in \(\mathbb F_2\); only elements of order 5 have trace \(\omega,\omega^2\). By (I3), \(T_p\bmod2\mathcal O\in\mathbb F_2\). The trace is 1 exactly when \(\mathrm{Frob}_p\) has order 3 in \(\mathrm{Gal}(H/\mathbb Q)\), i.e. when \(p\) splits in \(\mathbb Q(\sqrt{-23})\) into non-principal primes, i.e. \(p=2x^2+xy+3y^2\). ∎

**Corollaries.**
- (L1) \(v_p\) is even for every \(p\ne2,23\). A second, elementary proof: \(v_p=\sum_{m=1}^{p-1}b_mb_{p-m}\), and for odd \(p\) the terms pair off as \(m\leftrightarrow p-m\).
- (L2) \(a_p(f)\equiv p+1\pmod{\mathfrak P_{11}}\), with \(\mathfrak P_{11}=(11,\phi-4)\).
- (C) \(a_p(f)=0\Rightarrow p\equiv10\pmod{11}\). The only such primes up to \(6\cdot10^7\) are **43, 109, 1033**.

**Lemma R1 [P, given I3].** If \(p=x^2+xy+6y^2\) (\(\mathrm{Frob}_p\) trivial on \(J[2]\)), then \(T_p\equiv1+p\pmod{4\mathcal O}\).

*Proof.* \(\mathrm{Frob}=1+2X\) with \(X\in M_2(\mathcal O_2)\). Then \(T_p=2+2\operatorname{tr}X\) and \(p=1+2\operatorname{tr}X+4\det X\), so \(1+p-T_p=4\det X\). ∎

**Conjecture R2 [M, 0 exceptions].** If \(p=2x^2+xy+3y^2\), then \(v_p\equiv p-1\pmod 4\). This is **not** a formal consequence of the mod-2 data: writing \(\mathrm{Frob}=A+2\phi B\) leaves \(\operatorname{tr}B\bmod2\) free. It is a constraint on the 2-adic image itself, plausibly derivable from Schoof's analysis of \(J_0(23)[2]\). For inert \(p\), \(v_p\bmod4\) is uniform (0.4998 / 0.5002).

**Verification over all 3,562,113 primes \(\le6\cdot10^7\):** Theorem A (2,386 collisions), L1, L2, R1 (593,491 cases) and R2 (1,187,172 cases) all hold with zero exceptions. For R2 this is evidence only.

## 2. Mod-3 trace frequencies and an icosahedral model [M; image type SOURCE-REPORTED, not independently verified]

For \(p\equiv1\pmod3\), \(T_p\bmod3\in\mathbb F_9\) never equals \(\pm\phi^2\). For \(p\equiv2\pmod 3\) it never equals \(\pm1\). The observed frequencies are
- 0.25 for trace 0;
- 0.175 for each of \(\pm1\) (when \(p\equiv1\)) or \(\pm\phi^2\) (when \(p\equiv2\));
- 0.10 for each of the other four values.

The measured frequencies match the trace law computed for the explicitly constructed model group \(\langle i\rangle\cdot 2.A_5\subset GL_2(\mathbb F_9)\) with \(i=\phi^2\), of order 240. **This comparison does not determine the residual Galois image of the level-23 newform.** Literature: Billerey–Dieulefait, *Explicit large image theorems for modular forms* (J. LMS 89 (2014); arXiv:1210.5428), §5.3, treat this exact newform (weight 2, level 23, orbit 23.4a in Sage) at the prime above 3. They state that the projective mod-3 image is \(A_5\), with kernel field the splitting field of \(X^5+3X^3+6X^2+9\), and attribute the example to Ribet and to Kiming–Verrill (J. Number Theory 110 (2005), Ex. 3.2). An intermediate note in this dossier said the citation belonged to a weight-6, level-81 form. That was a misreading: §5.3 explicitly concerns \(S_2(\Gamma_0(23))\) and \(\ell=3\). We have not re-derived the image theoretically, and the 240-element model is our refinement, not theirs. **Independent consistency check** ([`scripts/mod3_quintic_check.py`](../scripts/mod3_quintic_check.py)): for all 17,981 primes \(p\le2\cdot10^5\) not dividing \(\operatorname{disc}=3^8\cdot23^2\), the cycle type of \(\mathrm{Frob}_p\) on the roots of \(X^5+3X^3+6X^2+9\) agrees with the projective class read off from \(r=t^2/\det\in\mathbb F_9\), \(t=T_p\bmod 3\): 5-cycle ⇔ \(r\notin\mathbb F_3\); (2,2,1) ⇔ \(t=0\); identity or 3-cycle ⇔ \(r=1\). There are **0 mismatches**, and the class counts 296 / 4,495 / 6,004 / 7,186 match the \(A_5\) proportions 1:15:20:24 out of 60. This is strong evidence, not a proof. Under the present *model assumption*, the finite-level local factors are \(F_3^{(1)}=51/40\) and \(F_3^{(2)}=53/40\), compared with the generic model value \(1.0969\). Recorded scaled divisibility frequencies: 1.2750, 1.3258, 1.3328, 1.3333 for \(k=1,2,3,4\), respectively. The extrapolation to a limiting \(F_3\) remains heuristic.

## 3. Density: a pre-registered failure, then a held-out success

Heuristic (Murty; Van Hirtum, arXiv:1508.07206):
\[
\Pr[v_p=0]\approx \tfrac{8\sqrt5}{3\pi^2}\hat F\,p^{-1/2},\qquad \#\{p\le x\}\approx\tfrac{8\sqrt5}{3\pi^2}\hat F\sum_{p\le x}p^{-1/2},
\]
with \(\hat F=\prod_\ell F_\ell\) and \(F_\ell=\lim_k\ell^k\Pr[\ell^k\mid v_p]\). Van Hirtum treated exceptional local factors empirically. Here the factors at several primes are computed from proposed image models, while the M2 factor at 2 remains **empirical**.

| ℓ | model | \(F_\ell\) | status |
|---|---|---|---|
| 2 | mod-2 image in \(GL_2(\mathbb F_2)\) (Theorem A) | 22/9 under full preimage; **2.554 empirical** | full-preimage model refuted at \(k\ge3\) (R2 structure) |
| 3 | binary icosahedral | 1.332 (model \(k=2\): 53/40) | fits data |
| 5 | ramified, full | 1.0347 | fits (1.0304, 1.0259 ± 0.01) |
| 11 | Eisenstein: \(\mathfrak P\)-component ≡ diag(1,*) mod 11 | 1.0992 | fits (1.0948, 1.1083 ± 0.01); a Borel model would give 1.0916 |

**M1** was frozen before any scan beyond \(p\le101\) (sha256 `f3d823ab…`) and used the generic \(F_3\). It is **refuted**: at \(2\cdot10^7\) it predicted 1,187.9 against 1,477 observed (z = +8.4). Its E-tests diagnosed exactly which \(\ell\) failed (3 and 2) and which passed (5, 11).

**M2** is exploratory and post hoc: the icosahedral \(F_3\) plus the empirical \(F_2\), which comes from local frequencies and not from collision counts. It was frozen before \((2\cdot10^7,6\cdot10^7]\) was computed (sha256 `c71a78e4…`):

| held-out range | observed | M2 prediction | z | M1 z |
|---|---|---|---|---|
| (2e7, 4e7] | 504 | 524.6 | −0.90 | +4.45 |
| (4e7, 6e7] | 405 | 390.3 | +0.75 | +5.55 |
| (2e7, 6e7] | 909 | 914.9 | **−0.19** | +7.00 |

\(\hat F_{M2}=4.026\), \(c=16\sqrt5\hat F/(3\pi^2)=4.864\). M2 **over-predicts at small x** (ratio 0.75 at \(10^5\), 0.90 at \(10^6\), 0.99 at \(10^7\)). The \(p^{-1/2}\) model ignores finite-\(p\) discreteness and edge effects; this is open.

**Comparison with Van Hirtum's finite-range normalization [INFERRED, not a new resolution].** For level 23 Van Hirtum reports a predicted \(\hat c=4.982\) and a fitted \(\tilde c=5.490\) (SOURCE-REPORTED). Our asymptotic \(c=4.864\) lies within 2.4% of his prediction. Because \(\sum_{p\le x}p^{-1/2}\big/(2\sqrt x/\log x)\) is still 1.146 at \(10^8\), a least-squares fit of \(c\sqrt x/\log x\) to M2-expected counts on \(x\le10^8\) returns **5.60**. That matches his fitted 5.49 to within 2%. The numerical comparison is consistent with the already-known finite-range correction in replacing \(\sum p^{-1/2}\) by \(2\sqrt x/\log x\); it does not independently rule out omitted local corrections. Caveat: Van Hirtum's exact fitting window was not re-inspected.

## 4. Claim ledger

| Claim | Status | Scope / dependency |
|---|---|---|
| \(v_p=[q^p](\eta\eta_{23})^2\), \(u_p=[q^p]\theta_2\eta\eta_{23}\) | ARTIFACT-VERIFIED, exact to \(n\le10^4\) vs PARI | identity in finite-dimensional \(M_2\) (Sturm bound 4) ⇒ PROVED |
| Theorem A, L1, L2, C | PROVED | conditional on Mazur, End ring, Eichler–Shimura (imported) |
| Lemma R1 | PROVED | Eichler–Shimura + free rank 2 |
| Conjecture R2 | MEASURED, 1.19M cases, 0 exceptions | open proof |
| mod-3 trace law matches \(\langle i\rangle 2.A_5\) model | MEASURED (frequency match); projective image \(A_5\) SOURCE-REPORTED (Billerey–Dieulefait §5.3; Ribet; Kiming–Verrill Ex. 3.2) | full 3-adic image still a model hypothesis |
| \(F_3, F_5, F_{11}\) values | DERIVED from image models | the model images are HYPOTHESES (full preimages) |
| \(F_2=2.554\) | MEASURED, empirical | 2-adic image not determined |
| M1 refuted; M2 held-out fit | MEASURED, pre-registered | M2 is post hoc; one held-out window |
| Van Hirtum finite-range normalization comparison | INFERRED, not new | known approximation effect; exact fit not re-checked |
| Novelty | UNKNOWN | brief search found no explicit statement of Theorem A or R2; Theorem A is a short corollary of known inputs and may be folklore |

**Not claimed:** infinitude of collisions (open, as for Lang–Trotter), any Lean build, or a new asymptotic theorem. Here a *collision* means equality of the Hecke eigenvalues of two Galois-conjugate newforms; it is not a term for any additional physical or musical structure.

## 5. Next discriminators

1. Prove or kill R2: determine the 2-adic image of \(J_0(23)\) (Schoof's \(J_0(23)[2]\) extension structure). The test is whether the full preimage of that image reproduces \(F_2=2.554\) exactly.
2. A finite-\(p\) correction to the archimedean factor: replace \(\int\mu^2\) by the exact lattice sum \(\sqrt5\sum_a\mu_p(a)^2\), restricted to the residue classes forced by Theorem A. Test it on \(p\le10^6\), where M2 is 10% high.
3. Decide whether the 3-adic image is the full preimage of \(\langle i\rangle 2.A_5\). The mod-3 projective image passed the quintic cycle-type check (§2); a theoretical re-derivation is still open.
4. Port Theorem A's two arithmetic steps to Lean: \(I\cap\mathbb Z=11\mathbb Z\) given \(\mathbb T=\mathcal O\), and the "order divides 6 in \(SL_2(\mathbb F_4)\) ⇒ trace in \(\mathbb F_2\)" lemma. Both are finite or algebraic.
