# J0-23-collision-primes

When the two Hecke voices of \(J_0(23)\) coincide: exact congruences (\(a_p \bmod 22\) is determined by \(p\)), an icosahedral mod-3 law, and a pre-registered Lang–Trotter density test with held-out confirmation. 3,562,113 primes up to \(6\cdot10^7\).

Spun out of the [Bolyai Research Atlas](https://github.com/Dojo-1/Bolyai-Research-Atlas), round R14 → R14b. Full mathematical dossier with proofs and claim ledger: [`docs/R14b_dossier.md`](docs/R14b_dossier.md).

## Object

\(f, f^\sigma \in S_2(\Gamma_0(23))\) have Hecke field \(\mathbb Q(\sqrt5)\). Write \(a_p(f)=u_p+v_p\phi\) with \(\phi^2=\phi+1\). A **collision prime** is \(p\ne 2,23\) with \(a_p(f)=a_p(f^\sigma)\), i.e. \(v_p=0\).

Exact engine:
\[
v_p=[q^p]\,\eta(z)^2\eta(23z)^2,\qquad u_p=[q^p]\,\theta_{2,1,3}\,\eta(z)\eta(23z).
\]
It is checked for exact equality against PARI `mfcoefs` (\(n\le10^4\)) and against independent PARI point counts (468 primes, 0 mismatches).

## Results

| Claim | Status |
|---|---|
| **Theorem A.** At a collision, \(a_p\equiv p+1 \pmod{11}\), and \(a_p\) is odd iff \(p=2x^2+xy+3y^2\). Hence \(a_p \bmod 22\) is fixed by \(p\), and \(a_p=0\Rightarrow p\equiv10 \pmod{11}\). | Proved, given Mazur's Eisenstein ideal, \(\operatorname{End}J_0(23)=\mathbb Z[\phi]\) and Eichler–Shimura. 0 exceptions in 2,386 collisions. |
| **L1.** \(v_p\) is even for all \(p\ne2,23\). | Proved (two independent proofs). |
| **R1.** If \(p=x^2+xy+6y^2\), then \(T_p\equiv 1+p \pmod{4\mathcal O}\). | Proved. 593,491 cases. |
| **R2.** If \(p=2x^2+xy+3y^2\), then \(v_p\equiv p-1 \pmod 4\). | **Conjecture.** 1,187,172 cases, 0 exceptions. |
| The mod-3 trace law equals that of \(\langle i\rangle\cdot 2.A_5\subset GL_2(\mathbb F_9)\) (icosahedral). | Measured, 14/14 frequencies. The image type is in the literature. |
| Density model **M1** (pre-registered, generic factor at 3). | **Refuted**, z = +8.4. |
| Density model **M2** (icosahedral \(F_3\), empirical \(F_2\)), frozen before held-out data. | Held-out window \((2\cdot10^7,6\cdot10^7]\): 909 observed vs 914.9 predicted (z = −0.19). |
| Van Hirtum's level-23 gap (fitted 5.490 vs predicted 4.982). | Inferred to be the second-order term of \(\sum p^{-1/2}\); M2 implies a fitted value of 5.60. |

**Not claimed:** infinitude of collisions, a Lean proof, or novelty. Theorem A is a short corollary of known inputs and may be folklore.

## Layout

```
experiments/collision_primes_r14b.py   self-contained module: coefficient engine, law audit, exact local factors, 2.A5 construction, density constants
tests/test_collision_primes_r14b.py    10 regression tests (test_10 needs cypari2; skipped otherwise)
results/R14b_PREREG_M1.md              frozen before any scan beyond p<=101   sha256 f3d823ab48929614090d8d97f3d8a8aca03b1614b503619eb9196192dc43ed49
results/R14b_PREREG_M2.md              frozen before computing (2e7,6e7]       sha256 c71a78e4594e252d00630b3195d28a3e9ffd54ab90363b525ea6ca1517862f12
results/*.json                         audit ledger, M1 + local-frequency tests, held-out scoring, mod-4 laws
data/ap_lean_60000000.npz              (p, u_p, v_p) for all primes 3<=p<=6e7, p!=23
scripts/                               provenance: the exact working scripts as run (paths are those of the original working directory)
docs/R14b_dossier.md                   proofs, claim ledger, next discriminators
```

## Reproduce

```bash
pip install -r requirements.txt
python3 -m unittest discover -s tests -v                 # about 2 s
python3 -m experiments.collision_primes_r14b 2000000     # law audit + M1/M2 predictions at 2e6
sha256sum results/R14b_PREREG_M*.md                      # must match the hashes above
```

A note on `scripts/`: these are the scripts exactly as they were run, kept unedited for provenance. They expect intermediate files (`ap_20000000.npz`, `primes_6e7.npy`, `G3.pkl`) in the working directory. Run `data.py 2e7` first to regenerate them; the 6e7 run needs about 5 GB of RAM.

No CI is configured in this repository by design.
