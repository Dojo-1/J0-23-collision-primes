# R14d: the full 3-adic image for J₀(23)

25 November 2026 · Addendum to [R14c](R14c_addendum.md). This note leaves the November 23 manuscript, preregistrations, and original forecasts intact.

**Result.** Let `O = Z₃[φ]`, where `φ²−φ−1=0`, and let `ρ₃` be the level-23 weight-two newform representation on a stable `O`-lattice. Its determinant is the 3-adic cyclotomic character. Write `H = im(ρ₃ mod 3)` and let `K` be its determinant-one subgroup. Then

```text
im(ρ₃) = { M ∈ GL₂(O) : det M ∈ Z₃× and M mod 3 ∈ H },
H = ⟨iI⟩ K,  i² = −1 in F₉,  K ≅ 2.A₅.
```

Here the equality uses the **actual** residual image on the chosen lattice, not an assumption that two particular embeddings are conjugate. The finite check below covers all residual placements relevant to the proof. The resulting local factor is **F₃ = 213/160** in the normalization of R14c.

## 1. Residual input

[Kiming–Verrill, Example 3.2](https://arxiv.org/pdf/math/0311082#page=8) identifies the candidate projective field as the splitting field of `x⁵+3x³+6x²+9`. Their example is stated relative to complex weight-one modularity. For this application one can instead use their Theorem 1(b), which starts from a modular residual representation with this projective field. The required bridge is as follows.

Take a finite-image two-dimensional Tate lift of the complex icosahedral projective representation, as in Kiming–Verrill §4, and reduce a stable lattice at 3. Reduction defines a homomorphism on the projective image: an element acting by a scalar before reduction still acts by a scalar after reduction. An element of projective order five retains its eigenvalue ratio of order five, so this homomorphism is nontrivial and is injective by simplicity of A₅. The representation is odd: complex conjugation is projectively nontrivial, and its lift has eigenvalues +1 and −1. It is absolutely irreducible because a finite triangular projective group is solvable. Thus Serre modularity ([Khare–Wintenberger](https://www.math.ucla.edu/~shekhar/papers/results.pdf)) supplies the residual modularity required in Theorem 1(b). Apply that theorem with the local data and characters of Example 3.2: weight two, level 23, trivial character. The single Galois orbit of level-23 weight-two forms identifies the required residual projective image as A₅. See also [Billerey–Dieulefait, §5.3](https://arxiv.org/pdf/1210.5428#page=24).

The projective field and its local computations at 3 and 23 are imported from Kiming–Verrill, not recomputed by the certificate. Complex weight-one modularity is not used.

The determinant-one subgroup K maps onto A₅, since the projective quotient by its image has order dividing two. A determinant-one lift of a nontrivial projective involution is traceless and squares to −I. Hence K contains −I and is the full inverse image of the projective A₅ in SL₂(F₉). For any determinant −1 element of H, a member of K has the same projective class, so their quotient is scalar with square −1. The cyclotomic determinant attains −1 modulo 3, yielding the claimed description of H.

## 2. Frobenius at 13

For the [modular curve model](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6DF1C7B794ECADB769EFF24D4E9796F7/S1461157015000157a.pdf/hyperelliptic-modular-curves-dollarx0ndollar-and-isogenies-of-elliptic-curves-over-quadratic-fields.pdf) `y²−(x³+x+1)y=−2x⁵−3x²+2x−2`, exact counts give 8 and 204 points over F₁₃ and F₁₆₉. The certificate counts both points at infinity using `u=1/x, v=y/x³`, checks good reduction, and obtains the same counts independently on the [curve with isogenous Jacobian](https://www.lmfdb.org/Genus2Curve/Q/529/a/529/1).

The counts give `a₁₃+a₁₃′=6` and `a₁₃²+(a₁₃′)²=18`, so both traces are 3. If `g=ρ₃(Frob₁₃)` for arithmetic Frobenius, Cayley–Hamilton gives `g²=3g−13I`, and thus `g⁴=−51g+52I`. In particular

```text
(g⁴−I)/3 ≡ ḡ+2I (mod 3),   tr(ḡ)=0,   det(ḡ)=1.
```

The trace of this first-order direction is 1 in F₃.

## 3. From modulo 9 to every level

In the group of matrices over O with rational 3-adic determinant, put `Uₙ={M : M ≡ I mod 3ⁿ}` for n≥1. The quotient `Uₙ/Uₙ₊₁` is identified via `I+3ⁿX ↦ X mod 3` with the seven-dimensional F₃-space `W={X ∈ M₂(F₉) : tr X ∈ F₃}`. The determinant expansion gives containment; for surjectivity, an initial lift has determinant congruent to a unit in Z₃ modulo 3ⁿ⁺¹. Correct its determinant to that unit by multiplying by a diagonal matrix congruent to I modulo 3ⁿ⁺¹.

Conjugating g⁴ by lifts from K and subtracting its first-order direction gives the vectors `h ḡ h⁻¹−ḡ` for h∈K. These span the six-dimensional traceless subspace. To see why the finite check is exhaustive, every trace-zero determinant-one matrix over F₉ has the same separable characteristic polynomial `T²+1` and is GL₂(F₉)-conjugate to `A=[[0,−1],[1,0]]`. Conjugation preserves rank and rational trace. In A₅, any involution can be completed to generators of respective orders 2 and 3 whose product has order 5: such a pair exists, and all involutions are conjugate. Thus, after fixing A as the first lift, it suffices to enumerate determinant-one B with projective orders `ord(B)=3` and `ord(AB)=5`. The certificate finds 64 choices, generating four distinct order-120 groups containing A, and checks the same ranks (6,7) in each. The seventh direction is `ḡ+2I`, whose trace is nonzero.

For n≥1, cubing sends a first-order direction at level n to the same direction at level n+1: `(I+3ⁿX)³ ≡ I+3ⁿ⁺¹X mod 3ⁿ⁺²`. Hence the Galois image surjects onto each graded quotient. Successive corrections approximate any element of U₁, and the image is closed by compactness. It therefore contains U₁. Together with the definition of H and the cyclotomic determinant, this proves the full-preimage statement.

## 4. Local factor and evidence boundary

The certificate counts 102 rational-trace elements in the 240-element H, exactly two of them scalar, for every placement examined. For each of the 100 nonscalar rational-trace residues, the nonrational trace condition cuts an extra factor of 1/3 at each higher level: its differential and that of the rational-determinant constraint are independent. For the two scalar residues, write a rational-trace lift as `tI+3X` with t rational and X traceless. The nonrational determinant condition reduces to the dot-product problem of [R14c §2](R14c_addendum.md#2-exact-local-factor-at-3-conditional-on-the-image-model). For independent uniform three-vectors over Z/3ᵐ, the scaled zero-dot-product probability tends to 13/12, by summing over the minimum valuation of the first vector. Thus the scaled local limit is

```text
3·100/240 + 9·(2/240)·13/12 = 213/160.
```

This removes the full-image hypothesis for the **single local factor at 3**. It does not determine the full 2-adic image, settle joint entanglement of local factors, or prove infinitude or a global collision asymptotic. The original 914.9 forecast and 909 observation stand as recorded; the 914.4 figure is retrospective.

## Reproduce

From the repository root, run `python3 review/R14d/verify.py` and compare with [expected.json](../review/R14d/expected.json). The certificate is standard-library Python: exact point enumeration, exact finite matrices, and row reduction over F₃. The p-adic group argument and residual arithmetic inputs are mathematical arguments/citations, not Lean formalizations.
