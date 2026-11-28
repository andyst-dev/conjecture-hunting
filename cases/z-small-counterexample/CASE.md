# Smaller counterexample to real-rootedness of the matroid Z-polynomial (after Eberhardt, arXiv:2609.20903 v1)

**Status:** no counterexample in the classes searched (listed under "Search coverage"); the only non-real-rooted case found is Eberhardt's own G24.
**Date:** 2026-11-28 · **Tools:** Python 3, python-flint, sympy (exact root counts), numpy, a small C enumerator · **AI-assisted (Claude), unrefereed**

## Statement
- Matroid M of rank r on E, lattice of flats L(M). The KL polynomial P_M is that of Elias–Proudfoot–Wakefield.
- The Z-polynomial is Z_M(t) = Σ_{F ∈ L(M)} t^{rk F} P_{M/F}(t) (Proudfoot–Xu–Young 2018). It is palindromic of degree r, and deg P_M < r/2.
- **Conjecture (PXY18, Conj. 5.1):** Z_M has only real (negative) roots for every matroid M.
- **Status:** refuted by J. N. Eberhardt, *Condensed configurations and valuative matroid invariants*, arXiv:2609.20903 v1 (17 Sep 2026). The matroid G24 of the extended binary Golay code, of rank 12 on 24 elements, has non-real-rooted Z and P. The paper calls it "the first counterexample in the literature".
- **Question pursued here:** is there a smaller counterexample, with lower rank or fewer elements? Eberhardt computes no other matroid with this in view and says nothing about minimality.

Conventions: polynomials are integer coefficient lists from low to high degree. "Real-rooted" is decided **exactly** by a sympy real-root count against the degree. Floating point is used only for screening, and every flagged case is re-checked exactly.

## Method
- **P without Möbius functions.** Z_M palindromic and deg P_M < r/2 give [tᶦ]P_M = [t^{r−i}]R − [tᶦ]R, where R = Σ_{F ≠ cl(∅)} t^{rk F} P_{M/F}. This uses integer arithmetic only (`zfast.py`).
- **Nested (Schubert) matroids** are given by their chain of cyclic flats, with labels (size, rank). They are computed on *profiles*, meaning intersection sizes with the layers of the chain. This is exact because Aut = ∏ S_{layer} acts transitively on flats with a given profile and on pairs F ≤ X with given profiles.
- **Eberhardt's Theorem 3.2:**
  - Φ(M) = Σ over block chains of ∏ (−(A⁻¹)_{B_{i−1}B_i}) · Φ(nested matroid with those labels), for any valuative invariant Φ.
  - A is the condensed incidence matrix of the lattice of cyclic flats.
  - `condense.py` builds A by **partition refinement**: it starts from the (size, rank) labels and splits blocks until every Y in a block C contains the same number of members of each block B.
  - It then sums over chains by DP on label sequences.
- **Cyclic flats of binary matroids** come from `cycflats.c`. Circuits are the dual codewords w with rank(supp w) = |supp w| − 1. The program runs a BFS from cl(∅), joining one circuit at a time.

## Prior status
- Eberhardt computes P and Z of G24 and the rank-3 Steiner 2-design matroids. In rank 3, Z = 1 + b(t + t²) + t³ is always real-rooted.
- Real-rootedness is proved for sparse paving matroids, which include all uniform matroids (arXiv:2609.07636), and for some other families (PXY18; fan, wheel and whirl matroids).
- The predecessor case [`z-polynomial-rank4`](../z-polynomial-rank4/CASE.md) derived the rank-4 criterion 4·a₂ ≤ W₃² + 8 and found no counterexample in random rank-4 and rank-5 linear matroids.

## Reproduction of the paper
| Object | Paper value | Our value |
|---|---|---|
| G24: P (from the paper's condensed matrix) | 1, 3311, 162656, 2419692, 6722716, 2923921 | identical |
| G24: Z (from the paper's condensed matrix) | 1, 3335, 205436, 4285820, 31050690, 95326605, 137174576, … (palindromic) | identical |
| G24 from scratch (Golay code → cycflats → refinement → Thm 3.2): number of cyclic flats | 39,516 | 39,516 |
| G24 from scratch: blocks and labels | 6: (0,0), (8,7), (12,10), (12,11), (16,11), (24,12) | identical |
| G24 from scratch: condensed matrix A | rows 1 1 1 1 1 1 / 0 1 3 0 30 759 / 0 0 1 0 140 35420 / 0 0 0 1 0 2576 / 0 0 0 0 1 759 | identical |
| G24 from scratch: P and Z | as above, both non-real-rooted | identical; exact root count confirms non-real-rooted |
| Golay weight distribution | 1, 759, 2576, 759, 1 (weights 0, 8, 12, 16, 24) | identical |
| Rank-3 Steiner formulas: Fano S(2,3,7), b = 7 | P = 1 + (b−v)t = 1, Z = 1 + 7t + 7t² + t³ | identical |
| Rank-3 Steiner formulas: AG(2,3) = S(2,3,9), b = 12 | P = 1 + 3t, Z = 1 + 12t + 12t² + t³ | identical |

Not recomputed: the number of flats of G24 (2,047,118). We never enumerate all flats.

**Cross-checks of the fast paths** (`test_zfast.py`, all pass), against `zpoly_ref.py`. That reference is the EPW recursion with Möbius functions, validated in the predecessor case.
- `kl_z_from_flats`: 60 random linear matroids over GF(2), GF(3) and GF(5).
- `nested_kl_z`: 60 random nested matroids, built explicitly from their rank function.
- Schubert expansion: 25 random paving matroids.
- The full condensed pipeline (`condense.py` + `cycflats.c`): 40 random binary matroids.

## Search coverage
| Class | Range | Method | Result |
|---|---|---|---|
| A. Paving matroids whose stressed hyperplanes all have one size h (Z = Z(U) + λ(Z(S_h) − Z(U))) | rank 4–10, n ≤ 60, every h, every λ ≤ the constant-weight-code bound U(n, h, r−2); integer λ when U ≤ 3000, a 600-point grid otherwise | float screen + exact check | no realisable failure known: see Findings 2 |
| B. Golay code shortened at s and punctured at p coordinates, G(s,p) | all s + p ≤ 5 (one code each, since M24 is 5-transitive): rank 7–12, n = 19–23 | exact (condensed expansion) | all 20 real-rooted; only G(0,0) = G24 fails |
| C. Extended QR [18,9,6], RM(1,4) [16,5,8], RM(2,4) [16,11,4], extended Hamming [8,4,4], and their duals; Golay as control | — | exact | all real-rooted except the Golay control |
| D. All binary cyclic codes of odd length n (one per divisor of xⁿ + 1), and their parity extensions | **complete for n ≤ 25 (lengths 7–26)**: 256 entries, every non-trivial code analysed. **Partial for n = 27, 29, 31** (lengths 27–32): 14 analysed, 35 skipped as too large (more than 60,000 cyclic flats, a nested chain with ∏(layer + 1) > 3000, or dual dimension > 26); length 31 was stopped after 28 of about 254 entries | exact | 202 matroids analysed (ranks 1–26); the only non-real-rooted ones are the two parity-extended Golay codes (n = 23, both generator polynomials), i.e. G24 again. Skips are listed with their reason in `data/cyclic.jsonl` |

## Findings
1. **G24 is fragile.** Removing any single element, by shortening or puncturing, already gives a real-rooted Z, and so do all G(s,p) with s + p ≤ 5. These derivatives are nonetheless the closest-to-failure matroids we have seen. In h(z), defined by Z(t) = t^m h(t + 1/t), two real roots nearly collide: the minimal relative gap is 2·10⁻⁴ to 1.5·10⁻³ in ranks 10–12. Random linear matroids in the predecessor case stayed above 0.86.
2. **Single-size paving matroids** (Scan A):
   - Non-real-rootedness occurs only for λ ≥ 96.5% of an upper bound U on packings P(r−1, h, n), and U is not attained except by Steiner systems.
   - Every known Steiner system in range gives real-rooted Z: S(3,5,17), S(3,6,22), S(3,6,26), S(3,8,50), S(4,5,11), S(4,7,23), S(5,6,12), S(5,7,28) and S(5,8,24).
   - The two Steiner-admissible failures in rank 4 would need S(3,7,22) or S(3,11,56). Neither exists: their derived designs S(2,6,21) (14 blocks < 21 points) and S(2,10,55) (33 < 55) violate Fisher's inequality.
   - So no realisable single-size paving counterexample is known in rank ≤ 10, n ≤ 60. This is not a proof: near-extremal packings are not ruled out in general.

3. **Cyclic codes** (Scan D): no new counterexample among the 202 binary cyclic and parity-extended cyclic codes analysed. Coverage is complete up to length 26; G24 is the only failure.

## Engineering notes
- **Nested engine.** P(M) comes from R = Σ_{F ≠ 0} #(flats of profile f)·t^{rk F}·P(M/F), where M/F is again nested and is computed recursively, cached by its chain.
  - *Lemma:* in a nested matroid, replacing an element by one from a lower layer never raises the rank, by rk(S) = min_i (r_i + |S \ X_i|). So a cyclic flat of M/F that contains a layer contains every lower layer. The chain of M/F therefore consists of prefixes of the remnant layers, and it is read off directly.
  - The first implementation enumerated all profiles and is kept in `test_zfast.py` as a reference. The new engine matches it on 80 random chains of up to 6 layers.
- **Bug found and fixed (2026-11-28).** sympy's Sturm `count_roots` counts *distinct* real roots, so a real-rooted polynomial with repeated roots, such as a direct sum, looked non-real-rooted.
  - It produced 20 false failures in an early run of Scan D.
  - `real_rooted` now uses the square-free part, and `test_zfast.py` has a regression test.
  - Scans A–C were unaffected: their "real-rooted" verdicts stand, and Scan A, rerun, is identical.
- **Speed limit.** Python nested computations with more than about 6 layers take minutes each, which is why Scan D has skips at lengths 27–32. A C version of the nested engine would remove them.

## Open leads
- Port the nested engine and the refinement to C, and finish Scan D at lengths 27–32.
- Self-dual and extremal binary codes of length 20–22 from the classification lists, and other 5-designs.
- Abstract condensed configurations near Golay's, then a search for realisations.
- Paving matroids with two or more stressed hyperplane sizes (Scan A covers one size only).
- Ternary codes. The ternary Golay matroid is sparse paving (hexads are circuit-hyperplanes), hence real-rooted.

## Files
- `zfast.py`: fast exact P and Z (flats engine, nested engine, Schubert expansion, exact real-rootedness)
- `zpoly_ref.py`: reference implementation (copy from `z-polynomial-rank4`)
- `cycflats.c`: cyclic flats of a binary matroid (`cc -O2 -o cycflats cycflats.c`)
- `condense.py`: partition refinement + Theorem 3.2 for binary matroids
- `codes.py`: Golay, cyclic, Reed–Muller codes; shortening and puncturing
- `packing_bounds.py`: constant-weight-code bounds A(n, d, w)
- `scan_paving.py`, `scan_golay.py`, `scan_codes.py`, `scan_cyclic.py`: scans A–D (`data/*.jsonl`)
- `test_zfast.py`: cross-checks of every fast path
