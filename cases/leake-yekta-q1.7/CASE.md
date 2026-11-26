# Leake–Mohammadi Yekta, Question 1.7 (arXiv:2609.39917, Question 1.7)

**Status:** no counterexample in the classes searched: 3,471 unimodular configurations of rank 2–5 with small domains, plus contingency tables up to 3×4 (float LP screen; any hit would be re-checked exactly)
**Date:** 2026-11-26 · **Tools:** Python 3, sympy (exact unimodularity, exact rational re-check), scipy/HiGHS (LP screen) · **AI-assisted (Claude), unrefereed**

## Statement
Source: J. Leake and M. Mohammadi Yekta, *Log-concavity and approximate counting for totally unimodular polytopes*, arXiv:2609.39917. The PDF read is dated "September 2026" and shows no version marker. It appeared as a replacement in the math.CO listing of 2026-11-25. Its abstract says "we conjecture a generalization of Barvinok's conjecture", and the body has this as Question 1.7.

- **Unimodular** (paper, §3): every non-singular rank(M) × rank(M) submatrix of M has determinant ±1. Total unimodularity is not required.
- N_M(b) := #{y ∈ Z^n_{≥0} : My = b}.
- **Theorem 1.4 (proved in the paper).** Let M be unimodular and b, b^±, b ∈ Z^m with b = (b⁻ + b⁺)/2, and let all fibers be finite. Then N(b)² ≥ N(b⁻)·N(b⁺).
- **Question 1.7.** Let M ∈ Z^{m×n} be unimodular and μ a probability distribution on Z^m (perhaps finitely supported) with mean b ∈ Z^m. Is it true that N(b) ≥ ∏_a N(a)^{μ(a)}?
- Barvinok's original question (2006 preprint "Brunn–Minkowski inequalities for contingency tables and integer flows", eq. 2.3) is the case of contingency tables, possibly weighted. There, all margins in the combination have the same total.
- The paper notes that over Z^m, convexity on lines does not imply global convexity. So Question 1.7 is not a consequence of Theorem 1.4.

**Our reading.** For a finite domain D of lattice points with N > 0, the inequality fails for some finitely supported μ on D exactly when log N(b) lies strictly below the upper concave envelope of {(a, log N(a)) : a ∈ D} at b. A counterexample with finite support therefore needs only a finite D. Collinear supports reduce to log-concavity on a line, which follows from Theorem 1.4. So counterexamples need m ≥ 2 and a support that is not collinear. We require finite fibers: there is an integer w with w·c_j ≥ 1 for every column c_j.

## Method
1. **Exact counts** (`lyq.counts`). Dynamic programming over columns gives N(b) for every b with w·b ≤ W. Every solution of My = b has Σ y_j (w·c_j) = w·b, so the count is complete. Exact Python integers.
2. **Screen** (`lyq.screen`). For each b ∈ D, an LP maximises Σ μ_a log N(a) subject to μ ≥ 0, Σ μ_a = 1 and Σ μ_a a = b. A positive gap means a candidate violation (float, HiGHS).
3. **Exact re-check** (`lyq.exact_check`). The LP support is turned into exact rational weights μ by solving the affine system. We then test N(b)^D < ∏ N(a)^{D μ_a} in integers. No hit reached this stage.

## Reproduction of the paper
The paper has no numerical tables. We reproduce classical counts and the paper's own theorem (`validate.py`).

| Object | Paper / known value | Our value |
|---|---|---|
| 2×2 tables, all margins with r1+r2 ≤ 14 | min(margins)+1 | identical |
| 3×3 tables, all margins n, n = 0..5 | 1, 6, 21, 55, 120, 231 | identical |
| 4×4 tables, all margins n, n = 0..3 (A001501) | 1, 24, 282, 2008 | identical |
| Theorem 1.4 (midpoint), 60 instances, all pairs in domain | no violation | 0 violations |

Cross-checks (`test_counts.py`, passing):
- The DP counter agrees with direct enumeration on 300 random pointed {−1,0,1} configurations (m ≤ 4).
- The screen and the exact checker both detect the known violation for the non-unimodular M = [1 2]: N(1)² = 1 < 2 = N(0)N(2), with gap log 2 / 2.

## Search coverage
| Class | instances | n (columns) | domain points (total / max) | max LP gap | float hits | exact violations |
|---|---|---|---|---|---|---|
| a3 m=3 | 2186 | 3–12 | 512918 / 819 | 7.4e-13 | 0 | 0 |
| ct2x2 | 1 | 4–4 | 3311 / 3311 | 1.5e-13 | 0 | 0 |
| ct2x3 | 1 | 6–6 | 1705 / 1705 | 6.7e-15 | 0 | 0 |
| ct2x4 | 1 | 8–8 | 1218 / 1218 | 3.6e-14 | 0 | 0 |
| ct2x5 | 1 | 10–10 | 1302 / 1302 | 2.7e-15 | 0 | 0 |
| ct3x3 | 1 | 9–9 | 812 / 812 | 6.4e-15 | 0 | 0 |
| ct3x4 | 1 | 12–12 | 798 / 798 | 0.0e+00 | 0 | 0 |
| rand m=2 | 133 | 2–7 | 7553 / 91 | 9.1e-14 | 0 | 0 |
| rand m=3 | 651 | 3–8 | 145396 / 819 | 4.3e-13 | 0 | 0 |
| rand m=4 | 439 | 4–9 | 259683 / 1365 | 7.7e-13 | 0 | 0 |
| rand m=5 | 62 | 5–9 | 36038 / 1287 | 1.9e-13 | 0 | 0 |

Classes:
- `rand`: random {−1,0,1} matrices that are unimodular (exact minors), have full row rank and have finite fibers. Repeated columns are allowed. Domain {w·b ≤ W} with W ∈ [5, 12], capped at 1500 points.
- `a3`: rank-3 sub-multisets of the signed A₃ configuration {e1, e2, e3, e1+e2, e2+e3, e1+e2+e3}, which represents M(K4). Multiplicities are up to 2 and either sign is allowed when the result is pointed. Every simple rank-3 regular matroid is a restriction of M(K4).
- `ct p×q`: contingency tables, which are Barvinok's original setting. The domain is all margins (r, c) with Σr ≤ W. Here W = 20, 9, 6, 5, 5, 4 for 2×2, 2×3, 2×4, 3×3, 2×5, 3×4.

All "max LP gap" values are at floating-point noise level (≤ 10⁻¹²). In those cases the envelope touches log N(b), usually through an exactly log-linear piece.

## Findings
- No violation of Question 1.7 was found. The screen is in floats. No candidate exceeded tolerance 10⁻⁹, so the exact re-check never ran on real data.
- This is evidence only for small domains near the apex of the cone. Effects at the boundary of the cone are included there.

## Open leads
- Larger domains: W ≥ 20 in rank 3–4, which needs a faster envelope computation (one convex hull instead of one LP per point).
- Rank ≥ 4 regular matroids that are not sub-configurations of a graphic one, for example cographic M*(K_{3,3}) and R10.
- The analogue of the Ferroni counterexample [Fer26] to the IDP Ehrhart conjecture ("dominos"): look at the polytopes used there.

## Files
- `lyq.py`: counting, unimodularity test, LP screen, exact checker
- `test_counts.py`: cross-checks
- `validate.py`: reproduction table
- `search.py <class> <seed> <seconds>`: searches; logs to `data/coverage_*.jsonl`, hits to `data/hits_*.jsonl`
- `summarize.py`: coverage table
