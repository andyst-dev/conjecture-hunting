# Proudfoot–Xu–Young, real-rootedness of the matroid Z-polynomial: rank 4 (exploratory)

**Status:** no counterexample in the classes searched (random, not exhaustive): random simple linear rank-4 matroids over GF(2), GF(3) and GF(5), and a hill-climb over connected rank-5 linear matroids over GF(2) (n ≤ 16) and GF(3) (n ≤ 14), 90 min each. Floating-point screen; the best rank-5 candidates were re-checked exactly. Also derived: a closed-form real-rootedness criterion for rank 4 (4·a2 ≤ W3² + 8); whether every rank-4 matroid satisfies it is not proved.
**Date:** 2026-11-25, long run 2026-11-26 · **AI-assisted (Claude), unrefereed**

## Why this target
It is a "domino" candidate. Real-rootedness for Chow, augmented Chow and chain polynomials (arXiv:2610.17135) and for dual Chow (Davis, Nov 2026) has been disproved, and Kazhdan–Lusztig polynomials of matroids need not be unimodal (arXiv:2607.24186). Real-rootedness of Z_M (PXY18, Conj. 5.1) is still open; it is proved for sparse paving matroids (arXiv:2609.07636).

## Rank-4 criterion (derived here)
For a simple matroid of rank 4, write W2 for the number of lines, W3 for the number of hyperplanes, and |X| for the number of points in a flat X. Then

  Z_M(t) = 1 + W3 t + a2 t² + W3 t³ + t⁴,  with a2 = W2 + Σ_H |H| − Σ_L |L|.

This uses Z = Σ_F t^{rk F} P_{M/F}, P = 1 in rank ≤ 2, and P = 1 + (W_{r−1} − W_1) t in ranks 3 and 4.

Z_M is real-rooted ⇔ **4·a2 ≤ W3² + 8** (the γ-polynomial 1 + (W3 − 4)x + (a2 − 2W3 + 2)x² has real roots) and a2 ≥ 2W3 − 2 (γ-positivity).

## Experiment (`zrank4.py`)
- Random simple linear rank-4 matroids over GF(2), GF(3) and GF(5), with n from 5 up to the full projective space. The GF(7) run ran out of memory.
- The largest value of 4a2/(W3² + 8) found is exactly **1**, at W3 = 6, a2 = 11: Z = (t² + 3t + 1)² = Z(U_{2,3})², i.e. the direct sum U_{2,3} ⊕ U_{2,3}, which has a double root.
- Every other sample is strictly below 1.

The coalescence families of arXiv:2610.17135 (M_q, G_{4,q}) do **not** threaten rank-4 Z. Their auxiliary points create about q⁸ three-point hyperplanes, so W3² dwarfs a2.

## General implementation (`zpoly.py`, `search_z.py`)
- KL polynomials of every interval of the lattice of flats come from the EPW identity; each one is checked against t^r P(1/t) − P(t) = R(t). Z follows from them, and real-rootedness is tested **exactly** with sympy root counting.
- `test_zpoly.py` checks it against known values: P(U_{3,4}) = 1 + 2t, P(M(K₄)) = 1 + t, P(M(K₅)) = 1 + 5t, palindromic Z, the rank-4 closed form, and Z(A ⊕ B) = Z(A)Z(B).
- `search_z.py` hill-climbs over **connected** linear matroids (direct sums are excluded, since they factor and give trivial double roots), minimising the γ-root gap. Short smoke runs (2 min) in rank 5 over GF(2) and GF(3) stay far from violation, with a best margin of about 0.86. Long runs are meant for a local machine (`RUN.md`).

## Long run, rank 5 (2026-11-26)
Two 90-minute runs, one per field, on a local Mac:

| Run | Command | Evaluations | Restarts logged | Best γ-root margin (min / median / max over restarts) | n range |
|---|---|---|---|---|---|
| GF(2) | `search_z.py --rank 5 --p 2 --maxn 16 --minutes 90 --seed 1` | 39,969 | 199 | 0.842 / 0.861 / 0.878 | 6–16 |
| GF(3) | `search_z.py --rank 5 --p 3 --maxn 14 --minutes 90 --seed 2` | 111,877 | 557 | 0.873 / 0.873 / 0.942 | 6–14 |

- **No candidate.** No evaluation had a negative margin (γ-polynomial non-real-rooted in floating point), so `data/z_candidates.jsonl` was never created. The margin is min gap / max |root| of the γ-roots, computed with numpy. It is a floating-point **screen**, not a proof for the matroids searched.
- The best Z found are Z = 1 + 22t + 95t² + 95t³ + 22t⁴ + t⁵ (GF(2), n = 13) and Z = 1 + 15t + 50t² + 50t³ + 15t⁴ + t⁵ (GF(3), n = 9). Both have 5 real roots by an **exact** sympy root count, and both factor as (1 + t)·(quartic). The vectors are in `data/run_r5_p*.log`.
- The margin never fell below about 0.84, the same as the 2-minute smoke runs. On GF(3) the median equals the minimum: the hill-climb keeps returning to the same optimum (Z = 1 + 15t + 50t² + …). The objective does not seem to push towards violation, so a different objective or move set would be needed for a more useful search.
- Coverage is random, not exhaustive: 756 restarts appended to `data/z_search.jsonl` (lines 76–831).

## Leads
- Prove 4a2 ≤ W3² + 8 for all rank-4 matroids. That would settle rank 4 and explain why equality comes from U_{2,3} ⊕ U_{2,3}.
- A counterexample, if any, is in rank ≥ 5, where the KL coefficients of contractions enter.
