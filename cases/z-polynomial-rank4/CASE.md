# Proudfoot–Xu–Young, real-rootedness of the matroid Z-polynomial: rank 4 (exploratory)

**Status:** exploratory. No counterexample; a closed-form criterion for rank 4.
**Date:** 2026-11-25 · **AI-assisted (Claude), unrefereed**

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

## Leads
- Prove 4a2 ≤ W3² + 8 for all rank-4 matroids. That would settle rank 4 and explain why equality comes from U_{2,3} ⊕ U_{2,3}.
- A counterexample, if any, is in rank ≥ 5, where the KL coefficients of contractions enter.
