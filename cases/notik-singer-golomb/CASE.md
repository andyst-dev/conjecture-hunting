# Notik: Singer rulers are unique shortest points of their chambers (arXiv:2611.22656 v1, §7)

**Status:** no counterexample up to m = 131. The conjecture is now verified for every prime power m ≤ 131 (up to 132 marks; the paper went to 104)
**Date:** 2026-11-25, extended 2026-11-26 · **Tools:** Python, HiGHS (scipy), python-flint, Fractions · **AI-assisted (Claude), unrefereed**

## Statement
Source: E. Notik, *Golomb rulers with real marks: the real and integer optima agree up to ten marks*, arXiv:2611.22656v1.

- A **real 1-separated ruler** has marks 0 = x₁ < … < x_n with |(x_j − x_i) − (x_l − x_k)| ≥ 1 for distinct pairs, and all differences ≥ 1.
- Gap coordinates: g_t = x_{t+1} − x_t, so g ∈ ℝ^k with k = n − 1.
- Fix the sign of every relation e·g, where e = 1_I or 1_I − 1_J over intervals I, J. Each choice defines a polyhedral **chamber** P_σ = {g : σ_e (e·g) ≥ 1}.
- **Singer ruler.** Take a perfect difference set D mod q = m² + m + 1, with m a prime power. For each multiplier t coprime to q, cut tD at its largest cyclic gap. This gives n = m + 1 marks.
- **Conjecture (§7, unnumbered).** For every prime power m, every Singer ruler is the only shortest point of its chamber.
- The paper verifies this for m ≤ 103: 37 prime powers, 24,450 distinct rulers, with exact certificates.

## Method
**Chain reduction (our observation).** Sort the N = n(n−1)/2 differences of the ruler, with interval forms f₁, …, f_N. Then

  P_σ = { g : f₁·g ≥ 1, (f_{r+1} − f_r)·g ≥ 1 for r < N }.

- "⊆": these rows are among the defining rows, with the correct signs.
- "⊇": chaining gives f_s·g − f_r·g ≥ s − r ≥ 1 for every r < s, and f_r·g ≥ 1.

So each chamber has N ≈ n²/2 rows instead of |R_n| ≈ n⁴/8. That makes n > 104 cheap.

**Certificate (as in the paper).** Find y ≥ 0 supported on the rows tight at the ruler, with Σ y_r a_r = (1, …, 1) and Σ y = length, such that the support has rank k. This proves the ruler is the unique minimiser.
- y is found with HiGHS, then certified in exact rationals (python-flint).
- On degenerate vertices, the certificate is built by averaging exact basic dual solutions until the support reaches rank k. If that cannot be done, the ruler is reported as "shortest, not unique".
- If no y exists, the chamber optimum is computed and certified exactly ("not shortest").

**Mirror reduction.** The reversal map g ↦ (g_k, …, g_1) sends a ruler's chamber onto its mirror's chamber and preserves length. For the extension, only one ruler per mirror pair is checked.

## Reproduction of the paper
| Object | Paper | Ours |
|---|---|---|
| Distinct Singer rulers (largest-gap cut, all t), m ≤ 103 | 24,450 | 24,450 |
| All cuts, m ≤ 17: not shortest / shortest but not unique | 396 / 18 | 396 / 18 |
| Prop. 2: X = {0,4,24,27,67,77,86,88,103,119,125,132} (12 marks, non-largest cut) | unique chamber optimum (0,4,22,25,123/2,141/2,157/2,161/2,94,217/2,114,121), length 121 | identical point, length 121 |
| Perfect difference sets for all prime powers used | n/a | checked exhaustively |

## Search coverage
| Range | Rulers | Result |
|---|---|---|
| m ≤ 17, all cuts | 3,412 | matches the paper (396 / 18 exceptions, all at non-largest cuts) |
| m = 107 (108 marks), largest-gap cut, all multipliers | 1,512 mirror classes (3,024 rulers) | all unique shortest points |
| m = 109 (110 marks) | 1,140 classes (2,280 rulers) | all unique |
| m = 113 (114 marks) | 1,980 classes (3,960 rulers) | all unique |
| m = 121 = 11² (122 marks) | 648 classes (1,296 rulers) | all unique |
| m = 125 = 5³ (126 marks) | 828 classes (1,656 rulers) | all unique |
| m = 127 (128 marks) | 1,806 classes (3,612 rulers) | all unique |
| m = 128 = 2⁷ (129 marks) | 336 classes (672 rulers) | all unique |
| m = 131 (132 marks) | 2,882 classes (5,764 rulers) | all unique |
| Independent `Fraction` re-check (`verify_certificates.py`) | 3 rulers each for m = 107, 113, 125, 127, 128 and 131, 2 for m = 121, plus every ruler for m = 13, 19 and 29 (the first 40) | all certified |

m = 125, 127, 128 and 131 were run on 2026-11-26 (`python run.py --list <m> --mirror`, one process per m). Wall time: 1,194 s, 2,623 s, 511 s and 4,089 s. Class counts agree with φ(q)/(6e) for m = pᵉ, since the multiplier group {pⁱ} has order 3e. Outputs: `data/results_<m>_mir.jsonl`, `logs/run_<m>_mir.log`, `logs/indep_<m>.txt`. Every prime power in 113 < m ≤ 131 (121, 125, 127, 128, 131) is now covered.

## Findings
1. **New range.** Every Singer ruler (largest-gap cut, every multiplier) for m = 107, 109, 113, 121, 125, 127, 128 and 131 is the unique shortest point of its chamber, certified in exact arithmetic. That is 11,132 mirror classes, or 22,264 rulers. Together with the paper, the conjecture now holds for all prime powers m ≤ 131.
2. **Chain reduction.** The chamber needs only the n(n−1)/2 consecutive-order rows (proof under Method). This cuts the LP from about n⁴/8 rows to about n²/2, which is what makes n > 104 cheap. The paper does not state this simplification, as far as we could read.
3. All of the paper's checkable data was reproduced exactly (see the table above).

## Files
- `singer.py`: Singer difference sets for any prime power, via GF(m) tables and a primitive cubic over GF(m)
- `chamber.py`: chain reduction, exact certificates and the exact better point
- `run.py`: driver (`python3 run.py 2 17 --allcuts`, or `python3 run.py 0 0 --list 107,109 --mirror`)
- `test_chain.py`: chain reduction vs the full chamber (all of R_n), 150 random rulers + all cuts for m = 3, 4, 5, agreement to about 1e-15
- `verify_certificates.py`: independent checker in Python `Fraction`s; re-checks the perfect difference set property and every certificate condition exactly (`python3 verify_certificates.py 107`)
- `data/`, `logs/`: outputs
