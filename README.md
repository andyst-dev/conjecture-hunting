# conjecture-hunting

These are computational attacks on recent, checkable conjectures, mostly from arXiv combinatorics preprints. Each one is done with exact arithmetic, after reproducing the paper's own values, with every attempt logged.

**Negative results are published too.** "No counterexample up to X" with a precise coverage table saves the next person from repeating the same search.

All work here is **AI-assisted and unrefereed**. The protocol every case follows is in [`CLAUDE.md`](CLAUDE.md).

## Cases

| Case | Source | Status |
|---|---|---|
| [Khan–Haidar, Conj. 7.6](cases/khan-haidar-7.6/CASE.md) | arXiv:2611.21670 | No counterexample in the classes searched; 2 resistance-regular graphs not among the paper's listed examples |
| [Costa–Luu–Wood–Yip, CDM](cases/cdm-alpha2/CASE.md) | arXiv:2512.17114 | No counterexample on ≤ 17 vertices, exhaustive (paper: ≤ 11); reduction to maximal triangle-free graphs |
| [Benedetti–Menara, Conj. 7.1](cases/benedetti-menara-7.1/CASE.md) | arXiv:2611.08071 | Proved for all m ≤ 8 by SAT, and r = 2 up to m = 12 exhaustively (paper: random tests only) |
| [Wild number of K_n, Conj. 6.1](cases/wild-number-kn/CASE.md) | arXiv:2610.22142 | No counterexample: exact SAT values reproduce n ≤ 11; sampled colourings up to n = 23 |
| [Z-polynomial, rank 4 (exploratory)](cases/z-polynomial-rank4/CASE.md) | PXY18 Conj. 5.1 | Closed-form rank-4 criterion; random linear matroids stay below the threshold (equality only for U_{2,3}⊕U_{2,3}) |
| [GMR inequality in small groups](cases/gmr-small-groups/CASE.md) | arXiv:2611.12776 (Q9.4) | No counterexample of order < 243 in the families searched (heuristic, calibrated on the paper's 243 witness) |
| [Notik, Singer rulers and chambers (§7)](cases/notik-singer-golomb/CASE.md) | arXiv:2611.22656 | Verified for m ≤ 113 and m = 121 (paper: m ≤ 103); chain reduction of the chamber LP |
| [Leake–Mohammadi Yekta, Question 1.7](cases/leake-yekta-q1.7/CASE.md) | arXiv:2609.39917 | No violation of the weighted log-concavity question on 3,471 small unimodular configurations (rank 2–5) and contingency tables up to 3×4 (float screen) |
| [Cherukupally, Conj. 3.4](cases/cherukupally-a-n/CASE.md) | arXiv:2603.17434 | No n with \|A(n)\| ≥ 4 up to 10^27, exhaustive (abstract: n < 10^15) |
| [Easlick, non-transitive dice, Conj. 1](cases/easlick-dice-conj1/CASE.md) | arXiv:2610.40321 | Proved for 7 ≤ n ≤ 60 and n = 70, 80, 90, 100 (explicit cycles of length n+1) under the reading the paper's evidence uses; the literal "a cycle of length n" reading fails at n = 7 |

## Layout
```
CLAUDE.md      protocol (read this first)
TEMPLATE.md    case write-up template
cases/<slug>/  one self-contained folder per conjecture
shared/        code reused by several cases (filled lazily)
```

## Setup
```
pip install -r requirements.txt
cd cases/khan-haidar-7.6 && python3 validate.py
```
