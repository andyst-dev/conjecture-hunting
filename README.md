# conjecture-hunting

These are computational attacks on recent, checkable conjectures, mostly from arXiv combinatorics preprints. Each one is done with exact arithmetic, after reproducing the paper's own values, with every attempt logged.

**Negative results are published too.** "No counterexample up to X" with a precise coverage table saves the next person from repeating the same search.

All work here is **AI-assisted and unrefereed**. The protocol every case follows is in [`CLAUDE.md`](CLAUDE.md).

## Cases

| Case | Source | Status |
|---|---|---|
| [Khan–Haidar, Conj. 7.6](cases/khan-haidar-7.6/CASE.md) | arXiv:2611.21670 | No counterexample in the classes searched; 2 resistance-regular graphs not among the paper's listed examples |
| [Costa–Luu–Wood–Yip, CDM](cases/cdm-alpha2/CASE.md) | arXiv:2512.17114 | No counterexample on ≤ 17 vertices, exhaustive (paper: ≤ 11); reduction to maximal triangle-free graphs |
| [GMR inequality in small groups](cases/gmr-small-groups/CASE.md) | arXiv:2611.12776 (Q9.4) | No counterexample of order < 243 in the families searched (heuristic, calibrated on the paper's 243 witness) |
| [Notik, Singer rulers and chambers (§7)](cases/notik-singer-golomb/CASE.md) | arXiv:2611.22656 | Verified for m ≤ 113 and m = 121 (paper: m ≤ 103); chain reduction of the chamber LP |

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
