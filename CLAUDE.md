# Protocol for agents working in this repo

This repo hunts recent, computationally checkable conjectures, mostly from arXiv preprints. Each case lives in `cases/<slug>/` and is documented in `cases/<slug>/CASE.md`, which follows `TEMPLATE.md`.

The goal is **trustworthy results, including negative ones**. A false counterexample is far worse than no result.

## Picking a target
- Prefer conjectures that are recent, niche, and decidable on finite objects: graphs, matroids, posets, polynomials, sequences.
- Before starting, search for the conjecture's current status (arXiv, the authors' pages, and repos that track conjectures). If someone already settled it, record that and stop.

## The non-negotiable steps
1. **State it exactly.** Copy every definition and the conjecture statement from the source into `CASE.md`. Record the version (for example arXiv v1) and the conventions used: normalisations, signs, pseudo-inverse versus inverse, and so on.
2. **Implement in exact arithmetic.** Use rationals (python-flint, sympy or Sage). Floating point is allowed only for screening.
3. **Reproduce the paper before searching.** Recompute the paper's own examples and stated values, and list them in `CASE.md` as "paper value / our value". If any of them disagrees, the implementation is wrong. Do not search until they all agree.
4. **Cross-check every fast path.** Any shortcut formula (quotients, Fourier decompositions, symmetry reductions) must be compared against the exact full computation on many random instances. Commit the test as `test_*.py`.
5. **Log coverage.** Record what was searched and how far: classes, sizes and parameters. Write every detection to `data/*.jsonl`.
6. **Claims need a standalone verifier.** A counterexample counts only once a single script rebuilds the object from scratch and verifies it in exact arithmetic, without importing the search code. Name it `verify_*.py`.
7. **Report honestly.** Give the case one status: `refuted`, `proved (partial: …)`, `no counterexample up to …` or `abandoned`. Separate what was verified exactly from what was only screened in floating point.

## Etiquette
- Contact the authors before announcing a refutation publicly.
- Credit the source paper and any prior work.
- Never present an AI-assisted result as refereed. State that it is unrefereed and AI-assisted.

## Code conventions
- Python 3 with numpy, networkx, python-flint and scipy (`requirements.txt`).
- Each case is self-contained. Move code to `shared/` only once a second case needs it.
- Scripts run from inside their case directory.
