# Benedetti–Menara, Conjecture 7.1 (arXiv:2611.08071 v1)

**Status:** proved for every m ≤ 8 (all r) by SAT, and for **r = 2 up to m = 12** (exhaustive over triangle-free graphs); the general case is open
**Date:** 2026-11-25 · **Tools:** PySAT (CaDiCaL 1.5.3) · **AI-assisted (Claude), unrefereed**

## Statement
Source: B. Benedetti and G. Menara, *Diameter bounds for complexes that are close to flag or close to normal*, arXiv:2611.08071.

- A family F of r-subsets of an m-set is **1-layered** if it is partitioned into nonempty layers L₀, …, L_N such that, for every nonempty set S, the indices t for which some member of L_t contains S form an interval.
- F is **flag** if, in the graph where two vertices are adjacent when they lie in a common member, no member is contained in an (r+1)-clique.
- **Conjecture 7.1.** A flag 1-layered family of r-subsets of an m-set has at most m − r + 1 layers.
- It would imply an interpolation between the Adiprasito–Benedetti and Kalai–Kleitman diameter bounds.
- The paper says the case r = 2 (triangle-free graphs) is already open. Its evidence is random testing only: tens of thousands of random families for r = 2, 3, 4.

## Method
SAT(m, r, L) asks whether a flag 1-layered family of r-subsets of [m] with exactly L nonempty layers exists (`sat.py`).

- Variables: member T in layer t; an edge variable for every pair.
- The interval property is enforced for every S with |S| < r via prefix and suffix indicator variables. |S| = r is a single member, and larger S lie in no member.
- Flagness: no member T plus an outside vertex w forms an (r+1)-clique.
- Merging two adjacent layers preserves 1-layeredness, so feasibility is monotone in L. **UNSAT at L = m − r + 2 proves the conjecture for (m, r).**

## Validation
- The bound itself is attained: L = m − r + 1 is SAT for (6,2), (6,3), (7,4).
- Flagness matters: without it, r = 2, m = 4 admits 4 > 3 layers (the solver finds the K₄ labelling, which matches a hand construction).
- `test_bruteforce.py` checks the maximum number of layers against a brute force from the definitions, with no SAT involved, for (m, r) = (3,2), (4,2), (5,2), (4,3), (5,3), with and without flagness where feasible. The maximum agrees in every case.

## Search coverage (all UNSAT, so the conjecture holds)
| r \ m | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|
| 2 | ✓ | ✓ | ✓ | ✓ | ✓ (1.4 s) | ✓ (25 s) |
| 3 | | ✓ | ✓ | ✓ | ✓ (2.3 s) | ✓ (168 s) |
| 4 | | | ✓ | ✓ | ✓ | ✓ (43 s) |
| 5 | | | | ✓ | ✓ | ✓ (3 s) |
| 6 | | | | | ✓ | ✓ |

r = 2, m = 9 did not finish within 5 minutes. Symmetry breaking is minimal (one clause).

## r = 2: graph-by-graph search (`r2_graphs.py`)
For r = 2, a flag 1-layered family is a triangle-free graph G with a surjective edge labelling in which every vertex sees an interval of labels. A counterexample needs m labels.

**Reductions.**
- G may be taken connected: each component's labels form an interval, so by induction a disconnected G has at most m − c labels.
- G needs ≥ m edges, so trees are excluded.

**Path bound.** The labels at v span at most deg(v) − 1. Summing deg − 1 along line-graph walks between two edges bounds their label difference. Graphs whose bound is ≤ m − 2 are skipped; every other graph gets an exact SAT check for m labels.

**Validation (`test_r2.py`).** For every connected triangle-free graph with m ≤ 8, the exact maximum (by SAT) never exceeds the path bound + 1. Its maximum over graphs is exactly m − 1, matching the global SAT results.

| m | graphs (connected, triangle-free, ≥ m edges) | pruned by bound | exact SAT checks | counterexamples |
|---|---|---|---|---|
| 9 | 1,333 | 1,158 | 175 | 0 |
| 10 | 9,726 | 8,640 | 1,086 | 0 |
| 11 | 90,607 | 83,367 | 7,240 | 0 |
| 12 | 1,143,510 (= 1,144,061 connected triangle-free − 551 trees) | 1,086,963 | 56,547 | 0 |

## Open leads
- r = 2, m = 13 (about 19M graphs) needs a C port of the path bound. Stronger symmetry breaking in `sat.py` would push r ≥ 3 past m = 8.
- For r = 2 the problem is purely about graphs: edge labellings of a triangle-free graph in which every vertex sees an interval of labels. Maximise the number of labels; the conjecture says it is at most m − 1.

## Files
`sat.py` (encoding and solver; `python3 sat.py r m`, `--L k`, `--noflag`), `test_bruteforce.py`, `data/sat_runs.jsonl`
