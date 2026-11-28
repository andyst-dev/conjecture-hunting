# Knor–Sedlar–Škrekovski, maximum interval index without pendant vertices (arXiv:2610.39966, Problem 17)

**Status:** no counterexample up to n ≤ 10 (exhaustive, exact integers); weak heuristic screen for 11 ≤ n ≤ 30 found nothing
**Date:** 2026-11-28 · **Tools:** Python 3, numpy, networkx, pynauty, scipy, C (own code) · **AI-assisted (Claude), unrefereed**

## Statement
Source: M. Knor, J. Sedlar, R. Škrekovski, *On the maximum interval index of a graph*, arXiv:2610.39966. The full text was read on 2026-11-28 through a fetcher; it shows no version line, so we treat it as v1. The definitions below are restated from that text.

- I(u,v) = {w : d(u,w) + d(w,v) = d(u,v)}, the vertices lying on some shortest u–v path.
- Int(G) = Σ_{{u,v}} (|I(u,v)| − 1).
- A triple of distinct vertices is *collinear* if one of them lies in the interval of the other two. noncol(G) counts the non-collinear triples.
- Lemma 4: Int(G) = C(n,2) + col(G) = C(n,2) + C(n,3) − noncol(G). So maximising Int means minimising noncol.
- D_n (n ≥ 8), the dumbbell: two disjoint C_4 joined by a path of length n − 7. Prop. 6: noncol(D_n) = 2(n − 4).
- Theorem 15 / Cor. 16: for n ≥ 31, every graph without pendant vertices has Int(G) ≤ (n³ − n)/6 − 2(n − 4), attained by D_n.

**Problem 17** (concluding remarks), for the remaining range 8 ≤ n ≤ 30:
- (a) decide whether noncol(G) ≥ 2(n − 4) for every 2-connected graph G of diameter ≥ 4 on 8 ≤ n ≤ 30 vertices;
- (b) determine the n ≥ 9 for which D_n is the unique graph without pendant vertices that maximises Int. The paper notes that n = 8 is excluded because K_{4,4} and C_8 also have noncol = 8.

"Without pendant vertices" is read as minimum degree ≥ 2, for connected graphs on n ≥ 3 vertices.

## Prior status
- The paper reports no computation at all. It declares that Claude was used to check and simplify the proofs.
- No follow-up was found (web search, 2026-11-28).

## Reproduction of the paper
`reproduce.py` prints 74 checks, all OK.

| Object | Paper value | Our value |
|---|---|---|
| Int(P_n) = (n³−n)/6, n = 1..12 | formula | identical |
| Int(P_4), Int(C_4) | 10, 10 | 10, 10 |
| noncol(K_{2,3}), K_{3,3}, K_{3,4}, K_{4,4} | 1, 2, 5, 8 | 1, 2, 5, 8 |
| noncol(K_{a,b}) = C(a,3)+C(b,3), a ≤ b ≤ 6 | formula | identical |
| noncol(C_8) | 8 | 8 |
| noncol(D_n) = 2(n−4), 8 ≤ n ≤ 40 | formula | identical |
| W(C_20); W(two triangles joined by a path, n = 20) | 1000; 1294 | 1000; 1294 |
| Lemma 4 on 3000 random graphs, plus noncol against brute force | identity | 0 failures |
| Theorem 5 (Int ≤ Int(P_n), equality iff P_n, or n = 4 and C_4), all connected graphs n ≤ 9 | holds | holds: unique maximiser P_n, plus C_4 at n = 4 (`exhaustive.py`) |

## Search coverage
| Class | Range | Method | Result |
|---|---|---|---|
| All connected graphs (Theorem 5, Problem 17 a/b) | n ≤ 9 (261,080 connected graphs at n = 9) | Exhaustive. Graphs from vertex augmentation with pynauty certificates; counts match OEIS A000088. Exact integer distances. `exhaustive.py`, `gen_graphs.py` | no violation |
| Connected graphs with δ ≥ 2 | n = 10 (113,666,555 labelled augmentations) | Exhaustive. Every connected graph has a non-cut vertex, so each δ ≥ 2 graph on 10 vertices is G + v for a connected 9-vertex G, with v joined to a set S of size ≥ 2 that contains every degree-1 vertex of G. noncol is counted exactly in C (`aug.c`), with early exit above 16. Validated: the same pipeline from n = 8 reproduces the Python histogram at n = 9. | D_10 is the unique minimiser (noncol 12). The next values are 14 (2 graphs) and 16 (4 graphs). No 2-connected graph of diameter ≥ 4 has noncol ≤ 16 |
| "Gadget–path–gadget" graphs H₁ –P– H₂ | all gadgets with ≤ 9 vertices (1,897,173 rooted gadgets), every path length, so every n | Exact formula noncol = I₁ + I₂ + p₁(n−k₁) + p₂(n−k₂), cross-checked on 3000 random instances (`test_gadget_formula.py`). Hence noncol − 2(n−4) = s₁ + s₂ + (p₁+p₂−2)(L−1) + (p₂−1)(k₁−4) + (p₁−1)(k₂−4), with s = I + 4p − k and p ≥ 1. `gadget_scan.py` computes s for all gadgets. `gadget_triangle_check.py` handles the only negative term, which comes from a triangle end (k = 3) | among gadgets with δ ≥ 2 and ≤ 9 vertices, only C_4 has s ≤ 0, and a triangle end always loses (3 + I + 3p − k ≥ 2). So D_n is the unique minimiser in this family, for every n |
| δ ≥ 2 graphs (mode A) and 2-connected graphs of diameter ≥ 4 (mode B) | 10 ≤ n ≤ 30 | Simulated annealing on noncol (`anneal.py`): for n ≤ 22, 4 restarts × 40,000 moves (seed 11); for 23 ≤ n ≤ 30, 3 restarts × 20,000 moves (seed 12). Half of the mode-A starts are relabelled dumbbells. **Heuristic and weak** | No hit (`data/anneal_hits_seed1*.jsonl` empty). Mode A's best equals 2(n−4) (the dumbbell) except at n = 23 and 24, where the search never reached it, which shows its weakness. Mode B's best is far above 2(n−4), e.g. 18 at n = 10 against an exhaustive minimum of ≥ 17, and 383 at n = 30. This says little about the true minimum |

Cross-checks: `test_gadget_formula.py` checks the gadget formula against the full computation. `test_aug_c.py` checks the C noncol and diameter against Python on 27,000 graphs.

## Findings
- **Problem 17(b), n = 9 and n = 10:** D_n is the unique maximiser of Int among graphs without pendant vertices. This is exhaustive.
- **Problem 17(a), n ≤ 10:** the bound holds. At n = 8 the minimum over 2-connected graphs of diameter ≥ 4 is 8 = 2(n−4), attained by three graphs (C_8 and two others). At n = 9 it is 13, against a bound of 10. At n = 10 it is > 16, against a bound of 12.
- **n = 8 has nine extremal graphs, not three.** Nine pairwise non-isomorphic graphs with δ ≥ 2 on 8 vertices have noncol = 8 (Int = 76): D_8, C_8, Q_3, K_{4,4}, and five others. All nine are bipartite with parts 4 + 4. The paper names only D_8, K_{4,4} and C_8, and it does not claim a full list, so this is a completion, not an error. Standalone check: `verify_n8_extremal.py`.
- **Structural fact behind the dumbbell** (proved by the formula, and checked): an end block attached at a cut vertex contributes p·(n − k) non-collinear triples, where p ≥ 1 is the number of pairs of its vertices that are incomparable with respect to the root. Among small blocks, only C_4 achieves the trade-off I + 4p = k.

## Open leads
- Make n = 11 exhaustive. This needs all 11,716,571 connected 10-vertex graphs as a list: write a C canonical-augmentation generator, or use nauty's `geng -c -d2 11` (not available in this sandbox).
- Mode B (2-connected, diameter ≥ 4) shows a large margin, so Problem 17(a) looks true. A proof for 11 ≤ n ≤ 30 would need the Theorem 15 argument made effective for small n.

## Files
- `interval.py`: Int, noncol, distances (definition-level)
- `families.py`: D_n and the two-triangle graph
- `reproduce.py`: paper values
- `gen_graphs.py`, `exhaustive.py`: exhaustive n ≤ 9 (`data/exhaustive_summary_n_le_9.json`, `data/exhaustive_n_le_9.jsonl`)
- `aug.c`, `postproc.py`: exhaustive n = 10 (`data/aug10_hits.txt`, `data/aug10_summary.txt`)
- `gadget.py`, `gadget_scan.py`, `gadget_triangle_check.py`, `test_gadget_formula.py`: gadget family (`data/gadgets_k_le_9.jsonl`)
- `anneal.py`: heuristic screen (`data/anneal_summary_seed1{1,2}.jsonl`, `data/anneal_hits_seed1{1,2}.jsonl`)
- `test_aug_c.py`: C against Python cross-check
- `verify_n8_extremal.py`: standalone check of the nine n = 8 extremal graphs
- `graphs/` (generated, not committed): run `python3 gen_graphs.py 9` first
