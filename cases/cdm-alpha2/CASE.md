# Costa–Luu–Wood–Yip, Conjecture CDM (arXiv:2512.17114 v2, Conjecture 1.12 / 10)

**Status:** no counterexample on up to 17 vertices (exhaustive; the paper verified up to 11)
**Date:** 2026-11-25 · **Tools:** C, nauty `geng` 2.8.8 with a custom PRUNE hook, Python/networkx for cross-checks · **AI-assisted (Claude), unrefereed**

## Statement
Source: J. Costa, E. Luu, D. R. Wood and J. H. Yip, *Verifying Hadwiger's Conjecture for Examples of Graphs with α(G) = 2*, arXiv:2512.17114.

- Disjoint vertex sets A and B are **adjacent** if some vertex of A is adjacent to some vertex of B.
- A matching M is **connected** if V(e) and V(f) are adjacent for every two edges e, f ∈ M.
- M is **dominating** if every vertex of V(G) ∖ V(M) is adjacent to V(e) for every e ∈ M.
- **Conjecture CDM.** Every connected graph G with α(G) = 2 has a non-empty connected dominating matching.
- The paper states CDM has been computationally verified for all graphs on up to 11 vertices. CDM implies SHC_{α=2}, which implies HC_{α=2}.

## Reduction (ours, elementary)
Write H for the complement of G. Then α(G) ≤ 2 if and only if H is triangle-free.

- A **dominating edge** xy of G (N(x) ∪ N(y) = V) is a CDM of size 1. In H, it is a non-edge xy with no common neighbour.
- So G has no dominating edge exactly when every non-edge of H has a common neighbour. That is, H is **maximal triangle-free** (MTF; for n ≥ 3 this means triangle-free with diameter ≤ 2).
- G is disconnected exactly when H = K_{a,b}. Indeed, a complete join of ≥ 2 parts inside a triangle-free graph forces 2 independent parts.
- **Counterexamples to CDM on n vertices are therefore complements of MTF graphs on n vertices other than K_{a,b}.** There are far fewer of these than triangle-free graphs (OEIS A216783 against A006785).

## Method
1. **Generation.** nauty `geng -tc` with a PRUNE hook (`mtfprune.c`). A deficient pair of an intermediate induced subgraph (a non-edge with no common neighbour) needs a later common neighbour. Pairs whose union is not independent need different later vertices. Pruning happens when a greedy count of pairwise-incompatible deficient pairs exceeds the number of vertices still to add. This is a necessary condition, so no MTF graph is lost; the counts are checked against OEIS below.
2. **Search.** For each MTF H ≠ K_{a,b}, `cdm.c` searches the complement G exhaustively for a CDM. Start from each edge e₁, maintain the set of required vertices ∪ C(e), where C(ab) = N_H(a) ∩ N_H(b), and branch on the G-edges covering the lowest uncovered required vertex. The search is complete: any valid M extending the partial matching must contain such an edge.

## Reproduction / validation
| Check | Expected | Ours |
|---|---|---|
| Triangle-free graphs from our `geng` build, n = 8..12 (A006785) | 410, 1897, 12172, 105071, 1262180 | identical |
| MTF graphs, n = 3..14 (A216783), via `geng -tc` + filter (n ≤ 13) and via the pruned generator | 1, 2, 3, 4, 6, 10, 16, 31, 61, 147, 392, 1274 | identical with both generators |
| MTF graphs, n = 15, 16, 17 (A216783) | 5036, 25617, 164796 | identical |
| CDM for n ≤ 11 (paper) | holds | holds |

Cross-checks in `test_cdm.py`, all passing:
1. The C search agrees with a brute force over all matchings on 600 random graphs, 244 of which have **no** CDM. So the search is shown to detect absence too.
2. A brute force with no reduction finds CDM for every connected graph with α ≤ 2 on n ≤ 9 vertices (12 + 35 + 104 + 406 + 1893 graphs for n = 5..9).
3. The reduction holds: for n ≤ 10, the number of connected α ≤ 2 graphs without a dominating edge equals the number of MTF graphs other than K_{a,b} (1, 1, 3, 6, 12, 26 for n = 5..10).

## Search coverage
| n | MTF graphs | checked (≠ K_{a,b}) | counterexamples |
|---|---|---|---|
| 5–11 | 3, 4, 6, 10, 16, 31, 61 | 1, 1, 3, 6, 12, 26, 56 | 0 (paper's range, reproduced) |
| 12 | 147 | 141 | 0 |
| 13 | 392 | 386 | 0 |
| 14 | 1,274 | 1,267 | 0 |
| 15 | 5,036 | 5,029 | 0 |
| 16 | 25,617 | 25,609 | 0 |
| 17 | 164,796 (= A216783; split 0/2 + 1/2 = 72,588 + 92,208) | 164,788 | 0 |

The coverage is exhaustive for each n listed (all graphs on n vertices).

## Findings
- CDM holds for **all** graphs on at most **17** vertices; the paper's range was 11.
- The size of the first CDM found grows with n. For example, at n = 12, 49 of the 141 graphs first yield a perfect matching. Small CDMs are not always available, so this is not a trivial case.

## Open leads
- n = 17 took about 37 min on 2 cores (generation dominates). n = 18 (1,337,848 MTF graphs) should take on the order of 10 CPU-hours with `res/mod` splitting. n = 19 (13.7M) needs about a day.
- A dedicated MTF generator (Brandt–Brinkmann–Harmuth) would reach n = 20–23.

## Files
- `cdm.c`: CDM checker; `--raw` prints has_cdm for every input (test mode)
- `mtfprune.c`: geng PRUNE hook. Build: see `build.sh`
- `test_cdm.py`: independent brute-force cross-checks
- `logs/`: per-n outputs
