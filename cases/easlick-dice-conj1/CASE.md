# Easlick, non-transitive dice with faces in [n] (arXiv:2610.40321, Conjecture 1)

**Status:** proved (partial: for every 7 ≤ n ≤ 60 and n ∈ {70, 80, 90, 100}, under the reading the paper's own evidence uses; see "Statement"). Under the other literal reading the conjecture already fails at n = 7.
**Date:** 2026-11-27 · **Tools:** Python 3, numpy/scipy integer arithmetic, `fractions.Fraction` in the verifier · **AI-assisted, unrefereed**

## Statement
Source: T. Easlick, *Computer-Assisted Results for n-Sided Non-Transitive Dice Taking Face Values in [n]*, arXiv:2610.40321 (document date 30 Oct 2026; version number not shown in the fetched PDF, so treat it as v1). The text was read through a web fetcher in fragments, so the quotes below are short.

- **Definition 1.** "We call a non-standard n-sided die to be a die taking face values in [n] such that face values may be repeated." Each face has probability 1/n, and a die is written as a non-decreasing n-tuple. The standard die (1, …, n) "will also occur in the set of non-standard dice", so there are binom(2n−1, n) dice.
- **Beats.** D_x beats D_y when an independent toss of both gives a strictly larger face on D_x, so P(D_x beats D_y) = #{(i,j) : x_i > y_j} / n². Ties count for neither die.
- **Definition 2.** "A non-transitive cycle of length k ≥ 3 is a sequence of distinct non-standard n-sided dice" D_1, …, D_k with P(D_i beats D_{i+1}) > 1/2 for every i, where D_{k+1} = D_1. The cycle has a *uniform* advantage if all these probabilities are equal. The advantage is *minimal* if it is "the smallest attainable probability greater than 1/2 for the given n".
- **Eq. (1).** The minimal value is 1/2 + 1/n² for even n and 1/2 + 1/(2n²) for odd n. In integers, each arc has count exactly ⌊n²/2⌋ + 1.
- **Graph G_n.** The vertices are the dice, with an arc X → Y iff the count is ⌊n²/2⌋ + 1. Uniform minimal cycles are exactly the simple directed cycles of G_n.
- **Theorem 1.** "n=6 is the smallest n≥3 such that the set of non-standard n-sided dice admits a non-transitive cycle of maximum length n with uniform and minimal winning advantage."
- **Conjecture 1.** "n = 6 is the largest n ≥ 3 such that […] admits a non-transitive cycle of maximum length n with uniform and minimal winning advantage."

**Ambiguity.** The paper never defines "maximum length n". There are two readings:
- **(A)** the longest uniform minimal cycle in G_n has length exactly n;
- **(B)** G_n has a uniform minimal cycle of length n.

The paper's evidence only makes sense under (A). The proof for n = 6 rules out cycles of length 7. Section 4 also says: "For n = 7, 8, 9 and 10, the corresponding directed graphs contain cycles of length strictly greater than n, providing computational evidence for Conjecture 1." Under (A), Conjecture 1 says that for every n ≥ 7, G_n has a cycle of length > n, or has no cycle at all.

## Prior status
The paper proves Theorem 1 by exhaustive computation for n = 3, …, 6. For n = 7, …, 10 it gives one cycle of length ≥ n + 1, found by a depth-first search that stops at the first such cycle. For n ≥ 11 it only "surmises" such cycles exist. The conclusion says the general case "remains open". Searches on 27 Nov 2026 found no follow-up work.

## Reproduction of the paper
| Object | Paper value | Our value |
|---|---|---|
| n=3: V, E, SCCs | 10, 3, 10 (all trivial) | 10, 3, 10 (all trivial) |
| n=4: V, E, SCCs | 35, 67, 35 | 35, 67, 35 |
| n=5: V, E, SCCs | 126, 620, 126 | 126, 620, 126 |
| n=6: V, E, SCCs | 462, 5491, 450 | 462, 5491, 450 |
| n=6: non-trivial SCCs | 2, each 7 vertices / 10 arcs | 2, each 7 vertices / 10 arcs |
| n=6: SCC vertex lists | as listed in §2 | identical (`data/reproduce.jsonl`) |
| n=6: simple cycles per SCC | 4: one of length 6, three of length 3 | lengths [3, 3, 3, 6] in each |
| n=7: V, E, SCCs, largest SCC | 1716, 89153, 636, 1081 V / 53389 E | 1716, 89153, 636, 1081 / 53389 |
| n=8: V, E, largest SCC | 6435, 956307, 4321 V / 674108 E | 6435, 956307, 4321 / 674108 (2115 SCCs) |
| n=9: V, E, SCCs, largest SCC | 24310, 12528131, 3860, 20451 V / 10933011 E | 24310, 12528131, 3860, 20451 / 10933011 |
| n=10: V, E, largest SCC | 92378, 157244770, 81791 V / 147620052 E | 92378, 157244770, 81791 / 147620052 (10588 SCCs) |

Every value agrees. The n = 6 row is Theorem 1 under reading (A): the longest cycle in G_6 has length 6, and G_3, G_4, G_5 have no cycles.

**Cross-checks** (`test_counts.py`, all pass):
- the fast arc count (matrix product of multiplicity and "faces below" vectors) equals the direct double sum on the full arc sets for n = 3..6 and on 3000 random pairs for each n = 3..12;
- the paper's constructions attain ⌊n²/2⌋ + 1 for n = 3..39.

## Search coverage
| Class | Range | Method | Result |
|---|---|---|---|
| Full G_n | n = 3..10 | exact integer arc computation, SCCs via scipy | statistics above |
| Simple cycles of length L, full graph | n = 7, L = 3..8; n = 8, L = 8, 9; n = 9, L = 9, 10 | DFS with distance-to-start pruning (`find_cycle.py`) | found for every L tried |
| Induced subgraph on 20,000 random dice (seed 0) | every n = 10..60, and n = 70, 80, 90, 100 | exact arcs on the sample, then DFS for cycles of length n and n+1 (`search_sampled.py`) | both lengths found for every n |

All cycles are in `data/cycles.jsonl` (full graphs) and `data/large_n.jsonl` (samples). A cycle in an induced subgraph is a cycle of G_n.

## Findings
1. **Reading (A): proved for 7 ≤ n ≤ 60 and n ∈ {70, 80, 90, 100}.** For each such n we give an explicit uniform minimal non-transitive cycle of length n + 1. So the longest such cycle is not n, and n is not a case of the property. Together with the paper's Theorem 1, the property "longest cycle has length exactly n" holds for n = 6 and for no other n ≤ 60. `python3 verify_cycles.py` re-checks every cycle from scratch with `Fraction` arithmetic and imports none of the search code.
2. **Reading (B): false at n = 7.** G_7 contains the uniform minimal 7-cycle
   (1,1,1,3,3,4,7) → (2,2,2,2,2,2,6) → (1,1,1,4,5,5,5) → (2,3,3,3,4,4,4) → (1,2,2,3,3,6,7) → (1,1,1,3,5,5,6) → (2,2,2,2,3,4,4) → back. Each arc has count 25 out of 49, so probability 1/2 + 1/98.
   Under (B) we also find cycles of length n for every n tested. Since the paper's evidence contradicts (B), this is almost certainly not the intended meaning. Recording it is useful only because the statement as printed is ambiguous. **It is not a refutation to announce.**
3. Cycles of every length 3..8 already exist in G_7.

## Open leads
- A uniform construction proving reading (A) for all n ≥ 7, for example by lifting a cycle of G_n to G_{n+2}. The random samples suggest that cycles of length > n are everywhere once n ≥ 7.
- Exact longest-cycle lengths of G_7 (SCC of 1081 vertices). This is a Hamiltonian-type question and may be within reach with a SAT solver.

## Files
- `dice.py`: dice, exact arc counts, CSR graph
- `reproduce.py`: paper statistics (`data/reproduce.jsonl`)
- `find_cycle.py`: exact-length cycle DFS on the full graph (`data/cycles.jsonl`)
- `search_sampled.py`: induced-subgraph search for large n (`data/large_n.jsonl`); `search_large_n.py` is an earlier variant for n ≤ 12
- `test_counts.py`: cross-checks of the fast path
- `verify_cycles.py`: standalone exact verifier
