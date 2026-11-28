# Easlick, non-transitive dice with faces in [n] (arXiv:2610.40321, Conjecture 1)

**Status:** proved (partial: in full for **every n ≥ 7** under reading (A), the reading the paper's own evidence uses; see "Statement" and "General construction"). Under the other literal reading (B) the conjecture already fails at n = 7.
**Date:** 2026-11-27, general construction 2026-11-28 · **Tools:** Python 3, numpy/scipy integer arithmetic, `fractions.Fraction` in the verifier · **AI-assisted, unrefereed**

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
0. **Reading (A): proved for every n ≥ 7** (2026-11-28). An explicit block construction gives a uniform minimal cycle of length about 3n in G_n for every n ≥ 49; see "General construction". Together with item 1, which covers 7 ≤ n ≤ 60, this settles Conjecture 1 under reading (A). Unrefereed, AI-assisted.
1. **Reading (A): proved for 7 ≤ n ≤ 60 and n ∈ {70, 80, 90, 100}** (2026-11-27, by sampling). For each such n we give an explicit uniform minimal non-transitive cycle of length n + 1. So the longest such cycle is not n, and n is not a case of the property. Together with the paper's Theorem 1, the property "longest cycle has length exactly n" holds for n = 6 and for no other n ≤ 60. `python3 verify_cycles.py` re-checks every cycle from scratch with `Fraction` arithmetic and imports none of the search code.
2. **Reading (B): false at n = 7.** G_7 contains the uniform minimal 7-cycle
   (1,1,1,3,3,4,7) → (2,2,2,2,2,2,6) → (1,1,1,4,5,5,5) → (2,3,3,3,4,4,4) → (1,2,2,3,3,6,7) → (1,1,1,3,5,5,6) → (2,2,2,2,3,4,4) → back. Each arc has count 25 out of 49, so probability 1/2 + 1/98.
   Under (B) we also find cycles of length n for every n tested. Since the paper's evidence contradicts (B), this is almost certainly not the intended meaning. Recording it is useful only because the statement as printed is ambiguous. **It is not a refutation to announce.**
3. Cycles of every length 3..8 already exist in G_7.

## General construction (2026-11-28): every n ≥ 49
Notation: g_n(X, Y) = 2·W(X, Y) − n², where W counts the pairs in which X shows the larger face. An arc of G_n is exactly g_n = δ(n), with δ(n) = 1 for odd n and 2 for even n, since ⌊n²/2⌋ + 1 = (n² + δ(n))/2.

**Lemma 1 (direct sums).** Let X = X₁ ⊕ … ⊕ X_s. Here block X_j is an m_j-die placed on the value range m₁+…+m_{j−1}+1 … m₁+…+m_j, and N = Σ m_j. Then g_N(X, Y) = Σ_j g_{m_j}(X_j, Y_j).
*Proof:* every face of a higher block beats every face of a lower block, so W(X, Y) = Σ_j W(X_j, Y_j) + Σ_{i<j} m_i m_j. Substitute into 2W − N².
For a block that stays fixed, g(Z, Z) = −Σ_v mult_Z(v)², because 2W(Z, Z) + #ties = m². So a fixed (1) costs 1, a fixed (1,1) costs 4 and a fixed (1,1,1) costs 9. Both formulas are cross-checked against direct counts in `test_additivity.py`.

**Ingredients** (`data/figure8_g7.json`, found by `find_figure8.py`):
- in G_7, a simple cycle A = (A₀ … A₂₃) of length a = 24 and a simple 3-cycle B = (B₀, B₁, B₂), with A₀ = B₀ = x = (1,2,2,4,7,7,7) and no other common vertex;
- in G_6, a 3-cycle C6 = (1,1,2,5,5,5) → (1,3,3,4,4,4) → (2,2,3,3,3,6).

**The cycle for N.** Let d = δ(N), K = ⌈(N + d)/8⌉ and e = 8K − N − d ∈ {0, 2, 4, 6}. The excess blocks E_e are: nothing for e = 0, one (1,1) for e = 2, two (1,1) for e = 4, one (1,1,1) for e = 6. Their size is s_e = 0, 2, 4, 3 and their total cost is s_e + e. Choose K₆ ∈ {0, 1, 2, 3}, the smallest that works. Set K₇ = K − K₆, r = N − 7K₇ − 6K₆ − s_e ≥ 0, and P = 24(K₇ − 1) + 3. Die number t, for t = 0 … P−1, is the direct sum of:
- K₇ *counter* 7-blocks. Block i (i = 0 … K₇−1) walks A^i B A^(K₇−1−i): at time t it is at A_{t mod 24} if t < 24i, at B_{t−24i} if 24i ≤ t < 24i + 3, and at A_{(t−24i−3) mod 24} otherwise;
- K₆ *passenger* 6-blocks, all at C6_{t mod 3} (3 divides P);
- the fixed blocks E_e and r fixed 1-blocks (1).

**Lemma 2 (arcs).** Consecutive dice, including die P−1 → die 0, form an arc of G_N.
*Proof:* each counter block moves along an arc of A or B, or along A₂₃ → x or B₂ → x, so it contributes δ(7) = 1. Each passenger contributes δ(6) = 2. The fixed blocks contribute −(s_e + e) − r. The total is K₇ + 2K₆ − s_e − e − r = 8K − N − e = d.

**Lemma 3 (distinct dice).** The P dice are pairwise distinct.
*Proof:* the block ranges are fixed, so a die determines its blocks. It is enough to recover t from the states of the counters. The windows [24i, 24i + 3) are disjoint.
- (i) Some counter i is at B₁ or B₂, which are not on A. Then t = 24i + (its index in B).
- (ii) Otherwise every counter is on A. Counter i shows A_ρ with ρ = t mod 24 while t ≤ 24i, and A_σ with σ = (t − 3) mod 24 once t ≥ 24i + 3. Since 3 ≢ 0 (mod 24), A_ρ ≠ A_σ.
  - So the state vector is (A_σ, …, A_σ, A_ρ, …, A_ρ), with q leading entries A_σ, where q = #{i : 24i + 3 ≤ t}.
  - The last counter always shows A_ρ, because t < P. It is constant only when q = 0, that is t = 0.
  - Otherwise q is the position where the value changes, and t lies in [24(q−1) + 3, 24q], an interval of 22 < 24 consecutive integers, so ρ fixes t.

**Lemma 4 (parameters).** For N ≥ 111, K₆ = 0 works: S := N − 7K ≥ (N − 7d − 49)/8 ≥ 6 ≥ s_e and P = 24(K − 1) + 3 ≥ 3N − 18 > N. For 49 ≤ N ≤ 110, the rule above also finds valid parameters; this finite check is run by `construct_all_n.py`, which reports no N ≥ 49 without a plan up to 100,000.

**Theorem.** For every N ≥ 49, G_N contains a uniform minimal non-transitive cycle of length P_N > N, with P_N ≈ 3N.
Together with the certified cycles of length n + 1 for 7 ≤ n ≤ 60 (`verify_cycles.py`), every G_n with n ≥ 7 has a cycle longer than n. Under reading (A), Conjecture 1 therefore holds: by Theorem 1, n = 6 is the only n ≥ 3 whose longest uniform minimal cycle has length exactly n.

**Exact verification** (`verify_construction.py`, standalone; it imports nothing from this case):
- It rebuilds the cycle for each N and checks the arc counts with the plain pair count, not with Lemma 1, together with distinctness, face ranges and P > N.
- Run on 2026-11-28 for every N in 49…400 and for N ∈ {500, 777, 1000, 1001, 2000, 2001}; all pass.
- An extended run covered every N in 49…1200 and N = 5000 (P = 15,003); all pass (`data/verify_construction_49_1200.log`).

## Open leads
- ~~A uniform construction proving reading (A) for all n ≥ 7~~: done (see "General construction"). The proof is computer-assisted only through finitely many arc counts in G_7 and G_6 and the finite range 49 ≤ N ≤ 110; it is unrefereed.
- Contact the author (T. Easlick) to ask which reading is intended before stating the result publicly.
- Exact longest-cycle lengths of G_7 (SCC of 1081 vertices). This is a Hamiltonian-type question and may be within reach with a SAT solver.

## Files
- `dice.py`: dice, exact arc counts, CSR graph
- `reproduce.py`: paper statistics (`data/reproduce.jsonl`)
- `find_cycle.py`: exact-length cycle DFS on the full graph (`data/cycles.jsonl`)
- `search_sampled.py`: induced-subgraph search for large n (`data/large_n.jsonl`); `search_large_n.py` is an earlier variant for n ≤ 12
- `test_counts.py`: cross-checks of the fast path
- `verify_cycles.py`: standalone exact verifier
- `find_figure8.py`: finds the G_7 figure-eight (`data/figure8_g7.json`)
- `construct_all_n.py`: parameters of the general construction (`data/construction.jsonl`, N = 7…2001)
- `test_additivity.py`: cross-check of Lemma 1 and the fixed-block cost
- `verify_construction.py`: standalone exact verifier of the general construction
