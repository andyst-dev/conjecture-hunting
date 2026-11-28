# Tamura–Palacios–Petojević, asymmetry term of trees: is A(T) ≤ A(P_n)? (arXiv:2610.05663, Sec. 7)

**Status:** no counterexample up to n ≤ 21 (all 3,490,526 trees with 4 ≤ n ≤ 21, exhaustive, exact integer arithmetic)
**Date:** 2026-11-28 · **Tools:** Python 3, networkx, python-flint · **AI-assisted (Claude), unrefereed**

## Statement
Source: S. Tamura, J. L. Palacios, A. Petojević, *A Decomposition of the Hitting Time Index Using the Kirchhoff Index and an Asymmetry Term*, arXiv:2610.05663. The full text was read on 2026-11-28 through a fetcher. No version line was shown, so we treat it as v1.

- Simple random walk, H_uv = E_u[T_v]. The hitting time index is HT(G) = Σ_{{u,v}} max{H_uv, H_vu}, a sum over unordered pairs.
- R_d(u) = Σ_x d(x) R(u,x), where R is effective resistance. The asymmetry term (Def. 3.1) is A(G) = ½ Σ_{{u,v}} |H_uv − H_vu| = ½ Σ |R_d(u) − R_d(v)|.
- Thm 3.1: HT = m·Kf + A. For a tree, Kf = W.
- **Question (Sec. 7):** "whether A(T) ≤ A(P_n) holds for every tree T on n vertices, with equality only for the path P_n." The paper reports no computation for it.

Related: Palacios, Tamura, Wiberg and Petojević (arXiv:2611.34295, 25 Nov 2026) disprove "the largest HT is attained by P_n" over all graphs (lollipops win, from n = 6). For trees, A(T) ≤ A(P_n) together with W(T) ≤ W(P_n) would give HT(T) ≤ HT(P_n). The omega-institute bot preregistered that tree statement on 2026-11-26 (trureturing issue #34417). It reports a Lean proof "awaiting repository publication" and still open registration. The asymmetry question itself is stronger, and we found no claim on it.

**Observation, checked for all trees n ≤ 12 and immediate from the formula:** for a tree, R_d(u) = 2·Tr(u) − (n − 1), where Tr(u) = Σ_v d(u,v). Hence A(T) = Σ_{{u,v}} |Tr(u) − Tr(v)|, the *total transmission irregularity* of T. The question may therefore already be covered by the transmission-irregularity literature. That has not been checked; see Open leads.

## Prior status
- The paper proves A(T) ≥ (n−1)(n−2), with equality iff T is the star (Thm 5.1), and A(T) ≤ (n−2)W(T) (Thm 5.2). It poses the path question without computation.

## Reproduction of the paper
`reproduce.py` prints 77 checks, all OK. Hitting times are computed exactly by solving the linear systems (python-flint `fmpq_mat`), independently of Tetali's formula.

| Object | Paper value | Our value |
|---|---|---|
| A(P_n) (Thm 6.1), n = 2..10 | n²(n²−4)/24 or (n²−1)(n²−3)/24 | identical |
| HT(P_n) (Cor 6.1), n = 2..10 | formula | identical |
| HT(S_n) (Rem 5.2), n = 3..9 | (n−1)(n²−n−1) | identical |
| A(S_n) (equality in Thm 5.1) | (n−1)(n−2) | identical |
| A(K_{p,q}) (Prop 6.1), Kf(K_{p,q}), HT = m·Kf + A | pq\|p−q\|; formula; Thm 3.1 | identical for 1 ≤ p ≤ 4, p ≤ q ≤ 5 |
| HT(P_6), HT(L(6,4)) (arXiv:2611.34295, as quoted in its certificate repo) | 223, 224 | 223, 224 |
| Tree fast path (integer R_d from distances) against exact hitting times | — | 300 random trees, 0 mismatches |

## Search coverage
| Class | Range | Method | Result |
|---|---|---|---|
| All trees (networkx `nonisomorphic_trees`; counts match OEIS A000055) | 4 ≤ n ≤ 21 (2,144,505 trees at n = 21) | exhaustive, exact integers (`search_trees.py`) | A(T) < A(P_n) for every non-path tree. No violation of Thm 5.1 or Thm 5.2 |

Gap at the second-best tree: A(P_n) − max_{T≠P_n} A(T) = 2, 2, 4, 6, 12, 16, 24, 32, 42, 50, 68, 76, 94, 110, 128, 146, 170, 188 for n = 4..21. The runner-up is the path with one leaf moved to the second vertex, a "broom" with a 2-fork at one end (see `argmax_nonpath` in the log). The margin grows, so a counterexample at larger n looks unlikely.

## Findings
- No counterexample for n ≤ 21.
- A(T) equals the total transmission irregularity of T (proof: R_d(u) = 2Tr(u) − (n−1) for trees).

## Open leads
- Search the literature on total transmission irregularity of trees (Dobrynin; Alizadeh–Klavžar; Xu et al.) for a theorem that the path maximises it. If one exists, the question is already answered.
- A proof sketch: moving a leaf toward the end of a longest path increases Σ|Tr(u) − Tr(v)|. This looks provable by a standard transformation argument.

## Files
- `asym.py`: exact A/HT (hitting-time systems, resistances) and the tree fast path
- `reproduce.py`: paper values
- `search_trees.py`: exhaustive tree search, log in `data/trees_summary.jsonl`
