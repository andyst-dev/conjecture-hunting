# Wild number of minimum proper colourings of K_n (arXiv:2610.22142, Conjecture 6.1)

**Status:** no counterexample among the colourings tested (sampling, not exhaustive over colourings)
**Date:** 2026-11-25 · **Tools:** PySAT (CaDiCaL), Python · **AI-assisted (Claude), unrefereed**

## Statement
- An edge colouring γ of K_n is fixed. A set W of edges is **colour-connecting** if G_i ∪ W is connected for every colour class G_i.
- The **wild number** wld_γ(K_n) is the minimum |W|.
- **Conjecture 6.1.** For every n ≥ 3 and **every** minimum proper edge colouring ρ of K_n (χ′ colours), wld_ρ(K_n) = ⌊4⌈n/2⌉/3⌋.
- The paper proves the lower bound for all ρ and equality for 3 ≤ n ≤ 12, by constructions (no computer enumeration). It also proves wld ≤ n − 2 (odd n ≥ 7) and ≤ n − 3 (even n ≥ 8).

## Method
- **Exact test** (`wild.py`). "wld ≤ k" is decided by SAT. Connectivity of each G_i ∪ W is encoded by BFS layers from vertex 0, and |W| ≤ k by a sequential counter. Every SAT witness is re-checked with union-find.
- **Colourings.**
  - the standard rotational colouring GK_n;
  - Kempe-switch randomisations of it (5 to 500 random (a,b)-chain swaps, which preserve properness and the colour count);
  - the affine colouring i ⊕ j of K_{2^k}, which has very few rainbow K₄'s, and its vertex-deleted version on K_{2^k − 1}.
- A failure of "wld ≤ ⌊4⌈n/2⌉/3⌋" would be a counterexample. Since the paper proves the lower bound, only the upper side can fail.

## Reproduction
| n | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|
| paper (Thm 3.23) | 2 | 2 | 4 | 4 | 5 | 5 | 6 | 6 | 8 |
| ours: exact wld on GK + 5 randomised colourings each | 2 | 2 | 4 | 4 | 5 | 5 | 6 | 6 | 8 |

Here "exact" means SAT at k and UNSAT at k − 1. This independently re-derives the paper's lower bound on these colourings (`reproduce.py`).

## Search coverage
| n | conjectured wld | colourings tested | result |
|---|---|---|---|
| 13, 14 | 9 | 50 each | all ≤ 9 |
| 15, 16 | 10 | 40 each + affine (XOR) colouring | all ≤ 10 |
| 17, 18 | 12 | 40 each | all ≤ 12 |
| 23 | 16 | 12 | all ≤ 16 |
| 4, 7, 8 | — | affine / vertex-deleted affine | match |

Each instance takes seconds for n ≤ 18 and about 40 s at n = 23. The affine colourings of K_31 and K_32 did not finish within 9 minutes.

## Findings
- No counterexample. Every tested colouring attains the conjectured value, including the structured affine colourings, which are the natural candidates for few rainbow K₄'s.
- The search samples colourings near GK (Kempe chains). It is not uniform over the astronomically many 1-factorisations, so this is evidence only.

## Files
`wild.py`, `reproduce.py`, `hunt.py` (`python3 hunt.py n trials seed`), `special.py`, `logs/`
