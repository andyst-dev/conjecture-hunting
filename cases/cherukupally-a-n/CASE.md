# Cherukupally, |A(n)| ≤ 3 (arXiv:2603.17434, Conjecture 3.4)

**Status:** no counterexample up to n ≤ 10^27 (exhaustive, exact integers; the abstract claims n < 10^15, the v1 full text n < 10^7)
**Date:** 2026-11-26 · **Tools:** Python 3, sympy · **AI-assisted (Claude), unrefereed**

## Statement
Source: Srikanth Cherukupally, *On the size of {a : 1 ≤ a < n, n | a²−1, a | n²−1} for number n: connections to Fibonacci-like sequences and twin primes*, arXiv:2603.17434. The full text read is v1 (2026-03-19). The paper appeared as a replacement in the math.NT listing of 2026-11-24. Its current abstract claims verification for n < 10^15, while the v1 text says n < 10^7.

- For n > 1, A(n) = {1 ≤ a < n : n | a² − 1 and a | n² − 1}. Here 1 and n − 1 always lie in A(n). This is OEIS A394398.
- G_0 = 1, G_1 = x, G_i = x·G_{i−1} − G_{i−2} (Lemma 2.1). The chain C_z is (G_i(z))_{i ≥ 2}.
- Lemma 3.2 / Cor. 3.3: for n > 2, |A(n)| − 2 = number of chains C_z (z ≥ 3) containing n.
- **Conjecture 3.4.** For every n > 1, |A(n)| ≤ 3. Equivalently, no integer lies in two chains C_{z1}, C_{z2} with z1 > z2 ≥ 3.

## Prior status
- The paper proves |A(n)| < log₂ n (v1, Lemma 3.6) and that the average of |A(n)| is 2 + o(1). It gives no explicit table.
- No later work was found (web search, 2026-11-26).

## Reproduction of the paper
| Object | Paper / OEIS value | Our value |
|---|---|---|
| A394398, n = 1..21 | 0,1,2,2,2,2,2,3,2,2,2,2,2,2,3,2,2,2,2,2,3 | identical |
| A(15), A(35) (twin-prime products) | size 3 (abstract) | {1,4,14}, {1,6,34} |
| Cor. 3.3: \|A(n)\| = 2 + #chains, 3 ≤ n < 20000 | holds | holds for all n (`test_chain.py`) |

## Search coverage
| Class | Range | Method | Result |
|---|---|---|---|
| All n | n ≤ 10^27 | Exhaustive collision search over all chain values G_i(z), z ≥ 3, i ≥ 2, exact integers (`search_chains.py`). The detector was checked by including the trivial chain z = 2, which it flags. | no n with \|A(n)\| ≥ 4 |

Logs: `data/run_1e24.txt`, `data/run_1e27.txt` (830 s for 10^27 on one core).

## Findings
- No counterexample up to 10^27. The full-text check was n < 10^7 and the abstract's was n < 10^15.
- Any counterexample is a solution of G_i(x) = G_j(y) with (i, x) ≠ (j, y) and x, y ≥ 3. The case i = 2, j = 3 is the elliptic curve y² = x³ − 2x + 1. The cases i = 2, j = 4 give Y² = (x² − 1)(x² − 2). Each pair (i, j) has finitely many solutions (Siegel), so a very large sporadic solution cannot be ruled out by search alone.

## Open leads
- Compute the integral points of y² = x³ − 2x + 1 rigorously (Sage/PARI `ellintegralpoints`; not available here). Treat the quartic case i = 2, j = 4 the same way. This would settle the cases most likely to hide a large solution.
- The Bilu–Tichy theory of equal values of polynomials may settle the cases with i, j ≥ 3.

## Files
- `test_chain.py`: definition against chain characterisation and OEIS terms
- `search_chains.py E`: exhaustive search up to 10^E
- `data/`: run logs
