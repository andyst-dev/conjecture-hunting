"""Standalone exact verifier (imports nothing from the search code).
For every cycle logged in data/cycles.jsonl and data/large_n.jsonl, check from
scratch that it is a non-transitive cycle in the sense of Easlick
(arXiv:2610.40321, Def. 2) with uniform and minimal winning advantage (Eq. (1)):
  * each die is a non-decreasing n-tuple with faces in {1..n};
  * the dice are pairwise distinct and the cycle has the logged length L;
  * P(D_i beats D_{i+1}) is, for every i (indices mod L), exactly
    1/2 + 1/n^2 (n even) or 1/2 + 1/(2 n^2) (n odd), as a Fraction.
Prints, per n, the cycle lengths certified."""
import json, sys
from fractions import Fraction
from collections import defaultdict

def p_beats(x, y):
    n = len(x)
    return Fraction(sum(1 for a in x for b in y if a > b), n * n)

def minimal(n):
    return Fraction(1, 2) + (Fraction(1, n * n) if n % 2 == 0 else Fraction(1, 2 * n * n))

def check(n, cyc):
    L = len(cyc)
    assert L >= 3
    for d in cyc:
        assert len(d) == n and list(d) == sorted(d) and all(1 <= v <= n for v in d)
    assert len({tuple(d) for d in cyc}) == L
    for i in range(L):
        assert p_beats(cyc[i], cyc[(i + 1) % L]) == minimal(n), (n, i)
    return L

cert = defaultdict(set)
for fn in ('data/cycles.jsonl', 'data/large_n.jsonl'):
    for line in open(fn):
        r = json.loads(line)
        n = r['n']
        cands = [r['cycle']] if 'cycle' in r else [r.get(f'cycle_len_{n}'), r.get(f'cycle_len_{n+1}')]
        for c in cands:
            if c:
                cert[n].add(check(n, c))
bad = []
for n in sorted(cert):
    has_n, has_more = n in cert[n], any(L > n for L in cert[n])
    print(f"n={n:3d}: certified cycle lengths {sorted(cert[n])}  length n: {has_n}  length > n: {has_more}")
    if not has_more: bad.append(n)
print("ALL CHECKS PASSED" if not bad else f"no cycle longer than n certified for n in {bad}")
