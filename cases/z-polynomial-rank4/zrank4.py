"""Rank-4 Z-polynomial (Proudfoot–Xu–Young Conj. 5.1: Z_M real-rooted) via flat counts.
For a simple rank-4 matroid: Z = 1 + W3 t + a2 t^2 + W3 t^3 + t^4,
a2 = W2 + sum_{H hyperplane} |H| - sum_{L line} |L|   (|X| = number of points).
(derivation: Z = sum_F t^{rk F} P_{M/F}, P of rank<=2 is 1, P of rank 3/4 = 1 + (W_{r-1} - W_1) t.)
Real-rooted  <=>  4*a2 <= W3^2 + 8  and  a2 >= 2*W3 - 2.
This script computes the ratio 4*a2/(W3^2+8) on random linear rank-4 matroids over GF(p)
and checks the formula against a direct computation of Z from KL polynomials on small cases."""
import itertools, random, sys
import numpy as np

def rank_mod_p(rows, p):
    M = [list(r) for r in rows]; r = 0; cols = len(M[0]) if M else 0
    for c in range(cols):
        piv = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
    return r

def flats_rank4(vecs, p):
    n = len(vecs); rk = {}
    def R(S): 
        S = tuple(sorted(S))
        if S not in rk: rk[S] = rank_mod_p([vecs[i] for i in S], p)
        return rk[S]
    def closure(S):
        r = R(S); return frozenset(i for i in range(n) if R(tuple(S) + (i,)) == r)
    lines = {closure((a, b)) for a, b in itertools.combinations(range(n), 2)}
    planes = set()
    for L in lines:
        a, b = sorted(L)[:2]
        for c in range(n):
            if c not in L: planes.add(closure((a, b, c)))
    return lines, planes

def stats(vecs, p):
    lines, planes = flats_rank4(vecs, p)
    W2, W3 = len(lines), len(planes)
    a2 = W2 + sum(len(H) for H in planes) - sum(len(L) for L in lines)
    return W3, a2

def random_simple(n, p, rng):
    """n pairwise non-parallel nonzero vectors in GF(p)^4 spanning rank 4"""
    while True:
        reps = set(); vecs = []
        while len(vecs) < n:
            v = tuple(rng.randrange(p) for _ in range(4))
            if not any(v): continue
            # normalise first nonzero to 1 for projective point
            f = next(x for x in v if x); inv = pow(f, p - 2, p); w = tuple((x * inv) % p for x in v)
            if w in reps: continue
            reps.add(w); vecs.append(w)
        if rank_mod_p(vecs, p) == 4: return vecs

if __name__ == "__main__":
    rng = random.Random(1); best = (0, None)
    for p in (2, 3, 5, 7):
        maxpts = (p**4 - 1) // (p - 1)
        for n in sorted(set([5, 6, 8, 10, 12, 15, 20, 30, min(40, maxpts), maxpts])):
            if n > maxpts: continue
            for _ in range(15 if n < maxpts else 1):
                vecs = random_simple(n, p, rng)
                W3, a2 = stats(vecs, p)
                ratio = 4 * a2 / (W3 * W3 + 8)
                if ratio > best[0]: best = (ratio, (p, n, W3, a2))
                assert a2 >= 2 * W3 - 2, (p, n, W3, a2)
        print(f"after GF({p}): best ratio 4a2/(W3^2+8) = {best[0]:.4f} at (p,n,W3,a2)={best[1]}", flush=True)
