"""Reproduce the two explicit witnesses of arXiv:2611.12776 before any search."""
import numpy as np
from gmr import Group, direct_product, evaluate, cyc

def perm(cycles, d):
    p = list(range(d))
    for c in cycles:
        for i in range(len(c)):
            p[c[i]] = c[(i + 1) % len(c)]
    return p

# --- A5 x A5, a=(345), b=(123) (1-based -> 0-based), Theorem 7.2 six sets
d = 5
a5gens = [perm([(0,1,2)], d), perm([(0,1,2,3,4)], d)]
A5 = Group(a5gens)
assert A5.order == 60
a = A5.el(perm([(2,3,4)], d)); b = A5.el(perm([(0,1,2)], d))
inv = A5.inv
G = direct_product(A5, A5); pr = G.pair
H1 = cyc(A5, a); K1 = cyc(A5, b)
sets = [[pr(h1, h2) for h1 in H1 for h2 in H1],
        [pr(inv[b], 0), pr(0, b)],
        [pr(b, 0), pr(0, inv[b])],
        [0, pr(inv[a], 0), pr(0, inv[a])],
        [0, pr(a, 0), pr(0, a)],
        [pr(k1, k2) for k1 in K1 for k2 in K1]]
S, N, R = evaluate(G.mult, sets)
print("A5xA5 Thm 7.2:  |S| =", S, " N =", N, " ratio =", R, float(R))
print("   paper:       |S| = 1053  N = [171, 567, 567, 405, 405, 135]  ratio = 371293/349125")

# --- SmallGroup(243,3) inside S18, n = 10
d = 18
g_a = perm([(6,7,8),(9,15,12),(10,16,13),(11,17,14)], d)
g_b = perm([(0,3,6),(1,4,7),(2,5,8),(15,16,17)], d)
G3 = Group([g_a, g_b])
print("order of <a,b> in S18:", G3.order)
A = G3.el(g_a); B = G3.el(g_b); m = G3.mult
a2 = int(m[A, A]); b2 = int(m[B, B])
sets = [cyc(G3, A), [A, b2], [0, B], [0, A], [0, a2], [B, b2], [B, b2], [0, a2], [0, A], cyc(G3, B)]
S, N, R = evaluate(m, sets)
print("SmallGroup(243,3) n=10: |S| =", S, " N =", N, " ratio =", R, float(R))
print("   paper:       |S| = 189  N1=N10=93, N2..N9=117, ratio = 21^9/(31^2 13^8) ~ 1.01322")
