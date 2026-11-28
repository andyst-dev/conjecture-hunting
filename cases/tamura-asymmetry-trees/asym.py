"""Asymmetry term A(G) and hitting time index HT(G) of Tamura-Palacios-Petojevic (arXiv:2610.05663).
Exact rational arithmetic (fractions) via resistance distances from the Laplacian, plus an integer fast path for trees."""
from fractions import Fraction
from itertools import combinations
import networkx as nx
import flint


def _q(x):
    return flint.fmpq(x.numerator, x.denominator)


def resistance_exact(G):
    """Exact effective resistances: R(u,v) = L+_uu + L+_vv - 2 L+_uv via (L + J/n)^-1 (Moore-Penrose for connected G)."""
    nodes = list(G.nodes()); n = len(nodes); idx = {v: i for i, v in enumerate(nodes)}
    L = [[Fraction(0)] * n for _ in range(n)]
    for u, v in G.edges():
        i, j = idx[u], idx[v]
        L[i][j] -= 1; L[j][i] -= 1; L[i][i] += 1; L[j][j] += 1
    Mx = flint.fmpq_mat([[_q(L[i][j] + Fraction(1, n)) for j in range(n)] for i in range(n)]).inv()
    q = lambda x: Fraction(int(x.p), int(x.q))
    M = [[q(Mx[i, j]) - Fraction(1, n) for j in range(n)] for i in range(n)]
    R = {}
    for a in nodes:
        for b in nodes:
            i, j = idx[a], idx[b]
            R[a, b] = M[i][i] + M[j][j] - 2 * M[i][j]
    return R


def hitting_times_exact(G):
    """H[u,v] = expected hitting time of v from u, by solving the linear system exactly (independent of Tetali)."""
    nodes = list(G.nodes()); n = len(nodes)
    H = {}
    for t in nodes:
        others = [v for v in nodes if v != t]; pos = {v: i for i, v in enumerate(others)}
        A = [[Fraction(0)] * (n - 1) for _ in range(n - 1)]
        for v in others:
            A[pos[v]][pos[v]] += 1
            d = G.degree(v)
            for w in G[v]:
                if w != t:
                    A[pos[v]][pos[w]] -= Fraction(1, d)
        x = flint.fmpq_mat([[_q(e) for e in row] for row in A]).solve(flint.fmpq_mat([[1] for _ in range(n - 1)]))
        for v in others:
            e = x[pos[v], 0]
            H[v, t] = Fraction(int(e.p), int(e.q))
        H[t, t] = Fraction(0)
    return H


def A_and_HT_exact(G):
    """Definition: HT = sum_{u<v} max(H_uv, H_vu); A = 1/2 sum_{u<v} |H_uv - H_vu|."""
    H = hitting_times_exact(G)
    HT = sum(max(H[u, v], H[v, u]) for u, v in combinations(G.nodes(), 2))
    A = sum(abs(H[u, v] - H[v, u]) for u, v in combinations(G.nodes(), 2)) / 2
    return A, HT


def A_tree(T):
    """Integer fast path for trees: A = 1/2 sum |R_d(u)-R_d(v)| with R_d(u) = sum_x deg(x) d(u,x)."""
    D = dict(nx.all_pairs_shortest_path_length(T)); deg = dict(T.degree())
    vals = sorted(sum(deg[x] * D[u][x] for x in T) for u in T)
    n = len(vals); s = sum((2 * i - n + 1) * v for i, v in enumerate(vals))   # sum_{i<j} (v_j - v_i)
    W = sum(D[u][v] for u in T for v in T) // 2
    return Fraction(s, 2), W


def A_path(n):
    return Fraction(n * n * (n * n - 4), 24) if n % 2 == 0 else Fraction((n * n - 1) * (n * n - 3), 24)
