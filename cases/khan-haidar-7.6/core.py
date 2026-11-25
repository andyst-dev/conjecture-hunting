"""Exact computations for resistance-regular graphs.

L+ = (L + J/n)^{-1} - J/n   (exact rationals via python-flint)
Omega(u,v) = L+_uu + L+_vv - 2 L+_uv
RR  <=> diag(L+) constant
S1(u) = sum_v Omega(u,v),  S2(u) = sum_v Omega(u,v)^2
"""
from flint import fmpq_mat, fmpq
import networkx as nx


def lplus_exact(G):
    nodes = list(G.nodes())
    idx = {v: i for i, v in enumerate(nodes)}
    n = len(nodes)
    M = fmpq_mat(n, n)
    for v in nodes:
        i = idx[v]
        M[i, i] = G.degree(v)
    for u, v in G.edges():
        M[idx[u], idx[v]] = -1
        M[idx[v], idx[u]] = -1
    Jn = fmpq(1, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = M[i, j] + Jn
    X = M.inv()
    for i in range(n):
        for j in range(n):
            X[i, j] = X[i, j] - Jn
    return nodes, X


def analyze(G):
    """Return dict with exact data: is_rr, per-vertex degree, S1, S2."""
    assert nx.is_connected(G)
    nodes, X = lplus_exact(G)
    n = len(nodes)
    diag = [X[i, i] for i in range(n)]
    S1, S2 = [], []
    for i in range(n):
        s1 = fmpq(0)
        s2 = fmpq(0)
        for j in range(n):
            om = X[i, i] + X[j, j] - 2 * X[i, j]
            s1 += om
            s2 += om * om
        S1.append(s1)
        S2.append(s2)
    deg = [G.degree(v) for v in nodes]
    return {
        "nodes": nodes,
        "deg": deg,
        "diag": diag,
        "S1": S1,
        "S2": S2,
        "is_rr": len(set(S1)) == 1,
        "X": X,
    }


def check_conjecture(res):
    """Conjecture: in an RR graph, S2 is a strictly increasing function of degree.
    Returns (holds, details). Checks: (a) S2 constant on each degree class,
    (b) deg u < deg w  =>  S2(u) < S2(w)."""
    by_deg = {}
    for d, s in zip(res["deg"], res["S2"]):
        by_deg.setdefault(d, set()).add(s)
    degs = sorted(by_deg)
    same_class_const = all(len(v) == 1 for v in by_deg.values())
    mono = all(max(by_deg[a]) < min(by_deg[b]) for a, b in zip(degs, degs[1:]))
    return mono, {d: sorted(by_deg[d]) for d in degs}, same_class_const
