"""Interval index Int(G) and non-collinear triples noncol(G) (Knor-Sedlar-Skrekovski, arXiv:2610.39966).

Definitions (from the paper):
  I(u,v) = {w : d(u,w) + d(w,v) = d(u,v)}
  Int(G) = sum over unordered pairs {u,v} of (|I(u,v)| - 1)
  a triple {u,v,w} is collinear if one of its vertices lies in the interval of the other two
  noncol(G) = number of non-collinear triples
All quantities are integers; computed exactly from integer BFS distances.
"""
import itertools
from collections import deque
import numpy as np


def distances(n, adj):
    """adj: list of neighbour lists. Returns integer distance matrix (-1 if disconnected)."""
    D = np.full((n, n), -1, dtype=np.int64)
    for s in range(n):
        D[s, s] = 0
        q = deque([s])
        while q:
            x = q.popleft()
            for y in adj[x]:
                if D[s, y] < 0:
                    D[s, y] = D[s, x] + 1
                    q.append(y)
    return D


def adj_from_edges(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v); adj[v].append(u)
    return adj


def between(D):
    """B[u,w,v] = True iff w in I(u,v)."""
    return (D[:, :, None] + D[None, :, :]) == D[:, None, :]


def interval_index_direct(D):
    """Int(G) straight from the definition (sum over unordered pairs of |I(u,v)|-1)."""
    n = D.shape[0]
    B = between(D)
    sizes = B.sum(axis=1)  # sizes[u,v] = |I(u,v)|
    iu = np.triu_indices(n, 1)
    return int((sizes[iu] - 1).sum())


def noncol(D):
    """Number of unordered triples {u,v,w} of distinct vertices with no vertex between the other two."""
    n = D.shape[0]
    if n < 3:
        return 0
    B = between(D)  # B[u,w,v]: w between u and v
    # collinear triple: w between u,v  or u between w,v or v between u,w
    i, j, k = np.array(list(itertools.combinations(range(n), 3))).T
    col = B[i, j, k] | B[j, i, k] | B[i, k, j]
    return int((~col).sum())


def noncol_bruteforce(D):
    n = D.shape[0]
    c = 0
    for a, b, e in itertools.combinations(range(n), 3):
        if D[a, b] + D[b, e] == D[a, e] or D[b, a] + D[a, e] == D[b, e] or D[a, e] + D[e, b] == D[a, b]:
            continue
        c += 1
    return c
