"""Completes the gadget argument. With L >= 1 and k1, k2 >= 4, every extra term in
noncol - 2(n-4) = s1 + s2 + (p1+p2-2)(L-1) + (p2-1)(k1-4) + (p1-1)(k2-4)
is >= 0. So only a triangle end (k=3, I=1, p=1, s=2) can contribute a negative term. With H1 = triangle:
noncol - 2(n-4) >= 3 + I2 + 3 p2 - k2, attained at L = 1. Check that this is > 0 for every delta>=2 gadget with <= KMAX vertices."""
import sys
import numpy as np, networkx as nx
from interval import distances, adj_from_edges, noncol
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 9
worst = None
for k in range(3, KMAX + 1):
    for line in open(f"graphs/conn{k}.g6"):
        H = nx.from_graph6_bytes(line.strip().encode())
        if min(d for _, d in H.degree()) < 2: continue
        D = distances(k, adj_from_edges(k, H.edges())); I = noncol(D)
        for r in range(k):
            dr = D[:, r]; M = (D + dr[None, :]) == dr[:, None]; C = M | M.T
            idx = [v for v in range(k) if v != r]
            p = int((~C[np.ix_(idx, idx)][np.triu_indices(k - 1, 1)]).sum())
            t = 3 + I + 3 * p - k
            if worst is None or t < worst[0]: worst = (t, line.strip(), r, I, p)
    print(k, "min 3+I+3p-k so far:", worst, flush=True)
assert worst[0] > 0
print("OK: a triangle end never beats the dumbbell (gadgets <= %d vertices)" % KMAX)
