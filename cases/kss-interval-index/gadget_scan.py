"""Scan all rooted gadgets (H, r): H connected on k <= KMAX vertices, every vertex except r has degree >= 2.
For G = H1 -path- H2 (n = k1+k2+L-1, L >= 1):  noncol(G) - 2(n-4) = s1 + s2 + (p1+p2-2)(L-1) + (p2-1)(k1-4) + (p1-1)(k2-4),
with s = I + 4p - k (derived from the gadget formula, cross-checked in test_gadget_formula.py).
The C4 gadget has s = 0, p = 1, k = 4 (giving the dumbbell). A gadget with s < 0 would beat the dumbbell."""
import sys, json, itertools
import numpy as np
import networkx as nx
from interval import distances, adj_from_edges, noncol
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 9
recs = []
best = {}
out = open(f"data/gadgets_k_le_{KMAX}.jsonl", "w")
for k in range(3, KMAX + 1):
    cnt = 0
    for line in open(f"graphs/conn{k}.g6"):
        code = line.strip()
        H = nx.from_graph6_bytes(code.encode())
        deg = np.array([H.degree(v) for v in range(k)])
        low = np.where(deg < 2)[0]
        if len(low) > 1:
            continue
        roots = list(low) if len(low) == 1 else list(range(k))
        D = distances(k, adj_from_edges(k, H.edges()))
        I = noncol(D)
        for r in roots:
            dr = D[:, r]
            M = (D + dr[None, :]) == dr[:, None]   # M[y,x]: x on a geodesic y -> r
            C = M | M.T
            idx = [v for v in range(k) if v != r]
            sub = C[np.ix_(idx, idx)]
            p = int((~sub[np.triu_indices(k - 1, 1)]).sum())
            s = I + 4 * p - k
            cnt += 1
            key = k
            if key not in best or s < best[key][0]:
                best[key] = (s, [])
            if s == best[key][0] and len(best[key][1]) < 10:
                best[key][1].append((code, r, I, p))
            if s <= 0:
                out.write(json.dumps({"k": k, "g6": code, "root": int(r), "I": I, "p": p, "s": s}) + "\n")
    print(k, "rooted gadgets", cnt, "min s", best[k][0], best[k][1][:5], flush=True)
