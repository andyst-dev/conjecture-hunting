"""For n >= 11 (graph too big to build whole): take a random set S of dice,
build the induced subgraph of G_n on S exactly, and look for simple cycles of
length n and n+1 in its largest SCC. A cycle in an induced subgraph is a cycle
in G_n. Each hit is asserted with the direct definition and logged."""
import sys, json, random, time
import numpy as np
import scipy.sparse.csgraph as cg
from math import comb
from dice import build_graph, beats_count, target, dice
from find_cycle import cycle_of_length

n = int(sys.argv[1]); size = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
t0 = time.time()
rng = np.random.default_rng(seed)
V = comb(2 * n - 1, n)
S = np.sort(rng.choice(V, size=min(size, V), replace=False))
D, B = build_graph(n, block=1000, rows=S, cols=S)
ncomp, lab = cg.connected_components(B, directed=True, connection='strong')
sizes = np.bincount(lab); big = np.argmax(sizes); idx = np.nonzero(lab == big)[0]
C = B[idx][:, idx].tocsr()
T = target(n)
rec = dict(n=n, V=V, sample=len(S), seed=seed, sub_E=int(B.nnz), largest_scc=int(sizes[big]))
for L in (n, n + 1):
    cyc = cycle_of_length(C, L, starts=list(range(0, len(idx), max(1, len(idx)//10)))) if len(idx) >= L else None
    if cyc:
        dd = [D[S[idx[c]]] for c in cyc]
        assert len(set(dd)) == L
        assert all(beats_count(dd[i], dd[(i + 1) % L]) == T for i in range(L))
        rec[f'cycle_len_{L}'] = dd
    else:
        rec[f'cycle_len_{L}'] = None
rec['secs'] = round(time.time() - t0, 1)
print(json.dumps(rec), flush=True)
with open('data/large_n.jsonl', 'a') as f: f.write(json.dumps(rec) + '\n')
