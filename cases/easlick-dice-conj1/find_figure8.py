import sys, json, random
# run from the case directory: python3 find_figure8.py 24 1
import numpy as np, scipy.sparse.csgraph as cg
from dice import build_graph, beats_count, target
from find_cycle import cycle_of_length
n = 7; D, A = build_graph(n); T = target(n)
nc, lab = cg.connected_components(A, directed=True, connection="strong")
big = np.argmax(np.bincount(lab)); idx = np.nonzero(lab == big)[0]
B0 = A[idx][:, idx].tocsr()
out = lambda u: B0.indices[B0.indptr[u]:B0.indptr[u+1]]
rng = random.Random(int(sys.argv[2]))
L = int(sys.argv[1])
for trial in range(50):
    x = rng.randrange(len(idx))
    tri = [(x, v, w) for v in out(x) for w in out(v) if w != x and x in set(out(w))]
    if not tri: continue
    Bc = list(tri[rng.randrange(len(tri))])
    keep = [u for u in range(len(idx)) if u not in (Bc[1], Bc[2])]
    pos = {u: i for i, u in enumerate(keep)}
    Sub = B0[keep][:, keep].tocsr()
    cyc = cycle_of_length(Sub, L, starts=[pos[x]], budget=3_000_000)
    if cyc:
        Ac = [keep[c] for c in cyc]
        for C in (Ac, Bc):
            dd = [D[idx[c]] for c in C]
            assert all(beats_count(dd[i], dd[(i+1) % len(dd)]) == T for i in range(len(dd)))
        assert Ac[0] == x == Bc[0] and not set(Ac[1:]) & set(Bc[1:])
        print(json.dumps({"A": [D[idx[c]] for c in Ac], "B": [D[idx[c]] for c in Bc]})); break
    print("trial", trial, "no", file=sys.stderr)
