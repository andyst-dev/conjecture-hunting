"""Search G_n for a simple directed cycle of exact length L (DFS with
distance-to-start pruning). Results are re-checked by verify_cycle.py."""
import sys, json, random
import numpy as np
import scipy.sparse.csgraph as cg
from dice import build_graph, beats_count, target

def cycle_of_length(A, L, starts, budget=2_000_000, rng=random.Random(1)):
    indptr, indices = A.indptr, A.indices
    AT = A.T.tocsr()
    for s in starts:
        # dist[v] = shortest distance v -> s
        dist = cg.shortest_path(AT, unweighted=True, indices=[s], directed=True)[0]
        path = [s]; onpath = {s}; steps = 0
        def nbrs(u):
            nb = list(indices[indptr[u]:indptr[u+1]]); rng.shuffle(nb); return nb
        stack = [iter(nbrs(s))]
        while stack:
            steps += 1
            if steps > budget: break
            k = len(path)  # vertices on path; need L-k more arcs after the next one
            advanced = False
            for v in stack[-1]:
                if k == L:
                    break
                if v == s:
                    continue
                if v in onpath:
                    continue
                if dist[v] > L - k:  # after moving to v we have L-k arcs left incl. return
                    continue
                path.append(v); onpath.add(v)
                if len(path) == L:
                    # need arc v -> s
                    if s in set(indices[indptr[v]:indptr[v+1]]):
                        return path
                    onpath.discard(path.pop())
                    continue
                stack.append(iter(nbrs(v))); advanced = True
                break
            if not advanced:
                stack.pop()
                if len(path) > 1 or not stack:
                    if path and path[-1] != s:
                        onpath.discard(path.pop())
    return None

if __name__ == '__main__':
    n = int(sys.argv[1]); Ls = [int(x) for x in sys.argv[2:]]
    D, A = build_graph(n)
    ncomp, lab = cg.connected_components(A, directed=True, connection='strong')
    big = np.argmax(np.bincount(lab)); idx = np.nonzero(lab == big)[0]
    B = A[idx][:, idx].tocsr()
    T = target(n)
    for L in Ls:
        cyc = cycle_of_length(B, L, starts=list(range(0, len(idx), max(1, len(idx)//20))))
        rec = dict(n=n, L=L, found=cyc is not None)
        if cyc:
            dd = [D[idx[c]] for c in cyc]
            assert len(set(dd)) == L
            assert all(beats_count(dd[i], dd[(i+1) % L]) == T for i in range(L))
            rec['cycle'] = dd
        print(json.dumps(rec), flush=True)
        with open('data/cycles.jsonl', 'a') as f: f.write(json.dumps(rec) + '\n')
