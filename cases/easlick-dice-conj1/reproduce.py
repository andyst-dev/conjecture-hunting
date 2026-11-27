"""Reproduce the paper's graph statistics (Theorem 1, Section 3, Appendices B-F)."""
import sys, json, time
import numpy as np
import scipy.sparse.csgraph as cg
import networkx as nx
from dice import build_graph

def stats(n):
    t0 = time.time()
    V = None
    blk = max(50, int(2e8 // (8 * 2 * 10**5)) if n >= 10 else 2000)
    D, A = build_graph(n, block=blk)
    V = A.shape[0]; E = A.nnz
    ncomp, lab = cg.connected_components(A, directed=True, connection='strong')
    sizes = np.bincount(lab)
    big = np.argmax(sizes)
    idx = np.nonzero(lab == big)[0]
    Ebig = A[idx][:, idx].nnz
    nontriv = [(int(s), int(A[np.nonzero(lab == c)[0]][:, np.nonzero(lab == c)[0]].nnz))
               for c, s in enumerate(sizes) if s > 1]
    out = dict(n=n, V=V, E=E, SCCs=int(ncomp), largest_scc_V=int(sizes[big]),
               largest_scc_E=int(Ebig), nontrivial=len(nontriv), secs=round(time.time()-t0,1))
    if n == 6:
        cyc = []
        for c, s in enumerate(sizes):
            if s > 1:
                ids = np.nonzero(lab == c)[0]
                G = nx.DiGraph()
                sub = A[ids][:, ids].tocoo()
                G.add_edges_from((D[ids[a]], D[ids[b]]) for a, b in zip(sub.row, sub.col))
                cyc.append(sorted(len(x) for x in nx.simple_cycles(G)))
                out.setdefault('scc_vertices', []).append(sorted(D[i] for i in ids))
        out['cycle_lengths_per_nontrivial_scc'] = cyc
        out['nontrivial_sizes_edges'] = nontriv
    return out

if __name__ == '__main__':
    for n in map(int, sys.argv[1:]):
        r = stats(n)
        print(json.dumps(r), flush=True)
        with open('data/reproduce.jsonl', 'a') as f:
            f.write(json.dumps(r) + '\n')
