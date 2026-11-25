"""r = 2 case of Benedetti–Menara 7.1, graph by graph.

A flag 1-layered family of 2-sets is a triangle-free graph G together with a surjective edge
labelling l: E -> {0..L-1} such that every vertex sees an interval of labels.
(The interval condition for S = a vertex is the only nontrivial one when r = 2.)
Conjecture: L <= m - 1. A counterexample needs L = m (merging adjacent labels is monotone).

Reductions:
 * connected G suffices: each component's labels form an interval, so a disconnected G gives
   L <= sum (m_k - 1) < m - 1 once the components satisfy the bound (induction on m);
 * L <= |E|, so G needs >= m edges.
Cheap bound (path bound): labels at v span at most deg(v) - 1, so for edges e, f joined by a walk
e = f_0, f_1, ..., f_k = f in the line graph through vertices v_1..v_k,
|l(e) - l(f)| <= sum (deg(v_i) - 1). If B(G) = max over e, f of these weighted distances is <= m - 2,
then L <= m - 1 and G is skipped. Otherwise a small SAT instance decides L = m exactly.
Usage: geng -tc m m:M | python3 r2_graphs.py m
"""
import sys, itertools, heapq, json
import networkx as nx
from pysat.formula import IDPool, CNF
from pysat.solvers import Cadical153


def path_bound(G):
    E = list(G.edges())
    idx = {frozenset(e): i for i, e in enumerate(E)}
    adj = [[] for _ in E]
    for v in G.nodes():
        inc = [idx[frozenset((v, u))] for u in G[v]]
        w = G.degree(v) - 1
        for a in inc:
            for b in inc:
                if a != b:
                    adj[a].append((b, w))
    best = 0
    for s in range(len(E)):
        dist = [None] * len(E); dist[s] = 0; pq = [(0, s)]
        while pq:
            d, a = heapq.heappop(pq)
            if d > dist[a]: continue
            for b, w in adj[a]:
                nd = d + w
                if dist[b] is None or nd < dist[b]:
                    dist[b] = nd; heapq.heappush(pq, (nd, b))
        best = max(best, max(dist))
    return best


def sat_labelling(G, L):
    pool = IDPool(); cnf = []
    E = list(G.edges())
    x = {(i, t): pool.id((i, t)) for i in range(len(E)) for t in range(L)}
    for i in range(len(E)):
        cnf.append([x[i, t] for t in range(L)])
        for a, b in itertools.combinations(range(L), 2):
            cnf.append([-x[i, a], -x[i, b]])
    for t in range(L):
        cnf.append([x[i, t] for i in range(len(E))])
    for v in G.nodes():
        inc = [i for i, e in enumerate(E) if v in e]
        y = [pool.id(("y", v, t)) for t in range(L)]
        for t in range(L):
            cnf.append([-y[t]] + [x[i, t] for i in inc])
            for i in inc: cnf.append([-x[i, t], y[t]])
        # no gap: y[t1] & y[t3] -> y[t2]
        for t1 in range(L):
            for t3 in range(t1 + 2, L):
                for t2 in range(t1 + 1, t3):
                    cnf.append([-y[t1], -y[t3], y[t2]])
    # symmetry: some edge has label 0 (always true for surjective) - fine as is
    with Cadical153(bootstrap_with=cnf) as s:
        if not s.solve(): return None
        mdl = set(l for l in s.get_model() if l > 0)
    return {f"{u}-{v}": t for i, (u, v) in enumerate(E) for t in range(L) if x[i, t] in mdl}


if __name__ == "__main__":
    m = int(sys.argv[1]); L = m
    n_in = n_pruned = n_sat = 0
    for line in sys.stdin:
        line = line.strip()
        if not line or line.startswith(">"): continue
        n_in += 1
        G = nx.from_graph6_bytes(line.encode())
        if G.number_of_edges() < L:
            n_pruned += 1; continue
        if path_bound(G) <= m - 2:
            n_pruned += 1; continue
        n_sat += 1
        lab = sat_labelling(G, L)
        if lab is not None:
            print("COUNTEREXAMPLE", line, json.dumps(lab), flush=True)
    print(f"m={m}: graphs={n_in} pruned_by_bound={n_pruned} sat_checked={n_sat}", flush=True)
