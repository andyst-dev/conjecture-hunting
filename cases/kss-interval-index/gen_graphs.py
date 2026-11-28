"""Generate all graphs on n <= NMAX vertices (vertex augmentation + pynauty certificates),
check counts against OEIS A000088, write connected ones as graph6 to graphs/conn<n>.g6."""
import sys, os, itertools
import networkx as nx
import pynauty
from interval import distances
NMAX = int(sys.argv[1])
A000088 = {1: 1, 2: 2, 3: 4, 4: 11, 5: 34, 6: 156, 7: 1044, 8: 12346, 9: 274668, 10: 12005168}
os.makedirs("graphs", exist_ok=True)
def cert(n, adj):
    return pynauty.certificate(pynauty.Graph(n, adjacency_dict={i: a for i, a in enumerate(adj)}))
level = {cert(1, [[]]): [[]]}
for n in range(1, NMAX + 1):
    if n > 1:
        new = {}
        for adj in level.values():
            m = n - 1
            for r in range(m + 1):
                for S in itertools.combinations(range(m), r):
                    a2 = [list(x) for x in adj] + [list(S)]
                    for s in S: a2[s].append(m)
                    c = cert(n, a2)
                    if c not in new: new[c] = a2
        level = new
    assert len(level) == A000088[n]
    with open(f"graphs/conn{n}.g6", "w") as f:
        for adj in level.values():
            if (distances(n, adj) < 0).any(): continue
            G = nx.Graph(); G.add_nodes_from(range(n))
            G.add_edges_from((u, v) for u in range(n) for v in adj[u] if u < v)
            f.write(nx.to_graph6_bytes(G, header=False).decode())
    print(n, len(level), flush=True)
