"""Standalone check (no import of the search code): the nine pairwise non-isomorphic connected graphs on 8 vertices
with minimum degree >= 2 found by the exhaustive search all have exactly 8 = 2(n-4) non-collinear triples,
i.e. Int(G) = C(8,2)+C(8,3)-8 = 76 (arXiv:2610.39966 names three of them: D_8, C_8 and K_{4,4})."""
from itertools import combinations
from collections import deque
import networkx as nx

G6 = ["G?r@d_", "G?qa`_", "G?ovD_", "G?ovf_", "G?zTb_", "G?zTf_", "G?zVf_", "G?zvf_", "G?~vf_"]

def dist(G):
    D = {}
    for s in G:
        D[s] = {s: 0}; q = deque([s])
        while q:
            x = q.popleft()
            for y in G[x]:
                if y not in D[s]:
                    D[s][y] = D[s][x] + 1; q.append(y)
    return D

def noncol(G):
    D = dist(G); c = 0
    for a, b, e in combinations(G.nodes(), 3):
        if D[a][b] + D[b][e] == D[a][e] or D[b][a] + D[a][e] == D[b][e] or D[a][e] + D[e][b] == D[a][b]:
            continue
        c += 1
    return c

def interval_index(G):
    D = dist(G); V = list(G.nodes()); s = 0
    for u, v in combinations(V, 2):
        s += sum(1 for w in V if D[u][w] + D[w][v] == D[u][v]) - 1
    return s

graphs = [nx.from_graph6_bytes(c.encode()) for c in G6]
for c, G in zip(G6, graphs):
    assert G.number_of_nodes() == 8 and nx.is_connected(G) and min(d for _, d in G.degree()) >= 2
    assert noncol(G) == 8 and interval_index(G) == 76, c
    print(c, "edges", G.number_of_edges(), "diam", nx.diameter(G), "noncol 8, Int 76")
for i, j in combinations(range(len(graphs)), 2):
    assert not nx.is_isomorphic(graphs[i], graphs[j])
# the three named in the paper
D8 = nx.Graph([(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4)])
named = {"D_8": D8, "C_8": nx.cycle_graph(8), "K_4,4": nx.complete_bipartite_graph(4, 4)}
for name, H in named.items():
    print(name, "is", [c for c, G in zip(G6, graphs) if nx.is_isomorphic(G, H)])
print("VERIFIED: 9 non-isomorphic extremal graphs at n = 8")
