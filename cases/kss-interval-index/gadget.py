"""Fast path for 'gadget - path - gadget' graphs.

G = H1 (root r1)  --- path of L>=1 edges ---  H2 (root r2),  n = k1 + k2 + L - 1.
Claim (checked against the full computation in test_gadget_formula.py):
  noncol(G) = I(H1) + I(H2) + p(H1,r1)*(n-k1) + p(H2,r2)*(n-k2)
where I(H) = noncol(H) and p(H,r) = #unordered pairs {x,y} in H\\{r} with neither x on a
geodesic y->r nor y on a geodesic x->r (distances inside H).
"""
import itertools
import networkx as nx
from interval import distances, adj_from_edges, noncol


def gadget_params(H, r):
    H = nx.relabel_nodes(H, {v: i for i, v in enumerate(H.nodes())}) if set(H.nodes()) != set(range(H.number_of_nodes())) else H
    k = H.number_of_nodes()
    D = distances(k, adj_from_edges(k, H.edges()))
    I = noncol(D)
    p = 0
    for x, y in itertools.combinations([v for v in range(k) if v != r], 2):
        if D[y, x] + D[x, r] == D[y, r] or D[x, y] + D[y, r] == D[x, r]:
            continue
        p += 1
    return k, I, p


def glue(H1, r1, H2, r2, L):
    G = nx.Graph()
    m1 = {v: ('a', v) for v in H1}; m2 = {v: ('b', v) for v in H2}
    G.add_edges_from((m1[u], m1[v]) for u, v in H1.edges())
    G.add_edges_from((m2[u], m2[v]) for u, v in H2.edges())
    G.add_nodes_from(m1.values()); G.add_nodes_from(m2.values())
    path = [m1[r1]] + [('p', i) for i in range(L - 1)] + [m2[r2]]
    nx.add_path(G, path)
    return nx.convert_node_labels_to_integers(G)
