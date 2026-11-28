"""Reproduce the values stated in arXiv:2610.39966 before searching."""
from math import comb
import networkx as nx
from interval import distances, adj_from_edges, interval_index_direct, noncol, noncol_bruteforce


def D_of(G):
    G = nx.convert_node_labels_to_integers(G)
    n = G.number_of_nodes()
    return distances(n, adj_from_edges(n, G.edges()))


from families import dumbbell, two_triangles


rows = []
def chk(name, paper, ours):
    ok = paper == ours
    rows.append((name, paper, ours, ok))
    print(f"{name:40s} paper={paper!s:>8} ours={ours!s:>8} {'OK' if ok else 'MISMATCH'}")

for n in range(1, 13):
    D = D_of(nx.path_graph(n))
    chk(f"Int(P_{n}) = (n^3-n)/6", (n**3 - n)//6, interval_index_direct(D))
chk("Int(P_4)", 10, interval_index_direct(D_of(nx.path_graph(4))))
chk("Int(C_4)", 10, interval_index_direct(D_of(nx.cycle_graph(4))))
for a, b, v in [(2, 3, 1), (3, 3, 2), (3, 4, 5), (4, 4, 8)]:
    chk(f"noncol(K_{a},{b})", v, noncol(D_of(nx.complete_bipartite_graph(a, b))))
for a in range(1, 7):
    for b in range(a, 7):
        chk(f"noncol(K_{a},{b}) = C(a,3)+C(b,3)", comb(a, 3) + comb(b, 3), noncol(D_of(nx.complete_bipartite_graph(a, b))))
chk("noncol(C_8)", 8, noncol(D_of(nx.cycle_graph(8))))
for n in range(8, 41):
    chk(f"noncol(D_{n}) = 2(n-4)", 2 * (n - 4), noncol(D_of(dumbbell(n))))
chk("W(C_20)", 1000, nx.wiener_index(nx.cycle_graph(20)))
chk("W(two triangles + path, n=20)", 1294, int(nx.wiener_index(two_triangles(20))))
# Lemma 4: Int = C(n,2) + col = C(n,2)+C(n,3)-noncol, on random connected graphs; also noncol vs brute force
import random
random.seed(1)
bad = 0
for t in range(3000):
    n = random.randint(3, 14)
    p = random.uniform(0.15, 0.9)
    G = nx.gnp_random_graph(n, p, seed=random.randint(0, 10**9))
    if not nx.is_connected(G):
        continue
    D = D_of(G)
    nc = noncol(D)
    if nc != noncol_bruteforce(D) or interval_index_direct(D) != comb(n, 2) + comb(n, 3) - nc:
        bad += 1
chk("Lemma 4 + noncol brute-force, random graphs (#bad)", 0, bad)
print("ALL OK" if all(r[3] for r in rows) else "SOME MISMATCH")
