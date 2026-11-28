"""Deduplicate aug.c hits up to isomorphism (pynauty certificates); histogram by noncol;
flag 2-connected graphs of diameter >= 4."""
import sys, json
from collections import Counter
import networkx as nx, pynauty
def cert(G):
    n = G.number_of_nodes()
    return pynauty.certificate(pynauty.Graph(n, adjacency_dict={v: list(G[v]) for v in G}))
def load(path):
    seen = {}
    for line in open(path):
        code, nc, dm = line.split()
        G = nx.from_graph6_bytes(code.encode())
        c = cert(G)
        if c not in seen:
            seen[c] = (code, int(nc), int(dm), nx.is_biconnected(G))
    return seen
if __name__ == "__main__":
    seen = load(sys.argv[1])
    h = Counter(v[1] for v in seen.values())
    h2 = Counter(v[1] for v in seen.values() if v[3] and v[2] >= 4)
    print("all delta>=2:", dict(sorted(h.items())))
    print("2-conn diam>=4:", dict(sorted(h2.items())))
    lo = min(h)
    print("minimisers:", [v for v in seen.values() if v[1] == lo])
    lo2 = min(h2) if h2 else None
    print("2c-diam4 minimisers:", [v for v in seen.values() if v[3] and v[2] >= 4 and v[1] == lo2])
