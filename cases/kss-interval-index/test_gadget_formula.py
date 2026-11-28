"""Cross-check of the gadget formula against the exact full computation on random instances."""
import random
import networkx as nx
from interval import distances, adj_from_edges, noncol
from gadget import gadget_params, glue

def random_gadget(rng):
    while True:
        k = rng.randint(2, 8)
        H = nx.gnp_random_graph(k, rng.uniform(0.2, 0.9), seed=rng.randint(0, 10**9))
        if nx.is_connected(H):
            return H, rng.randrange(k)

def test_formula(trials=3000, seed=7):
    rng = random.Random(seed)
    for _ in range(trials):
        H1, r1 = random_gadget(rng); H2, r2 = random_gadget(rng)
        L = rng.randint(1, 8)
        G = glue(H1, r1, H2, r2, L)
        n = G.number_of_nodes()
        k1, I1, p1 = gadget_params(H1, r1); k2, I2, p2 = gadget_params(H2, r2)
        pred = I1 + I2 + p1 * (n - k1) + p2 * (n - k2)
        actual = noncol(distances(n, adj_from_edges(n, G.edges())))
        assert pred == actual, (pred, actual, list(H1.edges()), r1, list(H2.edges()), r2, L)
    return trials

if __name__ == "__main__":
    print("formula agrees on", test_formula(), "random instances")
