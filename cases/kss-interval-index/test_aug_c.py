"""Cross-check the C noncol/diameter (aug.c, eval mode) against the Python implementation."""
import subprocess, random
import networkx as nx
from interval import distances, adj_from_edges, noncol

def test_eval():
    codes = []
    for k in range(3, 9):
        codes += [l.strip() for l in open(f"graphs/conn{k}.g6")]
    rng = random.Random(3)
    while len(codes) < 12000 + 15000:
        n = rng.randint(9, 16)
        G = nx.gnp_random_graph(n, rng.uniform(0.12, 0.7), seed=rng.randint(0, 10**9))
        if nx.is_connected(G):
            codes.append(nx.to_graph6_bytes(G, header=False).decode().strip())
    res = subprocess.run(["./aug", "eval"], input="\n".join(codes) + "\n", capture_output=True, text=True).stdout.split("\n")
    for code, r in zip(codes, res):
        G = nx.from_graph6_bytes(code.encode()); n = G.number_of_nodes()
        D = distances(n, adj_from_edges(n, G.edges()))
        assert r.split() == [str(noncol(D)), str(int(D.max()))], (code, r)
    return len(codes)

if __name__ == "__main__":
    print("C eval agrees with Python on", test_eval(), "graphs")
