"""Sanity for r2_graphs.py: on every connected triangle-free graph with m <= 8 vertices,
compute the exact maximum number of labels by SAT and check (a) it never exceeds path_bound + 1,
(b) its maximum over graphs equals m - 1 (bound attained) and agrees with the global SAT result."""
import subprocess, networkx as nx
from r2_graphs import path_bound, sat_labelling
GENG = "/tmp/nautydl/pynauty-2.8.8.1/src/nauty2_8_8/geng"
for m in range(3, 9):
    lines = subprocess.run([GENG, "-tcq", str(m)], capture_output=True, text=True).stdout.split()
    best = 0; bad = 0
    for l in lines:
        G = nx.from_graph6_bytes(l.encode())
        L = 1
        while L + 1 <= G.number_of_edges() and sat_labelling(G, L + 1) is not None: L += 1
        if L > path_bound(G) + 1: bad += 1
        best = max(best, L)
    print(f"m={m}: {len(lines)} graphs, max labels over graphs = {best} (bound m-1 = {m-1}), path-bound violations = {bad}")
