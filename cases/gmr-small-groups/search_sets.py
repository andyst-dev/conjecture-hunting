"""Unstructured search: simulated annealing over arbitrary subsets A_1..A_n of a small group
(moves: toggle one element in one set, keep sets nonempty). Exact integer ratio via gmr.evaluate.
Usage: python3 search_sets.py GROUP n restarts steps"""
import sys, json, math, random
import numpy as np
from gmr import evaluate
from groups import build_all

name, n, restarts, steps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
G = build_all()[name]
rng = random.Random(n * 1000 + restarts)
best_all = None
for r in range(restarts):
    k0 = rng.choice([2, 3, 4, 6])
    sets = [set(rng.sample(range(G.order), k0)) for _ in range(n)]
    for s in sets: s.add(0) if rng.random() < 0.5 else None
    S, N, R = evaluate(G.mult, [sorted(s) for s in sets]); lr = math.log(R)
    best = (R, [sorted(s) for s in sets], S, N)
    for t in range(steps):
        T = 0.08 * (1 - t / steps) + 1e-6
        i = rng.randrange(n); g = rng.randrange(G.order)
        s = sets[i]
        if g in s:
            if len(s) == 1: continue
            s.remove(g); undo = lambda: s.add(g)
        else:
            if len(s) >= 12: continue
            s.add(g); undo = lambda: s.remove(g)
        S2, N2, R2 = evaluate(G.mult, [sorted(x) for x in sets]); l2 = math.log(R2)
        if l2 >= lr or rng.random() < math.exp((l2 - lr) / T):
            lr = l2
            if R2 > best[0]: best = (R2, [sorted(x) for x in sets], S2, N2)
        else:
            undo()
    if best_all is None or best[0] > best_all[0]: best_all = best
R, sets, S, N = best_all
print(f"{name} n={n}: best ratio {float(R):.5f}  |S|={S} N={N}" + ("   <<< VIOLATION" if R > 1 else ""))
with open("data/sets.jsonl", "a") as f:
    f.write(json.dumps({"group": name, "n": n, "ratio": str(R), "sets": sets, "S": S, "N": N}) + "\n")
