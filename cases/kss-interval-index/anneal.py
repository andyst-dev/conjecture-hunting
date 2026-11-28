"""Simulated-annealing screen for small noncol (heuristic; exact integer arithmetic throughout).
mode A: connected graphs with min degree >= 2 (Problem 2 / 17(b)); target noncol < 2(n-4), or = 2(n-4) and not D_n.
mode B: 2-connected graphs of diameter >= 4 (Problem 17(a)); target noncol < 2(n-4).
Usage: python3 anneal.py NMIN NMAX STEPS RESTARTS SEED"""
import sys, json, math, random, time
import numpy as np
import networkx as nx
from scipy.sparse.csgraph import shortest_path
from interval import noncol
from families import dumbbell

def dist(A):
    D = shortest_path(A, unweighted=True, directed=False)
    if np.isinf(D).any():
        return None
    return D.astype(np.int64)

def valid(A, mode):
    if A.sum(axis=1).min() < 2:
        return None
    D = dist(A)
    if D is None:
        return None
    if mode == "B":
        if D.max() < 4 or not nx.is_biconnected(nx.from_numpy_array(A)):
            return None
    return D

def start(n, rng, mode):
    while True:
        A = np.zeros((n, n), dtype=np.int8)
        perm = list(range(n)); rng.shuffle(perm)
        for i in range(n):
            a, b = perm[i], perm[(i + 1) % n]; A[a, b] = A[b, a] = 1
        for _ in range(rng.randint(0, n // 3)):
            a, b = rng.sample(range(n), 2); A[a, b] = A[b, a] = 1
        if mode == "A" and rng.random() < 0.5:
            # start near the dumbbell with random relabelling
            G = dumbbell(n); P = rng.sample(range(n), n)
            A = np.zeros((n, n), dtype=np.int8)
            for u, v in G.edges(): A[P[u], P[v]] = A[P[v], P[u]] = 1
        D = valid(A, mode)
        if D is not None:
            return A, D

def run(n, mode, steps, restarts, rng, log):
    target = 2 * (n - 4)
    best = None
    Dn = dumbbell(n)
    for r in range(restarts):
        A, D = start(n, rng, mode)
        e = noncol(D)
        for t in range(steps):
            T = 3.0 * (0.02 / 3.0) ** (t / steps)
            B = A.copy()
            k = 1 if rng.random() < 0.7 else 2
            for _ in range(k):
                a, b = rng.sample(range(n), 2); B[a, b] = B[b, a] = 1 - B[a, b]
            D2 = valid(B, mode)
            if D2 is None:
                continue
            e2 = noncol(D2)
            if e2 <= e or rng.random() < math.exp((e - e2) / T):
                A, D, e = B, D2, e2
                if best is None or e < best[0]:
                    best = (e, nx.to_graph6_bytes(nx.from_numpy_array(A), header=False).decode().strip())
                if e < target or (e == target and mode == "B") or (e == target and mode == "A" and not nx.is_isomorphic(nx.from_numpy_array(A), Dn)):
                    log.write(json.dumps({"n": n, "mode": mode, "noncol": int(e), "2(n-4)": target,
                        "g6": nx.to_graph6_bytes(nx.from_numpy_array(A), header=False).decode().strip()}) + "\n"); log.flush()
    return best

if __name__ == "__main__":
    nmin, nmax, steps, restarts, seed = map(int, sys.argv[1:6])
    rng = random.Random(seed)
    log = open(f"data/anneal_hits_seed{seed}.jsonl", "a")
    summ = open(f"data/anneal_summary_seed{seed}.jsonl", "a")
    for n in range(nmin, nmax + 1):
        for mode in ("A", "B"):
            t0 = time.time()
            b = run(n, mode, steps, restarts, rng, log)
            rec = {"n": n, "mode": mode, "steps": steps, "restarts": restarts, "seed": seed,
                   "best_noncol": int(b[0]) if b else None, "2(n-4)": 2 * (n - 4), "best_g6": b[1] if b else None,
                   "secs": round(time.time() - t0, 1)}
            summ.write(json.dumps(rec) + "\n"); summ.flush()
            print(json.dumps(rec), flush=True)
