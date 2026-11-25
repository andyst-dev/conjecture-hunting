"""Hill-climbing over sequences of small 'word' sets built from two order-3 elements a, b
(the shape of the paper's order-243 witness), in small groups. Exact integer ratios."""
import sys, json, random, itertools
import numpy as np
from fractions import Fraction
from gmr import evaluate, cyc
from groups import build_all, order3_generating_pairs


def alphabet(G, a, b):
    m = G.mult
    a2, b2 = int(m[a, a]), int(m[b, b])
    T = {"e": 0, "a": a, "A": a2, "b": b, "B": b2}
    L = {}
    for r in (1, 2, 3):
        for combo in itertools.combinations(T, r):
            idx = sorted({T[c] for c in combo})
            L["{" + ",".join(combo) + "}"] = idx
    L["<a>"] = sorted(cyc(G, a)); L["<b>"] = sorted(cyc(G, b))
    # dedupe identical sets
    seen, out = set(), {}
    for k, v in L.items():
        if tuple(v) not in seen:
            seen.add(tuple(v)); out[k] = v
    return out


def climb(G, L, n, rng, steps, start=None, T0=0.05):
    """simulated annealing on log ratio; singletons are excluded from the alphabet to avoid the
    trivial plateau (ratio 1 for all-singleton sequences)"""
    import math
    keys = [k for k in L if len(L[k]) >= 2]
    seq = start[:] if start else [rng.choice(keys) for _ in range(n)]
    S, N, R = evaluate(G.mult, [L[k] for k in seq])
    lr = math.log(R)
    best = (R, seq[:], S, N)
    for t in range(steps):
        T = T0 * (1 - t / steps) + 1e-6
        i = rng.randrange(n); old = seq[i]; seq[i] = rng.choice(keys)
        S2, N2, R2 = evaluate(G.mult, [L[k] for k in seq])
        l2 = math.log(R2)
        if l2 >= lr or rng.random() < math.exp((l2 - lr) / T):
            lr = l2
            if R2 > best[0]:
                best = (R2, seq[:], S2, N2)
        else:
            seq[i] = old
    return best


if __name__ == "__main__":
    names = None if len(sys.argv) < 2 or sys.argv[1] == "all" else sys.argv[1].split(";")
    nmin, nmax = (int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "4-12").split("-"))
    restarts = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    steps = int(sys.argv[4]) if len(sys.argv) > 4 else 200
    rng = random.Random(0)
    groups = build_all()
    out = open("data/words.jsonl", "a")
    for name, G in groups.items():
        if names and name not in names: continue
        from groups import order_of, gen_closure
        pairs = []
        if order_of(G, 1) == 3 and order_of(G, 2) == 3 and len(gen_closure(G, [1, 2])) == G.order:
            pairs.append((1, 2))
        pairs += [p for p in order3_generating_pairs(G, limit=2) if p not in pairs]
        if not pairs:
            print(name, "no generating pair of order-3 elements"); continue
        for (a, b) in pairs:
            L = alphabet(G, a, b)
            for n in range(nmin, nmax + 1):
                best = None
                for r in range(restarts):
                    res = climb(G, L, n, rng, steps)
                    if best is None or res[0] > best[0]: best = res
                R, seq, S, N = best
                rec = {"group": name, "order": G.order, "a": a, "b": b, "n": n, "ratio": str(R), "ratio_f": float(R), "seq": seq, "S": S, "N": N}
                out.write(json.dumps(rec) + "\n"); out.flush()
                flag = "   <<< VIOLATION" if R > 1 else ""
                print(f"{name:22s} |G|={G.order:3d} pair=({a},{b}) n={n:2d} best ratio={float(R):.5f}{flag}", flush=True)
