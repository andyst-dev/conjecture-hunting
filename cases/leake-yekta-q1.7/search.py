"""Screen Question 1.7 on families of unimodular configurations.
usage: python3 search.py <class> <seed> <seconds>
classes: ct (contingency tables), rand (random {-1,0,1} unimodular), a3 (signed
sub-multisets of the A3 root configuration, rank 3)."""
import sys, json, time, random
from lyq import *
from validate import ct_matrix  # noqa (validate prints its table on import)

A3 = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,1,1),(1,1,1)]

def run_instance(M, w, W, cls, cov, hitf, maxdom=1500):
    N = counts(M, w, W)
    if len(N) > maxdom:
        return None
    hits, maxgap, npts = screen(N, 1e-9)
    rec = {"class": cls, "M": M, "w": list(w), "W": W, "domain": npts, "maxgap": maxgap,
           "float_hits": len(hits)}
    exact = 0
    for gap, b, sup in hits:
        viol, mu = exact_check(N, b, sup)
        h = {"class": cls, "M": M, "w": list(w), "W": W, "b": list(b), "gap": gap,
             "support": [list(a) for a in sup],
             "mu": None if mu is None else [str(x) for x in mu], "exact_violation": viol}
        hitf.write(json.dumps(h) + "\n"); hitf.flush()
        exact += bool(viol)
    rec["exact_hits"] = exact
    cov.write(json.dumps(rec) + "\n"); cov.flush()
    return rec

def main():
    cls, seed, T = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    rng = random.Random(seed); t0 = time.time(); done = 0; seen = set()
    cov = open("data/coverage_%s_%d.jsonl" % (cls, seed), "a")
    hitf = open("data/hits_%s_%d.jsonl" % (cls, seed), "a")
    if cls == "ct":
        for p, q, W in [(2,2,20),(2,3,9),(2,4,6),(3,3,5),(2,5,5),(3,4,4)]:
            M, w = ct_matrix(p, q)
            r = run_instance(M, w, W, "ct%dx%d" % (p, q), cov, hitf, maxdom=4000)
            print(p, q, W, r and {k: r[k] for k in ("domain", "maxgap", "exact_hits")}, flush=True)
        return
    while time.time() - t0 < T:
        if cls == "rand":
            m = rng.choice([2, 3, 3, 4, 4, 5]); n = rng.randint(m, m + 5)
            M = [[rng.choice([-1, 0, 0, 1, 1]) for _ in range(n)] for _ in range(m)]
        elif cls == "a3":
            cols = []
            for v in A3:
                k = rng.choice([0, 1, 1, 2, -1, -1, -2])
                s = 1 if k > 0 else -1
                cols += [tuple(s * x for x in v)] * abs(k)
            if not cols: continue
            M = [[c[i] for c in cols] for i in range(3)]; m = 3
        key = str(sorted(zip(*M)))
        if key in seen: continue
        seen.add(key)
        if Matrix(M).rank() < m or not is_unimodular(M): continue
        w = pointing_functional(M, 2)
        if w is None: continue
        W = rng.randint(5, 12)
        r = run_instance(M, w, W, cls, cov, hitf)
        if r: done += 1
        if r and r["exact_hits"]:
            print("EXACT HIT", r, flush=True)
    print(cls, seed, "instances screened:", done, flush=True)

if __name__ == "__main__":
    main()
