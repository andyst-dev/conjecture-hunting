"""Benedetti–Menara Conjecture 7.1 (arXiv:2611.08071): a flag 1-layered family of r-subsets of an
m-set has at most m - r + 1 layers.

1-layered: the family is partitioned into nonempty layers L_0..L_N such that for every nonempty
set S, the indices t with some member of L_t containing S form an interval.
Flag: in the graph where u~v iff u,v lie in a common member, no member is contained in an
(r+1)-clique.

SAT(m, r, L): does a flag 1-layered family of r-subsets of [m] with exactly L layers exist?
Merging two adjacent layers keeps the interval property, so L-layer feasibility implies
(L-1)-layer feasibility. Hence UNSAT at L = m - r + 2 proves the conjecture for (m, r).
Usage: python3 sat.py r m [--noflag]
"""
import sys, itertools, time, json
from pysat.formula import IDPool, CNF
from pysat.solvers import Cadical153


def build(m, r, L, flag=True):
    pool = IDPool(); cnf = CNF()
    V = range(m)
    T = list(itertools.combinations(V, r))
    x = {(Tt, t): pool.id(("x", Tt, t)) for Tt in T for t in range(L)}
    mem = {Tt: pool.id(("mem", Tt)) for Tt in T}
    for Tt in T:
        lits = [x[Tt, t] for t in range(L)]
        # at most one layer
        for a, b in itertools.combinations(lits, 2):
            cnf.append([-a, -b])
        # mem <-> OR layers
        cnf.append([-mem[Tt]] + lits)
        for l in lits: cnf.append([-l, mem[Tt]])
    # every layer nonempty
    for t in range(L):
        cnf.append([x[Tt, t] for Tt in T])
    # interval property for every S with 1 <= |S| <= r-1 (|S| = r is a single member: trivial)
    for k in range(1, r):
        for S in itertools.combinations(V, k):
            sup = [Tt for Tt in T if set(S) <= set(Tt)]
            y = [pool.id(("y", S, t)) for t in range(L)]
            for t in range(L):
                lits = [x[Tt, t] for Tt in sup]
                cnf.append([-y[t]] + lits)
                for l in lits: cnf.append([-l, y[t]])
            # prefix p[t] = OR_{t'<=t} y ; suffix q[t] = OR_{t'>=t} y ; p[t] & q[t] -> y[t]
            p = [pool.id(("p", S, t)) for t in range(L)]
            q = [pool.id(("q", S, t)) for t in range(L)]
            for t in range(L):
                # p[t] <-> y[t] or p[t-1]
                prev = [p[t - 1]] if t > 0 else []
                cnf.append([-p[t], y[t]] + prev)
                cnf.append([-y[t], p[t]])
                if t > 0: cnf.append([-p[t - 1], p[t]])
                nxt = [q[t + 1]] if t < L - 1 else []
                cnf.append([-q[t], y[t]] + nxt)
                cnf.append([-y[t], q[t]])
                if t < L - 1: cnf.append([-q[t + 1], q[t]])
                cnf.append([-p[t], -q[t], y[t]])
    if flag:
        # edge variables
        E = {}
        for u, v in itertools.combinations(V, 2):
            e = pool.id(("e", u, v)); E[u, v] = E[v, u] = e
            sup = [Tt for Tt in T if u in Tt and v in Tt]
            cnf.append([-e] + [mem[Tt] for Tt in sup])
            for Tt in sup: cnf.append([-mem[Tt], e])
        # no member T together with an outside vertex w forms an (r+1)-clique
        for Tt in T:
            for w in V:
                if w in Tt: continue
                cnf.append([-mem[Tt]] + [-E[w, u] for u in Tt])
    # symmetry breaking (light): some member in layer 0 contains vertex 0
    cnf.append([x[Tt, 0] for Tt in T if 0 in Tt])
    return cnf, pool, T


def solve(m, r, L, flag=True):
    cnf, pool, T = build(m, r, L, flag)
    with Cadical153(bootstrap_with=cnf.clauses) as s:
        ok = s.solve()
        model = set(l for l in s.get_model() if l > 0) if ok else None
    if not ok:
        return None
    layers = {}
    for Tt in T:
        for t in range(L):
            if pool.id(("x", Tt, t)) in model:
                layers.setdefault(t, []).append(list(Tt))
    return [layers[t] for t in range(L)]


if __name__ == "__main__":
    r, m = int(sys.argv[1]), int(sys.argv[2])
    flag = "--noflag" not in sys.argv
    L = m - r + 2
    if "--L" in sys.argv: L = int(sys.argv[sys.argv.index("--L") + 1])
    t0 = time.time()
    sol = solve(m, r, L, flag)
    dt = time.time() - t0
    res = {"r": r, "m": m, "layers_asked": L, "flag": flag, "sat": sol is not None, "secs": round(dt, 2), "family": sol}
    print(json.dumps({k: v for k, v in res.items() if k != "family"}) + ("  <<< COUNTEREXAMPLE" if sol and flag and L > m - r + 1 else ""))
    if sol:
        print(json.dumps(sol))
    with open("data/sat_runs.jsonl", "a") as f:
        f.write(json.dumps(res) + "\n")
