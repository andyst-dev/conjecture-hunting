"""Symmetric blow-ups: module sizes & types constant on Aut(H)-orbits.
#equations = #orbits - 1, #unknowns = #orbits -> families of solutions, allows large sizes."""
import itertools, sys, json
import numpy as np
import networkx as nx
from networkx.generators.atlas import graph_atlas_g
from networkx.algorithms.isomorphism import GraphMatcher
from search_gen import tables, quotient


def orbits(H):
    k = H.number_of_nodes()
    parent = list(range(k))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for m in GraphMatcher(H, H).isomorphisms_iter():
        for a, b in m.items():
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    reps = sorted(set(find(x) for x in range(k)))
    return [[x for x in range(k) if find(x) == r] for r in reps]


def orbit_partitions(H, maxcells):
    k = H.number_of_nodes()
    auts = [tuple(m[i] for i in range(k)) for m in GraphMatcher(H, H).isomorphisms_iter()]
    def orb(gens):
        parent = list(range(k))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for g in gens:
            for a in range(k):
                ra, rb = find(a), find(g[a])
                if ra != rb: parent[ra] = rb
        cells = {}
        for x in range(k): cells.setdefault(find(x), []).append(x)
        return tuple(sorted(tuple(c) for c in cells.values()))
    parts = set()
    gensets = [(a,) for a in auts]
    if len(auts) <= 120:
        gensets += list(itertools.combinations(auts, 2))
    for gs in gensets:
        p = orb(gs)
        if 2 <= len(p) <= maxcells: parts.add(p)
    return [list(map(list, p)) for p in parts]


def run(kmax, maxorb, S, types, out, tol=1e-10):
    emax = (kmax - 1) * S
    W, W2, R, VALID = tables(types, S, emax)
    hits = viol = 0
    for H in graph_atlas_g():
        k = H.number_of_nodes()
        if k < 2 or k > kmax or not nx.is_connected(H):
            continue
        A = nx.to_numpy_array(H, nodelist=range(k))
        for orbs in orbit_partitions(H, maxorb):
          no = len(orbs)
          if True:
            osz = np.array(list(itertools.product(range(1, S + 1), repeat=no)), dtype=np.int64)
            N = np.zeros((len(osz), k), dtype=np.int64)
            for oi, ob in enumerate(orbs):
                N[:, ob] = osz[:, [oi]]
            e, qd, qd2 = quotient(A, N)
            ei = e.round().astype(int)
            rep = [ob[0] for ob in orbs]
            for ts in itertools.product(types, repeat=no):
                ok = np.ones(len(osz), bool)
                c = qd[:, rep].copy(); c2 = qd2[:, rep].copy(); deg = e[:, rep].copy()
                for oi, t in enumerate(ts):
                    n = osz[:, oi]; i = rep[oi]
                    ok &= VALID[t][n] & ((n > 1) | (t == "I"))
                    c[:, oi] += W[t][n, ei[:, i]]
                    c2[:, oi] += W2[t][n, ei[:, i]]
                    deg[:, oi] += R[t][n]
                idx = np.nonzero(ok & (c.max(1) - c.min(1) < tol * np.abs(c.mean(1))) & (deg.max(1) - deg.min(1) > .5))[0]
                for j in idx:
                    hits += 1
                    d, q = deg[j], c2[j]
                    bad = any(d[a] < d[b] - .5 and not q[a] < q[b] - 1e-12 * abs(q[b]) for a in range(no) for b in range(no))
                    rec = dict(H=[list(x) for x in H.edges()], k=k, orbits=orbs, types=list(ts), sizes=osz[j].tolist(),
                               deg=d.tolist(), c=float(c[j, 0]), c2=q.tolist(), viol=bool(bad))
                    out.write(json.dumps(rec) + "\n"); out.flush()
                    if bad:
                        viol += 1
                        print("VIOLATION?", rec, flush=True)
    return hits, viol


if __name__ == "__main__":
    kmax, maxorb, S = map(int, sys.argv[1:4])
    types = sys.argv[4].split(",") if len(sys.argv) > 4 else ["K", "I"]
    tag = "-".join(types)
    with open(f"orb_k{kmax}_o{maxorb}_S{S}_{tag}.jsonl", "w") as f:
        h, v = run(kmax, maxorb, S, types, f)
    print(f"orbit search kmax={kmax} maxorb={maxorb} S={S} types={types}: hits={h} violations={v}")
