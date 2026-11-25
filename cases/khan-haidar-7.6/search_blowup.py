"""Fast float search over module blow-ups H[W_1..W_k], W_i clique (K) or independent (I).

Quotient formulas (derived; validated against exact full-graph computation):
  e_i   = sum_{j ~ i in H} n_j
  lam_i = n_i + e_i (K) | e_i (I)      eigenvalue on within-module zero-sum vectors
  deg_i = n_i - 1 + e_i (K) | e_i (I)
  M     = diag(e) - A_H * sqrt(n_i n_j)  (symmetrised quotient Laplacian)
  diag L+   at module i = (1-1/n_i)/lam_i   + (M+)_ii / n_i
  diag L+^2 at module i = (1-1/n_i)/lam_i^2 + (M+^2)_ii / n_i
RR <=> diag L+ equal across modules.  Conjecture <=> diag L+^2 strictly increasing in degree.
"""
import itertools, sys, json
import numpy as np
import networkx as nx
from networkx.generators.atlas import graph_atlas_g


def eval_batch(A, types, N):
    """A: (k,k) 0/1, types: bool array (k,) True=clique, N: (B,k) sizes. Returns c, c2, deg."""
    B, k = N.shape
    Nf = N.astype(float)
    e = Nf @ A
    isK = types.astype(float)
    lam = e + isK * Nf
    deg = e + isK * (Nf - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where(N > 1, (1 - 1 / Nf) / lam, 0.0)
        w2 = np.where(N > 1, (1 - 1 / Nf) / lam**2, 0.0)
    s = np.sqrt(Nf)
    M = -A[None, :, :] * (s[:, :, None] * s[:, None, :])
    M[:, np.arange(k), np.arange(k)] += e
    u = s / np.sqrt(Nf.sum(1, keepdims=True))
    uu = u[:, :, None] * u[:, None, :]
    P = np.linalg.inv(M + uu) - uu
    P2 = P @ P
    dP = np.diagonal(P, axis1=1, axis2=2)
    dP2 = np.diagonal(P2, axis1=1, axis2=2)
    c = w + dP / Nf
    c2 = w2 + dP2 / Nf
    return c, c2, deg


def search(k, S, out, chunk=200000, tol=1e-10):
    graphs = [g for g in graph_atlas_g() if g.number_of_nodes() == k and nx.is_connected(g)]
    hits = viol = 0
    sizes = np.array(list(itertools.product(range(1, S + 1), repeat=k)), dtype=np.int64)
    for gi, H in enumerate(graphs):
        A = nx.to_numpy_array(H, nodelist=range(k))
        for types in itertools.product([0, 1], repeat=k):
            t = np.array(types, dtype=bool)
            for st in range(0, len(sizes), chunk):
                N = sizes[st:st + chunk]
                # skip ambiguity: size-1 modules canonical as I (types irrelevant) -> require t=0 where N==1
                keep = ~np.any((N == 1) & t[None, :], axis=1)
                N = N[keep]
                if len(N) == 0:
                    continue
                c, c2, deg = eval_batch(A, t, N)
                spread = c.max(1) - c.min(1)
                rr = spread < tol * c.mean(1)
                nonreg = deg.max(1) - deg.min(1) > 0.5
                idx = np.nonzero(rr & nonreg)[0]
                for i in idx:
                    hits += 1
                    d, q = deg[i], c2[i]
                    bad = False
                    for a in range(k):
                        for b in range(k):
                            if d[a] < d[b] - 0.5 and not (q[a] < q[b] - 1e-12 * abs(q[b])):
                                bad = True
                    rec = dict(k=k, H=list(map(list, H.edges())), types="".join("K" if x else "I" for x in types),
                               sizes=N[i].tolist(), deg=d.tolist(), c=float(c[i, 0]), c2=q.tolist(), viol=bad)
                    out.write(json.dumps(rec) + "\n")
                    if bad:
                        viol += 1
                        print("VIOLATION?", rec, flush=True)
    return hits, viol


if __name__ == "__main__":
    k, S = int(sys.argv[1]), int(sys.argv[2])
    with open(f"hits_k{k}_S{S}.jsonl", "w") as f:
        h, v = search(k, S, f)
    print(f"k={k} S={S}: RR non-regular hits={h}, float violations={v}")
