"""Generalised blow-up search: module W_i is a vertex-transitive regular graph from a menu.
Only (n_i, r_i, nontrivial Laplacian spectrum of W_i) matter:
  deg_i = r_i + e_i
  diag L+   = (1/n_i) sum_mu 1/(mu+e_i)    + (M+)_ii/n_i
  diag L+^2 = (1/n_i) sum_mu 1/(mu+e_i)^2  + (M+^2)_ii/n_i
Quotient M identical to plain blow-ups (internal regular graph cancels on module-constant vectors).
"""
import itertools, sys, json
import numpy as np
import networkx as nx
from networkx.generators.atlas import graph_atlas_g


def spectrum(t, n):
    """nontrivial Laplacian eigenvalues (list of (mu, mult)) and degree r, or None if invalid."""
    if t == "K":
        return [(n, n - 1)], n - 1
    if t == "I":
        return [(0, n - 1)], 0
    if n < 4 or n % 2:
        if t in ("CP", "PM", "KB"):
            return None
    if t == "CP":   # cocktail party K_n minus perfect matching
        return [(n, n // 2 - 1), (n - 2, n // 2)], n - 2
    if t == "PM":   # perfect matching
        return [(0, n // 2 - 1), (2, n // 2)], 1
    if t == "KB":   # K_{n/2,n/2}
        a = n // 2
        return [(a, 2 * a - 2), (2 * a, 1)], a
    if t == "C":
        if n < 5:   # C3=K3, C4=KB(4)
            return None
        return [(2 - 2 * np.cos(2 * np.pi * j / n), 1) for j in range(1, n)], 2
    if t == "CC":   # complement of cycle
        if n < 5:
            return None
        return [(n - (2 - 2 * np.cos(2 * np.pi * j / n)), 1) for j in range(1, n)], n - 3
    raise ValueError(t)


def tables(types, S, emax):
    W, W2, R, VALID = {}, {}, {}, {}
    for t in types:
        W[t] = np.zeros((S + 1, emax + 1)); W2[t] = np.zeros((S + 1, emax + 1))
        R[t] = np.zeros(S + 1); VALID[t] = np.zeros(S + 1, bool)
        for n in range(1, S + 1):
            sp = spectrum(t, n) if n > 1 else ([], 0)
            if sp is None:
                continue
            mus, r = sp
            VALID[t][n] = True; R[t][n] = r
            e = np.arange(emax + 1, dtype=float)
            with np.errstate(divide="ignore"):
                W[t][n] = sum(m / (mu + e) for mu, m in mus) / n if mus else 0
                W2[t][n] = sum(m / (mu + e) ** 2 for mu, m in mus) / n if mus else 0
    return W, W2, R, VALID


def quotient(A, N):
    B, k = N.shape
    Nf = N.astype(float)
    e = Nf @ A
    s = np.sqrt(Nf)
    M = -A[None] * (s[:, :, None] * s[:, None, :])
    M[:, np.arange(k), np.arange(k)] += e
    u = s / np.sqrt(Nf.sum(1, keepdims=True))
    uu = u[:, :, None] * u[:, None, :]
    P = np.linalg.inv(M + uu) - uu
    return e, np.diagonal(P, 0, 1, 2) / Nf, np.diagonal(P @ P, 0, 1, 2) / Nf


def search(k, S, types, out, tol=1e-10, need_nontrivial=True):
    emax = (k - 1) * S
    W, W2, R, VALID = tables(types, S, emax)
    graphs = [g for g in graph_atlas_g() if g.number_of_nodes() == k and nx.is_connected(g)]
    sizes = np.array(list(itertools.product(range(2, S + 1), repeat=k)), dtype=np.int64)
    hits = viol = 0
    for H in graphs:
        A = nx.to_numpy_array(H, nodelist=range(k))
        e, qd, qd2 = quotient(A, sizes)          # quotient part independent of module types
        ei = e.round().astype(int)
        for ts in itertools.product(types, repeat=k):
            if need_nontrivial and all(t in ("K", "I") for t in ts):
                continue
            ok = np.ones(len(sizes), bool)
            c = qd.copy(); c2 = qd2.copy(); deg = e.copy()
            for i, t in enumerate(ts):
                n = sizes[:, i]
                ok &= VALID[t][n]
                c[:, i] += W[t][n, ei[:, i]]
                c2[:, i] += W2[t][n, ei[:, i]]
                deg[:, i] += R[t][n]
            spread = c.max(1) - c.min(1)
            idx = np.nonzero(ok & (spread < tol * np.abs(c.mean(1))) & (deg.max(1) - deg.min(1) > 0.5))[0]
            for j in idx:
                hits += 1
                d, q = deg[j], c2[j]
                bad = any(d[a] < d[b] - .5 and not q[a] < q[b] - 1e-12 * abs(q[b]) for a in range(k) for b in range(k))
                rec = dict(k=k, H=[list(x) for x in H.edges()], types=list(ts), sizes=sizes[j].tolist(),
                           deg=d.tolist(), c=float(c[j, 0]), c2=q.tolist(), viol=bool(bad))
                out.write(json.dumps(rec) + "\n")
                if bad:
                    viol += 1
                    print("VIOLATION?", rec, flush=True)
    return hits, viol


if __name__ == "__main__":
    k, S = int(sys.argv[1]), int(sys.argv[2])
    types = sys.argv[3].split(",") if len(sys.argv) > 3 else ["K", "I", "CP", "PM", "KB", "C", "CC"]
    with open(f"gen_k{k}_S{S}.jsonl", "w") as f:
        h, v = search(k, S, types, f)
    print(f"gen k={k} S={S} types={types}: hits={h} violations={v}")
