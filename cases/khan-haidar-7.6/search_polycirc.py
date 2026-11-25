"""k-orbit polycirculant graphs on Z_N x {0..k-1}: (a,i) ~ (b,j) iff b-a in S_ij.
S_ii symmetric subset of Z_N\\{0}, S_ji = -S_ij. Fourier: L = (+)_chi L_chi with
L_chi[i][j] = delta_ij deg_i - sum_{s in S_ij} w^(chi s)  (w = e^{2 pi i/N})
diag L+   (orbit i) = (1/N) sum_chi (L_chi^+)_ii
diag L+^2 (orbit i) = (1/N) sum_chi ((L_chi^+)^2)_ii
Paper searched 2-orbit (bi-Cayley); here k=3 (tri-circulant) is the new territory.
Usage: python3 search_polycirc.py k N [mode=exh|rand] [samples]
"""
import itertools, sys, json
import numpy as np


GROUP = None  # (dims tuple)


def set_group(dims):
    """abelian group Z_d1 x ... ; element index = mixed radix."""
    global GROUP, ELEMS, NEG, CT
    GROUP = tuple(dims)
    ELEMS = list(itertools.product(*[range(d) for d in dims]))
    idx = {e: t for t, e in enumerate(ELEMS)}
    NEG = [idx[tuple((-x) % d for x, d in zip(e, dims))] for e in ELEMS]
    CT = np.array([[np.exp(2j * np.pi * sum(c * x / d for c, x, d in zip(chi, e, dims))) for e in ELEMS] for chi in ELEMS])
    return len(ELEMS)


def subsets(N):
    return [frozenset(c) for r in range(N + 1) for c in itertools.combinations(range(N), r)]


def sym_subsets(N):
    out = []
    for S in subsets(N):
        if 0 in S: continue
        if all(NEG[s] in S for s in S): out.append(S)
    return out


def charsum(S, N):
    return CT[:, sorted(S)].sum(1) if S else np.zeros(N, complex)


def evaluate(k, N, Sdiag, Soff):
    """Sdiag: list (len k) of arrays (B,N) of char sums for S_ii; Soff: dict (i,j)->(B,N) char sums for S_ij (i<j);
    degrees: arrays. Returns c (B,k), c2 (B,k), valid mask."""
    raise NotImplementedError


def batch_eval(k, N, cs_diag, deg_diag, cs_off, deg_off):
    """cs_diag[i]: (B,N) ; deg_diag[i]: (B,) ; cs_off[(i,j)]: (B,N); deg_off[(i,j)]: (B,)"""
    B = cs_diag[0].shape[0]
    deg = np.zeros((B, k))
    for i in range(k):
        deg[:, i] += deg_diag[i]
    for (i, j), d in deg_off.items():
        deg[:, i] += d; deg[:, j] += d
    L = np.zeros((B, N, k, k), complex)
    for i in range(k):
        L[:, :, i, i] = deg[:, i][:, None] - cs_diag[i]
    for (i, j), cs in cs_off.items():
        L[:, :, i, j] = -cs
        # S_ji = -S_ij -> char sum is conjugate
        L[:, :, j, i] = -np.conj(cs)
    J = np.full((k, k), 1.0 / k)
    L[:, 0] += J
    # singular nontrivial blocks => disconnected; detect via det
    det = np.linalg.det(L)
    valid = np.all(np.abs(det) > 1e-9, axis=1)
    L[~valid] = np.eye(k)
    P = np.linalg.inv(L)
    P[:, 0] -= J
    dP = np.real(np.diagonal(P, 0, 2, 3)).sum(1) / N
    P2 = P @ P
    dP2 = np.real(np.diagonal(P2, 0, 2, 3)).sum(1) / N
    return dP, dP2, deg, valid


def run(k, N, mode, samples, out, tol=1e-10, seed=0):
    rng = np.random.default_rng(seed)
    SS = sym_subsets(N)
    AS = subsets(N)
    csS = np.array([charsum(S, N) for S in SS]); dS = np.array([len(S) for S in SS])
    csA = np.array([charsum(S, N) for S in AS]); dA = np.array([len(S) for S in AS])
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    # normalisation: for pairs (0,j) we may translate orbit j so S_0j contains 0 (or is empty)
    norm_idx = np.array([t for t, S in enumerate(AS) if (not S) or 0 in S])
    all_idx = np.arange(len(AS))
    hits = viol = 0
    if mode == "exh":
        choices = [range(len(SS))] * k + [norm_idx if i == 0 else all_idx for (i, j) in pairs]
        total = np.prod([len(c) for c in choices])
        print("configs:", total, flush=True)
        it = itertools.product(*choices)
        def batches():
            while True:
                chunk = list(itertools.islice(it, 300000))
                if not chunk: return
                yield np.array(chunk)
    else:
        def batches():
            for _ in range(samples // 200000):
                cols = [rng.integers(0, len(SS), 200000) for _ in range(k)]
                cols += [rng.choice(norm_idx if i == 0 else all_idx, 200000) for (i, j) in pairs]
                yield np.stack(cols, 1)
    for X in batches():
        cs_diag = [csS[X[:, i]] for i in range(k)]
        deg_diag = [dS[X[:, i]] for i in range(k)]
        cs_off = {p: csA[X[:, k + t]] for t, p in enumerate(pairs)}
        deg_off = {p: dA[X[:, k + t]] for t, p in enumerate(pairs)}
        c, c2, deg, valid = batch_eval(k, N, cs_diag, deg_diag, cs_off, deg_off)
        m = valid & (c.max(1) - c.min(1) < tol * np.abs(c.mean(1))) & (deg.max(1) - deg.min(1) > .5)
        for b in np.nonzero(m)[0]:
            hits += 1
            d, q = deg[b], c2[b]
            bad = any(d[a] < d[e] - .5 and not q[a] < q[e] - 1e-12 * abs(q[e]) for a in range(k) for e in range(k))
            rec = dict(k=k, G=list(GROUP), Sdiag=[sorted(SS[X[b, i]]) for i in range(k)],
                       Soff={f"{i}{j}": sorted(AS[X[b, k + t]]) for t, (i, j) in enumerate(pairs)},
                       deg=d.tolist(), c=float(c[b, 0]), c2=q.tolist(), viol=bool(bad))
            out.write(json.dumps(rec) + "\n"); out.flush()
            if bad:
                viol += 1; print("VIOLATION?", rec, flush=True)
    return hits, viol


if __name__ == "__main__":
    k = int(sys.argv[1]); dims = [int(x) for x in sys.argv[2].split("x")]
    N = set_group(dims)
    mode = sys.argv[3] if len(sys.argv) > 3 else "exh"
    samples = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    with open(f"pc_k{k}_G{sys.argv[2]}_{mode}_{seed}.jsonl", "w") as f:
        h, v = run(k, N, mode, samples, f, seed=seed)
    print(f"polycirculant k={k} N={N} {mode}: hits={h} violations={v}")
