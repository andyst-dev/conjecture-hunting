"""Like search_large_n.py, but samples dice directly (uniform multisets via
stars and bars), so it never lists all binom(2n-1, n) dice. Usage:
    python3 search_sampled.py n sample_size seed"""
import sys, json, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.csgraph as cg
from dice import beats_count, target
from find_cycle import cycle_of_length

n, size, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
t0 = time.time()
rng = np.random.default_rng(seed)
seen = set()
while len(seen) < size:
    pos = np.sort(rng.choice(2 * n - 1, size=n, replace=False))
    seen.add(tuple(int(p - k + 1) for k, p in enumerate(pos)))  # multiset of [n]
D = sorted(seen)
M = np.zeros((len(D), n), dtype=np.int64)
for k, t in enumerate(D):
    for v in t: M[k, v - 1] += 1
C = np.zeros_like(M); C[:, 1:] = np.cumsum(M, axis=1)[:, :-1]
T = target(n); Ct = C.T.copy(); src, dst = [], []
for s in range(0, len(D), 1000):
    i, j = np.nonzero(M[s:s + 1000] @ Ct == T)
    src.append((i + s).astype(np.int32)); dst.append(j.astype(np.int32))
src = np.concatenate(src); dst = np.concatenate(dst)
B = sp.csr_matrix((np.ones(len(src), dtype=np.int8), (src, dst)), shape=(len(D), len(D)))
ncomp, lab = cg.connected_components(B, directed=True, connection='strong')
sizes = np.bincount(lab); big = np.argmax(sizes); idx = np.nonzero(lab == big)[0]
Cg = B[idx][:, idx].tocsr()
rec = dict(n=n, sample=len(D), seed=seed, sub_E=int(B.nnz), largest_scc=int(sizes[big]))
for L in (n, n + 1):
    cyc = cycle_of_length(Cg, L, starts=list(range(0, len(idx), max(1, len(idx)//10)))) if len(idx) >= L else None
    if cyc:
        dd = [D[idx[c]] for c in cyc]
        assert len(set(dd)) == L
        assert all(beats_count(dd[i], dd[(i + 1) % L]) == T for i in range(L))
    rec[f'cycle_len_{L}'] = [list(x) for x in dd] if cyc else None
rec['secs'] = round(time.time() - t0, 1)
print(json.dumps(rec), flush=True)
with open('data/large_n.jsonl', 'a') as f: f.write(json.dumps(rec) + '\n')
