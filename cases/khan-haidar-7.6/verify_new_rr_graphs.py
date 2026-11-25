"""Standalone exact verifier (imports no search code).
Builds two C6 blow-ups with non-trivial modules and checks, in exact rationals:
  (1) resistance-regularity: diag(L+) constant,
  (2) non-regularity,
  (3) Conjecture 7.6: deg u < deg v  =>  S2(u) < S2(v),  S2(v) = sum_w r_vw^2.
Module types: K = clique, CP = cocktail party (K_n minus a perfect matching),
PM = perfect matching. Consecutive modules around C6 are completely joined.
"""
import itertools
from flint import fmpq_mat, fmpq


def module_edges(t, n):
    if t == "K":
        return list(itertools.combinations(range(n), 2))
    if t == "PM":
        return [(2 * i, 2 * i + 1) for i in range(n // 2)]
    if t == "CP":
        return [(a, b) for a, b in itertools.combinations(range(n), 2) if b != a + 1 or a % 2]
    raise ValueError(t)


def c6_blowup(types, sizes):
    off = [sum(sizes[:i]) for i in range(6)]
    n = sum(sizes)
    E = set()
    for i, (t, s) in enumerate(zip(types, sizes)):
        for a, b in module_edges(t, s):
            E.add((off[i] + a, off[i] + b))
    for i in range(6):
        j = (i + 1) % 6
        for a in range(sizes[i]):
            for b in range(sizes[j]):
                u, v = off[i] + a, off[j] + b
                E.add((min(u, v), max(u, v)))
    return n, sorted(E)


def check(types, sizes):
    n, E = c6_blowup(types, sizes)
    deg = [0] * n
    M = fmpq_mat(n, n)
    for u, v in E:
        deg[u] += 1; deg[v] += 1
        M[u, v] = -1; M[v, u] = -1
    for i in range(n):
        M[i, i] = deg[i]
    J = fmpq(1, n)
    for i in range(n):
        for j in range(n):
            M[i, j] += J
    X = M.inv()
    for i in range(n):
        for j in range(n):
            X[i, j] -= J
    diag = {X[i, i] for i in range(n)}
    S2 = []
    for i in range(n):
        S2.append(sum((X[i, i] + X[j, j] - 2 * X[i, j]) ** 2 for j in range(n)))
    by = {}
    for d, s in zip(deg, S2):
        by.setdefault(d, set()).add(s)
    ds = sorted(by)
    conj = all(max(by[a]) < min(by[b]) for a, b in zip(ds, ds[1:]))
    print(f"{types} {sizes}: n={n}, |E|={len(E)}, degrees={ds}")
    print(f"  resistance-regular: {len(diag) == 1}  (L+_vv = {next(iter(diag)) if len(diag) == 1 else diag})")
    print(f"  S2 per degree: { {d: [str(x) for x in by[d]] for d in ds} }")
    print(f"  Conjecture 7.6 holds: {conj}")
    assert len(diag) == 1 and len(ds) > 1 and conj


if __name__ == "__main__":
    check(["PM", "CP", "CP", "PM", "CP", "CP"], [4, 8, 8, 4, 8, 8])
    check(["CP", "K", "K", "CP", "K", "K"], [20, 36, 36, 20, 36, 36])
