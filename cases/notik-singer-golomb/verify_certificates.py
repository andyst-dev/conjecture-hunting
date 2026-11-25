"""Independent re-check (Python Fractions only; imports no search code) for one prime power m.

For every Singer ruler (all multipliers t coprime to q, largest-gap cut, one per mirror pair):
  1. checks that D is a perfect difference set mod q = m^2+m+1 (the construction itself is not trusted);
  2. builds the chain description of the chamber (sorted differences, consecutive forms);
  3. gets a dual vector from HiGHS and checks EXACTLY, in Fractions:
     y >= 0, y supported on rows tight at the ruler, sum_r y_r a_r = (1,...,1), sum y = length,
     and rank(support rows) = k  =>  the ruler is the unique shortest point of its chamber.
A candidate y is accepted only if all of these exact checks pass.
Usage: python3 verify_certificates.py m [max_rulers]
"""
import sys, math
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog


def singer_from_file_or_build(m):
    from singer import singer_set   # construction only; property is re-checked below
    return singer_set(m)


def perfect(D, q):
    s = set()
    for a in D:
        for b in D:
            if a != b:
                s.add((a - b) % q)
    return len(s) == q - 1 and len(D) * (len(D) - 1) == q - 1


def rows_of(marks):
    n = len(marks); k = n - 1
    diffs = sorted((marks[j] - marks[i], i, j) for i in range(n) for j in range(i + 1, n))
    def f(i, j):
        v = [0] * k
        for t in range(i, j): v[t] = 1
        return v
    rows = [f(diffs[0][1], diffs[0][2])]; vals = [diffs[0][0]]
    for r in range(len(diffs) - 1):
        a = f(diffs[r + 1][1], diffs[r + 1][2]); b = f(diffs[r][1], diffs[r][2])
        rows.append([x - y for x, y in zip(a, b)]); vals.append(diffs[r + 1][0] - diffs[r][0])
    return k, rows, vals


def exact_rank(M):
    M = [[Fraction(x) for x in row] for row in M]
    r = 0; cols = len(M[0]) if M else 0
    for c in range(cols):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                fac = M[i][c] / M[r][c]
                M[i] = [x - fac * y for x, y in zip(M[i], M[r])]
        r += 1
    return r


def dual_vector(rows, vals, k, delta=1e-3):
    """Maximise sum_r min(y_r, delta) over {y >= 0 on tight rows, A^T y = 1}: every tight row that
    can carry positive weight in some dual optimum gets weight, so the support has maximal rank."""
    T = [i for i, v in enumerate(vals) if v == 1]
    nT = len(T)
    AT = np.array([rows[i] for i in T], dtype=float).T
    from scipy.sparse import hstack, eye, csr_matrix
    A_eq = hstack([csr_matrix(AT), csr_matrix((k, nT))])
    A_ub = hstack([-eye(nT), eye(nT)])          # z - y <= 0
    c = np.concatenate([np.zeros(nT), -np.ones(nT)])
    bounds = [(0, None)] * nT + [(0, delta)] * nT
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(nT), A_eq=A_eq, b_eq=np.ones(k), bounds=bounds, method="highs")
    if res.status != 0:
        return None, T
    return res.x[:nT], T


def rational_y(rows, T, y, k):
    """Round the float y on its support to an exact solution: solve least-norm system exactly on an
    independent subset, then average with the float support pattern is not exact -> we instead solve
    exactly for each of up to k basic solutions. Simple robust route: exact solve on support via
    Fraction Gaussian elimination of A_S^T y = 1 with free variables set to their rounded values."""
    S = [T[i] for i in range(len(T)) if y[i] > 1e-7]
    yv = {T[i]: Fraction(y[i]).limit_denominator(10**6) for i in range(len(T)) if y[i] > 1e-7}
    # pick k independent rows of S as pivots, fix the others at their rounded values, solve exactly
    cols = S
    M = [[Fraction(rows[c][t]) for c in cols] for t in range(k)]
    piv_cols, r = [], 0
    A = [row[:] for row in M]
    for c in range(len(cols)):
        p = next((i for i in range(r, k) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        for i in range(k):
            if i != r and A[i][c] != 0:
                fac = A[i][c] / A[r][c]
                A[i] = [x - fac * z for x, z in zip(A[i], A[r])]
        piv_cols.append(c); r += 1
        if r == k: break
    if r < k:
        return None
    free = [c for c in range(len(cols)) if c not in piv_cols]
    sol = {cols[c]: yv[cols[c]] for c in free}
    # solve M[:, piv] y_piv = 1 - M[:, free] y_free exactly
    rhs = [Fraction(1) - sum(M[t][c] * sol[cols[c]] for c in free) for t in range(k)]
    B = [[M[t][c] for c in piv_cols] + [rhs[t]] for t in range(k)]
    for c in range(k):
        p = next(i for i in range(c, k) if B[i][c] != 0)
        B[c], B[p] = B[p], B[c]
        for i in range(k):
            if i != c and B[i][c] != 0:
                fac = B[i][c] / B[c][c]
                B[i] = [x - fac * z for x, z in zip(B[i], B[c])]
    for c in range(k):
        sol[cols[piv_cols[c]]] = B[c][k] / B[c][c]
    return sol


def check_ruler(marks):
    k, rows, vals = rows_of(marks)
    y, T = dual_vector(rows, vals, k)
    if y is None:
        return False, "no dual"
    sol = rational_y(rows, T, y, k)
    if sol is None:
        return False, "support rank < k"
    if any(v < 0 for v in sol.values()):
        return False, "negative after rationalisation"
    if any(vals[r] != 1 for r in sol):
        return False, "non-tight row in support"
    for t in range(k):
        if sum(v * rows[r][t] for r, v in sol.items()) != 1:
            return False, "A^T y != 1"
    if sum(sol.values()) != marks[-1]:
        return False, "sum y != length"
    S = [r for r, v in sol.items() if v > 0]
    if exact_rank([rows[r] for r in S]) != k:
        return False, "positive support rank < k"
    return True, None


if __name__ == "__main__":
    m = int(sys.argv[1]); cap = int(sys.argv[2]) if len(sys.argv) > 2 else 10**9
    q = m * m + m + 1
    D = singer_from_file_or_build(m)
    assert perfect(D, q), "not a perfect difference set"
    seen, ok, bad = set(), 0, []
    for t in range(1, q):
        if math.gcd(t, q) != 1: continue
        S = sorted((t * d) % q for d in D)
        G, i = max(((S[(i + 1) % len(S)] - S[i]) % q, i) for i in range(len(S)))
        st = S[(i + 1) % len(S)]
        mk = tuple(sorted((s - st) % q for s in S))
        key = min(mk, tuple(sorted(mk[-1] - x for x in mk)))
        if key in seen: continue
        seen.add(key)
        good, why = check_ruler(list(mk))
        print(f"  t={t}: {'certified' if good else why}", flush=True)
        if good: ok += 1
        else: bad.append((t, why))
        if len(seen) >= cap: break
    print(f"m={m}: {ok} rulers certified exactly, {len(bad)} not certified by this checker {bad[:5]}")
