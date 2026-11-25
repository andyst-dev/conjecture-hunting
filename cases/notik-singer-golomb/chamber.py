"""Chamber optimality and uniqueness for an integer Golomb ruler (Notik, arXiv:2611.22656).

Gap coordinates g in R^k (k = n-1). Each difference x_j - x_i is the interval form 1_I . g.
The chamber P_sigma of a ruler is fixed by the sign of every relation e.g, where e is 1_I or 1_I - 1_J.

Key simplification (our observation): sort the N = n(n-1)/2 differences of the ruler,
d_(1) < ... < d_(N), with interval forms f_1..f_N. Then
    P_sigma = { g : f_1.g >= 1,  (f_{r+1} - f_r).g >= 1 for r = 1..N-1 }.
"⊆": these are among the defining rows, with the right signs.
"⊇": the chain gives f_s.g - f_r.g >= s - r >= 1 for every r < s, and f_r.g >= 1.
So the chamber has only N rows instead of |R_n| ~ n^4/8.

The ruler g* is the unique shortest point of P_sigma iff some y >= 0 on the rows tight at g*
satisfies sum y_r a_r = (1,...,1), and the support of y has rank k.
We find y with HiGHS, then certify it in exact rationals (python-flint) on a k x k basis.
If no y >= 0 exists, g* is not the chamber optimum; a strictly shorter chamber point is then
computed and certified exactly (see `better_point`).
"""
import numpy as np
from scipy.optimize import linprog
from flint import fmpq_mat, fmpq


class Span:
    """Incremental orthonormal basis (modified Gram-Schmidt) for fast 'is r outside span' tests."""
    def __init__(self, k):
        self.Q = np.zeros((0, k)); self.k = k
    def residual(self, v):
        v = v.astype(float)
        if len(self.Q):
            v = v - self.Q.T @ (self.Q @ v)
            v = v - self.Q.T @ (self.Q @ v)
        return v
    def outside(self, v, tol=1e-8):
        return np.linalg.norm(self.residual(v)) > tol * max(1.0, np.linalg.norm(v))
    def add(self, v, tol=1e-8):
        r = self.residual(v); nr = np.linalg.norm(r)
        if nr > tol * max(1.0, np.linalg.norm(v)):
            self.Q = np.vstack([self.Q, r / nr]); return True
        return False
    @property
    def dim(self):
        return self.Q.shape[0]


def chain_rows(marks):
    """marks sorted, start 0. Returns k, list of rows (dense int vectors) with b=1, values at g*, g*."""
    n = len(marks); k = n - 1
    g = np.diff(marks)
    diffs = []
    for i in range(n):
        for j in range(i + 1, n):
            diffs.append((marks[j] - marks[i], i, j))
    diffs.sort()
    vals = [d for d, _, _ in diffs]
    assert len(set(vals)) == len(vals), "not a Golomb ruler"

    def form(i, j):
        v = np.zeros(k, dtype=np.int64); v[i:j] = 1; return v
    rows = [form(diffs[0][1], diffs[0][2])]
    rhs_val = [vals[0]]
    for r in range(len(diffs) - 1):
        _, i1, j1 = diffs[r]; _, i2, j2 = diffs[r + 1]
        rows.append(form(i2, j2) - form(i1, j1))
        rhs_val.append(vals[r + 1] - vals[r])
    return k, np.array(rows), np.array(rhs_val), g


def certify(marks):
    """Returns ('unique', None) | ('shortest_not_unique', info) | ('not_shortest', info) | ('undecided', info)."""
    k, A, val, g = chain_rows(marks)
    T = np.nonzero(val == 1)[0]          # tight rows (b = 1)
    if len(T) == 0:
        return "not_shortest", better_point(A, k, marks)
    AT = A[T].T.astype(float)            # k x |T|
    res = linprog(np.zeros(len(T)), A_eq=AT, b_eq=np.ones(k), bounds=(0, None), method="highs")
    if res.status == 2:                  # infeasible: g* not optimal in its chamber
        return "not_shortest", better_point(A, k, marks)
    if res.status != 0:
        return "undecided", f"lp status {res.status}"
    y = res.x
    supp = T[y > 1e-9]
    cert = exact_basis_certificate(A, supp, k)
    if cert is not None:
        return "unique", cert
    # degenerate vertex: build an exact dual solution whose support has maximal rank
    return union_certificate(A, T, k)


def exact_basic_solution(A, supp, k):
    """Basic feasible dual solution from an LP vertex: its support is linearly independent.
    Re-solve exactly on that support (least squares not needed: independent rows, consistent)."""
    rows = []
    sp = Span(k)
    for r in supp:
        if sp.add(A[r]):
            rows.append(int(r))
    s = len(rows)
    # solve A_rows^T y = 1 exactly: k equations, s unknowns, full column rank -> use normal equations
    Bt = fmpq_mat(k, s)
    for a, r in enumerate(rows):
        for t in range(k):
            Bt[t, a] = int(A[r, t])
    N = Bt.transpose() * Bt
    rhs = Bt.transpose() * fmpq_mat(k, 1, [1] * k)
    y = N.solve(rhs)
    ys = [y[i, 0] for i in range(s)]
    # check exact consistency A^T y = 1 and y >= 0
    chk = Bt * y
    if any(chk[t, 0] != 1 for t in range(k)) or any(v < 0 for v in ys):
        return None
    return dict(zip(rows, ys))


def union_certificate(A, T, k):
    """Collect exact basic dual solutions maximising y_r for uncovered tight rows r, until the union
    of supports has rank k (=> unique shortest point) or no tight row outside the span can be
    made positive (=> optimal face is not a single point)."""
    AT = A[T].T.astype(float)
    sols = []
    sp = Span(k)
    candidates = list(range(len(T)))
    tried = set()
    while sp.dim < k:
        # pick a tight row outside the current span that has not been tried
        pick = None
        for c in candidates:
            if c in tried: continue
            if sp.outside(A[T[c]]):
                pick = c; break
        if pick is None:
            return "shortest_not_unique", f"dual support rank {sp.dim} < k={k}"
        tried.add(pick)
        obj = np.zeros(len(T)); obj[pick] = -1
        res = linprog(obj, A_eq=AT, b_eq=np.ones(k), bounds=(0, 1e7), method="highs")
        if res.status != 0 or res.x[pick] <= 1e-9:
            continue        # row pick can never be positive
        sol = exact_basic_solution(A, T[res.x > 1e-9], k)
        if sol is None:
            return "undecided", "exact basic solution failed"
        sols.append(sol)
        for r in sol:
            sp.add(A[r])
    # average of exact solutions: y >= 0, A^T y = 1, support = union, rank k
    ysum = {}
    for sol in sols:
        for r, v in sol.items():
            ysum[r] = ysum.get(r, fmpq(0)) + v / len(sols)
    # exact rank check of the support
    S = [r for r, v in ysum.items() if v > 0]
    Mq = fmpq_mat(len(S), k, [int(A[r, t]) for r in S for t in range(k)])
    if Mq.rank() != k:
        return "undecided", "exact rank < k"
    return "unique", {"support": len(S), "sum_y": str(sum(ysum.values(), fmpq(0)))}


def exact_basis_certificate(A, supp, k):
    """Pick k independent rows from supp (greedy, float rank), solve B^T y = 1 exactly, require y > 0.
    If y_B > 0 on a full-rank basis, g* is the unique shortest point (the support has rank k)."""
    rows = []
    sp = Span(k)
    for r in supp:
        if sp.add(A[r]):
            rows.append(r)
            if len(rows) == k:
                break
    if len(rows) < k:
        return None
    B = fmpq_mat(k, k)
    for a, r in enumerate(rows):
        for t in range(k):
            B[t, a] = int(A[r, t])       # B = A_rows^T
    ones = fmpq_mat(k, 1, [1] * k)
    try:
        y = B.solve(ones)
    except ZeroDivisionError:
        return None
    ys = [y[i, 0] for i in range(k)]
    if all(v > 0 for v in ys):
        return {"rows": [int(r) for r in rows], "sum_y": str(sum(ys, fmpq(0))), "min_y": str(min(ys))}
    return None


def better_point(A, k, marks):
    """Exact chamber optimum when g* is not optimal: LP vertex, then re-solve its tight basis exactly."""
    res = linprog(np.ones(k), A_ub=-A.astype(float), b_ub=-np.ones(len(A)), bounds=(None, None), method="highs")
    gx = res.x
    tight = np.nonzero(np.abs(A @ gx - 1) < 1e-7)[0]
    rows, M = [], np.zeros((0, k))
    for r in tight:
        cand = np.vstack([M, A[r]])
        if np.linalg.matrix_rank(cand) > M.shape[0]:
            M = cand; rows.append(r)
            if len(rows) == k: break
    B = fmpq_mat(k, k)
    for a, r in enumerate(rows):
        for t in range(k):
            B[a, t] = int(A[r, t])
    g = B.solve(fmpq_mat(k, 1, [1] * k))
    gq = [g[i, 0] for i in range(k)]
    # exact feasibility check of the better point
    for r in range(len(A)):
        if sum((int(A[r, t]) * gq[t] for t in range(k)), fmpq(0)) < 1:
            return {"error": "rationalised point infeasible"}
    L = sum(gq, fmpq(0))
    pts = [fmpq(0)]
    for v in gq: pts.append(pts[-1] + v)
    return {"length": str(L), "ruler_length": int(marks[-1]), "marks": [str(p) for p in pts]}
