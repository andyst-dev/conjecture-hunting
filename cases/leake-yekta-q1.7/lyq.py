"""Fiber counts N(b) = #{y in Z^n_{>=0} : My = b} and a screen for Question 1.7
of Leake--Mohammadi Yekta (arXiv:2609.39917).

Exact integer arithmetic for all counts; the LP screen uses floats, and any hit is
re-checked exactly (rational weights, integer power comparison)."""
import itertools, math
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
from sympy import Matrix

def is_unimodular(M):
    """Paper, Sec. 3: all non-singular rank(M) x rank(M) submatrices have det +-1."""
    A = Matrix(M); r = A.rank(); m, n = A.shape
    for rows in itertools.combinations(range(m), r):
        for cols in itertools.combinations(range(n), r):
            d = A.extract(list(rows), list(cols)).det()
            if d not in (0, 1, -1):
                return False
    return True

def pointing_functional(M, R=3):
    """Integer w with w.c_j >= 1 for every column c_j (exists iff all fibers are finite
    and there is no zero column).  Small search, then None."""
    m, n = len(M), len(M[0])
    cols = [tuple(M[i][j] for i in range(m)) for j in range(n)]
    best = None
    for w in itertools.product(range(-R, R + 1), repeat=m):
        if all(sum(a * b for a, b in zip(w, c)) >= 1 for c in cols):
            s = sum(sum(a * b for a, b in zip(w, c)) for c in cols)
            if best is None or s < best[0]:
                best = (s, w)
    return None if best is None else best[1]

def counts(M, w, W):
    """Exact N(b) for every b with w.b <= W (b reachable).  Requires w.c_j >= 1.
    Every solution y of My=b has sum_j y_j (w.c_j) = w.b <= W, so the DP sees all."""
    m, n = len(M), len(M[0])
    lev = lambda b: sum(a * x for a, x in zip(w, b))
    N = {(0,) * m: 1}
    for j in range(n):
        c = tuple(M[i][j] for i in range(m))
        new = dict(N)
        # process in increasing level; new[s+c] += new[s]
        order = sorted(new.keys(), key=lev)
        # states can be created on the fly; use a level-bucketed queue
        buckets = {}
        for s in order:
            buckets.setdefault(lev(s), []).append(s)
        L = 0
        while L <= W:
            for s in buckets.get(L, []):
                t = tuple(a + b for a, b in zip(s, c))
                lt = lev(t)
                if lt <= W:
                    if t not in new:
                        new[t] = 0
                        buckets.setdefault(lt, []).append(t)
                    new[t] += new[s]
            L += 1
        N = new
    return N

def counts_brute(M, w, W):
    """Independent enumeration of all y with sum_j y_j (w.c_j) <= W."""
    m, n = len(M), len(M[0])
    cols = [tuple(M[i][j] for i in range(m)) for j in range(n)]
    wc = [sum(a * b for a, b in zip(w, c)) for c in cols]
    out = {}
    def rec(j, budget, b):
        if j == n:
            out[b] = out.get(b, 0) + 1
            return
        t = 0
        while t * wc[j] <= budget:
            rec(j + 1, budget - t * wc[j], tuple(x + t * y for x, y in zip(b, cols[j])))
            t += 1
    rec(0, W, (0,) * m)
    return out

def exact_check(N, b, support):
    """support: list of lattice points a_i.  Solve sum mu_i a_i = b, sum mu_i = 1 exactly
    (least-squares-free: requires affinely independent support).  Returns
    (is_violation, mu) where violation means N(b) < prod N(a_i)^{mu_i}."""
    k = len(support); m = len(b)
    A = Matrix([[Fraction(a[r]) for a in support] for r in range(m)] + [[1] * k])
    rhs = Matrix([Fraction(x) for x in b] + [1])
    if A.rank() < k:
        return None, None
    sol = A.solve_least_squares(rhs) if A.shape[0] != k else A.solve(rhs)
    if A * sol != rhs:
        return None, None
    mu = [Fraction(int(x.p), int(x.q)) for x in sol]
    if any(x < 0 for x in mu):
        return None, mu
    D = 1
    for x in mu:
        D = D * x.denominator // math.gcd(D, x.denominator)
    lhs = N[tuple(b)] ** D
    rhs_v = 1
    for a, x in zip(support, mu):
        rhs_v *= N[tuple(a)] ** int(x * D)
    return lhs < rhs_v, mu

def screen(N, tol=1e-9):
    """For each b in the domain, maximise sum mu_a log N(a) over mu >= 0 supported on
    the domain with mean b.  Returns list of (gap, b, support) with gap > tol."""
    pts = [b for b, v in N.items() if v > 0]
    P = np.array(pts, dtype=float).T
    L = np.array([math.log(N[b]) for b in pts])
    A_eq = np.vstack([P, np.ones(len(pts))])
    hits = []; maxgap = -1.0
    for idx, b in enumerate(pts):
        res = linprog(-L, A_eq=A_eq, b_eq=np.append(np.array(b, float), 1.0),
                      bounds=(0, None), method='highs')
        if res.status != 0:
            continue
        gap = -res.fun - L[idx]
        maxgap = max(maxgap, gap)
        if gap > tol:
            sup = [pts[i] for i in range(len(pts)) if res.x[i] > 1e-9]
            hits.append((gap, b, sup))
    return hits, maxgap, len(pts)

def midpoint_violations(N):
    """Theorem 1.4 of the paper (exact): N(b)^2 >= N(b-d) N(b+d)."""
    pts = list(N.keys()); S = set(pts); bad = []
    for p, q in itertools.combinations(pts, 2):
        s = tuple(x + y for x, y in zip(p, q))
        if all(x % 2 == 0 for x in s):
            b = tuple(x // 2 for x in s)
            if b in S and N[b] ** 2 < N[p] * N[q]:
                bad.append((b, p, q))
    return bad
