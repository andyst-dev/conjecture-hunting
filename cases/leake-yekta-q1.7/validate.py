"""Reproduction / sanity values.  The paper (arXiv:2609.39917) gives no numerical
tables, so we reproduce (i) classical contingency-table counts and (ii) the paper's
Theorem 1.4 (midpoint log-concavity on lines), which must hold on every instance."""
import random
from lyq import *

def ct_matrix(p, q):
    """Transportation matrix without the last column-sum row (full row rank).
    Variables x_{ij} in row-major order; rows: r_1..r_p, c_1..c_{q-1}."""
    M = [[1 if a == i else 0 for a in range(p) for _ in range(q)] for i in range(p)]
    M += [[1 if b == j else 0 for _ in range(p) for b in range(q)] for j in range(q - 1)]
    w = (1,) * p + (0,) * (q - 1)
    return M, w

if __name__ == '__main__':
    rows = []
    # 2x2: N = min(r1,r2,c1,c2)+1
    M, w = ct_matrix(2, 2); N = counts(M, w, 14)
    ok = all(v == min(b[0], b[1], b[2], b[0] + b[1] - b[2]) + 1 for b, v in N.items())
    rows.append(("2x2 tables, all margins with r1+r2<=14", "min(margins)+1", "identical" if ok else "MISMATCH"))
    # 3x3 magic: (n+1)(n+2)(n^2+3n+4)/8
    M, w = ct_matrix(3, 3); N = counts(M, w, 15)
    ours = [N[(n, n, n, n, n)] for n in range(6)]
    paper = [(n + 1) * (n + 2) * (n * n + 3 * n + 4) // 8 for n in range(6)]
    rows.append(("3x3 tables, all margins n, n=0..5", str(paper), str(ours)))
    # 4x4 with all margins n (OEIS A001501): 1, 24, 282, 2008
    M, w = ct_matrix(4, 4); N = counts(M, w, 12)
    ours = [N[(n,) * 7] for n in range(4)]
    rows.append(("4x4 tables, all margins n, n=0..3 (A001501)", "[1, 24, 282, 2008]", str(ours)))

    # Theorem 1.4 on contingency tables and random unimodular configurations
    inst = 0; bad = 0
    for (p, q, W) in [(2, 3, 9), (3, 3, 7), (2, 4, 7)]:
        M, w = ct_matrix(p, q); N = counts(M, w, W); inst += 1
        bad += len(midpoint_violations(N))
    rng = random.Random(1)
    while inst < 60:
        m = rng.choice([2, 3, 4]); n = rng.randint(m, m + 4)
        M = [[rng.choice([-1, 0, 0, 1, 1]) for _ in range(n)] for _ in range(m)]
        if Matrix(M).rank() < m or not is_unimodular(M): continue
        w = pointing_functional(M, 2)
        if w is None: continue
        N = counts(M, w, 8); inst += 1
        bad += len(midpoint_violations(N))
    rows.append(("Theorem 1.4 (midpoint), %d instances, all pairs in domain" % inst, "no violation", "%d violations" % bad))

    print("| Object | Paper / known value | Our value |\n|---|---|---|")
    for r in rows: print("| %s | %s | %s |" % r)
