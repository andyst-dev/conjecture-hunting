"""Cross-check the chain reduction against the full chamber of the paper.
Full chamber: all relation vectors e in R_n (1_I and 1_I - 1_J, normalised, deduplicated),
with signs taken from the ruler. Chain chamber: f_1 and the consecutive differences
f_{r+1} - f_r of the sorted difference list. We compare LP optima (exact rational vertex via
basis re-solve) on random integer Golomb rulers and on Singer rulers for small m.
"""
import itertools, random
import numpy as np
from scipy.optimize import linprog
from chamber import chain_rows


def full_rows(marks):
    n = len(marks); k = n - 1
    ints = [(i, j) for i in range(n) for j in range(i + 1, n)]
    def f(i, j):
        v = np.zeros(k, dtype=np.int64); v[i:j] = 1; return v
    forms = [f(i, j) for i, j in ints]
    vecs = {}
    for a in range(len(forms)):
        vecs[tuple(forms[a])] = forms[a]
        for b in range(a + 1, len(forms)):
            e = forms[a] - forms[b]
            nz = np.nonzero(e)[0]
            if len(nz) == 0: continue
            if e[nz[0]] < 0: e = -e
            vecs[tuple(e)] = e
    g = np.diff(marks)
    rows = []
    for e in vecs.values():
        s = int(np.sign(e @ g)); assert s != 0
        rows.append(s * e)
    return np.array(rows)


def lp_opt(A, k):
    r = linprog(np.ones(k), A_ub=-A.astype(float), b_ub=-np.ones(len(A)), bounds=(None, None), method="highs")
    assert r.status == 0
    return r.fun, r.x


def random_golomb(n, L, rng):
    while True:
        marks = sorted(rng.sample(range(1, L), n - 2)) + [L]
        marks = [0] + marks
        d = [marks[j] - marks[i] for i in range(n) for j in range(i + 1, n)]
        if len(set(d)) == len(d):
            return marks


if __name__ == "__main__":
    rng = random.Random(7)
    worst = 0.0; cnt = 0
    for _ in range(150):
        n = rng.randint(4, 8)
        marks = random_golomb(n, rng.randint(n * n // 2 + 5, 3 * n * n), rng)
        k = n - 1
        _, Ac, _, _ = chain_rows(marks)
        Af = full_rows(marks)
        fc, xc = lp_opt(Ac, k); ff, xf = lp_opt(Af, k)
        # same optimum value, and each optimal point satisfies the other system
        worst = max(worst, abs(fc - ff), max(0, 1 - (Af @ xc).min()), max(0, 1 - (Ac @ xf).min()))
        cnt += 1
    print(f"{cnt} random rulers (n=4..8): max discrepancy chain vs full chamber = {worst:.2e}")
    from run import rulers_for
    for m in (3, 4, 5):
        for t, rank, marks in rulers_for(m, allcuts=True):
            k = len(marks) - 1
            _, Ac, _, _ = chain_rows(list(marks)); Af = full_rows(list(marks))
            fc, xc = lp_opt(Ac, k); ff, xf = lp_opt(Af, k)
            worst = max(worst, abs(fc - ff), max(0, 1 - (Af @ xc).min()))
    print(f"incl. all cuts of Singer rulers m=3,4,5: max discrepancy = {worst:.2e}")
