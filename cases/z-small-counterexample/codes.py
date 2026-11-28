"""Binary codes used as matroid sources. A code is given by a k x n generator matrix; the matroid is
the column matroid. Columns are returned as integers (bit j = row j)."""
import itertools


def rref_gf2(rows, n):
    rows = [r for r in rows]
    out = []; pivs = []
    for c in range(n):
        p = next((i for i, r in enumerate(rows) if r >> c & 1), None)
        if p is None: continue
        pr = rows.pop(p)
        rows = [r ^ pr if r >> c & 1 else r for r in rows]
        out = [r ^ pr if r >> c & 1 else r for r in out]
        out.append(pr); pivs.append(c)
    return out


def cyclic_code(n, g):
    """generator polynomial g (int, bit i = coeff of x^i); rows = x^i g(x)."""
    deg = g.bit_length() - 1
    return [g << i for i in range(n - deg)]


def golay23():
    return cyclic_code(23, 0b110001110101)   # x^11+x^10+x^6+x^5+x^4+x^2+1


def extend_parity(rows, n):
    return [r | ((bin(r).count("1") & 1) << n) for r in rows]


def golay24():
    return extend_parity(golay23(), 23)


def columns(rows, n):
    rows = rref_gf2(rows, n)
    k = len(rows)
    cols = []
    for e in range(n):
        c = 0
        for j, r in enumerate(rows):
            if r >> e & 1: c |= 1 << j
        cols.append(c)
    return cols, k


def codewords(rows):
    words = [0]
    for r in rows:
        words += [w ^ r for w in words]
    return words


def weight_distribution(rows, n):
    from collections import Counter
    return dict(sorted(Counter(bin(w).count("1") for w in codewords(rows)).items()))


def puncture(rows, n, positions):
    """delete coordinates in positions (column deletion of the generator matrix)."""
    keep = [e for e in range(n) if e not in positions]
    new = []
    for r in rows:
        v = 0
        for i, e in enumerate(keep):
            if r >> e & 1: v |= 1 << i
        new.append(v)
    return [r for r in rref_gf2(new, len(keep)) if r], len(keep)


def shorten(rows, n, positions):
    """codewords vanishing on positions, then delete those coordinates."""
    words = [w for w in codewords(rows) if all(not (w >> e & 1) for e in positions)]
    basis = rref_gf2([w for w in words if w], n)
    return puncture(basis, n, positions)


def reed_muller(r, m):
    n = 2 ** m
    pts = list(itertools.product([0, 1], repeat=m))
    rows = []
    for deg in range(r + 1):
        for S in itertools.combinations(range(m), deg):
            v = 0
            for i, p in enumerate(pts):
                if all(p[j] for j in S): v |= 1 << i
            rows.append(v)
    return rows, n
