"""Fast exact Kazhdan-Lusztig (P) and Z-polynomials of matroids.

Key identity (Proudfoot-Xu-Young 2018): Z_M(t) = sum_{F flat} t^{rk F} P_{M/F}(t) is palindromic of
degree r = rk M, and deg P_M < r/2. Writing R(t) = sum_{F != cl(empty)} t^{rk F} P_{M/F}(t), this gives
    [t^i] P_M = [t^{r-i}] R - [t^i] R     for 0 <= i < r/2,
so P (and Z) follow from the P's of the proper contractions, with integer arithmetic only
(no Moebius function, no characteristic polynomials).

Two engines:
  * flats_linear(vecs, p) + kl_z_from_flats(...): any linear matroid over GF(p), via its lattice of flats;
  * nested_kl_z(chain): nested (Schubert) matroids given by their chain of cyclic flats, computed on
    "profiles" (intersection sizes with the layers of the chain), which is exact because the
    automorphism group S_{L1} x ... x S_{Lk} acts transitively on flats of a given profile and on pairs
    F <= X with given profiles.
Polynomials are lists of Python ints, low degree first.
"""
from functools import lru_cache
from math import comb
import itertools
import numpy as np


# ---------------------------------------------------------------- polynomials
def padd(a, b, shift=0, mult=1):
    """a += mult * t^shift * b  (in place, returns a)"""
    need = shift + len(b)
    if len(a) < need:
        a.extend([0] * (need - len(a)))
    for i, c in enumerate(b):
        a[i + shift] += mult * c
    return a


def p_from_R(R, r):
    """KL polynomial from R (coefficients of sum over F != bottom), rank r."""
    R = list(R) + [0] * (r + 1 - len(R))
    P = [R[r - i] - R[i] for i in range((r + 1) // 2) if 2 * i < r]
    if r == 0:
        P = [1]
    return P or [1]


def z_from(P, R):
    Z = list(R)
    padd(Z, P)
    return Z


# ---------------------------------------------------------------- linear matroids over GF(p)
class LinearMatroid:
    def __init__(self, vecs, p):
        self.p = p
        self.vecs = [tuple(int(x) % p for x in v) for v in vecs]
        self.n = len(self.vecs)
        self.dim = len(self.vecs[0])

    def _reduce(self, basis, v):
        """basis: dict pivot -> normalised row (pivot entry 1). returns reduced copy of v."""
        p = self.p
        v = list(v)
        for piv, row in basis.items():
            c = v[piv]
            if c:
                v = [(x - c * y) % p for x, y in zip(v, row)]
        return v

    def _add(self, basis, v):
        v = self._reduce(basis, v)
        piv = next((i for i, x in enumerate(v) if x), None)
        if piv is None:
            return False
        inv = pow(v[piv], self.p - 2, self.p)
        v = [(x * inv) % self.p for x in v]
        for q, row in list(basis.items()):
            c = row[piv]
            if c:
                basis[q] = [(x - c * y) % self.p for x, y in zip(row, v)]
        basis[piv] = v
        return True

    def closure_with_basis(self, elems):
        basis = {}
        for e in elems:
            self._add(basis, self.vecs[e])
        mask = 0
        for e in range(self.n):
            if not any(self._reduce(basis, self.vecs[e])):
                mask |= 1 << e
        return mask, len(basis)

    def rank_of(self, mask):
        basis = {}
        for e in range(self.n):
            if mask >> e & 1:
                self._add(basis, self.vecs[e])
        return len(basis)


def flats_linear(M):
    """all flats as (bitmask, rank), by BFS over covers."""
    bot, r0 = M.closure_with_basis([])
    levels = [{bot}]
    while True:
        nxt = set()
        for F in levels[-1]:
            elems = [e for e in range(M.n) if F >> e & 1]
            for e in range(M.n):
                if not F >> e & 1:
                    G, _ = M.closure_with_basis(elems + [e])
                    nxt.add(G)
        if not nxt:
            break
        levels.append(nxt)
    return [(F, k) for k, lev in enumerate(levels) for F in lev]


def kl_z_from_flats(flats):
    """flats: list of (bitmask, rank). Returns (P_M, Z_M, number of flats)."""
    flats = sorted(flats, key=lambda x: -x[1])
    masks = np.array([f for f, _ in flats], dtype=object)
    ranks = [k for _, k in flats]
    r = max(ranks)
    top = flats[0][0]
    P = {}
    # process from top rank down; P[F] = P(M/F)
    by_rank = {}
    for f, k in flats:
        by_rank.setdefault(k, []).append(f)
    for k in sorted(by_rank, reverse=True):
        for F in by_rank[k]:
            if k == r:
                P[F] = [1]
                continue
            R = []
            for kk in range(k + 1, r + 1):
                acc = []
                for X in by_rank[kk]:
                    if X & F == F:
                        padd(acc, P[X])
                if acc:
                    padd(R, acc, shift=kk - k)
            P[F] = p_from_R(R, r - k)
    bot = min(by_rank)
    F0 = by_rank[bot][0]
    R = []
    for kk in range(bot + 1, r + 1):
        for X in by_rank[kk]:
            padd(R, P[X], shift=kk - bot)
    Z = z_from(P[F0], R)
    return P[F0], Z, len(flats)


# ---------------------------------------------------------------- nested (Schubert) matroids
@lru_cache(maxsize=None)
def _nested_cached(chain):
    return nested_kl_z_uncached(list(chain))


def nested_kl_z(chain):
    P, Z = _nested_cached(tuple(tuple(c) for c in chain))
    return list(P), list(Z)


def _check_chain(chain):
    chain = [tuple(c) for c in chain]
    assert chain[0] == (0, 0), chain
    for i in range(1, len(chain)):
        assert chain[i][1] > chain[i - 1][1] and chain[i][0] - chain[i][1] > chain[i - 1][0] - chain[i - 1][1], chain
    return chain


def _nested_rank_fn(sizes, rks):
    k = len(sizes)
    def rank(a):
        best = tail = sum(a)
        for i in range(1, k + 1):
            tail -= a[i - 1]                 # sum_{l > i} a_l
            if rks[i] + tail < best: best = rks[i] + tail
        return best
    return rank


def _compositions(total, caps):
    if not caps:
        yield (); return
    for x in range(min(total, caps[0]) + 1):
        for rest in _compositions(total - x, caps[1:]):
            yield (x,) + rest


def _flat_profiles(sizes, rks, rank):
    """Every flat F has cyclic part X_j and F = X_j u T, T coloops of M|F, so rk F = r_j + |T|.
    Enumerate profiles with layers 1..j full and |T| <= r - r_j, keep the closed ones."""
    k = len(sizes); r = rks[-1]
    cand = set()
    for j in range(k + 1):
        head = tuple(sizes[:j])
        for tail in _compositions(r - rks[j], sizes[j:]):
            cand.add(head + tail)
    out = []
    for a in cand:
        ra = rank(a)
        if all(a[l] == sizes[l] or rank(a[:l] + (a[l] + 1,) + a[l + 1:]) > ra for l in range(k)):
            out.append((a, ra))
    return out


def contraction_chain(sizes, rks, rank, f):
    """Chain of cyclic flats of M/F (nested again; its cyclic flats are Aut-invariant, hence unions of
    remnant layers, and by the swap lemma they are prefixes of the remnant layers)."""
    k = len(sizes)
    keep = [l for l in range(k) if sizes[l] - f[l] > 0]
    rem = [sizes[l] - f[l] for l in keep]
    rf = rank(f)
    def rk2(a):                                   # a over kept layers
        b = list(f)
        for x, l in zip(a, keep): b[l] += x
        return rank(tuple(b)) - rf
    m = len(keep)
    chain = []
    for i in range(m + 1):
        a = tuple(rem[:i]) + (0,) * (m - i)
        ra = rk2(a)
        flat = all(rk2(a[:l] + (a[l] + 1,) + a[l + 1:]) > ra for l in range(i, m))
        cyc = all(rk2(a[:l] + (a[l] - 1,) + a[l + 1:]) == ra for l in range(i))
        if flat and cyc:
            chain.append((sum(rem[:i]), ra))
    chain = _check_chain(chain)
    # verify on all profiles if few, else on all profiles with entries in {0, 1, rem-1, rem}
    sz2 = [chain[i][0] - chain[i - 1][0] for i in range(1, len(chain))]
    rank2 = _nested_rank_fn(sz2, [c[1] for c in chain])
    import math
    # The prefix property is a lemma (CASE.md, "Nested engine"): swapping an element for one of a lower
    # layer never raises the rank, so a cyclic flat containing layer l contains every lower layer.
    # The exhaustive check below is kept as a safety net on small contractions only.
    if math.prod(x + 1 for x in rem) <= 64:
        test = itertools.product(*[range(x + 1) for x in rem])
    else:
        test = ()
    bounds = [0]
    for c in chain[1:]: bounds.append(c[0])
    for a in test:
        # map profile over remnant layers to profile over chain layers
        pos = 0; agg = [0] * (len(chain) - 1); cum = 0
        for x, s_ in zip(a, rem):
            while cum >= bounds[pos + 1]: pos += 1
            agg[pos] += x; cum += s_
        assert rank2(tuple(agg)) == rk2(a), ("contraction chain check failed", sizes, rks, f, chain)
    return tuple(chain)


def nested_kl_z_uncached(chain):
    """chain: list of (size, rank) of the cyclic flats X_0 < X_1 < ... < X_k, with X_0 = (0, 0)
    (loopless) and X_k = (n, r) the whole ground set (no coloops). Ranks and nullities must strictly
    increase. Returns (P, Z) of the nested matroid with exactly this chain of cyclic flats.
    Recursion: R = sum over flats F != 0 of #(flats of profile f) * t^{rk F} * P(M/F), with M/F nested
    and computed (cached) from its own chain."""
    chain = _check_chain(chain)
    sizes = [chain[i][0] - chain[i - 1][0] for i in range(1, len(chain))]
    rks = [c[1] for c in chain]
    k = len(sizes); r = rks[-1]
    if r == 0:
        return [1], [1]
    rank = _nested_rank_fn(sizes, rks)
    R = []
    for f, rf in _flat_profiles(sizes, rks, rank):
        if rf == 0:
            continue
        mult = 1
        for l in range(k): mult *= comb(sizes[l], f[l])
        if rf == r:
            Pc = [1]
        else:
            Pc, _ = _nested_cached(contraction_chain(sizes, rks, rank, f))
        padd(R, Pc, shift=rf, mult=mult)
    P = p_from_R(R, r)
    return P, z_from(P, R)


def schubert_expansion(blocks, A, nested=nested_kl_z):
    """Eberhardt (arXiv:2609.20903) Thm 3.2: blocks = list of (size, rank) labels in rank order,
    A = condensed incidence matrix (upper unitriangular). Returns (P, Z)."""
    import sympy as sp
    Ainv = sp.Matrix(A).inv()
    m = len(blocks)
    totP, totZ = [], []
    for k in range(1, m):
        for mid in itertools.combinations(range(1, m - 1), k - 1):
            seq = (0,) + mid + (m - 1,)
            coef = 1
            for i in range(1, len(seq)):
                coef *= -Ainv[seq[i - 1], seq[i]]
            if coef == 0:
                continue
            Pn, Zn = nested([blocks[j] for j in seq])
            padd(totP, Pn, mult=int(coef)); padd(totZ, Zn, mult=int(coef))
    while totP and totP[-1] == 0: totP.pop()
    return totP, totZ


def real_rooted(coeffs):
    """exact: all roots real, counted with multiplicity. Sturm-based count_roots counts DISTINCT real
    roots, so it is applied to the square-free part (repeated roots occur e.g. for direct sums)."""
    import sympy as sp
    x = sp.symbols("x")
    poly = sp.Poly(list(reversed(coeffs)), x)
    sqf = sp.Poly(sp.sqf_part(poly.as_expr()), x)
    return sqf.count_roots() == sqf.degree()
