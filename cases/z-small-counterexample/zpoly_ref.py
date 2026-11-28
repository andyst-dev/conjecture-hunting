"""General Kazhdan–Lusztig and Z-polynomials of matroids, and an exact real-rootedness test.

Matroid input: a list of vectors over GF(p) (linear matroid), optionally truncated or summed.
Lattice of flats is built rank by rank with bitmasks.
KL polynomials of all intervals [F, G] via the defining identity (Elias–Proudfoot–Wakefield):
   t^{rk[F,G]} P_{[F,G]}(1/t) = sum_{F <= X <= G} chi_{[F,X]}(t) P_{[X,G]}(t),   deg P < rk/2.
Z_M(t) = sum_F t^{rk F} P_{[F,E]}(t)  (Proudfoot–Xu–Young). Real-rootedness: exact root count
with sympy (Sturm / real root isolation), not floating point.
"""
import itertools
from functools import lru_cache
import sympy as sp

t = sp.symbols("t")


# ---------- linear matroids over GF(p) ----------
def rank_mod_p(rows, p):
    M = [list(r) for r in rows]; r = 0
    if not M: return 0
    for c in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p); M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
    return r


class Matroid:
    """Matroid given by a rank oracle on frozensets of 0..n-1 (cached)."""
    def __init__(self, n, rank_fn):
        self.n = n; self._rank = rank_fn; self._cache = {}
        self.r = self.rank(frozenset(range(n)))

    def rank(self, S):
        S = frozenset(S)
        if S not in self._cache: self._cache[S] = self._rank(S)
        return self._cache[S]

    def closure(self, S):
        S = frozenset(S); r = self.rank(S)
        return frozenset(list(S) + [e for e in range(self.n) if e not in S and self.rank(S | {e}) == r])

    @staticmethod
    def linear(vecs, p):
        return Matroid(len(vecs), lambda S: rank_mod_p([vecs[i] for i in sorted(S)], p))

    def truncate(self, k):
        return Matroid(self.n, lambda S: min(self.rank(S), k))

    @staticmethod
    def direct_sum(A, B):
        n = A.n + B.n
        return Matroid(n, lambda S: A.rank({i for i in S if i < A.n}) + B.rank({i - A.n for i in S if i >= A.n}))


def lattice_of_flats(M):
    """returns list of flats by rank: levels[k] = set of frozensets"""
    levels = [{M.closure(frozenset())}]
    for k in range(M.r):
        nxt = set()
        for F in levels[k]:
            for e in range(M.n):
                if e not in F:
                    nxt.add(M.closure(F | {e}))
        levels.append(nxt)
    return levels


# ---------- KL and Z polynomials on the lattice ----------
class FlatLattice:
    def __init__(self, levels):
        self.flats = [F for lev in levels for F in lev]
        self.rk = {F: k for k, lev in enumerate(levels) for F in lev}
        self.idx = {F: i for i, F in enumerate(self.flats)}
        self.top = max(self.flats, key=len); self.bot = min(self.flats, key=len)
        self.up = {F: [G for G in self.flats if F <= G] for F in self.flats}
        self._mu = {}; self._P = {}

    def mu(self, F, G):
        key = (F, G)
        if key in self._mu: return self._mu[key]
        if F == G: v = 1
        else: v = -sum(self.mu(F, X) for X in self.up[F] if X <= G and X != G)
        self._mu[key] = v; return v

    def chi(self, F, G):
        r = self.rk[G] - self.rk[F]
        return sp.expand(sum(self.mu(F, X) * t ** (self.rk[G] - self.rk[X]) for X in self.up[F] if X <= G))

    def P(self, F, G):
        key = (F, G)
        if key in self._P: return self._P[key]
        r = self.rk[G] - self.rk[F]
        if r == 0:
            v = sp.Integer(1)
        else:
            # R(t) = sum_{F < X <= G} chi_[F,X](t) P_[X,G](t);  t^r P(1/t) - P(t) = R(t)
            R = sp.expand(sum(self.chi(F, X) * self.P(X, G) for X in self.up[F] if X <= G and X != F))
            # P has degree < r/2: coefficients of P are minus the low-degree coefficients of R
            Rp = sp.Poly(R, t)
            coeffs = {m[0]: c for m, c in zip(Rp.monoms(), Rp.coeffs())}
            v = sp.Integer(0)
            for i in range(0, (r + 1) // 2):
                if 2 * i < r:
                    v += -coeffs.get(i, 0) * t ** i
            v = sp.expand(v)
            # sanity: t^r P(1/t) - P(t) == R
            chk = sp.expand(t ** r * v.subs(t, 1 / t) - v - R)
            assert chk == 0, (r, v, R)
        self._P[key] = v; return v

    def Z(self):
        return sp.expand(sum(t ** self.rk[F] * self.P(F, self.top) for F in self.flats))

    def KL(self):
        return self.P(self.bot, self.top)


def real_rooted(poly):
    """exact: number of real roots (with multiplicity) equals degree"""
    Pp = sp.Poly(poly, t)
    deg = Pp.degree()
    return Pp.count_roots() == deg, deg


def analyse(M):
    L = FlatLattice(lattice_of_flats(M))
    Z = L.Z(); P = L.KL()
    rr, deg = real_rooted(Z)
    return {"rank": M.r, "n": M.n, "flats": len(L.flats), "KL": str(P), "Z": str(sp.Poly(Z, t).all_coeffs()[::-1]),
            "Z_real_rooted": rr, "KL_real_rooted": real_rooted(P)[0] if sp.Poly(P, t).degree() > 0 else True}
