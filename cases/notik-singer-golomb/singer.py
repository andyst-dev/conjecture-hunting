"""Singer perfect difference sets mod q = m^2+m+1 for any prime power m.

GF(m) is built as GF(p)[y]/(h), with h irreducible of degree e (m = p^e), using log/antilog tables.
GF(m^3) is built as GF(m)[x]/(f), with f a primitive monic cubic over GF(m).
D = { i in [0,q) : x^i lies in span_GF(m){1, x} }, i.e. the x^2-coefficient of x^i is 0.
This is a projective line in PG(2,m), so |D| = m+1.
"""
import itertools
from functools import lru_cache


def factor(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def prime_power(m):
    f = factor(m)
    if len(f) != 1:
        return None
    (p, e), = f.items()
    return p, e


class GFm:
    """GF(p^e) elements as ints 0..m-1 (base-p digits = poly coefficients)."""
    def __init__(self, p, e):
        self.p, self.e, self.m = p, e, p ** e
        m = self.m
        if e == 1:
            g = next(g for g in range(2, p) if all(pow(g, (p - 1) // r, p) != 1 for r in factor(p - 1))) if p > 2 else 1
            elems = [pow(g, i, p) for i in range(p - 1)]
        else:
            # find primitive poly of degree e over GF(p): brute force
            for coeffs in itertools.product(range(p), repeat=e):
                h = list(coeffs) + [1]  # monic, low->high
                if h[0] == 0:
                    continue
                elems = self._powers(h)
                if elems is not None:
                    break
        self.exp = elems + elems
        self.log = {v: i for i, v in enumerate(elems)}
        self.add_t = [[self._add(a, b) for b in range(m)] for a in range(m)]
        self.neg = [self._neg(a) for a in range(m)]

    def _digits(self, a):
        return [(a // self.p ** i) % self.p for i in range(self.e)]

    def _num(self, d):
        return sum(c * self.p ** i for i, c in enumerate(d))

    def _add(self, a, b):
        return self._num([(x + y) % self.p for x, y in zip(self._digits(a), self._digits(b))])

    def _neg(self, a):
        return self._num([(-x) % self.p for x in self._digits(a)])

    def _powers(self, h):
        p, e = self.p, self.e
        cur = [1] + [0] * (e - 1)
        out, seen = [], set()
        for _ in range(p ** e - 1):
            v = sum(c * p ** i for i, c in enumerate(cur))
            if v in seen:
                return None
            seen.add(v); out.append(v)
            # multiply by y
            nxt = [0] + cur[:-1]
            top = cur[-1]
            nxt = [(nxt[i] - top * h[i]) % p for i in range(e)]
            cur = nxt
        return out if len(seen) == p ** e - 1 else None

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def add(self, a, b):
        return self.add_t[a][b]

    def sub(self, a, b):
        return self.add_t[a][self.neg[b]]


def singer_set(m):
    pe = prime_power(m)
    assert pe, m
    F = GFm(*pe)
    q = m * m + m + 1
    order = m ** 3 - 1
    primes = list(factor(order))
    # search primitive monic cubic f = x^3 + c2 x^2 + c1 x + c0
    for c0 in range(1, m):
        for c1 in range(m):
            for c2 in range(m):
                f = (c0, c1, c2)
                D = _try(F, f, m, q, order, primes)
                if D is not None:
                    return sorted(D)
    raise RuntimeError


def _mulx(F, v, f):
    a0, a1, a2 = v
    c0, c1, c2 = f
    # x*(a0 + a1 x + a2 x^2) = a2 x^3 + a1 x^2 + a0 x ; x^3 = -(c2 x^2 + c1 x + c0)
    return (F.sub(0, F.mul(a2, c0)), F.sub(a0, F.mul(a2, c1)), F.sub(a1, F.mul(a2, c2)))


def _try(F, f, m, q, order, primes):
    # quick irreducibility: no root in GF(m)
    for r in range(m):
        r2 = F.mul(r, r); r3 = F.mul(r2, r)
        val = F.add(F.add(r3, F.mul(f[2], r2)), F.add(F.mul(f[1], r), f[0]))
        if val == 0:
            return None
    # walk powers of x up to q collecting D, then check primitivity (order of x == m^3-1)
    v = (1, 0, 0)
    D = []
    for i in range(q):
        if v[2] == 0:
            D.append(i)
        v = _mulx(F, v, f)
    # x^q must be a nonzero scalar for D to be valid; primitivity check via powering
    if len(D) != m + 1:
        return None
    if not _is_primitive(F, f, order, primes):
        return None
    return D


def _polymul(F, a, b, f):
    # multiply two elements of GF(m)[x]/(f)
    res = [0] * 5
    for i in range(3):
        if a[i] == 0: continue
        for j in range(3):
            if b[j] == 0: continue
            res[i + j] = F.add(res[i + j], F.mul(a[i], b[j]))
    c0, c1, c2 = f
    for d in (4, 3):
        t = res[d]
        if t:
            res[d] = 0
            res[d - 1] = F.sub(res[d - 1], F.mul(t, c2))
            res[d - 2] = F.sub(res[d - 2], F.mul(t, c1))
            res[d - 3] = F.sub(res[d - 3], F.mul(t, c0))
    return tuple(res[:3])


def _pow(F, a, n, f):
    r = (1, 0, 0)
    while n:
        if n & 1:
            r = _polymul(F, r, a, f)
        a = _polymul(F, a, a, f); n >>= 1
    return r


def _is_primitive(F, f, order, primes):
    x = (0, 1, 0)
    return all(_pow(F, x, order // r, f) != (1, 0, 0) for r in primes)


def is_perfect_difference_set(D, q):
    seen = set()
    for a in D:
        for b in D:
            if a != b:
                d = (a - b) % q
                if d in seen:
                    return False
                seen.add(d)
    return len(seen) == q - 1


if __name__ == "__main__":
    for m in [2, 3, 4, 5, 7, 8, 9, 11, 16, 25, 27, 32]:
        D = singer_set(m)
        q = m * m + m + 1
        print(m, q, len(D), is_perfect_difference_set(D, q))
