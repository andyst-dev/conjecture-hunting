"""Cross-checks of the fast engines in zfast.py against the reference zpoly_ref.py (EPW recursion with
Moebius functions and sympy, validated in cases/z-polynomial-rank4/test_zpoly.py).
  1. kl_z_from_flats on random linear matroids over GF(2), GF(3), GF(5), ranks 2-5;
  2. nested_kl_z against the same matroid built explicitly from its rank function;
  3. schubert_expansion on random paving matroids (cyclic flats = stressed hyperplanes) against the
     explicit computation.
Run: python3 test_zfast.py"""
import itertools, random
import sympy as sp
from zpoly_ref import Matroid, FlatLattice, lattice_of_flats, t
from zfast import LinearMatroid, flats_linear, kl_z_from_flats, nested_kl_z, schubert_expansion


def ref(M):
    L = FlatLattice(lattice_of_flats(M))
    P = sp.Poly(sp.expand(L.KL()), t).all_coeffs()[::-1]
    Z = sp.Poly(sp.expand(L.Z()), t).all_coeffs()[::-1]
    return [int(c) for c in P], [int(c) for c in Z]


def strip(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0: a.pop()
    return a


def rand_vecs(rng, p, dim, n):
    pts = set()
    while len(pts) < n:
        v = [rng.randrange(p) for _ in range(dim)]
        if any(v):
            f = next(x for x in v if x); inv = pow(f, p - 2, p)
            pts.add(tuple((x * inv) % p for x in v))
    return list(pts)


def test_linear(trials=60, seed=0):
    rng = random.Random(seed)
    done = 0
    while done < trials:
        p = rng.choice([2, 3, 5]); dim = rng.randint(2, 5)
        maxn = {2: 2 ** dim - 1, 3: (3 ** dim - 1) // 2, 5: (5 ** dim - 1) // 4}[p]
        n = rng.randint(dim, min(maxn, 10))
        vecs = rand_vecs(rng, p, dim, n)
        LM = LinearMatroid(vecs, p)
        P, Z, _ = kl_z_from_flats(flats_linear(LM))
        Pr, Zr = ref(Matroid.linear(vecs, p))
        assert strip(P) == strip(Pr) and strip(Z) == strip(Zr), (vecs, p, P, Pr, Z, Zr)
        done += 1
    return done


def nested_rank_fn(chain):
    """explicit nested matroid on 0..n-1: layer l = X_l \\ X_{l-1}, rk(S) = min_i r_i + |S \\ X_i|."""
    n = chain[-1][0]
    X = [frozenset(range(c[0])) for c in chain]
    return Matroid(n, lambda S: min(c[1] + len(set(S) - Xi) for c, Xi in zip(chain, X)))


def rand_chain(rng):
    while True:
        k = rng.randint(1, 3)
        ch = [(0, 0)]
        ok = True
        for _ in range(k):
            s, r = ch[-1]
            ds = rng.randint(1, 4); dr = rng.randint(1, 3)
            if dr >= ds: ok = False; break      # nullity must increase: ds - dr >= 1
            ch.append((s + ds, r + dr))
        if ok and ch[-1][0] <= 10:
            return ch


def test_nested(trials=60, seed=1):
    rng = random.Random(seed)
    for _ in range(trials):
        ch = rand_chain(rng)
        P, Z = nested_kl_z(ch)
        Pr, Zr = ref(nested_rank_fn(ch))
        assert strip(P) == strip(Pr) and strip(Z) == strip(Zr), (ch, P, Pr, Z, Zr)
    return trials


def test_paving_expansion(trials=25, seed=2):
    """random paving matroids of rank r on n points from a family of pairwise 'r-2 intersecting'
    hyperplanes; compare explicit computation with the Schubert expansion."""
    rng = random.Random(seed)
    done = 0
    while done < trials:
        r = rng.randint(3, 4); n = rng.randint(r + 3, 9)
        hyps = []
        cand = [frozenset(c) for h in range(r, n - 1) for c in itertools.combinations(range(n), h)]
        rng.shuffle(cand)
        for H in cand[:400]:
            if all(len(H & G) <= r - 2 for G in hyps):
                hyps.append(H)
            if len(hyps) >= rng.randint(1, 6): break
        if not hyps: continue
        def rk(S, hyps=hyps, r=r):
            S = set(S)
            if len(S) < r: return len(S)
            return r - 1 if any(S <= H for H in hyps) else r
        Pr, Zr = ref(Matroid(n, rk))
        sizes = sorted({len(H) for H in hyps})
        blocks = [(0, 0)] + [(h, r - 1) for h in sizes] + [(n, r)]
        m = len(blocks)
        A = [[0] * m for _ in range(m)]
        for i in range(m): A[i][i] = 1
        for j in range(1, m): A[0][j] = 1
        for i, h in enumerate(sizes, start=1): A[i][m - 1] = sum(1 for H in hyps if len(H) == h)
        P, Z = schubert_expansion(blocks, A)
        assert strip(P) == strip(Pr) and strip(Z) == strip(Zr), (r, n, hyps, P, Pr, Z, Zr)
        done += 1
    return done


def test_condense(trials=40, seed=3):
    """cycflats + partition refinement + Theorem 3.2 (condense.py) against the explicit computation, on
    random binary matroids without loops or coloops."""
    from condense import kl_z_binary
    rng = random.Random(seed)
    done = 0
    while done < trials:
        k = rng.randint(2, 5); n = rng.randint(k + 2, 11)
        cols = [rng.randrange(1, 2 ** k) for _ in range(n)]
        vecs = [[(c >> j) & 1 for j in range(k)] for c in cols]
        LM = LinearMatroid(vecs, 2)
        if LM.rank_of((1 << n) - 1) != k: continue
        if any(LM.rank_of(((1 << n) - 1) ^ (1 << e)) < k for e in range(n)): continue   # coloop
        P, Z, _ = kl_z_from_flats(flats_linear(LM))
        res = kl_z_binary(cols, k, verbose=False)
        assert strip(res["P"]) == strip(P) and strip(res["Z"]) == strip(Z), (cols, k, res, P, Z)
        done += 1
    return done


def nested_full_profiles(chain):
    """independent slow version: enumerate ALL profiles (the first implementation of nested_kl_z)."""
    from math import comb
    from zfast import p_from_R, padd, z_from
    sizes = [chain[i][0] - chain[i - 1][0] for i in range(1, len(chain))]
    rks = [c[1] for c in chain]; k = len(sizes); r = chain[-1][1]
    def rank(a):
        return min([sum(a)] + [rks[i] + sum(a[i:]) for i in range(1, k + 1)])
    prof = list(itertools.product(*[range(s + 1) for s in sizes]))
    rk = {a: rank(a) for a in prof}
    flats = [a for a in prof if all(a[l] == sizes[l] or rk[a[:l] + (a[l] + 1,) + a[l + 1:]] > rk[a] for l in range(k))]
    flats.sort(key=lambda a: -rk[a]); P = {}
    for a in flats:
        if rk[a] == r: P[a] = [1]; continue
        R = []
        for x in flats:
            if rk[x] > rk[a] and all(xi >= ai for xi, ai in zip(x, a)):
                m = 1
                for l in range(k): m *= comb(sizes[l] - a[l], x[l] - a[l])
                padd(R, P[x], shift=rk[x] - rk[a], mult=m)
        P[a] = p_from_R(R, r - rk[a])
    a0 = tuple([0] * k); R = []
    for x in flats:
        if x != a0:
            m = 1
            for l in range(k): m *= comb(sizes[l], x[l])
            padd(R, P[x], shift=rk[x], mult=m)
    return P[a0], z_from(P[a0], R)


def test_nested_long(trials=80, seed=4):
    """restricted-profile nested engine vs full profile enumeration, chains of up to 6 layers, n <= 22."""
    rng = random.Random(seed)
    done = 0
    while done < trials:
        ch = [(0, 0)]
        for _ in range(rng.randint(1, 6)):
            s, r = ch[-1]; ds = rng.randint(1, 5); dr = rng.randint(1, ds - 1) if ds > 1 else None
            if dr is None: break
            ch.append((s + ds, r + dr))
        if len(ch) < 2 or ch[-1][0] > 22: continue
        assert [strip(x) for x in nested_kl_z(ch)] == [strip(x) for x in nested_full_profiles(ch)], ch
        done += 1
    return done


def test_real_rooted():
    """repeated roots must count as real (bug found 2026-11-28: Sturm counts distinct roots)."""
    from zfast import real_rooted, nested_kl_z
    assert real_rooted([1, 2, 1]) and real_rooted([1, 3, 3, 1])
    assert not real_rooted([1, 0, 1]) and not real_rooted([1, 2, 3, 2, 1])
    _, Z = nested_kl_z([(0, 0), (9, 8)])                     # U_{8,9}
    Z3 = [1]
    for _ in range(3):
        Z3 = [sum(Z3[i] * Z[j - i] for i in range(len(Z3)) if 0 <= j - i < len(Z)) for j in range(len(Z3) + len(Z) - 1)]
    assert real_rooted(Z3)                                   # Z of U89 + U89 + U89: repeated roots
    return True


if __name__ == "__main__":
    print("real-rootedness with multiplicities:", test_real_rooted())
    print("nested, long chains vs full profiles:", test_nested_long(), "agree")
    print("condensed Schubert expansions (binary):", test_condense(), "agree")
    print("linear matroids:", test_linear(), "agree")
    print("nested matroids:", test_nested(), "agree")
    print("paving Schubert expansions:", test_paving_expansion(), "agree")
