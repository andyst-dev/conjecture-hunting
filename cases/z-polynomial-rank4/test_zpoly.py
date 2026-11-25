"""Validation of zpoly.py against known values:
 - U_{n-1,n} (uniform corank 1): KL coefficients known (e.g. P_{U_{3,4}} = 1 + 2t ... ) -> checked via PXY palindromicity
 - braid matroid M(K4) (rank 3): P = 1 + t ; M(K5) (rank 4): P = 1 + 5t ... (EPW: P_{K4}=1+t? computed below)
 - Z palindromic of degree rk M for all tests
 - rank-4 closed form Z = 1 + W3 t + a2 t^2 + W3 t^3 + t^4 matches the general computation
 - Z(M1 + M2) = Z(M1) Z(M2)"""
import itertools, random
import sympy as sp
from zpoly import Matroid, FlatLattice, lattice_of_flats, t, real_rooted

def uniform(r, n):
    return Matroid(n, lambda S: min(len(S), r))

def graphic(edges, nv):
    def rk(S):
        parent = list(range(nv)); c = 0
        def f(x):
            while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for i in S:
            a, b = f(edges[i][0]), f(edges[i][1])
            if a != b: parent[a] = b; c += 1
        return c
    return Matroid(len(edges), rk)

def info(M):
    L = FlatLattice(lattice_of_flats(M))
    return L, sp.expand(L.KL()), sp.expand(L.Z())

K4 = graphic(list(itertools.combinations(range(4), 2)), 4)
K5 = graphic(list(itertools.combinations(range(5), 2)), 5)
for name, M in [("U(3,4)", uniform(3, 4)), ("U(3,5)", uniform(3, 5)), ("U(4,6)", uniform(4, 6)), ("M(K4)", K4), ("M(K5)", K5)]:
    L, P, Z = info(M)
    Zp = sp.Poly(Z, t).all_coeffs()
    print(f"{name}: rank {M.r}, KL = {P}, Z = {Z}, palindromic={Zp == Zp[::-1]}, Z real-rooted={real_rooted(Z)[0]}")
# known: P_{M(K4)} = 1 + t? (EPW Prop: braid K4 rank 3: coefficient = W2 - W1 = 7 - 6 = 1) ; P_{M(K5)} = 1 + 5t (W3 - W1 = 15 - 10 = 5)
# rank-4 closed form on random GF(3) matroids
from zpoly import rank_mod_p
rng = random.Random(0)
for _ in range(5):
    vecs = set()
    while len(vecs) < 9:
        v = tuple(rng.randrange(3) for _ in range(4))
        if any(v):
            f = next(x for x in v if x); inv = pow(f, 1, 3) if f == 1 else 2
            vecs.add(tuple((x * inv) % 3 for x in v))
    vecs = list(vecs)
    M = Matroid.linear(vecs, 3)
    if M.r != 4: continue
    levels = lattice_of_flats(M)
    W2, W3 = len(levels[2]), len(levels[3])
    a2 = W2 + sum(len(H) for H in levels[3]) - sum(len(Lf) for Lf in levels[2])
    L = FlatLattice(levels); Z = sp.expand(L.Z())
    print("rank-4 closed form matches:", Z == sp.expand(1 + W3*t + a2*t**2 + W3*t**3 + t**4))
# multiplicativity
A, B = uniform(2, 3), uniform(3, 4)
ZA, ZB, ZAB = info(A)[2], info(B)[2], info(Matroid.direct_sum(A, B))[2]
print("Z(A+B) == Z(A)Z(B):", sp.expand(ZA * ZB - ZAB) == 0)
