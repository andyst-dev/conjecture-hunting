"""Long-running search for a matroid whose Z-polynomial is not real-rooted (PXY18, Conj. 5.1).

Strategy: hill-climbing over linear matroids (sets of projective points of GF(p)^r, plus optional
truncation to rank k < r). Objective: the gamma-polynomial of the palindromic Z,
  Z(t) = sum_i gamma_i t^i (1+t)^(d-2i);
Z is real-rooted iff gamma(x) = sum gamma_i x^i has only real roots (gamma_i >= 0 is known).
Margin = smallest |difference| between two roots of gamma, normalised by the largest root
(roots about to collide and turn complex give a small margin). Moves: add / remove / replace a point.
Any non-real-rooted Z is written to data/z_candidates.jsonl and must then be checked with
verify_z.py (exact, independent).
Usage: python3 search_z.py --rank 5 [--trunc 5] --p 3 --minutes 60 --seed 1
"""
import argparse, json, math, random, time
import numpy as np
import sympy as sp
from zpoly import Matroid, FlatLattice, lattice_of_flats, t, real_rooted


def gamma_coeffs(Zc):
    """Zc: coefficient list low->high of a palindromic polynomial of degree d"""
    d = len(Zc) - 1
    rem = list(Zc); g = []
    for i in range(d // 2 + 1):
        gi = rem[i]; g.append(gi)
        # subtract gi * t^i (1+t)^(d-2i)
        for j in range(d - 2 * i + 1):
            rem[i + j] -= gi * math.comb(d - 2 * i, j)
    assert all(abs(x) == 0 for x in rem), rem
    return g


def margin(g):
    g = [float(x) for x in g]
    while len(g) > 1 and g[-1] == 0: g.pop()
    if len(g) <= 2: return float("inf")
    roots = np.roots(g[::-1])
    if np.max(np.abs(roots.imag)) > 1e-9 * max(1, np.max(np.abs(roots))): return -1.0
    rr = np.sort(roots.real)
    return float(np.min(np.diff(rr)) / max(1e-12, np.max(np.abs(rr))))


def projective_points(r, p):
    pts = []
    for v in np.ndindex(*([p] * r)):
        if any(v):
            f = next(x for x in v if x)
            if f == 1: pts.append(tuple(int(x) for x in v))
    return pts


def connected_linear(vecs, p):
    """matroid connectivity of the column matroid: RREF, basis = pivot columns; element x outside the
    basis is joined to every basis element in its fundamental circuit; components = graph components."""
    r = len(vecs[0]); n = len(vecs)
    A = [[vecs[j][i] % p for j in range(n)] for i in range(r)]
    piv = []; row = 0
    for c in range(n):
        pr = next((i for i in range(row, r) if A[i][c]), None)
        if pr is None: continue
        A[row], A[pr] = A[pr], A[row]
        inv = pow(A[row][c], p - 2, p); A[row] = [(x * inv) % p for x in A[row]]
        for i in range(r):
            if i != row and A[i][c]:
                f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[row])]
        piv.append(c); row += 1
        if row == r: break
    parent = list(range(n))
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for c in range(n):
        if c in piv: continue
        for k, pc in enumerate(piv):
            if A[k][c]: parent[find(c)] = find(pc)
    return len({find(x) for x in range(n)}) == 1


def evaluate(vecs, p, trunc):
    if not trunc and not connected_linear(vecs, p):
        return -1, None          # skip direct sums (their Z factorises; double roots are trivial)
    M = Matroid.linear(vecs, p)
    if trunc and M.r > trunc: M = M.truncate(trunc)
    L = FlatLattice(lattice_of_flats(M))
    Z = sp.Poly(sp.expand(L.Z()), t).all_coeffs()[::-1]
    Z = [int(c) for c in Z]
    return M.r, Z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank", type=int, default=5); ap.add_argument("--trunc", type=int, default=0)
    ap.add_argument("--p", type=int, default=3); ap.add_argument("--minutes", type=float, default=30)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--maxn", type=int, default=16)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    allpts = projective_points(a.rank, a.p)
    target_rank = a.trunc or a.rank
    end = time.time() + 60 * a.minutes
    best = (float("inf"), None); evals = 0
    log = open("data/z_search.jsonl", "a")
    while time.time() < end:
        # restart
        n = rng.randint(a.rank + 1, a.maxn)
        cur = rng.sample(allpts, n)
        try:
            rk, Z = evaluate(cur, a.p, a.trunc)
        except AssertionError:
            continue
        if rk != target_rank: continue
        cm = margin(gamma_coeffs(Z)); evals += 1
        for step in range(200):
            if time.time() > end: break
            cand = list(cur)
            mv = rng.random()
            if mv < 0.34 and len(cand) < a.maxn:
                cand.append(rng.choice([q for q in allpts if q not in cand]))
            elif mv < 0.67 and len(cand) > a.rank + 1:
                cand.pop(rng.randrange(len(cand)))
            else:
                cand[rng.randrange(len(cand))] = rng.choice([q for q in allpts if q not in cand])
            rk2, Z2 = evaluate(cand, a.p, a.trunc); evals += 1
            if rk2 != target_rank: continue
            m2 = margin(gamma_coeffs(Z2))
            if m2 < 0:
                rec = {"p": a.p, "rank": a.rank, "trunc": a.trunc, "vectors": cand, "Z": Z2}
                exact_rr = real_rooted(sum(c * t ** i for i, c in enumerate(Z2)))[0]
                rec["exact_real_rooted"] = exact_rr
                with open("data/z_candidates.jsonl", "a") as f: f.write(json.dumps(rec) + "\n")
                print("CANDIDATE (float says non-real; exact says real_rooted=%s)" % exact_rr, rec, flush=True)
            if m2 <= cm:
                cur, cm = cand, m2
            if cm < best[0]:
                best = (cm, {"n": len(cur), "Z": Z2 if cur is cand else Z, "vectors": cur})
        log.write(json.dumps({"p": a.p, "rank": a.rank, "trunc": a.trunc, "restart_best_margin": cm, "n": len(cur)}) + "\n"); log.flush()
    print(f"done: {evals} evaluations, best margin {best[0]:.4g}", json.dumps(best[1])[:400], flush=True)


if __name__ == "__main__":
    main()
