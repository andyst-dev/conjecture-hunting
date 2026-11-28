"""Scan A: paving matroids with stressed hyperplanes of one size h.
For a paving matroid of rank r on n elements whose stressed hyperplanes (rank r-1, size >= r) all have
size h, with lam of them, the Schubert expansion gives (Eberhardt Thm 3.2, cross-checked in test_zfast.py)
    Z = Z(U_{r,n}) + lam * (Z(S_h) - Z(U_{r,n})),   S_h = nested chain (0,0) < (h, r-1) < (n, r).
Every (r-1)-set lies in exactly one hyperplane, so lam * C(h, r-1) <= C(n, r-1), with equality iff the
hyperplanes form a Steiner system S(r-1, h, n).
For each (r, n, h) we screen lam = 1..U, U = packing_upper(n, h, r-2) <= lam_max, (all integers if lam_max <= 3000, else a grid) in floating
point and confirm every hit exactly (sympy root count). Output: data/scan_paving.jsonl
Usage: python3 scan_paving.py RMIN RMAX NMAX"""
import sys, json, math
import numpy as np
from zfast import nested_kl_z, real_rooted
from packing_bounds import packing_upper


def float_nonreal(Z):
    roots = np.roots(Z[::-1])
    return np.max(np.abs(roots.imag)) > 1e-7 * max(1.0, np.max(np.abs(roots)))


def lam_grid(lmax):
    if lmax <= 3000:
        return list(range(1, lmax + 1))
    g = sorted({int(round(x)) for x in np.geomspace(1, lmax, 600)} | {lmax})
    return g


if __name__ == "__main__":
    rmin, rmax, nmax = map(int, sys.argv[1:4])
    out = open("data/scan_paving.jsonl", "a")
    for r in range(rmin, rmax + 1):
        for n in range(r + 3, nmax + 1):
            _, ZU = nested_kl_z([(0, 0), (n, r)])
            for h in range(r, n - 1):
                _, ZS = nested_kl_z([(0, 0), (h, r - 1), (n, r)])
                D = [a - b for a, b in zip(ZS, ZU)]
                lmax = math.comb(n, r - 1) // math.comb(h, r - 1)
                U = packing_upper(n, h, r - 2)          # constant-weight-code bound (tighter)
                first_bad = None
                for lam in lam_grid(U):
                    Z = [u + lam * d for u, d in zip(ZU, D)]
                    if float_nonreal(np.array(Z, dtype=float)):
                        if not real_rooted(Z):
                            first_bad = lam; break
                # behaviour far beyond the packing bound (is D itself the culprit?)
                big = [u + 10**12 * d for u, d in zip(ZU, D)]
                steiner = (math.comb(n, r - 1) % math.comb(h, r - 1) == 0) and U == lmax
                rec = {"r": r, "n": n, "h": h, "lam_max": lmax, "U": U, "steiner_bound_attainable": steiner, "first_bad_lam": first_bad,
                       "bad_at_1e12": (not real_rooted(big))}
                out.write(json.dumps(rec) + "\n"); out.flush()
                if first_bad is not None:
                    print("NON-REAL-ROOTED within packing bound:", rec, flush=True)
        print(f"r={r} done up to n={nmax}", flush=True)
