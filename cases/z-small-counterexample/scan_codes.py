"""Scan C: other dense, highly symmetric binary codes of length < 24 (and their duals).
Output: data/codes.jsonl.  Usage: python3 scan_codes.py"""
import json, time
import sympy as sp
from codes import cyclic_code, extend_parity, reed_muller, columns, weight_distribution, codewords, rref_gf2
from condense import kl_z_binary
from zfast import real_rooted

def qr_extended(p):
    import flint
    f = flint.nmod_poly([1] + [0] * (p - 1) + [1], 2)          # x^p + 1 over GF(2)
    facs = [g for g, _ in f.factor()[1] if g.degree() == (p - 1) // 2]
    g = sum(int(c) << i for i, c in enumerate(facs[0].coeffs()))
    return extend_parity(cyclic_code(p, g), p), p + 1

def dual(rows, n):
    words = codewords(rows)
    k = len(rref_gf2(rows, n))
    # kernel of generator matrix: vectors orthogonal to all rows
    basis = []
    for v in range(1, 1 << n) if n <= 20 else []:
        pass
    # linear algebra: solve via rref
    R = rref_gf2(rows, n)
    piv = [next(c for c in range(n) if r >> c & 1) for r in R]
    free = [c for c in range(n) if c not in piv]
    out = []
    for f in free:
        v = 1 << f
        for r, pc in zip(R, piv):
            if r >> f & 1: v |= 1 << pc
        out.append(v)
    assert all(bin(a & b).count('1') % 2 == 0 for a in out for b in rows)
    return out, n

cands = []
q, n = qr_extended(17); cands.append(("extQR[18,9,6]", q, n))
r1, n16 = reed_muller(1, 4); cands.append(("RM(1,4)[16,5,8]", r1, n16))
r2, _ = reed_muller(2, 4); cands.append(("RM(2,4)[16,11,4]", r2, n16))
q7, n8 = qr_extended(7); cands.append(("extHamming[8,4,4]", q7, n8))
q23, n24 = qr_extended(23); cands.append(("extQR23=Golay[24,12,8] (control)", q23, n24))
out = open("data/codes.jsonl", "a")
for name, rows, n in cands:
    for which, (R, nn) in (("code", (rows, n)), ("dual", dual(rows, n))):
        cols, k = columns(R, nn)
        t0 = time.time()
        try:
            res = kl_z_binary(cols, k, verbose=False)
        except AssertionError as e:
            print(name, which, "skipped:", e, flush=True); continue
        zr = real_rooted(res["Z"]); pr = real_rooted(res["P"]) if len(res["P"]) > 2 else True
        rec = {"name": name, "which": which, "n": nn, "rank": k, "weights": weight_distribution(R, nn),
               "cyclic_flats": res["cyclic_flats"], "blocks": res["blocks"], "P": res["P"], "Z": res["Z"],
               "Z_real_rooted": zr, "P_real_rooted": pr}
        out.write(json.dumps(rec) + "\n"); out.flush()
        print(f"{name} {which}: rank {k} n {nn} cyclic flats {res['cyclic_flats']} blocks {res['blocks']} "
              f"Z real-rooted {zr} P real-rooted {pr} ({time.time()-t0:.1f}s)", flush=True)
