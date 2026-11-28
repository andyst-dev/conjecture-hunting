"""Scan B: matroids of the extended binary Golay code shortened at s and punctured at p coordinates.
M24 is 5-transitive, so for s + p <= 5 each (s, p) gives one code up to equivalence. The matroid
G(s, p) has rank 12 - s on 24 - s - p elements; G(s, p)* = G(p, s) (the Golay code is self-dual).
Output: data/golay_family.jsonl.  Usage: python3 scan_golay.py MAXSUM"""
import sys, json, time
from codes import golay24, shorten, puncture, columns, weight_distribution
from condense import kl_z_binary
from zfast import real_rooted

maxsum = int(sys.argv[1]) if len(sys.argv) > 1 else 5
out = open("data/golay_family.jsonl", "a")
for tot in range(0, maxsum + 1):
    for s in range(tot, -1, -1):
        p = tot - s
        rows, n = golay24(), 24
        if s:
            rows, n = shorten(rows, n, list(range(s)))
        if p:
            # after shortening, coordinates were renumbered 0..n-1; puncture the first p of them
            rows, n = puncture(rows, n, list(range(p)))
        cols, k = columns(rows, n)
        t0 = time.time()
        try:
            res = kl_z_binary(cols, k, verbose=False)
        except AssertionError as e:
            print(f"G({s},{p}): skipped ({e})", flush=True); continue
        zr, pr = real_rooted(res["Z"]), (real_rooted(res["P"]) if len(res["P"]) > 2 else True)
        rec = {"s": s, "p": p, "n": n, "rank": k, "cyclic_flats": res["cyclic_flats"], "blocks": res["blocks"],
               "labels": res["labels"], "P": res["P"], "Z": res["Z"], "Z_real_rooted": zr, "P_real_rooted": pr,
               "weights": weight_distribution(rows, n), "secs": round(time.time() - t0, 1)}
        out.write(json.dumps(rec) + "\n"); out.flush()
        print(f"G({s},{p}): rank {k}, n {n}, cyclic flats {res['cyclic_flats']}, blocks {res['blocks']}, "
              f"Z real-rooted {zr}, P real-rooted {pr}  ({rec['secs']}s)", flush=True)
