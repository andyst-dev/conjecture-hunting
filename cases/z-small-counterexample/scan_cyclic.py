"""Scan D: every binary cyclic code of odd length n in [NMIN, NMAX] (one per divisor g of x^n + 1 over
GF(2), excluding the zero and full codes), and its extension by an overall parity bit (length n + 1).
The column matroid is analysed when it has no loops and no coloops and at most MAXCF cyclic flats.
Output: data/cyclic.jsonl.  Usage: python3 scan_cyclic.py NMIN NMAX [MAXCF]"""
import sys, json, time, itertools
import flint
from codes import cyclic_code, extend_parity, columns, weight_distribution
from condense import cyclic_flats_binary, condense, expansion
from zfast import real_rooted

nmin, nmax = int(sys.argv[1]), int(sys.argv[2])
maxcf = int(sys.argv[3]) if len(sys.argv) > 3 else 60000
out = open("data/cyclic.jsonl", "a")


def analyse(rows, n, tag):
    cols, k = columns(rows, n)
    if k == 0 or k == n or any(c == 0 for c in cols):
        return {"tag": tag, "n": n, "rank": k, "skip": "trivial or loops"}
    flats, info = cyclic_flats_binary(cols, k)
    if len(flats) > maxcf:
        return {"tag": tag, "n": n, "rank": k, "skip": f"{len(flats)} cyclic flats > {maxcf}"}
    members, labels, A = condense(flats)
    if labels[-1][0] != n:
        return {"tag": tag, "n": n, "rank": k, "skip": "coloops"}
    P, Z, nseq = expansion(labels, A)
    return {"tag": tag, "n": n, "rank": k, "cyclic_flats": len(flats), "blocks": len(labels),
            "P": P, "Z": Z, "Z_real_rooted": real_rooted(Z),
            "P_real_rooted": real_rooted(P) if len(P) > 2 else True}


for n in range(nmin | 1, nmax + 1, 2):
    f = flint.nmod_poly([1] + [0] * (n - 1) + [1], 2)
    facs = [g for g, e in f.factor()[1] for _ in range(e)]
    seen = set()
    for mask in range(1, 2 ** len(facs) - 1):
        g = flint.nmod_poly([1], 2)
        for i, h in enumerate(facs):
            if mask >> i & 1: g *= h
        gi = sum(int(c) << i for i, c in enumerate(g.coeffs()))
        if gi in seen: continue
        seen.add(gi)
        rows = cyclic_code(n, gi)
        for ext in (False, True):
            R, N = (extend_parity(rows, n), n + 1) if ext else (rows, n)
            tag = f"n={n} g={bin(gi)}" + (" +parity" if ext else "")
            t0 = time.time()
            rec = analyse(R, N, tag)
            rec["secs"] = round(time.time() - t0, 1)
            out.write(json.dumps(rec) + "\n"); out.flush()
            if "skip" not in rec:
                flag = "" if rec["Z_real_rooted"] else "   <-- Z NOT REAL-ROOTED"
                print(f"{tag}: rank {rec['rank']} n {N} cf {rec['cyclic_flats']} Z rr {rec['Z_real_rooted']} "
                      f"P rr {rec['P_real_rooted']}{flag}", flush=True)
    print(f"length {n} done", flush=True)
