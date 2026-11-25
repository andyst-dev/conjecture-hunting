"""For each prime power m in a range: build every Singer ruler (multiplier t, largest-gap cut)
and certify that it is the unique shortest point of its chamber.
Optionally (--allcuts) also test every other cyclic cut.
Usage: python3 run.py MMIN MMAX [--allcuts] [--out file]"""
import sys, json, math, time
from singer import singer_set, prime_power
from chamber import certify


def rulers_for(m, allcuts=False, mirror=False):
    """mirror=True keeps one ruler per {ruler, reversed ruler} pair: the reversal map
    g -> (g_k,...,g_1) sends the chamber of a ruler onto the chamber of its mirror and preserves
    length, so 'unique shortest point' holds for both or neither."""
    q = m * m + m + 1
    D = singer_set(m)
    seen = set()
    out = []
    for t in range(1, q):
        if math.gcd(t, q) != 1:
            continue
        S = sorted((t * d) % q for d in D)
        gaps = [((S[(i + 1) % len(S)] - S[i]) % q, i) for i in range(len(S))]
        cuts = sorted(gaps, reverse=True)
        if not allcuts:
            cuts = cuts[:1]
        for rank, (G, i) in enumerate(cuts):
            start = S[(i + 1) % len(S)]
            marks = sorted((s - start) % q for s in S)
            key = tuple(marks)
            if mirror:
                L = key[-1]
                key = min(key, tuple(sorted(L - x for x in key)))
            if key in seen:
                continue
            seen.add(key)
            out.append((t, rank, key))
    return out


if __name__ == "__main__":
    if "--list" in sys.argv:
        ms = [int(x) for x in sys.argv[sys.argv.index("--list") + 1].split(",")]
        mmin, mmax = min(ms), max(ms)
    else:
        mmin, mmax = int(sys.argv[1]), int(sys.argv[2]); ms = list(range(mmin, mmax + 1))
    allcuts = "--allcuts" in sys.argv
    mirror = "--mirror" in sys.argv
    outf = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else f"results_{mmin}_{mmax}{'_all' if allcuts else ''}{'_mir' if mirror else ''}.jsonl"
    tot = {}
    with open(outf, "w") as f:
        for m in ms:
            if not prime_power(m):
                continue
            t0 = time.time()
            R = rulers_for(m, allcuts, mirror)
            cnt = {}
            for t, rank, marks in R:
                status, info = certify(list(marks))
                if status == "unique" and info.get("sum_y") != str(marks[-1]):
                    status = "CHECK_FAILED"   # dual objective must equal the ruler length exactly
                cnt[status] = cnt.get(status, 0) + 1
                if status != "unique" or rank > 0:
                    f.write(json.dumps({"m": m, "t": t, "cut_rank": rank, "marks": marks, "status": status,
                                        "info": info if status != "unique" else None}) + "\n")
                if status != "unique" and rank == 0:
                    print("LARGEST-GAP CUT NOT UNIQUE SHORTEST:", m, t, status, info, flush=True)
            f.write(json.dumps({"m": m, "n": m + 1, "rulers": len(R), "counts": cnt, "secs": round(time.time() - t0, 1)}) + "\n")
            f.flush()
            for s, c in cnt.items():
                tot[s] = tot.get(s, 0) + c
            print(f"m={m} n={m+1} rulers={len(R)} {cnt} {time.time()-t0:.1f}s", flush=True)
    print("TOTAL", tot)
