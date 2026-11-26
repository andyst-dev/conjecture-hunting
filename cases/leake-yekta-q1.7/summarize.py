"""Coverage table from data/coverage_*.jsonl."""
import json, glob, collections
rows = collections.defaultdict(lambda: [0, 0, 0, 0.0, 0, 0, set()])
for f in sorted(glob.glob("data/coverage_*.jsonl")):
    for l in open(f):
        r = json.loads(l); m = len(r["M"]); n = len(r["M"][0])
        cls = r["class"] if r["class"].startswith("ct") else "%s m=%d" % (r["class"], m)
        x = rows[cls]; x[0] += 1; x[1] += r["domain"]; x[2] = max(x[2], r["domain"])
        x[3] = max(x[3], r["maxgap"]); x[4] += r["float_hits"]; x[5] += r["exact_hits"]; x[6].add(n)
print("| Class | instances | n (columns) | domain points (total / max) | max LP gap | float hits | exact violations |")
print("|---|---|---|---|---|---|---|")
for k in sorted(rows):
    x = rows[k]
    print("| %s | %d | %d–%d | %d / %d | %.1e | %d | %d |" % (k, x[0], min(x[6]), max(x[6]), x[1], x[2], x[3], x[4], x[5]))
