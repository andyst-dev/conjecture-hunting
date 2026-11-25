"""For n in range, test many minimum proper colourings (GK + Kempe-randomised) for wld <= conj(n).
An UNSAT answer means wld > conj(n): a counterexample to Conjecture 6.1 (then re-verified)."""
import sys, random, time, json
from wild import gk_colouring, kempe_shuffle, check_proper, wld_at_most, verify_W, conj_value
n, trials, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(seed); k = conj_value(n); t0 = time.time(); ok = 0
for trial in range(trials):
    col = gk_colouring(n) if trial == 0 else kempe_shuffle(gk_colouring(n), n, rng.choice([5, 50, 500]), rng)
    assert check_proper(col, n)
    W = wld_at_most(col, n, k)
    if W is None:
        print(f"n={n} trial={trial}: NO colour-connecting set of size {k}  <<< CANDIDATE COUNTEREXAMPLE", flush=True)
        with open(f"data/candidates_n{n}.jsonl", "a") as f:
            f.write(json.dumps({"n": n, "k": k, "colouring": [[a, b, c] for (a, b), c in col.items()]}) + "\n")
    else:
        assert verify_W(col, n, W) and len(W) <= k
        ok += 1
print(f"n={n}: {ok}/{trials} colourings have wld <= {k} ({time.time()-t0:.0f}s)", flush=True)
