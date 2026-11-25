"""Reproduce Theorem 3.23 values (n = 3..12) on GK and on Kempe-randomised colourings."""
import random, time
from wild import gk_colouring, kempe_shuffle, check_proper, wld_exact, verify_W, conj_value
paper = {3: 2, 4: 2, 5: 4, 6: 4, 7: 5, 8: 5, 9: 6, 10: 6, 11: 8, 12: 8}
rng = random.Random(0)
for n in range(3, 13):
    vals = []
    t = time.time()
    for trial in range(6):
        col = gk_colouring(n) if trial == 0 else kempe_shuffle(gk_colouring(n), n, 200, rng)
        assert check_proper(col, n)
        k, W = wld_exact(col, n)
        assert verify_W(col, n, W)
        vals.append(k)
    print(f"n={n}: paper {paper[n]}, formula {conj_value(n)}, ours {vals}  ({time.time()-t:.1f}s)", flush=True)
