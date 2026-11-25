"""Highly structured colourings: affine (XOR) 1-factorisation of K_{2^k}: colour(i,j) = i xor j,
where every K_4 {a,b,c,d} with a^b^c^d = 0 carries only 3 colours (few rainbow K_4s), and
its vertex-deleted near-1-factorisation of K_{2^k - 1}."""
import sys, itertools, time
from wild import check_proper, wld_at_most, verify_W, conj_value

def xor_col(n):
    return {(i, j): i ^ j for i, j in itertools.combinations(range(n), 2)}

def xor_deleted(n):   # n = 2^k - 1: delete vertex 0 of K_{2^k}
    return {(i - 1, j - 1): i ^ j for i, j in itertools.combinations(range(1, n + 1), 2)}

for n, col in [(4, xor_col(4)), (8, xor_col(8)), (16, xor_col(16)), (7, xor_deleted(7)), (15, xor_deleted(15))] + \
              ([(32, xor_col(32)), (31, xor_deleted(31))] if "--big" in sys.argv else []):
    assert check_proper(col, n), n
    t = time.time(); k = conj_value(n)
    W = wld_at_most(col, n, k)
    res = "wld <= %d" % k if W is not None else "NO set of size %d  <<< CANDIDATE" % k
    if W is not None: assert verify_W(col, n, W)
    print(f"XOR n={n}: {res} ({time.time()-t:.1f}s)", flush=True)
