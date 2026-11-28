"""Upper bounds on the number of h-subsets of [n] pairwise meeting in <= s points
(= constant-weight code size A(n, d, h) with d = 2(h - s)).
Uses A(n,d,w) = A(n,d,n-w), Johnson A(n,d,w) <= floor(n/w * A(n-1,d,w-1)), and the trivial cases."""
from functools import lru_cache

@lru_cache(None)
def A_up(n, d, w):
    if w < 0 or w > n: return 0
    if d > 2 * min(w, n - w): return 1          # any two distinct words would be too close
    if d <= 2 and True: pass
    if w == 0 or w == n: return 1
    b1 = (n * A_up(n - 1, d, w - 1)) // w
    b2 = (n * A_up(n - 1, d, n - w - 1)) // (n - w)   # complement: weight n-w
    return min(b1, b2)

def packing_upper(n, h, s):
    return A_up(n, 2 * (h - s), h)
