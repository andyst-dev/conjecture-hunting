"""Cross-check the DP fiber counter against direct enumeration, and the exact
violation checker on a known non-unimodular violation."""
import random
from lyq import *

def random_config(rng, m, n):
    while True:
        M = [[rng.choice([-1, 0, 0, 1, 1]) for _ in range(n)] for _ in range(m)]
        w = pointing_functional(M, 2)
        if w is not None and Matrix(M).rank() == m:
            return M, w

def test_dp_vs_brute(trials=300, seed=0):
    rng = random.Random(seed)
    for _ in range(trials):
        m = rng.choice([1, 2, 3, 4]); n = rng.randint(m, m + 4)
        M, w = random_config(rng, m, n)
        W = rng.randint(3, 9)
        assert counts(M, w, W) == counts_brute(M, w, W), (M, w, W)
    return trials

def test_exact_checker():
    # M = [1 2] is not unimodular; N = 1,1,2,2,3,... and N(1)^2 < N(0) N(2).
    M = [[1, 2]]; N = counts(M, (1,), 8)
    viol, mu = exact_check(N, (1,), [(0,), (2,)])
    assert viol and mu == [Fraction(1, 2), Fraction(1, 2)]
    hits, gap, _ = screen(N)
    assert hits and abs(gap - math.log(2) / 2) < 1e-9
    # and a unimodular case that must pass the midpoint test (Theorem 1.4)
    M = [[1, 0, 1], [0, 1, 1]]; N = counts(M, (1, 1), 12)
    assert not midpoint_violations(N)

if __name__ == "__main__":
    print("dp vs brute:", test_dp_vs_brute(), "random instances agree")
    test_exact_checker(); print("exact checker / screen detect the known violation: ok")
