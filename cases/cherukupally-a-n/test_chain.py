"""Cross-check the chain characterisation (paper, Lemma 3.2 / Cor. 3.3) against the
definition of A(n), and reproduce OEIS A394398."""
from sympy import divisors
from collections import Counter

def A_def(n):
    """A(n) = {1 <= a < n : n | a^2-1, a | n^2-1}, straight from the definition."""
    return [a for a in range(1, n) if (a * a - 1) % n == 0 and (n * n - 1) % a == 0]

def A_size_div(n):
    return sum(1 for a in divisors(n * n - 1) if a < n and (a * a - 1) % n == 0)

def chain_counts(X):
    """#{(i,z): i >= 2, z >= 3, G_i(z) = v} for v <= X, G_0=1, G_1=z, G_i = z G_{i-1} - G_{i-2}."""
    c = Counter(); z = 3
    while z * z - 1 <= X:
        g0, g1 = 1, z
        while True:
            g0, g1 = g1, z * g1 - g0
            if g1 > X: break
            c[g1] += 1
        z += 1
    return c

if __name__ == "__main__":
    oeis = [0, 1, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 3]  # A394398, n=1..21
    ours = [len(A_def(n)) for n in range(1, 22)]
    assert ours == oeis, ours
    print("A394398 n=1..21 reproduced:", ours)
    assert A_def(15) == [1, 4, 14] and A_def(35) == [1, 6, 34] and A_def(21) == [1, 8, 20]
    for n in range(2, 400):
        assert len(A_def(n)) == A_size_div(n)
    N = 20000; c = chain_counts(N); dist = Counter()
    for n in range(3, N):
        s = A_size_div(n)
        assert s == 2 + c.get(n, 0), n
        dist[s] += 1
    print("|A(n)| = 2 + #chain occurrences for all 3 <= n < %d; distribution %s" % (N, dict(dist)))
