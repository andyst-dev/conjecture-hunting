"""Cross-check of the direct-sum lemma used by the general construction.

For dice X = X_1 (+) ... (+) X_s (block j an m_j-die shifted above the previous blocks),
g_N(X, Y) := 2 W(X, Y) - N^2 equals sum_j (2 W(X_j, Y_j) - m_j^2), and g(Z, Z) = -sum_v mult_Z(v)^2.
Compared against the direct pair count on random block decompositions."""
import random
from collections import Counter


def W(x, y):
    return sum(1 for f in x for g in y if f > g)


def rand_die(m, rng):
    return sorted(rng.randint(1, m) for _ in range(m))


def direct_sum(blocks):
    die, off = [], 0
    for blk in blocks:
        die += [v + off for v in blk]
        off += len(blk)
    return die


def test_additivity(trials=3000, seed=0):
    rng = random.Random(seed)
    for _ in range(trials):
        sizes = [rng.randint(1, 9) for _ in range(rng.randint(1, 6))]
        X = [rand_die(m, rng) for m in sizes]
        Y = [rand_die(m, rng) for m in sizes]
        N = sum(sizes)
        lhs = 2 * W(direct_sum(X), direct_sum(Y)) - N * N
        rhs = sum(2 * W(x, y) - len(x) ** 2 for x, y in zip(X, Y))
        assert lhs == rhs


def test_fixed_block(trials=2000, seed=1):
    rng = random.Random(seed)
    for _ in range(trials):
        m = rng.randint(1, 12)
        z = rand_die(m, rng)
        assert 2 * W(z, z) - m * m == -sum(c * c for c in Counter(z).values())


if __name__ == "__main__":
    test_additivity(); test_fixed_block()
    print("additivity and fixed-block formulas agree with direct counts")
