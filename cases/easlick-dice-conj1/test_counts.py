"""Cross-check the matrix fast path (M @ C^T) against the direct double sum."""
import random
import numpy as np
from dice import dice, mult_matrix, below_matrix, beats_count, build_graph, target

def test_full_graph_small():
    for n in range(3, 7):
        D, A = build_graph(n)
        T = target(n)
        brute = {(i, j) for i, x in enumerate(D) for j, y in enumerate(D) if beats_count(x, y) == T}
        coo = A.tocoo()
        assert set(zip(coo.row.tolist(), coo.col.tolist())) == brute, n

def test_random_pairs():
    rng = random.Random(7)
    for n in range(3, 13):
        D = dice(n); M = mult_matrix(D, n); C = below_matrix(M)
        for _ in range(3000):
            a, b = rng.randrange(len(D)), rng.randrange(len(D))
            assert int(M[a] @ C[b]) == beats_count(D[a], D[b])

def test_target_is_minimal_above_half():
    # counts are integers in [0, n^2]; the paper's constructions attain floor(n^2/2)+1
    for n in range(3, 40):
        k = n // 2
        if n % 2 == 0:
            x = [1]*(k-1) + [2] + [3]*k; y = [1] + [2]*(2*k-1)
        else:
            x = [1]*k + [2] + [3]*k; y = [1]*(k+1) + [2]*k
        assert beats_count(x, y) == target(n) and 2*target(n) > n*n and 2*(target(n)-1) <= n*n

if __name__ == '__main__':
    test_full_graph_small(); test_random_pairs(); test_target_is_minimal_above_half(); print('all tests passed')
