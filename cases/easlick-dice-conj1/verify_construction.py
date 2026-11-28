"""Standalone exact verifier for the block construction (imports no other code of this case).

For each N it rebuilds, from data/figure8_g7.json and the rules in CASE.md ("General construction"),
the sequence of P dice and checks from the definitions alone:
  * every die is a non-decreasing N-tuple with faces in {1..N};
  * the P dice are pairwise distinct and P > N;
  * for every i (indices mod P), #{(f, g) : f in D_i, g in D_{i+1}, f > g} = floor(N^2/2) + 1,
    i.e. P(D_i beats D_{i+1}) = 1/2 + 1/N^2 (N even) or 1/2 + 1/(2N^2) (N odd).
The count is the plain pair count (via a running "faces below" table), not the block additivity.
For N <= 60 a few arcs are also recounted with the naive double loop.
Usage: python3 verify_construction.py [NMIN NMAX] [extra N ...]   default: 49 400 and some large N
"""
import json, sys

fig = json.load(open("data/figure8_g7.json"))
A, B, C6 = fig["A"], fig["B"], fig["C6"]


def rules(N):
    d = 1 if N % 2 else 2
    K = (N + d + 7) // 8
    e = 8 * K - N - d
    extra = {0: [], 2: [[1, 1]], 4: [[1, 1], [1, 1]], 6: [[1, 1, 1]]}[e]
    for K6 in range(0, min(3, K - 1) + 1):
        K7 = K - K6
        r = N - 7 * K7 - 6 * K6 - sum(map(len, extra))
        P = len(A) * (K7 - 1) + len(B)
        if r >= 0 and P > N:
            return K7, K6, extra, r, P
    return None


def state(i, t):
    a, b = len(A), len(B)
    if t < a * i:
        return A[t % a]
    if t < a * i + b:
        return B[t - a * i]
    return A[(t - a * i - b) % a]


def dice_for(N):
    K7, K6, extra, r, P = rules(N)
    out = []
    for t in range(P):
        blocks = [state(i, t) for i in range(K7)] + [C6[t % 3]] * K6 + extra + [[1]] * r
        die, off = [], 0
        for blk in blocks:
            die.extend(v + off for v in blk)
            off += len(blk)
        out.append(sorted(die))
    return out


def wins(x, y, N):
    below = [0] * (N + 2)          # below[v] = #{faces of y < v}
    for g in y:
        below[g + 1] += 1
    for v in range(1, N + 2):
        below[v] += below[v - 1]
    return sum(below[f] for f in x)


def check(N):
    cyc = dice_for(N)
    P, T = len(cyc), N * N // 2 + 1
    assert P > N
    for d in cyc:
        assert len(d) == N and all(1 <= v <= N for v in d) and d == sorted(d)
    assert len({tuple(d) for d in cyc}) == P
    for i in range(P):
        assert wins(cyc[i], cyc[(i + 1) % P], N) == T, (N, i)
    if N <= 60:
        for i in range(0, P, max(1, P // 7)):
            x, y = cyc[i], cyc[(i + 1) % P]
            assert sum(1 for f in x for g in y if f > g) == T
    return P


def check_ingredients():
    """The facts about A, B, C6 that Lemmas 2 and 3 use for every N (not only the N run below)."""
    for cyc, m in ((A, 7), (B, 7), (C6, 6)):
        T = m * m // 2 + 1
        assert all(len(d) == m and d == sorted(d) and all(1 <= v <= m for v in d) for d in cyc)
        assert len({tuple(d) for d in cyc}) == len(cyc)                  # simple cycle
        for i in range(len(cyc)):
            assert sum(1 for f in cyc[i] for g in cyc[(i + 1) % len(cyc)] if f > g) == T
    assert len(A) == 24 and len(B) == 3 and len(C6) == 3
    assert A[0] == B[0]                                                  # common vertex x
    assert not ({tuple(d) for d in A} & {tuple(d) for d in B[1:]})       # B1, B2 not on A
    print("ingredients: A (24-cycle) and B (3-cycle) in G_7 meet only at x =", A[0], "; C6 is a 3-cycle of G_6")


if __name__ == "__main__":
    check_ingredients()
    args = [int(x) for x in sys.argv[1:]]
    lo, hi = (args[0], args[1]) if len(args) >= 2 else (49, 400)
    extra_N = args[2:] if len(args) > 2 else ([] if len(args) >= 2 else [500, 777, 1000, 1001, 2000, 2001])
    Ns = list(range(lo, hi + 1)) + extra_N
    for N in Ns:
        P = check(N)
        if N <= 60 or N % 50 == 0 or N in extra_N:
            print(f"N={N}: cycle of length {P} certified", flush=True)
    print(f"ALL CHECKS PASSED for {len(Ns)} values of N ({lo}..{hi}" + (f" and {extra_N})" if extra_N else ")"))
