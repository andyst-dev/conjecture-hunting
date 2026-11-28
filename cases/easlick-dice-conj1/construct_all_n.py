"""Block construction of a uniform minimal non-transitive cycle of length P > N in G_N, for N >= 49.

Ingredients (data/figure8_g7.json): in G_7 a simple 24-cycle A and a simple 3-cycle B with
A[0] = B[0] = x and no other common vertex; in G_6 a 3-cycle C6.

Direct sum: an N-die X = X_1 (+) ... (+) X_s, where block X_j is an m_j-die placed on the value
range just above the previous blocks. With g_n(X, Y) = 2 W(X, Y) - n^2 one has
g_N(X, Y) = sum_j g_{m_j}(X_j, Y_j), and an arc of G_n is exactly g_n = delta(n)
(delta = 1 for odd n, 2 for even n). A fixed block Z contributes g(Z, Z) = -sum_v mult_Z(v)^2.

Blocks used for N:
  K7 counter 7-blocks, block i walks A^i B A^(K7-1-i)       (period P = 24 (K7-1) + 3, +1 each)
  K6 passenger 6-blocks walking C6 repeatedly (3 | P)         (+2 each)
  excess blocks, chosen by e = 8(K7+K6) - N - delta(N) in {0,2,4,6}:
      e = 2: one (1,1); e = 4: two (1,1); e = 6: one (1,1,1)    ((1,1): -4, (1,1,1): -9)
  r fixed 1-blocks (1)                                        (-1 each)
The excess blocks have size s_e (0, 2, 4, 3) and cost s_e + e, so the total is
K7 + 2 K6 - (s_e + e) - r = delta(N) with 7 K7 + 6 K6 + s_e + r = N.
Usage: python3 construct_all_n.py NMIN NMAX   (writes data/construction.jsonl)
"""
import json, sys

FIG = json.load(open("data/figure8_g7.json"))
A, B, C6 = FIG["A"], FIG["B"], FIG["C6"]
a, b = len(A), len(B)


def delta(n):
    return 1 if n % 2 else 2


EXCESS = {0: [], 2: [[1, 1]], 4: [[1, 1], [1, 1]], 6: [[1, 1, 1]]}


def plan(N):
    d = delta(N)
    K = -(-(N + d) // 8)
    e = 8 * K - N - d                      # in {0, 2, 4, 6}
    ex = EXCESS[e]
    se = sum(len(z) for z in ex)
    for K6 in range(0, min(3, K - 1) + 1):
        K7 = K - K6
        S = N - 7 * K7 - 6 * K6            # faces left for fixed blocks
        if S >= se:
            P = a * (K7 - 1) + b
            if P > N:
                return dict(K7=K7, K6=K6, e=e, excess=ex, r=S - se, P=P)
    return None


def counter_state(i, t):
    if t < a * i:
        return A[t % a]
    if t < a * i + b:
        return B[t - a * i]
    return A[(t - a * i - b) % a]


def build(N):
    p = plan(N)
    if p is None:
        return None, None
    K7, K6, ex, r, P = p["K7"], p["K6"], p["excess"], p["r"], p["P"]
    cyc = []
    for t in range(P):
        blocks = [counter_state(i, t) for i in range(K7)] + [C6[t % 3]] * K6 + ex + [[1]] * r
        die, off = [], 0
        for blk in blocks:
            die += [v + off for v in blk]
            off += len(blk)
        assert off == N
        cyc.append(tuple(sorted(die)))
    return p, cyc


if __name__ == "__main__":
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    with open("data/construction.jsonl", "w") as f:
        for N in range(lo, hi + 1):
            p = plan(N)
            f.write(json.dumps({"N": N, "plan": p}) + "\n")
    fails = [N for N in range(lo, hi + 1) if plan(N) is None]
    print("no plan for N in", fails)
