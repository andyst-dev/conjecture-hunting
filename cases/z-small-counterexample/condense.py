"""Valuative invariants (KL polynomial P, Z-polynomial) of a binary matroid from its cyclic flats.

Pipeline:
  1. cycflats (C) enumerates all cyclic flats with rank and size;
  2. partition refinement: start from blocks = (size, rank) labels and split until, for all blocks
     B, C, every Y in C contains the same number A(B, C) of members of B (Eberhardt's condensed
     configuration, arXiv:2609.20903 Sec. 3);
  3. Theorem 3.2: Phi(M) = sum over block chains 0 = B0 < ... < Bk = 1 of prod(-(A^-1)_{B_{i-1} B_i}) *
     Phi(nested matroid with chain labels), computed by DP over chains grouped by label sequence.
Usage as a module: kl_z_binary(columns, n, k) -> dict.
"""
import subprocess, os
from fractions import Fraction
import numpy as np
from zfast import nested_kl_z, padd

HERE = os.path.dirname(os.path.abspath(__file__))


def cyclic_flats_binary(cols, k):
    exe = os.path.join(HERE, "cycflats")
    if not os.path.exists(exe):
        subprocess.check_call(["cc", "-O2", "-o", exe, os.path.join(HERE, "cycflats.c")])
    inp = f"{len(cols)} {k}\n" + " ".join(str(c) for c in cols) + "\n"
    res = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True)
    out = [tuple(map(int, line.split())) for line in res.stdout.split("\n") if line.strip()]
    return out, res.stderr.strip()


def condense(flats):
    """flats: list of (mask, rank, size). Returns blocks (list of index lists), labels, A (ints)."""
    N = len(flats)
    masks = np.array([f[0] for f in flats], dtype=np.uint64)
    lab = [(f[2], f[1]) for f in flats]
    # block id from labels
    keys = sorted(set(lab), key=lambda x: (x[1], x[0]))
    kpos = {l: i for i, l in enumerate(keys)}
    bid = np.array([kpos[l] for l in lab], dtype=np.int64)
    order = np.argsort([f[1] for f in flats], kind="stable")
    while True:
        nb = int(bid.max()) + 1
        sig = []
        for y in range(N):
            Y = masks[y]
            sub = (masks & Y) == masks
            sub[y] = False
            sig.append((int(bid[y]), tuple(np.bincount(bid[sub], minlength=nb).tolist())))
        first = {}
        for y, sg in enumerate(sig):
            first.setdefault(sg, y)
        newkeys = sorted(first, key=lambda sg: (flats[first[sg]][1], sg))
        pos = {sg: i for i, sg in enumerate(newkeys)}
        new_bid = np.array([pos[sg] for sg in sig], dtype=np.int64)
        if len(newkeys) == nb:
            break
        bid = new_bid
    nb = int(bid.max()) + 1
    members = [np.nonzero(bid == b)[0] for b in range(nb)]
    # sort blocks by rank
    brank = [flats[m[0]][1] for m in members]
    perm = sorted(range(nb), key=lambda b: (brank[b], flats[members[b][0]][2]))
    members = [members[b] for b in perm]
    inv = {b: i for i, b in enumerate(perm)}
    bid = np.array([inv[b] for b in bid])
    A = [[0] * nb for _ in range(nb)]
    for c in range(nb):
        y = members[c][0]
        Y = masks[y]
        sub = (masks & Y) == masks
        cnt = np.bincount(bid[sub], minlength=nb)
        for b in range(nb):
            A[b][c] = int(cnt[b])
    labels = [(flats[m[0]][2], flats[m[0]][1]) for m in members]
    return members, labels, A


def unitri_inverse(A):
    n = len(A)
    Inv = [[Fraction(0)] * n for _ in range(n)]
    for j in range(n):
        Inv[j][j] = Fraction(1)
        for i in range(j - 1, -1, -1):
            s = sum(Fraction(A[i][l]) * Inv[l][j] for l in range(i + 1, j + 1))
            Inv[i][j] = -s
    return Inv


class TooExpensive(Exception):
    pass


def chain_cost(seq):
    import math
    return math.prod(seq[i][0] - seq[i - 1][0] + 1 for i in range(1, len(seq)))


def expansion(labels, A, max_cost=None):
    """Theorem 3.2 with DP over chains; returns (P, Z, number of label sequences).
    max_cost: refuse (TooExpensive) if some nested chain has prod(layer size + 1) > max_cost."""
    m = len(labels)
    Inv = unitri_inverse(A)
    assert labels[0] == (0, 0), "matroid has loops"
    # dp[c] = {label sequence (tuple of labels incl. labels[c]) : coefficient}
    dp = [dict() for _ in range(m)]
    dp[0][(labels[0],)] = Fraction(1)
    for c in range(1, m):
        acc = {}
        for b in range(c):
            w = -Inv[b][c]
            if w == 0 or not dp[b]: continue
            for seq, coef in dp[b].items():
                key = seq + (labels[c],)
                acc[key] = acc.get(key, 0) + coef * w
        dp[c] = {s: v for s, v in acc.items() if v != 0}
    if max_cost is not None:
        worst = max(chain_cost(s) for s in dp[m - 1])
        if worst > max_cost:
            raise TooExpensive(f"nested chain cost {worst} > {max_cost}")
    P, Z = [], []
    cache = {}
    for seq, coef in dp[m - 1].items():
        assert coef.denominator == 1
        if seq not in cache: cache[seq] = nested_kl_z(list(seq))
        Pn, Zn = cache[seq]
        padd(P, Pn, mult=int(coef)); padd(Z, Zn, mult=int(coef))
    while len(P) > 1 and P[-1] == 0: P.pop()
    return P, Z, len(dp[m - 1])


def kl_z_binary(cols, k, verbose=True):
    flats, info = cyclic_flats_binary(cols, k)
    if verbose: print(info, flush=True)
    members, labels, A = condense(flats)
    top = labels[-1]
    n = len(cols)
    assert top[0] == n, f"matroid has coloops (top cyclic flat has size {top[0]} < {n})"
    P, Z, nseq = expansion(labels, A)
    return {"cyclic_flats": len(flats), "blocks": len(labels), "labels": labels, "A": A,
            "label_sequences": nseq, "P": P, "Z": Z}
