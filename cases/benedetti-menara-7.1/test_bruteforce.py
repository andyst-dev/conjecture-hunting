"""Independent brute force of the maximum number of layers (no SAT), compared with the SAT encoding.
For small (m, r): enumerate every family of r-subsets and every surjective layering, check the
interval property and (optionally) flagness directly from the definitions."""
import itertools, sys
from sat import solve

def is_flag(fam, m, r):
    adj = set()
    for T in fam:
        for u, v in itertools.combinations(T, 2): adj.add((u, v)); adj.add((v, u))
    for T in fam:
        for w in range(m):
            if w not in T and all((w, u) in adj for u in T): return False
    return True

def layered_ok(fam, lab, r):
    for k in range(1, r):
        S_all = {S for T in fam for S in itertools.combinations(T, k)}
        for S in S_all:
            idx = sorted({lab[i] for i, T in enumerate(fam) if set(S) <= set(T)})
            if idx[-1] - idx[0] + 1 != len(idx): return False
    return True

def brute_max(m, r, flag):
    allT = list(itertools.combinations(range(m), r)); best = 0
    for size in range(1, len(allT) + 1):
        if size <= best: continue
        for fam in itertools.combinations(allT, size):
            if flag and not is_flag(fam, m, r): continue
            for L in range(size, best, -1):
                found = False
                for lab in itertools.product(range(L), repeat=size):
                    if len(set(lab)) != L or lab[0] != 0 and False: continue
                    if layered_ok(fam, lab, r): found = True; break
                if found: best = max(best, L); break
    return best

def sat_max(m, r, flag):
    L = 1
    while solve(m, r, L + 1, flag) is not None: L += 1
    return L

for (m, r, flags) in [(3, 2, (True, False)), (4, 2, (True, False)), (5, 2, (True,)), (4, 3, (True, False)), (5, 3, (True,))]:
    for flag in flags:
        b, s = brute_max(m, r, flag), sat_max(m, r, flag)
        print(f"m={m} r={r} flag={flag}: brute force max layers={b}, SAT max layers={s}, bound m-r+1={m-r+1}", "OK" if b == s else "MISMATCH")
