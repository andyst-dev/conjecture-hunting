"""Exhaustive check over all trees on n vertices (networkx.nonisomorphic_trees; counts vs OEIS A000055):
  Question (Sec. 7): A(T) <= A(P_n), equality only for P_n.
  Thm 5.1: A(T) >= (n-1)(n-2), equality iff star.   Thm 5.2: A(T) <= (n-2) W(T).
Integer arithmetic only (A(T) is half an integer sum; compared exactly as Fractions)."""
import sys, json, time
from fractions import Fraction
import networkx as nx
from asym import A_tree, A_path
A000055 = [1, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159, 7741, 19320, 48629, 123867, 317955, 823065, 2144505]
nmin, nmax = int(sys.argv[1]), int(sys.argv[2])
log = open("data/trees_summary.jsonl", "a")
for n in range(nmin, nmax + 1):
    t0 = time.time(); ap = A_path(n); cnt = 0; viol = []; best2 = None; th51 = th52 = 0
    for T in nx.nonisomorphic_trees(n):
        cnt += 1
        a, W = A_tree(T)
        degs = sorted(d for _, d in T.degree())
        is_path = degs[-1] <= 2
        is_star = degs[-1] == n - 1
        if not is_path:
            if a >= ap: viol.append(nx.to_graph6_bytes(T, header=False).decode().strip())
            if best2 is None or a > best2[0]: best2 = (a, nx.to_graph6_bytes(T, header=False).decode().strip())
        if a < (n - 1) * (n - 2) or (a == (n - 1) * (n - 2) and not is_star): th51 += 1
        if a > (n - 2) * W: th52 += 1
    assert cnt == A000055[n], (n, cnt)
    rec = {"n": n, "trees": cnt, "A(P_n)": str(ap), "max_A_nonpath": str(best2[0]) if best2 else None,
           "argmax_nonpath": best2[1] if best2 else None, "violations_question": viol,
           "violations_thm5.1": th51, "violations_thm5.2": th52, "secs": round(time.time() - t0, 1)}
    log.write(json.dumps(rec) + "\n"); log.flush(); print(json.dumps(rec), flush=True)
