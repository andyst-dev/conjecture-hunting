"""Reproduce the stated values of arXiv:2610.05663 (and HT(P_6)=223 < 224 = HT(L(6,4)) from arXiv:2611.34295 as quoted in its certificate repo)."""
import random
from fractions import Fraction
import networkx as nx
from asym import A_and_HT_exact, A_tree, A_path, resistance_exact

ok = True
def chk(name, paper, ours):
    global ok
    good = Fraction(paper) == Fraction(ours); ok &= good
    print(f"{name:45s} paper={str(paper):>10} ours={str(ours):>10} {'OK' if good else 'MISMATCH'}")

for n in range(2, 11):
    A, HT = A_and_HT_exact(nx.path_graph(n))
    chk(f"A(P_{n}) Thm 6.1", A_path(n), A)
    hp = Fraction(n * (n - 1) * (n * n - 1), 6) + A_path(n)
    chk(f"HT(P_{n}) Cor 6.1", hp, HT)
for n in range(3, 10):
    A, HT = A_and_HT_exact(nx.star_graph(n - 1))
    chk(f"HT(S_{n}) Rem 5.2", (n - 1) * (n * n - n - 1), HT)
    chk(f"A(S_{n}) = (n-1)(n-2) Thm 5.1 equality", (n - 1) * (n - 2), A)
for p in range(1, 5):
    for q in range(p, 6):
        G = nx.complete_bipartite_graph(p, q)
        A, HT = A_and_HT_exact(G)
        chk(f"A(K_{p},{q}) = pq|p-q|", p * q * abs(p - q), A)
        R = resistance_exact(G)
        Kf = sum(R[u, v] for u in G for v in G if u < v)
        chk(f"Kf(K_{p},{q})", Fraction(p * (p - 1), q) + Fraction(q * (q - 1), p) + p + q - 1, Kf)
        chk(f"HT = m Kf + A (Thm 3.1), K_{p},{q}", G.number_of_edges() * Kf + A, HT)
def lollipop(n, k):
    G = nx.complete_graph(k); prev = 0
    for v in range(k, n): G.add_edge(prev, v); prev = v
    return G
chk("HT(P_6) [2611.34295]", 223, A_and_HT_exact(nx.path_graph(6))[1])
chk("HT(L(6,4)) [2611.34295]", 224, A_and_HT_exact(lollipop(6, 4))[1])
# cross-check the tree fast path against the exact hitting-time computation
rng = random.Random(5); bad = 0
for t in range(300):
    n = rng.randint(3, 11)
    T = nx.random_labeled_tree(n, seed=rng.randint(0, 10**9)) if hasattr(nx, "random_labeled_tree") else nx.random_tree(n, seed=rng.randint(0, 10**9))
    A, HT = A_and_HT_exact(T); a2, W = A_tree(T)
    if A != a2 or HT != (n - 1) * W + A: bad += 1
chk("tree fast path vs exact hitting times (#bad)", 0, bad)
print("ALL OK" if ok else "SOME MISMATCH")
