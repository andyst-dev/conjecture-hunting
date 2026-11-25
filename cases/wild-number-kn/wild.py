"""Wild number of edge-coloured K_n (arXiv:2610.22142, Conjecture 6.1).

gamma: proper edge colouring of K_n with chi'(K_n) colours (n-1 for even n, n for odd n).
W subset of E is colour-connecting if G_i + W is connected for every colour class G_i.
wld = min |W|. Conjecture 6.1: for EVERY minimum proper colouring, wld = floor(4*ceil(n/2)/3).
(The paper proves wld >= this value for every rho, and equality for n <= 12.)

Exact decision 'wld <= k' by SAT. Connectivity of G_i + W uses BFS layers from vertex 0:
reach[i][v][d] -> reach[i][v][d-1] or (some neighbour u in G_i with reach[i][u][d-1])
                   or (some u with w_uv and reach[i][u][d-1]); every v must be reached at depth n-1.
Cardinality sum w <= k uses a sequential counter.
"""
import itertools, random
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153


def conj_value(n):
    return (4 * ((n + 1) // 2)) // 3


def gk_colouring(n):
    """standard rotational colouring: odd n: colour(i,j) = (i+j) mod n (n colours, proper);
    even n: GK_n (vertex inf = n-1, colour classes = near-1-factors of K_{n-1} + edge to inf)."""
    col = {}
    if n % 2 == 1:
        for i, j in itertools.combinations(range(n), 2):
            col[i, j] = (i + j) % n
    else:
        m = n - 1
        for i, j in itertools.combinations(range(m), 2):
            col[i, j] = (i + j) % m
        for i in range(m):
            col[i, m] = (2 * i) % m     # the colour missing at vertex i in K_m
    return col


def check_proper(col, n):
    for v in range(n):
        seen = [col[min(v, u), max(v, u)] for u in range(n) if u != v]
        if len(seen) != len(set(seen)): return False
    ncol = len(set(col.values()))
    return ncol == (n if n % 2 else n - 1)


def kempe_shuffle(col, n, steps, rng):
    """random Kempe switches: pick colours a,b and a vertex, swap a/b on its 2-coloured component.
    Preserves properness and the number of colours."""
    col = dict(col)
    colours = sorted(set(col.values()))
    def c(u, v): return col[min(u, v), max(u, v)]
    for _ in range(steps):
        a, b = rng.sample(colours, 2); v0 = rng.randrange(n)
        # walk the (a,b)-path/cycle through v0
        comp = set(); stack = [v0]; edges = set()
        while stack:
            v = stack.pop()
            if v in comp: continue
            comp.add(v)
            for u in range(n):
                if u != v and c(u, v) in (a, b):
                    edges.add((min(u, v), max(u, v))); stack.append(u)
        for e in edges:
            col[e] = b if col[e] == a else a
    return col


def wld_at_most(col, n, k, timeout=None):
    pool = IDPool(); cl = []
    E = list(itertools.combinations(range(n), 2))
    w = {e: pool.id(("w", e)) for e in E}
    W = lambda u, v: w[min(u, v), max(u, v)]
    colours = sorted(set(col.values()))
    cls = {c: set() for c in colours}
    for e in E: cls[col[e]].add(e)
    D = n - 1
    for c in colours:
        adjc = {v: set() for v in range(n)}
        for (u, v) in cls[c]: adjc[u].add(v); adjc[v].add(u)
        R = {(v, d): pool.id(("r", c, v, d)) for v in range(n) for d in range(D + 1)}
        cl.append([R[0, 0]])
        for v in range(1, n): cl.append([-R[v, 0]])
        for d in range(1, D + 1):
            for v in range(n):
                lits = [R[v, d - 1]] + [R[u, d - 1] for u in adjc[v]]
                for u in range(n):
                    if u == v or u in adjc[v]: continue
                    a = pool.id(("a", c, u, v, d))      # a -> (w_uv and R[u,d-1])
                    cl.append([-a, W(u, v)]); cl.append([-a, R[u, d - 1]])
                    lits.append(a)
                cl.append([-R[v, d]] + lits)
        for v in range(n): cl.append([R[v, D]])
    card = CardEnc.atmost(lits=list(w.values()), bound=k, vpool=pool, encoding=EncType.seqcounter)
    cl.extend(card.clauses)
    with Cadical153(bootstrap_with=cl) as s:
        ok = s.solve()
        if ok:
            mdl = set(l for l in s.get_model() if l > 0)
            return [e for e in E if w[e] in mdl]
        return None


def wld_exact(col, n, lo=None):
    """returns (wld, witness W). Starts at the proven lower bound conj_value(n) - 1 check."""
    k = conj_value(n)
    sol = wld_at_most(col, n, k)
    if sol is not None:
        below = wld_at_most(col, n, k - 1)
        if below is not None:
            return k - 1, below          # would contradict the paper's lower bound
        return k, sol
    # wld > conjectured value: find it
    while sol is None:
        k += 1; sol = wld_at_most(col, n, k)
    return k, sol


def verify_W(col, n, W):
    """independent check (union-find) that W is colour-connecting"""
    colours = set(col.values())
    for c in colours:
        parent = list(range(n))
        def f(x):
            while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for e, cc in col.items():
            if cc == c: parent[f(e[0])] = f(e[1])
        for (u, v) in W: parent[f(u)] = f(v)
        if len({f(x) for x in range(n)}) != 1: return False
    return True
