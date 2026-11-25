"""Cross-checks for the CDM pipeline (independent brute force from the definitions, no reductions).
1. The C search (`cdm --raw`) vs brute force over ALL matchings, on random graphs where CDM can fail.
   Here H is arbitrary and G = complement(H) is arbitrary too, so the search's completeness is tested
   on graphs with and without a CDM.
2. CDM by brute force on every connected graph G with alpha(G) <= 2 for n <= 9
   (complements of all triangle-free graphs from geng -t). No MTF reduction is used.
3. Reduction check: for n <= 10, graphs with alpha<=2, connected, without a dominating edge are
   exactly the complements of MTF graphs other than K_{a,b}.
"""
import itertools, random, subprocess, networkx as nx

GENG = "/tmp/nautydl/pynauty-2.8.8.1/src/nauty2_8_8/geng"


def brute_cdm(G):
    V = set(G.nodes()); E = [tuple(e) for e in G.edges()]
    adjset = lambda A, B: any(G.has_edge(a, b) for a in A for b in B)
    # enumerate matchings by increasing size
    def matchings(start, used, cur):
        yield cur
        for i in range(start, len(E)):
            a, b = E[i]
            if a in used or b in used: continue
            yield from matchings(i + 1, used | {a, b}, cur + [E[i]])
    for M in matchings(0, frozenset(), []):
        if not M: continue
        if not all(adjset(e, f) for e, f in itertools.combinations(M, 2)): continue
        VM = {v for e in M for v in e}
        if all(adjset({z}, e) for e in M for z in V - VM):
            return True
    return False


def g6(H):
    return nx.to_graph6_bytes(H, header=False).decode().strip()


def run_c(graphs):
    inp = "\n".join(g6(H) for H in graphs) + "\n"
    out = subprocess.run(["./cdm_test", "--raw"], input=inp, capture_output=True, text=True).stdout.split("\n")
    return [int(l.split()[1]) for l in out if l.strip() and not l.startswith(("n=", "size"))]


if __name__ == "__main__":
    rng = random.Random(1)
    # 1. random graphs
    Hs = [nx.gnp_random_graph(rng.randint(4, 9), rng.choice([0.3, 0.5, 0.7]), seed=rng.randint(0, 10**9)) for _ in range(600)]
    Hs = [nx.convert_node_labels_to_integers(H) for H in Hs]
    c_res = run_c(Hs)
    b_res = [int(brute_cdm(nx.complement(H))) for H in Hs]
    agree = sum(a == b for a, b in zip(c_res, b_res))
    print(f"1. C search vs brute force on 600 random graphs: agree {agree}/600 "
          f"(graphs without CDM among them: {b_res.count(0)})")
    # 2 & 3
    for n in range(3, 11):
        lines = subprocess.run([GENG, "-tq", str(n)], capture_output=True, text=True).stdout.split()
        nconn = fail = no_dom = mtf_nonkab = 0
        for l in lines:
            H = nx.from_graph6_bytes(l.encode()); G = nx.complement(H)
            if not nx.is_connected(G): continue
            nconn += 1
            dom = any(set(G[x]) | set(G[y]) | {x, y} == set(G) for x, y in G.edges())
            if not dom: no_dom += 1
            is_mtf = nx.is_connected(H) and nx.diameter(H) <= 2
            kab = nx.is_bipartite(H) and H.number_of_edges() == (lambda a: a[0] * a[1])(
                [len(s) for s in nx.bipartite.sets(H)]) if nx.is_connected(H) else False
            if is_mtf and not kab: mtf_nonkab += 1
            if n <= 9 and not brute_cdm(G): fail += 1
        msg = f"   n={n}: connected alpha<=2 graphs={nconn}, without dominating edge={no_dom}, MTF non-K_ab={mtf_nonkab}"
        if n <= 9: msg += f", CDM failures (brute force)={fail}"
        print(("2/3." if n == 3 else "    ") + msg)
