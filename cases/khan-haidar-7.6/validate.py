import networkx as nx, itertools, sys
sys.path.insert(0, '.')
from core import analyze, check_conjecture

def F(p):
    h = p//2
    G = nx.Graph()
    Q = {b: [('q',b,i) for i in range(p)] for b in (0,1)}
    X = {b: [('x',b,c,t) for c in (0,1) for t in range(h)] for b in (0,1)}
    for b in (0,1):
        G.add_edges_from(itertools.combinations(Q[b],2))
        G.add_edges_from((q,x) for q in Q[b] for x in X[b])
    for c in (0,1):
        Y = [('x',b,c,t) for b in (0,1) for t in range(h)]
        for u,v in itertools.combinations(Y,2):
            if u[3]==v[3] and u[1]!=v[1]: continue
            G.add_edge(u,v)
    return G

def blowup(base_edges, mods):
    """mods: list of (size, 'K' or 'I')"""
    G = nx.Graph(); M=[]
    for i,(s,t) in enumerate(mods):
        vs=[(i,j) for j in range(s)]; M.append(vs); G.add_nodes_from(vs)
        if t=='K': G.add_edges_from(itertools.combinations(vs,2))
    for i,j in base_edges:
        G.add_edges_from((u,v) for u in M[i] for v in M[j])
    return G

def report(name, G):
    r = analyze(G)
    ok, cls, const = check_conjecture(r)
    print(name, 'n=',G.number_of_nodes(),'m=',G.number_of_edges(),'RR=',r['is_rr'],
          'diagL+=',set(r['diag']) if r['is_rr'] else '-', 'conj holds=',ok,'S2 const per deg=',const)
    print('   S2 by degree:', {d:[str(x) for x in v] for d,v in cls.items()})

if __name__=='__main__':
  report('F(4)', F(4))
  report('F(6)', F(6))
  C6 = [(i,(i+1)%6) for i in range(6)]
  report('C6(2,4)', blowup(C6, [(2,'I'),(4,'K'),(4,'K'),(2,'I'),(4,'K'),(4,'K')]))
  report('C6(3,9)', blowup(C6, [(3,'I'),(9,'K'),(9,'K'),(3,'I'),(9,'K'),(9,'K')]))
  P8 = [(i,i+1) for i in range(7)]
  report('Gamma40 guess', blowup(P8, [(9,'K')]+[(s,'I') for s in (6,3,2,2,3,6,9)]))
