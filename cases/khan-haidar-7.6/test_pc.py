import numpy as np, networkx as nx, random, sys
sys.path.insert(0,'.')
import search_polycirc as P
from search_polycirc import batch_eval, charsum, sym_subsets, subsets
from core import lplus_exact
def add(a,s):
    E=P.ELEMS; e=tuple((x+y)%d for x,y,d in zip(E[a],E[s],P.GROUP)); return E.index(e)
def build(k,N,Sd,So):
    G=nx.Graph(); G.add_nodes_from((a,i) for a in range(N) for i in range(k))
    for i in range(k):
        for a in range(N):
            for s in Sd[i]: G.add_edge((a,i),(add(a,s),i))
    for (i,j),S in So.items():
        for a in range(N):
            for s in S: G.add_edge((a,i),(add(a,s),j))
    return G
random.seed(5); err=0; tested=0
for _ in range(200):
    k=random.choice([2,3]); dims=random.choice([[4],[5],[6],[7],[2,2],[2,4],[2,3],[3,3]]); N=P.set_group(dims)
    SS=sym_subsets(N); AS=subsets(N)
    Sd=[random.choice(SS) for _ in range(k)]
    So={(i,j):random.choice(AS) for i in range(k) for j in range(i+1,k)}
    G=build(k,N,Sd,So)
    if not nx.is_connected(G): continue
    c,c2,deg,valid=batch_eval(k,N,[charsum(S,N)[None] for S in Sd],[np.array([len(S)]) for S in Sd],
        {p:charsum(S,N)[None] for p,S in So.items()},{p:np.array([len(S)]) for p,S in So.items()})
    assert valid[0]
    nodes,X=lplus_exact(G); Xf=np.array([[float(X[a,b]) for b in range(len(nodes))] for a in range(len(nodes))]); X2=Xf@Xf
    for i in range(k):
        for a in range(N):
            t=nodes.index((a,i)); err=max(err,abs(Xf[t,t]-c[0,i]),abs(X2[t,t]-c2[0,i]),abs(G.degree((a,i))-deg[0,i]))
    tested+=1
print('tested',tested,'max err',err)
