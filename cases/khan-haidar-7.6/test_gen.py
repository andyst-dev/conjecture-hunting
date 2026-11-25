import numpy as np, networkx as nx, random, itertools, sys
sys.path.insert(0,'.')
from search_gen import tables, quotient
from core import lplus_exact
def internal(t,n):
    if t=='K': return nx.complete_graph(n)
    if t=='I': return nx.empty_graph(n)
    if t=='CP': return nx.complement(nx.from_edgelist([(2*i,2*i+1) for i in range(n//2)]))
    if t=='PM': return nx.from_edgelist([(2*i,2*i+1) for i in range(n//2)])
    if t=='KB': return nx.complete_bipartite_graph(n//2,n//2)
    if t=='C': return nx.cycle_graph(n)
    if t=='CC': return nx.complement(nx.cycle_graph(n))
def build(H,ts,ns):
    G=nx.Graph(); M=[]
    for i,(t,n) in enumerate(zip(ts,ns)):
        W=internal(t,n); W.add_nodes_from(range(n))
        G.add_nodes_from((i,v) for v in range(n)); G.add_edges_from(((i,a),(i,b)) for a,b in W.edges()); M.append([(i,v) for v in range(n)])
    for i,j in H.edges(): G.add_edges_from((u,v) for u in M[i] for v in M[j])
    return G
if __name__=='__main__':
    random.seed(3); err=0
    types=['K','I','CP','PM','KB','C','CC']
    W,W2,R,V=tables(types,8,40)
    for _ in range(40):
        k=random.randint(2,5)
        while True:
            H=nx.gnp_random_graph(k,.6,seed=random.randint(0,10**6))
            if nx.is_connected(H): break
        ts=[random.choice(types) for _ in range(k)]
        ns=[random.choice([n for n in range(2,9) if V[t][n]]) for t in ts]
        A=nx.to_numpy_array(H,nodelist=range(k)); N=np.array([ns])
        e,qd,qd2=quotient(A,N)
        G=build(H,ts,ns); nodes,X=lplus_exact(G)
        Xf=np.array([[float(X[i,j]) for j in range(len(nodes))] for i in range(len(nodes))]); X2=Xf@Xf
        for i in range(k):
            ei=int(round(e[0,i])); c=qd[0,i]+W[ts[i]][ns[i],ei]; c2=qd2[0,i]+W2[ts[i]][ns[i],ei]
            for v in range(ns[i]):
                j=nodes.index((i,v)); err=max(err,abs(Xf[j,j]-c),abs(X2[j,j]-c2),abs(G.degree((i,v))-(R[ts[i]][ns[i]]+e[0,i])))
    print('gen formula max err',err)
