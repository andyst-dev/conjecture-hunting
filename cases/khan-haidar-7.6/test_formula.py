import numpy as np, networkx as nx, random, sys
sys.path.insert(0,'.')
from search_blowup import eval_batch
from validate import blowup
from core import lplus_exact
random.seed(1)
maxerr=0
for trial in range(40):
    k=random.randint(2,6)
    while True:
        H=nx.gnp_random_graph(k,0.6,seed=random.randint(0,10**6))
        if nx.is_connected(H): break
    t=np.array([random.random()<0.5 for _ in range(k)])
    n=np.array([random.randint(1,5) for _ in range(k)])
    A=nx.to_numpy_array(H,nodelist=range(k))
    c,c2,deg=eval_batch(A,t,n[None,:])
    G=blowup(list(H.edges()),[(int(n[i]),'K' if t[i] else 'I') for i in range(k)])
    nodes,X=lplus_exact(G)
    Xf=np.array([[float(X[i,j]) for j in range(len(nodes))] for i in range(len(nodes))])
    X2=Xf@Xf
    for i in range(k):
        j=nodes.index((i,0))
        maxerr=max(maxerr,abs(Xf[j,j]-c[0,i]),abs(X2[j,j]-c2[0,i]),abs(G.degree((i,0))-deg[0,i]))
print('max abs error formula vs exact:',maxerr)
