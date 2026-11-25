import numpy as np, networkx as nx, sys, json
sys.path.insert(0,'.')
from validate import F, blowup
from test_gen import build
def info(name,G):
    L=nx.laplacian_matrix(G).toarray().astype(float); w=np.linalg.eigvalsh(L)
    mus=sorted(set(np.round(w[w>1e-8],8)))
    X=np.linalg.pinv(L); dg=np.diag(X); X2=X@X
    deg=np.array([d for _,d in G.degree()])
    e1=sum(mus); 
    pairs=[(a,b) for a in set(deg) for b in set(deg) if a<b]
    print(f"{name}: n={len(deg)} m={len(mus)} mus={[round(x,4) for x in mus]} RR={np.ptp(dg)<1e-9} degs={sorted(set(deg))}",
          " e1-1-du-dv:",[round(e1-1-a-b,3) for a,b in pairs])
C6=nx.cycle_graph(6)
info('F(4)',F(4)); info('F(6)',F(6))
info('C6(2,4)',blowup(list(C6.edges()),[(2,'I'),(4,'K'),(4,'K'),(2,'I'),(4,'K'),(4,'K')]))
info('C6(3,9)',blowup(list(C6.edges()),[(3,'I'),(9,'K'),(9,'K'),(3,'I'),(9,'K'),(9,'K')]))
info('C6[PM4,CP8]',build(C6,['PM','CP','CP','PM','CP','CP'],[4,8,8,4,8,8]))
