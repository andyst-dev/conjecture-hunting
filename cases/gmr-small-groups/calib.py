import random, math, time
from groups import build_all
from gmr import evaluate
from search_words import alphabet, climb
G=build_all()['SG(243,3)']; L=alphabet(G,1,2)
for T0 in (0.02, 0.1, 0.3):
    rng=random.Random(1); t=time.time(); best=0
    for r in range(6):
        R=climb(G,L,10,rng,6000,T0=T0)[0]; best=max(best,R)
    print(T0, float(best), round(time.time()-t,1),'s')
