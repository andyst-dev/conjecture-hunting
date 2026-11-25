import itertools, sys; sys.path.insert(0,'.')
from validate import blowup, report
from core import analyze
P8=[(i,i+1) for i in range(7)]; sizes=(9,6,3,2,2,3,6,9)
for ts in itertools.product('KI', repeat=8):
    G=blowup(P8, list(zip(sizes,ts)))
    r=analyze(G)
    if r['is_rr']:
        print(''.join(ts), set(r['diag']), sorted(set(r['deg']))); report('G40',G)
