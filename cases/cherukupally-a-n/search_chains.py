"""Exhaustive search for n <= 10^E lying in two chains C_z (z >= 3), i.e. |A(n)| >= 4.
By the chain characterisation (test_chain.py) this is equivalent to a counterexample to
Conjecture 3.4 with n <= 10^E.  Values with i >= 4 are stored; i = 3 values are streamed;
the i = 2 values z^2-1 are tested by a perfect-square check.  Exact integers throughout.
usage: python3 search_chains.py E"""
import sys, time
from math import isqrt
X=10**int(sys.argv[1])
t=time.time()
hi={}  # values with i>=4
z=3
while z**4 <= 2*X:
    g0,g1=z,z*z-1; i=2
    while True:
        g0,g1=g1,z*g1-g0; i+=1
        if g1>X: break
        if i>=4: hi.setdefault(g1,[]).append((i,z))
    z+=1
hits=[]
for v,l in hi.items():
    if len(l)>1: hits.append((v,l))
    r=isqrt(v+1)
    if r*r==v+1: hits.append((v,l+[(2,r)]))
z=3
while True:
    v=z*z*z-2*z
    if v>X: break
    l=[(3,z)]
    if v in hi: l+=hi[v]
    r=isqrt(v+1)
    if r*r==v+1: l.append((2,r))
    if len(l)>1: hits.append((v,l))
    z+=1
print("X=%g"%X,"hits",hits,"time %.1fs"%(time.time()-t))
