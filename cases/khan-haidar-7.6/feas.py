"""Can Delta = s1 - 1 - a - b <= 0 for an RR graph with m=4, degrees a<b?
Necessary vertex-level conditions (both degree classes w in {a,b}):
  p_w(k) = mu_k/prod_{j!=k}(mu_k-mu_j) * [d^2+d - s1^(k) d + s2^(k) s - s3^(k) c] >= 0
  N_w = s1(d^2+d) - s2 d + s3 s - s4 c - d^3 - 2d^2  with  d <= N_w <= d*b
  mu4 >= b+1, mu1 <= a (alg. connectivity <= min degree), mu4 <= n, s = 1-1/n.
Maximise -Delta by random multistart SLSQP over (mu1..mu4, a, b, c, n) continuous."""
import numpy as np, itertools
from scipy.optimize import minimize
rng=np.random.default_rng(0)
def esym(m):
    s1=sum(m); s2=sum(x*y for x,y in itertools.combinations(m,2)); s3=sum(x*y*z for x,y,z in itertools.combinations(m,3)); s4=np.prod(m)
    return s1,s2,s3,s4
def cons_vals(x):
    mu=np.sort(x[:4]); a,b,c,n=x[4:]; s=1-1/n
    s1,s2,s3,s4=esym(mu); out=[]
    for d in (a,b):
        for k in range(4):
            o=[mu[j] for j in range(4) if j!=k]; t1,t2,t3=sum(o),o[0]*o[1]+o[0]*o[2]+o[1]*o[2],o[0]*o[1]*o[2]
            pref=mu[k]/np.prod([mu[k]-mu[j] for j in range(4) if j!=k])
            out.append(pref*(d*d+d-t1*d+t2*s-t3*c))
        N=s1*(d*d+d)-s2*d+s3*s-s4*c-d**3-2*d*d
        out+= [N-d, d*b-N]
    out+=[mu[3]-b-1, a-mu[0], n-mu[3], b-a-1, a-1, n-b-1]
    return np.array(out)
best=None
for t in range(4000):
    n=rng.uniform(8,200); b=rng.uniform(2,n-1); a=rng.uniform(1,b-1)
    x0=np.array(sorted(rng.uniform(0.1,n,4))+[a,b,rng.uniform(1e-3,1),n])
    r=minimize(lambda x: (np.sum(np.sort(x[:4]))-1-x[4]-x[5]), x0, constraints=[{'type':'ineq','fun':cons_vals}],
               bounds=[(1e-3,None)]*4+[(1,None),(1,None),(1e-6,None),(4,None)], method='SLSQP', options={'maxiter':500})
    if r.success and np.all(cons_vals(r.x)>-1e-7):
        if best is None or r.fun<best.fun: best=r
print('min Delta found:',best.fun if best else None)
if best is not None:
    x=best.x; print('mu',np.sort(x[:4]),'a,b',x[4],x[5],'c',x[6],'n',x[7])
