import numpy as np
from fractions import Fraction
def comp(l,k,a_ok,e_ok,b_ok,c_ok):
    m=l**k; r=np.arange(m)
    b=r[b_ok(r)]; c=r[c_ok(r)]
    M=np.bincount((np.outer(b,c)%m).ravel(),minlength=m).astype(np.int64)
    a=r[a_ok(r)]; e=r[e_ok(r)]
    A,E=np.meshgrid(a,e,indexing='ij'); prod=(A*E)%m; tr=(A+E)%m
    N=np.zeros((m,m),dtype=np.int64)  # N[t,D]
    for D in range(m):
        if D%l==0: continue
        cnt=M[(prod-D)%m]
        N[:,D]=np.bincount(tr.ravel(),weights=cnt.ravel(),minlength=m).astype(np.int64)
    return N
def F(l,k,g,h):
    Ng=comp(l,k,*g); Nh=comp(l,k,*h)
    num=int((Ng*Nh).sum()); den=int((Ng.sum(0)*Nh.sum(0)).sum())
    return Fraction(l**k*num,den)
T=lambda r: np.ones_like(r,bool)
full=(T,T,T,T)
l=11
diag=(lambda r:r%l==1, T, lambda r:r%l==0, lambda r:r%l==0)
borel=(lambda r:r%l==1, T, T, lambda r:r%l==0)
for k in (1,2,3):
    print(k, 'diag-model %.7f'%float(F(l,k,diag,full)), 'borel-model %.7f'%float(F(l,k,borel,full)), 'full %.7f'%float(F(l,k,full,full)),flush=True)
