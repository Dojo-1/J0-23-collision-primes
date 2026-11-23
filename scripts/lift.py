import numpy as np, pickle, itertools
def mulv(x1,y1,x2,y2,m): return ((x1*x2+y1*y2)%m,(x1*y2+x2*y1+y1*y2)%m)
def stats(base, l, k, cond_mod=None):
    """base: list of 2x2 matrices over O/l (tuples (x,y)); full preimage in GL2(O/l^k) with det in Z/l^k.
    returns l^k*P(tr in Z), optionally per det class mod cond_mod."""
    m=l**k; L=l**(k-1)
    offs=np.array(list(itertools.product(range(L),repeat=8)))  # 8 coordinates (x,y per entry) lifts
    out={}
    num=0;den=0;numc={};denc={}
    for A in base:
        e=[A[0][0],A[0][1],A[1][0],A[1][1]]
        X=[ (e[j][0]+l*offs[:,2*j])%m for j in range(4)]
        Y=[ (e[j][1]+l*offs[:,2*j+1])%m for j in range(4)]
        ax,ay=mulv(X[0],Y[0],X[3],Y[3],m); bx,by=mulv(X[1],Y[1],X[2],Y[2],m)
        dx=(ax-bx)%m; dy=(ay-by)%m
        ok=(dy==0); tro=((Y[0]+Y[3])%m==0)
        num+=int((ok&tro).sum()); den+=int(ok.sum())
        if cond_mod:
            for r in range(cond_mod):
                s=ok&(dx%cond_mod==r)
                if s.any():
                    numc[r]=numc.get(r,0)+int((s&tro).sum()); denc[r]=denc.get(r,0)+int(s.sum())
    res=m*num/den
    return res,{r:numc[r]/denc[r] for r in denc}
G3=pickle.load(open('G3.pkl','rb'))
for k in (1,2):
    r,c=stats(G3,3,k,cond_mod=3); print('l=3 icosahedral model k=%d: F=%.5f  P(tr in Z|det mod 3)='%(k,r),{a:round(b,4) for a,b in c.items()})
# l=2: Gamma2 = matrices over F4 with entries in F2 (y=0) and invertible
F2m=[((a,b),(c,d)) for a in [(0,0),(1,0)] for b in [(0,0),(1,0)] for c in [(0,0),(1,0)] for d in [(0,0),(1,0)] if (a[0]*d[0]-b[0]*c[0])%2==1]
for k in (2,3):
    r,c=stats(F2m,2,k,cond_mod=4); print('l=2 Gamma2 k=%d: F=%.5f P(2^k|v | det mod 4)='%(k,r),{a:round(b,4) for a,b in c.items()})
