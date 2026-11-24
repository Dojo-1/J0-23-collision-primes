#!/usr/bin/env python3
"""Independent point counts on the manuscript's genus-two curve."""
from verify import eta_coeffs,theta,convolution
import json
P=[3,5,7,11,13,19,29,43,97,109,223]
N=max(P);B=eta_coeffs(N);u=convolution(theta(2,1,3,N),B,N);v=convolution(B,B,N)
rows=[]
for p in P:
    qr={i*i%p for i in range(1,p)}
    chi=lambda z:0 if z%p==0 else (1 if z%p in qr else -1)
    d=next(z for z in range(1,p) if chi(z)==-1)
    # Completed square F=(x^3+x+1)^2-4*x^5. Monic, two points at infinity.
    n1=p+2+sum(chi((x**3+x+1)**2-4*x**5) for x in range(p))
    def mul(a,b):return ((a[0]*b[0]+d*a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0])%p)
    def add(a,b):return ((a[0]+b[0])%p,(a[1]+b[1])%p)
    def sub(a,b):return ((a[0]-b[0])%p,(a[1]-b[1])%p)
    total=0
    for a in range(p):
        for b in range(p):
            x=(a,b);x2=mul(x,x);x3=mul(x,x2);x5=mul(x2,x3)
            h=add(add(x3,x),(1,0));f=sub(mul(h,h),mul((4,0),x5))
            total+=chi(f[0]*f[0]-d*f[1]*f[1])
    n2=p*p+2+total
    T=p+1-n1;T2=p*p+1-n2
    assert T==2*u[p]+v[p]
    assert 2*(T2+4*p)==T*T+5*v[p]*v[p]
    rows.append(dict(p=p,N1=n1,N2=n2,u=u[p],v=v[p],T=T,T2=T2))
print(json.dumps(rows,indent=2));open('pointcounts.json','w').write(json.dumps(rows,indent=2)+'\n')
