"""Recompute expected counts independently of the author coefficient dataset."""
from math import isqrt, fsum, sqrt
import json
N=60_000_000
s=bytearray(b'\1')*(N+1); s[:2]=b'\0\0'
for p in range(2,isqrt(N)+1):
    if s[p]:s[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
s[2]=s[23]=0
K1=1.916948;K2=2.432194
out=[]
for a,b in [(0,100_000),(0,1_000_000),(0,20_000_000),(20_000_000,40_000_000),(40_000_000,60_000_000),(20_000_000,60_000_000)]:
    w=fsum(1/sqrt(p) for p in range(max(a+1,3),b+1,2) if s[p]) if a%2==0 else fsum(1/sqrt(p) for p in range(a+2,b+1,2) if s[p])
    out.append(dict(a=a,b=b,weight=w,M1=K1*w,M2=K2*w,corrected_M2=K2*w*(213/160)/1.332))
r=dict(primes=s.count(1),expected=out)
print(json.dumps(r,indent=2));open('predictions.json','w').write(json.dumps(r,indent=2)+'\n')
