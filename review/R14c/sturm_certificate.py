#!/usr/bin/env python3
"""Standalone, exact certificate for the mod-4 congruence. Python 3."""
from math import isqrt
from hashlib import sha256

N = 16 * 23**2
M = N * 3 // 2 * 24 // 23 // 6
assert (N, M) == (8464, 2208)
# Quadratic theta series, independently of the eta-product algorithm.
def theta(a,b,c):
    d=4*a*c-b*b; t=[0]*(M+1)
    for y in range(-isqrt(4*a*M//d),isqrt(4*a*M//d)+1):
        s=isqrt(4*a*M-d*y*y)
        for x in range(-((s+b*y)//(2*a)),(s-b*y)//(2*a)+1):
            t[a*x*x+b*x*y+c*y*y]+=1
    return t

def conv(a,b):
    return [sum(a[j]*b[n-j] for j in range(n+1)) for n in range(M+1)]

def v2(a):
    return [a[n//2] if n%2==0 else 0 for n in range(M+1)]

t1,t2=theta(1,1,6),theta(2,1,3)
assert all((a-b)%2==0 for a,b in zip(t1,t2))
B=[(a-b)//2 for a,b in zip(t1,t2)]
v=conv(B,B); E=conv(v2(t1),B); W=conv(v2(t1),v2(B))
assert all((a-b)%2==0 for a,b in zip(v,W))
G=[(a-b)//2 for a,b in zip(v,W)]
qr={r*r%23 for r in range(1,23)}
D=[G[n]-(n%4==3)*E[n] if n%2 else 0 for n in range(M+1)]
P=[D[n] if n%23 in qr else 0 for n in range(M+1)]
assert all(p%2==0 for p in P)
assert (G[3],E[3],P[3]%2)==(-1,1,0)
# Negative control: switch the sign of the quadratic mod-4 projection.
assert (G[3]-(3%4==1)*E[3])%2==1
print('level =',N,'Sturm bound =',M,'coefficients checked =',M+1)
print('nonzero coefficients mod 2 =',sum(p%2 for p in P))
print('mod-2 vector SHA256 =',sha256(bytes(p%2 for p in P)).hexdigest())
print('negative control at p=3: detected')
