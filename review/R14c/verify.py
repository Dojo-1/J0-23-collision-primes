#!/usr/bin/env python3
"""Exact checks; standard library only. No imported author data."""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt, lcm
from collections import Counter
from itertools import product
import json

# B=q*Euler(q)*Euler(q^23), using the pentagonal number theorem.
def eta_coeffs(N):
    e = {j*(3*j-1)//2: (-1 if j%2 else 1)
         for j in range(-isqrt(N)-2,isqrt(N)+3)
         if j*(3*j-1)//2 < N}
    B=[0]*(N+1)
    for n,a in e.items():
        for m,b in e.items():
            k=1+n+23*m
            if k<=N:B[k]+=a*b
    return B

def theta(a,b,c,N):
    # Complete square: (2*a*x+b*y)^2 + (4*a*c-b*b)*y^2 = 4*a*Q.
    d=4*a*c-b*b; Y=isqrt(4*a*N//d); out=[0]*(N+1)
    for y in range(-Y,Y+1):
        s=isqrt(4*a*N-d*y*y)
        lo=-((s+b*y)//(2*a)); hi=(s-b*y)//(2*a)
        for x in range(lo,hi+1):out[a*x*x+b*x*y+c*y*y]+=1
    return out

def convolution(a,b,N):
    # Exact signed integer convolution by nonnegative split and carry-free packing.
    # Each block exceeds the worst coefficient bound; no floating point.
    mx=max(map(abs,a))*max(map(abs,b))*(N+1)
    w=max(1,(mx.bit_length()+7)//8)
    def pack(x,sgn):
        return int.from_bytes(b''.join(max(0,sgn*i).to_bytes(w,'little') for i in x),'little')
    A,Am=pack(a,1),pack(a,-1); B,Bm=pack(b,1),pack(b,-1)
    def unpack(z):
        z=(z % (1<<(8*w*(N+1)))).to_bytes(w*(N+1),'little')
        return [int.from_bytes(z[i:i+w],'little') for i in range(0,len(z),w)]
    pp=unpack(A*B); mm=unpack(Am*Bm); pm=unpack(A*Bm); mp=unpack(Am*B)
    return [pp[i]+mm[i]-pm[i]-mp[i] for i in range(N+1)]

def sturm_check(N):
    B=eta_coeffs(N); t1=theta(1,1,6,N); t2=theta(2,1,3,N)
    assert all(2*b==x-y for b,x,y in zip(B,t1,t2))
    v=convolution(B,B,N)
    qr={i*i%23 for i in range(1,23)}
    # At odd n, G_n=v_n/2; E_n=B_n (mod 2).
    assert all(v[n]%2==0 for n in range(1,N+1,2))
    parity=[0]*(N+1)
    for n in range(1,N+1,2):
        if n%23 in qr:parity[n]=(v[n]//2-(n%4==3)*B[n])%2
    assert not any(parity)
    # Negative control: reverse the p mod 4 mask; p=3 gives failure.
    bad=[n for n in range(1,N+1,2) if n%23 in qr and (v[n]//2-(n%4==1)*B[n])%2]
    assert 3 in bad and bad[0]==1
    return dict(bound=N,nonzero_checked=sum(parity),split_odd_indices=sum(n%23 in qr for n in range(1,N+1,2)),
                parity_sha256=sha256(bytes(parity)).hexdigest(),negative_control_first=bad[0],negative_control_prime=3),B,t1,t2,v

def sieve(N):
    a=bytearray(b'\1')*(N+1);a[:2]=b'\0\0'
    for p in range(2,isqrt(N)+1):
        if a[p]:a[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
    return [p for p in range(3,N+1) if a[p] and p!=23]

# Field operations in F_p^2; phi^2=phi+1.
def model(l):
    R=list(product(range(l),repeat=2)); zero=(0,0); one=(1,0)
    def add(x,y):return ((x[0]+y[0])%l,(x[1]+y[1])%l)
    def neg(x):return ((-x[0])%l,(-x[1])%l)
    def mul(x,y):return ((x[0]*y[0]+x[1]*y[1])%l,(x[0]*y[1]+x[1]*y[0]+x[1]*y[1])%l)
    def det(A):return add(mul(A[0],A[3]),neg(mul(A[1],A[2])))
    def tr(A):return add(A[0],A[3])
    def mm(A,B):return tuple(add(mul(A[2*i],B[j]),mul(A[2*i+1],B[j+2])) for i in range(2) for j in range(2))
    I=(one,zero,zero,one)
    def gen(gs):
        S={I}; todo=[I]
        for a in todo:
            for g in gs:
                b=mm(a,g)
                if b not in S:S.add(b);todo.append(b)
        return S
    if l==2:
        H={A for A in product([zero,one],repeat=4) if det(A)==one}
    else:
        # (2,3,5) projective generators: X^2=-I, Y^3=-I, (XY)^5=-I.
        X=(zero,neg(one),one,zero)
        Y=(one,(1,1),(0,1),zero) # replace by deterministic search
        SL=[A for A in product(R,repeat=4) if det(A)==one]
        for Y in SL:
            if tr(Y)!=one:continue
            S=gen([X,Y])
            if len(S)==120:break
        else:raise AssertionError('no subgroup')
        i=(1,1);H=S|{tuple(mul(i,e) for e in A) for A in S}
    assert len(H)==(6 if l==2 else 240)
    assert all(mm(A,B) in H for A in H for B in H)
    good=[A for A in H if tr(A)[1]==0]
    scalar=[A for A in good if A[1]==zero and A[2]==zero and A[0]==A[3]]
    # Brute lift all matrices one layer, no use of the analytic recurrence.
    m=l*l
    def mulm(x,y):return ((x[0]*y[0]+x[1]*y[1])%m,(x[0]*y[1]+x[1]*y[0]+x[1]*y[1])%m)
    total=hit=0
    off=list(product(range(l),repeat=8))
    for A in H:
        for z in off:
            a,b,c,d=[(e[0]+l*z[2*j],e[1]+l*z[2*j+1]) for j,e in enumerate(A)]
            if (mulm(a,d)[1]-mulm(b,c)[1])%m:continue
            total+=1;hit+=(a[1]+d[1])%m==0
    f2=F(m*hit,total)
    # Inner-product formula for scalar residue classes.
    def inner(m):return sum(F(l**j*(l**3-1),l**(3*(j+1))) for j in range(m))+F(l**m,l**(3*m))
    ns=len(good)-len(scalar); s=len(scalar); h=len(H)
    factors=[F(l*len(good),h)] + [F(l*ns,h)+F(l*l*s,h)*inner(k-2) for k in range(2,9)]
    limit=F(l*ns,h)+F(l*l*s,h)*F(l**3-1,l*(l*l-1))
    assert f2==factors[1]
    return dict(l=l,group=h,rational=len(good),scalar=s,trace_histogram={str(k):v for k,v in sorted(Counter(map(tr,H)).items())},
                brute_mod_l2_total=total,brute_mod_l2_hits=hit,factors=list(map(str,factors)),limit=str(limit))

def main():
    # Use the conservative level N=46*16*23^2; also record sharper level 16*23^2.
    level=46*16*23**2; bound=level*3//2*24//23//6
    assert level==389344 and bound==101568
    cert,B,t1,t2,v=sturm_check(2208)
    # Validate coefficient formula and theorem tests in an independently generated prefix.
    u=convolution(t2,B,2208); P=sieve(2208); qr={i*i%23 for i in range(1,23)}
    collisions=[p for p in P if v[p]==0]
    for p in P:
        assert v[p]%2==0 and (u[p]+4*v[p]-1-p)%11==0
        if p%23 in qr and u[p]%2==0:assert (u[p]-1-p)%4==0 and v[p]%4==0
        if p%23 in qr and u[p]%2:assert (v[p]-p+1)%4==0
    # Exact quadratic-form check, independent of parity selector.
    for p in P:
        assert bool(t2[p])==(p%23 in qr and u[p]%2==1)
    out=dict(sturm=cert,conservative_level=level,sharp_level=16*23**2,sharp_sturm=2208,
             prefix=dict(primes=len(P),collisions=len(collisions),first=collisions[:12],zero_trace=[p for p in collisions if u[p]==0],
                         bound=2208),models=[model(2),model(3)])
    assert out["prefix"]["first"][:10]==[13,19,29,43,109,223,229,271,463,673]
    print(json.dumps(out,indent=2));open('results.json','w').write(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
