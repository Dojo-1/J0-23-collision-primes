import itertools, numpy as np, random
from collections import Counter
# F9 = F3[phi]/(phi^2-phi-1); element (x,y)=x+y*phi ; encode as 3x+y? use tuples
def mul(a,b,m=3):
    x1,y1=a;x2,y2=b; return ((x1*x2+y1*y2)%m,(x1*y2+x2*y1+y1*y2)%m)
def add(a,b,m=3): return ((a[0]+b[0])%m,(a[1]+b[1])%m)
def neg(a,m=3): return ((-a[0])%m,(-a[1])%m)
F9=[(x,y) for x in range(3) for y in range(3)]
def mm(A,B): return tuple(tuple(add(mul(A[i][0],B[0][j]),mul(A[i][1],B[1][j])) for j in range(2)) for i in range(2))
def det(A): return add(mul(A[0][0],A[1][1]),neg(mul(A[0][1],A[1][0])))
def tr(A): return add(A[0][0],A[1][1])
I=(((1,0),(0,0)),((0,0),(1,0)))
SL=[ ((a,b),(c,d)) for a in F9 for b in F9 for c in F9 for d in F9 if det(((a,b),(c,d)))==(1,0)]
def order(A):
    B=A;k=1
    while B!=I: B=mm(B,A);k+=1
    return k
def gen(gs):
    S={I};fr=[I]
    while fr:
        nf=[]
        for x in fr:
            for g in gs:
                y=mm(x,g)
                if y not in S: S.add(y);nf.append(y)
        fr=nf
        if len(S)>720: break
    return S
random.seed(0)
o4=[A for A in SL if order(A)==4]; o3=[A for A in SL if order(A)==3]
target={(0,0):30,(1,0):21,(2,0):21,(1,1):12,(2,2):12,(0,1):12,(0,2):12}
found=None
for _ in range(4000):
    S=gen([random.choice(o4),random.choice(o3)])
    if len(S)==120:
        found=S;break
C=Counter(tr(A) for A in found); print(len(found),dict(C))
# trace set of 2.A5 in F9 (det 1) vs empirical p≡1 mod 3 pattern: missing should be (1,1),(2,2)?
i=((1,1)) # phi^2 = 1+phi
G3=set(found)|{tuple(tuple(mul(i,e) for e in row) for row in A) for A in found}
print('|G3|',len(G3), 'dets',Counter(det(A) for A in G3))
for dd in [(1,0),(2,0)]:
    c=Counter(tr(A) for A in G3 if det(A)==dd); n=sum(c.values()); print(dd,{k:round(v/n,4) for k,v in sorted(c.items())})
import pickle; pickle.dump(sorted(G3),open('G3.pkl','wb'))
