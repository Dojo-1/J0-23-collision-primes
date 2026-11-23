import numpy as np, json
from math import sqrt, log
d=np.load('ap_lean_60000000.npz'); p,u,v=d['p'],d['u'],d['v']
K2=2.432194; K1=1.916948
QR=set((x*x)%23 for x in range(1,23)); eps=np.where(np.isin(p%23,list(QR)),1,-1)
col=v==0
def rep2(P):  # p = 2x^2+xy+3y^2 ?
    ymax=int(sqrt(8*P/23))+1
    for y in range(-ymax,ymax+1):
        D=y*y-8*(3*y*y-P)
        if D<0: continue
        s=int(sqrt(D))
        for t in (s-1,s,s+1):
            if t>=0 and t*t==D and (-y+t)%4==0: return True
    return False
pc,uc=p[col],u[col]
R1=(eps==1)&(u%2==0); R2=(eps==1)&(u%2==1)
out=dict(primes=int(len(p)), collisions_total=int(col.sum()),
 T1_viol=int((v%2!=0).sum()), T2_viol=int((((u+4*v)-1-p)%11!=0).sum()),
 T3_viol=int((((uc-pc-1)%11)!=0).sum()),
 T4_viol=int(sum((int(a)%2==1)!=rep2(int(q)) for q,a in zip(pc,uc))),
 R1_n=int(R1.sum()), R1_viol=int((((u[R1]-1-p[R1])%4!=0)|(v[R1]%4!=0)).sum()),
 R2_n=int(R2.sum()), R2_viol=int(((v[R2]-p[R2]+1)%4!=0).sum()),
 mod3_forbidden_viol=int((((p%3==1)&(((u%3==1)&(v%3==1))|((u%3==2)&(v%3==2)))) | ((p%3==2)&(v%3==0)&(u%3!=0))).sum()),
 supersingular_a0=[int(q) for q in pc[uc==0]][:12], a0_count=int((uc==0).sum()),
 a0_p_mod11_ok=bool(np.all(pc[uc==0]%11==10)))
ho=[]
for a,b in ((2e7,4e7),(4e7,6e7),(2e7,6e7)):
    m=(p>a)&(p<=b); S=float((p[m]**-0.5).sum()); n=int((col&m).sum()); pr=K2*S
    ho.append(dict(range=f"({a:.0e},{b:.0e}]",obs=n,pred_M2=round(pr,1),z_M2=round((n-pr)/sqrt(pr),2),pred_M1=round(K1*S,1),z_M1=round((n-K1*S)/sqrt(K1*S),2)))
out['heldout']=ho
cum=[]
for x in (1e5,1e6,1e7,2e7,4e7,6e7):
    m=p<=x; S=float((p[m]**-0.5).sum()); n=int((col&m).sum())
    cum.append(dict(x=x,obs=n,pred_M2=round(K2*S,1),ratio=round(n/(K2*S),4)))
out['cumulative']=cum
print(json.dumps(out,indent=1)); json.dump(out,open('heldout.json','w'),indent=1)
