import numpy as np, json, cypari2, random
from math import sqrt, pi, log
d=np.load('ap_20000000.npz'); p,u,v=d['p'],d['u'],d['v']
good=(p!=2)&(p!=23); p,u,v=p[good],u[good],v[good]
a=u  # at collision a_p = u
res={}
# independent point-count cross-check (PARI hyperellcharpoly on y^2=(x^3+x+1)^2-4x^5)
pari=cypari2.Pari(); pari.allocatemem(10**9); F=pari('(x^3+x+1)^2-4*x^5')
random.seed(1); idx=[i for i in range(len(p)) if p[i]<3000]+random.sample([i for i in range(len(p)) if 3000<=p[i]<12000],40)
bad=0
for i in idx:
    cp=pari.hyperellcharpoly(F*pari.Mod(1,int(p[i]))); P=int(p[i])
    # (x^2-a1 x+p)(x^2-a2 x+p), a1=u+v*phi, a2=u+v*phibar: t=2u+v, a1*a2=u^2+uv-v^2
    t=2*int(u[i])+int(v[i]); n=int(u[i])**2+int(u[i])*int(v[i])-int(v[i])**2
    want=[1,-t,n+2*P,-P*t,P*P]
    got=[int(cp.polcoef(4-j)) for j in range(5)]
    bad+= (want!=got)
res['pointcount_crosscheck']={'primes':len(idx),'mismatches':bad}
# T1
res['T1_v_even_violations']=int((v%2!=0).sum())
# T2
ok={r:bool(np.all(((u+v*r)-1-p)%11==0)) for r in (4,8)}
res['T2_root_holds_for_all_p']=ok
# collisions
c=v==0; pc,ac=p[c],a[c]
res['T3_violations']=int((((ac-pc-1)%11)!=0).sum())
# T4: p represented by 2x^2+xy+3y^2  <=> (-23/p)=1 and p not rep by x^2+xy+6y^2
def rep(P,A,B,C):
    ymax=int(sqrt(4*A*P/(4*A*C-B*B)))+1
    for y in range(-ymax,ymax+1):
        D=(B*y)**2-4*A*(C*y*y-P)
        if D<0: continue
        s=int(round(sqrt(D)))
        for ss in (s-1,s,s+1):
            if ss>=0 and ss*ss==D and (-B*y+ss)%(2*A)==0: return True
    return False
t4=sum(((int(x)%2)==1)!=rep(int(P),2,1,3) for P,x in zip(pc,ac))
res['T4_violations']=t4
res['collisions_first']=[(int(P),int(x)) for P,x in zip(pc[:20],ac[:20])]
# counts vs prereg
K=1.916948; out=[]
for X in (10**5,10**6,10**7,2*10**7):
    m=p<=X; S=float((p[m]**-0.5).sum() - (2**-0.5 if False else 0)); 
    S=float((p[m&(p>=3)]**-0.5).sum())
    obs=int((c&m).sum()); pred=K*S
    out.append(dict(x=X,obs=obs,pred=round(pred,1),z=round((obs-pred)/sqrt(pred),2),ratio=round(obs/pred,4),
                    fitted_c_sqrtx_logx=round(obs/(sqrt(X)/log(X)),4)))
res['counts']=out
# E tests
E={}
n=len(v)
for l,ks,pred in ((2,(4,8,16,32),(0.583333,0.302083,0.152344,None)),(3,(3,9,27),(0.3625,0.121759,None)),(5,(5,25,125),(0.206667,0.041387,None)),(11,(11,121),(0.099167,0.009084)),(7,(7,49),(None,None)),(13,(13,169),(None,None))):
    for kk,pr in zip(ks,pred):
        f=float((v%kk==0).mean()); se=sqrt(f*(1-f)/n)
        E[f'{kk}|v']=dict(freq=round(f,6),scaled=round(f*kk,5),model=pr,z=(round((f-pr)/sqrt(pr*(1-pr)/n),2) if pr else None))
res['E']=E
json.dump(res,open('results.json','w'),indent=1); print(json.dumps(res,indent=1))
