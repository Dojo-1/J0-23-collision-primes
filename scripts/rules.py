import numpy as np, json
d=np.load('ap_20000000.npz'); p,u,v=d['p'],d['u'],d['v']
g=(p>3)&(p!=23); p,u,v=p[g],u[g],v[g]
# Kronecker (-23/p) = (p/23) by reciprocity for odd p != 23
QR=set((x*x)%23 for x in range(1,23)); eps=np.where(np.isin(p%23,list(QR)),1,-1)
R1=(eps==1)&(u%2==0); R2=(eps==1)&(u%2==1)
out=dict(
 n=int(len(p)),
 R1_count=int(R1.sum()), R1_viol=int((((u[R1]-1-p[R1])%4!=0)|(v[R1]%4!=0)).sum()),
 R2_count=int(R2.sum()), R2_viol=int(((v[R2]-p[R2]+1)%4!=0).sum()),
 inert_v_mod4_dist={str(k):round(float((v[eps==-1]%4==k).mean()),4) for k in (0,2)},
 F2_emp={str(2**k):round(float((v%2**k==0).mean()*2**k),4) for k in range(1,8)},
 F3_emp={str(3**k):round(float((v%3**k==0).mean()*3**k),4) for k in range(1,5)},
)
print(json.dumps(out,indent=1)); json.dump(out,open('rules.json','w'),indent=1)
