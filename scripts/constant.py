from math import pi, sqrt, log
from sympy import primerange
def vh_inert(l): return l**3*(l*l+l+1)/((l+1)*(l**4-1))
def vh_split(l): return l*l*(l**3+l*l-l-2)/((l*l-1)**2*(l+1))
exc={2:22/9, 5:1.03472, 11:1.09924}   # model factors (2: Gamma2 limit; 5: ramified full, extrapolated; 11: Eisenstein diag model)
P=1.0; Pgen=1.0
for l in primerange(2,2*10**6):
    g = vh_split(l) if l%5 in (1,4) else (vh_inert(l) if l!=5 else 1.0)
    Pgen*= g if l not in exc else 1.0
    P*= exc.get(l,g)
F=P; c=16*sqrt(5)*F/(3*pi**2)
print("F_hat(model)=%.5f  c_model=%.4f"%(F,c))
print("generic part (all l except 2,5,11)=%.5f"%Pgen)
print("F with generic l=2 (56/45), l=11 split, l=5 ->1: %.5f"%(Pgen*56/45*vh_split(11)))
for x in (1e5,1e6,1e7,1e8):
    S=sum(p**-0.5 for p in primerange(3,int(x)+1) if p!=23) if x<=1e7 else None
    print(x, 'ratio sum p^-1/2 / (2 sqrt x/log x) =', (S/(2*sqrt(x)/log(x))) if S else '')
