"""Independent consistency check of the source-reported projective mod-3 image A5 (Billerey-Dieulefait §5.3).
For p not dividing disc(X^5+3X^3+6X^2+9) = 3^8 * 23^2, compare the cycle type of Frob_p on the quintic's roots
with r = t^2/det in F9, t = T_p mod 3 (u_p + v_p*phi), det = p:
  5-cycle <-> r not in F3;  (2,2,1) <-> t = 0;  1^5 or (3,1,1) <-> r = 1.
Run from repo root: python3 scripts/mod3_quintic_check.py [bound]   (needs data/ap_lean_60000000.npz, sympy)"""
import sys, numpy as np, sympy as sp
from collections import Counter
B = int(float(sys.argv[1])) if len(sys.argv) > 1 else 200000
d = np.load('data/ap_lean_60000000.npz'); P, U, V = d['p'], d['u'], d['v']
x = sp.symbols('x'); Q = x**5 + 3*x**3 + 6*x**2 + 9
mul = lambda a, b: ((a[0]*b[0]+a[1]*b[1]) % 3, (a[0]*b[1]+a[1]*b[0]+a[1]*b[1]) % 3)
inv = {(1, 0): (1, 0), (2, 0): (2, 0)}
tally, mism = Counter(), 0
for p, u, v in zip(P, U, V):
    p = int(p)
    if p > B: break
    if p == 3: continue
    t = (int(u) % 3, int(v) % 3); r = mul(mul(t, t), inv[(p % 3, 0)])
    cls = 't0' if t == (0, 0) else ('r1' if r == (1, 0) else ('out' if r[1] else 'other'))
    ct = tuple(sorted(sum([[f.degree()]*e for f, e in sp.Poly(Q, x, modulus=p).factor_list()[1]], []), reverse=True))
    exp = {(5,): 'out', (2, 2, 1): 't0', (1, 1, 1, 1, 1): 'r1', (3, 1, 1): 'r1'}.get(ct, 'BAD')
    tally[ct] += 1; mism += exp != cls
print('primes', sum(tally.values()), 'mismatches', mism, dict(tally))
