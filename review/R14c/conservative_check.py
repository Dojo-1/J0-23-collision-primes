#!/usr/bin/env python3
"""Independent larger-level check, via an exact carry-free modular square."""
from verify import eta_coeffs
from hashlib import sha256
N=101568; B=eta_coeffs(N)
w=3; assert 9*(N+1)<256**w
z=int.from_bytes(b''.join((b%4).to_bytes(w,'little') for b in B),'little')
raw=(z*z % (1<<(8*w*(N+1)))).to_bytes(w*(N+1),'little')
qr={r*r%23 for r in range(1,23)}; par=bytearray(N+1)
for n in range(1,N+1,2):
    v=int.from_bytes(raw[w*n:w*(n+1)],'little')%4
    assert v%2==0
    if n%23 in qr: par[n]=(v//2-(n%4==3)*B[n])%2
assert not any(par)
print('bound:',N,'nonzero coefficients:',sum(par),'SHA256:',sha256(par).hexdigest())
