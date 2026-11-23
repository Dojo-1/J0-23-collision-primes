"""Memory-lean variant: v_p = [q^p] B^2, u_p = [q^p] th2*B, with B = (th1-th2)/2 = eta(z)eta(23z).
(Identities established by exact agreement with PARI in data.py.)  Saves values at primes only."""
import sys, gc, numpy as np, scipy.fft as sf
from data import theta

N = int(float(sys.argv[1]))
t1 = theta(1, 1, 6, N).astype(np.int16); t2 = theta(2, 1, 3, N).astype(np.int16)
B = (t1 - t2); del t1; gc.collect()
assert np.all(B % 2 == 0); B //= 2
L = 1 << int(np.ceil(np.log2(2 * N + 2)))
buf = np.zeros(L); buf[:N+1] = B; del B
FB = sf.rfft(buf, overwrite_x=True); del buf; gc.collect()
P = np.load('primes_6e7.npy'); P = P[P <= N]
def coeffs(spec):
    r = sf.irfft(spec, L, overwrite_x=True)
    vals = r[P]; err = np.abs(vals - np.rint(vals)).max()
    del r
    return np.rint(vals).astype(np.int64), err
S = FB * FB
v, e1 = coeffs(S); del S; gc.collect()
buf = np.zeros(L); buf[:N+1] = t2; del t2
F2 = sf.rfft(buf, overwrite_x=True); del buf; gc.collect()
F2 *= FB; del FB; gc.collect()
u, e2 = coeffs(F2); del F2
print('rounding residual at primes', e1, e2, flush=True)
# overlap check with the PARI-validated 2e7 file
d = np.load('ap_20000000.npz'); m = np.isin(d['p'], P)
idx = np.searchsorted(P, d['p'][m])
assert np.array_equal(v[idx], d['v'][m]) and np.array_equal(u[idx], d['u'][m]), 'overlap mismatch'
print('overlap with validated 2e7 data: identical on', int(m.sum()), 'primes', flush=True)
np.savez_compressed(f'ap_lean_{N}.npz', p=P, u=u, v=v)
