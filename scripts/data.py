"""Hecke eigenvalues a_p(f) = u_p + v_p*phi of the level-23 newform for all p <= N.
Method: M_2(Gamma0(23)) is spanned by theta products th1^2, th1*th2, th2^2, where
th1 = sum q^{x^2+xy+6y^2}, th2 = sum q^{2x^2+xy+3y^2} (weight 1, character (-23/.)).
f = g0 + phi*g1 with g0, g1 integral; combination solved from PARI mfeigenbasis and
checked EXACTLY against PARI up to n = CHECK."""
import sys, numpy as np, cypari2
from fractions import Fraction

N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**6
CHECK = 10**4


def theta(a, b, c, N):
    """r_Q(n) for Q=a x^2+bxy+cy^2 (positive definite), n<=N."""
    out = np.zeros(N + 1, dtype=np.int64)
    disc = 4 * a * c - b * b
    ymax = int(np.sqrt(4 * a * N / disc)) + 1
    chunks = []
    for y in range(-ymax, ymax + 1):
        # a x^2 + b y x + c y^2 <= N
        A, B, C = a, b * y, c * y * y - N
        D = B * B - 4 * A * C
        if D < 0:
            continue
        s = np.sqrt(D)
        lo = int(np.floor((-B - s) / (2 * A))) - 1
        hi = int(np.ceil((-B + s) / (2 * A))) + 1
        x = np.arange(lo, hi + 1, dtype=np.int64)
        v = a * x * x + b * x * y + c * y * y
        chunks.append(v[(v >= 0) & (v <= N)])
        if sum(len(c) for c in chunks) > 2*10**7:
            out += np.bincount(np.concatenate(chunks), minlength=N + 1); chunks.clear()
    if chunks:
        out += np.bincount(np.concatenate(chunks), minlength=N + 1)
    return out


def conv(u, w, N):
    L = 1 << int(np.ceil(np.log2(2 * N + 2)))
    U = np.fft.rfft(u.astype(np.float64), L)
    W = U if w is u else np.fft.rfft(w.astype(np.float64), L)
    r = np.fft.irfft(U * W, L)[: N + 1]
    del U, W
    ri = np.rint(r)
    err = np.abs(r - ri).max()
    assert err < 0.2, err
    return ri.astype(np.int64), err


def main():
    t1 = theta(1, 1, 6, N); t2 = theta(2, 1, 3, N)
    P11, e1 = conv(t1, t1, N); P12, e2 = conv(t1, t2, N); P22, e3 = conv(t2, t2, N)
    del t1, t2
    print("fft max rounding residual", max(e1, e2, e3), flush=True)
    pari = cypari2.Pari(); pari.allocatemem(10**9)
    pari('mf=mfinit([23,2],0); L=mfeigenbasis(mf)')
    co = pari(f'mfcoefs(L[1],{CHECK})')
    g0 = [int(pari.lift(c).polcoef(0)) for c in co]
    g1 = [int(pari.lift(c).polcoef(1)) for c in co]
    # solve g = x*P11 + y*P12 + z*P22 on n=0..10 (exact rationals)
    import sympy as sp
    Mx = sp.Matrix([[int(P11[n]), int(P12[n]), int(P22[n])] for n in range(12)])
    res = []
    for g in (g0, g1):
        sol, params = Mx.gauss_jordan_solve(sp.Matrix(g[:12]))
        assert len(params) == 0
        res.append([sp.Rational(s) for s in sol])
    print("combinations:", res, flush=True)
    G = []
    for (x, y, z), g in zip(res, (g0, g1)):
        den = sp.ilcm(*[sp.Rational(t).q for t in (x, y, z)])
        X, Y, Z = int(x * den), int(y * den), int(z * den)
        val = X * P11 + Y * P12 + Z * P22
        assert np.all(val % den == 0)
        val //= den
        assert all(int(val[n]) == g[n] for n in range(CHECK + 1)), "mismatch vs PARI"
        G.append(val)
    print(f"exact agreement with PARI mfcoefs for n<={CHECK}: OK", flush=True)
    from sympy import primerange
    ps = np.array(list(primerange(2, N + 1)), dtype=np.int64)
    np.savez_compressed(f"ap_{N}.npz", p=ps, u=G[0][ps], v=G[1][ps])
    print("saved", len(ps), "primes")


if __name__ == "__main__":
    main()
