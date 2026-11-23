"""R14b — collision primes of J_0(23): exact congruences, icosahedral mod-3 law, density constant.

Collision prime: p != 2, 23 with a_p(f) = a_p(f^sigma) for the level-23 newform pair, i.e. v_p = 0 where
a_p(f) = u_p + v_p*phi,  phi^2 = phi + 1  (PARI's mfeigenbasis coefficient field y^2 - y - 1).

Coefficient engine (exact, numpy FFT, validated against PARI mfcoefs to n = 10^4 in R14b runs):
    th1 = sum q^{x^2+xy+6y^2},  th2 = sum q^{2x^2+xy+3y^2},  B = (th1 - th2)/2 = eta(z) eta(23z)
    v_p = [q^p] B^2,            u_p = [q^p] th2 * B.
Requires numpy (scipy optional for big N). Run:  python3 -m experiments.collision_primes_r14b 200000
"""
from __future__ import annotations
import itertools, json, math, sys
from fractions import Fraction
import numpy as np

# ---------------------------------------------------------------- coefficients
def theta(a: int, b: int, c: int, N: int) -> np.ndarray:
    out = np.zeros(N + 1, dtype=np.int64)
    ymax = int(math.sqrt(4 * a * N / (4 * a * c - b * b))) + 1
    chunks = []
    for y in range(-ymax, ymax + 1):
        D = (b * y) ** 2 - 4 * a * (c * y * y - N)
        if D < 0:
            continue
        s = math.sqrt(D)
        x = np.arange(int((-b * y - s) // (2 * a)) - 1, int((-b * y + s) // (2 * a)) + 2, dtype=np.int64)
        val = a * x * x + b * x * y + c * y * y
        chunks.append(val[(val >= 0) & (val <= N)])
    out += np.bincount(np.concatenate(chunks), minlength=N + 1)
    return out


def primes_upto(N: int) -> np.ndarray:
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def eigen_uv(N: int):
    """Return (p, u_p, v_p) for primes 3 <= p <= N, p != 23 (exact integers)."""
    t1 = theta(1, 1, 6, N); t2 = theta(2, 1, 3, N)
    B = t1 - t2
    assert np.all(B % 2 == 0)
    B //= 2
    L = 1 << int(math.ceil(math.log2(2 * N + 2)))
    FB = np.fft.rfft(B.astype(np.float64), L)
    F2 = np.fft.rfft(t2.astype(np.float64), L)
    v = np.fft.irfft(FB * FB, L)[: N + 1]
    u = np.fft.irfft(FB * F2, L)[: N + 1]
    P = primes_upto(N); P = P[(P >= 3) & (P != 23)]
    for arr in (u, v):
        assert np.abs(arr[P] - np.rint(arr[P])).max() < 1e-3
    return P, np.rint(u[P]).astype(np.int64), np.rint(v[P]).astype(np.int64), B


# ---------------------------------------------------------------- arithmetic laws
QR23 = sorted({(x * x) % 23 for x in range(1, 23)})


def kronecker_m23(p: np.ndarray) -> np.ndarray:
    """(-23/p) = (p/23) for odd p != 23."""
    return np.where(np.isin(p % 23, QR23), 1, -1)


def rep_by(P: int, a: int, b: int, c: int) -> bool:
    ymax = int(math.sqrt(4 * a * P / (4 * a * c - b * b))) + 1
    for y in range(-ymax, ymax + 1):
        D = (b * y) ** 2 - 4 * a * (c * y * y - P)
        if D < 0:
            continue
        s = math.isqrt(D)
        if s * s == D and ((-b * y + s) % (2 * a) == 0 or (-b * y - s) % (2 * a) == 0):
            return True
    return False


def law_audit(p, u, v) -> dict:
    eps = kronecker_m23(p)
    col = v == 0
    pc, ac = p[col], u[col]
    R1 = (eps == 1) & (u % 2 == 0)
    R2 = (eps == 1) & (u % 2 == 1)
    forb3 = ((p % 3 == 1) & (((u % 3 == 1) & (v % 3 == 1)) | ((u % 3 == 2) & (v % 3 == 2)))) | \
            ((p % 3 == 2) & (v % 3 == 0) & (u % 3 != 0))
    return dict(
        primes=int(len(p)), collisions=int(col.sum()),
        L1_v_even_violations=int((v % 2).astype(bool).sum()),
        L2_eisenstein_P11_violations=int(((u + 4 * v - 1 - p) % 11 != 0).sum()),   # P11 = (11, phi - 4)
        T_mod11_violations=int(((ac - pc - 1) % 11 != 0).sum()),
        T_parity_violations=int(sum((int(a) % 2 == 1) != rep_by(int(q), 2, 1, 3) for q, a in zip(pc, ac))),
        R1_cases=int(R1.sum()),
        R1_violations=int((((u[R1] - 1 - p[R1]) % 4 != 0) | (v[R1] % 4 != 0)).sum()),
        R2_cases=int(R2.sum()), R2_violations=int(((v[R2] - p[R2] + 1) % 4 != 0).sum()),
        icosahedral_mod3_violations=int(forb3.sum()),
        a_p_zero_primes=[int(q) for q in pc[ac == 0]],
    )


# ---------------------------------------------------------------- local factors (exact counting)
def _mul(a, b, m):
    return ((a[0] * b[0] + a[1] * b[1]) % m, (a[0] * b[1] + a[1] * b[0] + a[1] * b[1]) % m)


def preimage_factor(base, l: int, k: int, cond_mod: int | None = None):
    """l^k * Prob[ tr g in Z/l^k ] for g Haar on the FULL preimage in GL2(O/l^k) (det in (Z/l^k)^x)
    of a set `base` of 2x2 matrices over O/l, entries (x, y) = x + y*phi."""
    m, Lk = l ** k, l ** (k - 1)
    offs = np.array(list(itertools.product(range(Lk), repeat=8)), dtype=np.int64)
    num = den = 0; numc, denc = {}, {}
    for A in base:
        e = [A[0][0], A[0][1], A[1][0], A[1][1]]
        X = [(e[j][0] + l * offs[:, 2 * j]) % m for j in range(4)]
        Y = [(e[j][1] + l * offs[:, 2 * j + 1]) % m for j in range(4)]
        ax, ay = _mul((X[0], Y[0]), (X[3], Y[3]), m); bx, by = _mul((X[1], Y[1]), (X[2], Y[2]), m)
        ok = ((ay - by) % m == 0); tro = ((Y[0] + Y[3]) % m == 0)
        num += int((ok & tro).sum()); den += int(ok.sum())
        if cond_mod:
            dx = (ax - bx) % m
            for r in range(cond_mod):
                s = ok & (dx % cond_mod == r)
                if s.any():
                    numc[r] = numc.get(r, 0) + int((s & tro).sum()); denc[r] = denc.get(r, 0) + int(s.sum())
    return Fraction(m * num, den), {r: Fraction(numc[r], denc[r]) for r in denc}


def gl2_over_O_mod_l(l: int, pred=lambda e: True):
    F = [(x, y) for x in range(l) for y in range(l) if pred((x, y))]
    out = []
    for a, b, c, d in itertools.product(F, repeat=4):
        D = ((_mul(a, d, l)[0] - _mul(b, c, l)[0]) % l, (_mul(a, d, l)[1] - _mul(b, c, l)[1]) % l)
        if D[1] == 0 and D[0] != 0:
            out.append(((a, b), (c, d)))
    return out


def binary_icosahedral_G3(seed: int = 0):
    """Subgroup <i> * 2.A5 of GL2(F9), F9 = F3[phi], i = phi^2 (i^2 = -1). Order 240."""
    import random
    F9 = [(x, y) for x in range(3) for y in range(3)]
    I = (((1, 0), (0, 0)), ((0, 0), (1, 0)))
    def mm(A, B):
        return tuple(tuple(tuple((s + t) % 3 for s, t in zip(_mul(A[i][0], B[0][j], 3), _mul(A[i][1], B[1][j], 3)))
                           for j in range(2)) for i in range(2))
    def det(A):
        p, q = _mul(A[0][0], A[1][1], 3), _mul(A[0][1], A[1][0], 3)
        return ((p[0] - q[0]) % 3, (p[1] - q[1]) % 3)
    SL = [((a, b), (c, d)) for a in F9 for b in F9 for c in F9 for d in F9 if det(((a, b), (c, d))) == (1, 0)]
    def order(A):
        B, k = A, 1
        while B != I:
            B, k = mm(B, A), k + 1
        return k
    o4 = [A for A in SL if order(A) == 4]; o3 = [A for A in SL if order(A) == 3]
    rnd = random.Random(seed)
    while True:
        gens = [rnd.choice(o4), rnd.choice(o3)]
        S, fr = {I}, [I]
        while fr and len(S) <= 120:
            nf = []
            for x in fr:
                for g in gens:
                    y = mm(x, g)
                    if y not in S:
                        S.add(y); nf.append(y)
            fr = nf
        if len(S) == 120:
            break
    i = (1, 1)
    return sorted(S | {tuple(tuple(_mul(i, e, 3) for e in row) for row in A) for A in S})


def vh_inert(l):
    return l ** 3 * (l * l + l + 1) / ((l + 1) * (l ** 4 - 1))


def vh_split(l):
    return l * l * (l ** 3 + l * l - l - 2) / ((l * l - 1) ** 2 * (l + 1))


EXCEPTIONAL = {
    "M1": {2: 22 / 9, 3: None, 5: 1.03472, 11: 1.09924},      # 3 -> generic (wrong; see dossier)
    "M2": {2: 2.554, 3: 1.332, 5: 1.03472, 11: 1.09924},      # 2 empirical, 3 icosahedral model
}


def density_constant(model: str, lmax: int = 2 * 10 ** 6):
    F = 1.0
    for l in primes_upto(lmax):
        l = int(l)
        g = 1.0 if l == 5 else (vh_split(l) if l % 5 in (1, 4) else vh_inert(l))
        x = EXCEPTIONAL[model].get(l)
        F *= g if x is None else x
    K = 8 * math.sqrt(5) * F / (3 * math.pi ** 2)     # Prob(collision at p) ~ K p^{-1/2}
    return F, K, 16 * math.sqrt(5) * F / (3 * math.pi ** 2)


if __name__ == "__main__":
    N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 200000
    p, u, v, _ = eigen_uv(N)
    print(json.dumps(law_audit(p, u, v), indent=1))
    for mdl in ("M1", "M2"):
        F, K, c = density_constant(mdl, 10 ** 5)
        print(mdl, "F_hat=%.4f K=%.4f c=%.4f" % (F, K, c), "pred N(%d)=%.1f" % (N, K * (p.astype(float) ** -.5).sum()),
              "obs", int((v == 0).sum()))
