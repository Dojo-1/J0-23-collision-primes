"""Exact local factors F_l^(k) = l^k * Prob_{g in G_l mod l^k}[ trace(g) in Z/l^k ]
for model l-adic images G_l of the Galois representation on T_l J_0(23), O = Z[phi], phi^2 = phi + 1.

Elements of O/l^k are pairs (x, y) = x + y*phi.  'trace in Z' <=> y == 0 (mod l^k).
G_l = { g in GL_2(O/l^k) : det g in (Z/l^k)^x , entries mod l in allowed residue sets }.
Counting: sum over diagonal (a, e); number of off-diagonal (b, c) with bc = ae - D via table M.
"""
import itertools, sys
from fractions import Fraction


def ring(l, k):
    m = l ** k
    els = [(x, y) for x in range(m) for y in range(m)]
    idx = {e: i for i, e in enumerate(els)}
    def mul(u, v):
        x1, y1 = u; x2, y2 = v
        return ((x1 * x2 + y1 * y2) % m, (x1 * y2 + x2 * y1 + y1 * y2) % m)
    return m, els, idx, mul


def factor(l, k, allowed=lambda e: True, diag_allowed=None, offdiag_allowed=None):
    """allowed(e): predicate on ring element e (applied to all four entries) unless the
    diag/offdiag predicates are given."""
    da = diag_allowed or allowed
    oa = offdiag_allowed or allowed
    m, els, idx, mul = ring(l, k)
    n = len(els)
    # M[x] = #{(b,c) allowed : b*c = x}
    M = [0] * n
    offs = [e for e in els if oa(e)]
    for b in offs:
        for c in offs:
            M[idx[mul(b, c)]] += 1
    diags = [e for e in els if da(e)]
    units = [d for d in range(m) if d % l != 0]
    tot = 0; good = 0
    for a in diags:
        for e in diags:
            ae = mul(a, e)
            t_int = ((a[1] + e[1]) % m == 0)
            s = 0
            for D in units:
                x = ((ae[0] - D) % m, ae[1])
                s += M[idx[x]]
            tot += s
            if t_int:
                good += s
    return Fraction(l ** k * good, tot)


def split_component_counts(l, k, diag_pred, off_pred):
    """Return N[(t, D)] for matrices over Z/l^k with entry constraints, det D unit."""
    m = l ** k
    M = [0] * m
    offs = [x for x in range(m) if off_pred(x)]
    for b in offs:
        for c in offs:
            M[(b * c) % m] += 1
    N = {}
    for a in range(m):
        if not diag_pred(a, 0):
            continue
        for e in range(m):
            if not diag_pred(e, 1):
                continue
            for D in range(m):
                if D % l == 0:
                    continue
                cnt = M[(a * e - D) % m]
                if cnt:
                    N[(a + e) % m, D] = N.get(((a + e) % m, D), 0) + cnt
    return N


def split_factor(l, k, g_diag, g_off, h_diag=lambda x, i: True, h_off=lambda x: True):
    m = l ** k
    Ng = split_component_counts(l, k, g_diag, g_off)
    Nh = split_component_counts(l, k, h_diag, h_off)
    num = 0; den = 0
    for D in range(m):
        if D % l == 0:
            continue
        g_tot = sum(v for (t, d), v in Ng.items() if d == D)
        h_tot = sum(v for (t, d), v in Nh.items() if d == D)
        den += g_tot * h_tot
        num += sum(Ng.get((t, D), 0) * Nh.get((t, D), 0) for t in range(m))
    return Fraction(l ** k * num, den)


def vh_inert(l):
    return Fraction(l**3 * (l*l + l + 1), (l + 1) * (l**4 - 1))


def vh_split(l):
    return Fraction(l*l * (l**3 + l*l - l - 2), (l*l - 1)**2 * (l + 1))


if __name__ == "__main__":
    out = {}
    # --- generic checks of Van Hirtum closed forms (independent enumeration) ---
    for k in (1, 2):
        f = factor(3, k)
        print("l=3 inert full, k=%d: %s = %.6f   VH closed form %.6f" % (k, f, float(f), float(vh_inert(3))), flush=True)
    for k in (1, 2):
        f = split_factor(3, k, lambda x, i: True, lambda x: True)  # l=3 is not split in Q(sqrt5); pure check of split formula
        print("split-model l=3 full, k=%d: %.6f  VH %.6f" % (k, float(f), float(vh_split(3))), flush=True)
    # --- l = 2: mod-2 image inside GL2(F2) (entries with y even) ---
    F2 = lambda e: e[1] % 2 == 0
    for k in (1, 2, 3, 4):
        f = factor(2, k, allowed=F2)
        print("l=2 model Gamma2 (g mod 2 in GL2(F2)), k=%d: %s = %.6f" % (k, f, float(f)), flush=True)
    for k in (1, 2, 3):
        f = factor(2, k)
        print("l=2 FULL image (control), k=%d: %.6f" % (k, float(f)), flush=True)
    # --- l = 5 ramified, full image ---
    for k in (1, 2):
        f = factor(5, k)
        print("l=5 ramified full, k=%d: %s = %.6f" % (k, f, float(f)), flush=True)
    # --- l = 11 split, Eisenstein at one prime ---
    l = 11
    diagE = lambda x, i: (x % l == 1) if i == 0 else True   # g = diag(1,*) mod 11 (+ off-diagonal 0 mod 11)
    offE = lambda x: x % l == 0
    borelE_off_b = None
    for k in (1, 2):
        f = split_factor(l, k, diagE, offE)
        print("l=11 split, p-component = preimage of diag(1,*), k=%d: %.6f" % (k, float(f)), flush=True)
        g = split_factor(l, k, lambda x, i: True, lambda x: True)
        print("l=11 split FULL (control), k=%d: %.6f   VH %.6f" % (k, float(g), float(vh_split(l))), flush=True)
