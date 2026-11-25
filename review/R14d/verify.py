#!/usr/bin/env python3
"""Exact R14d finite checks; Python standard library, no downloads or input data."""
from fractions import Fraction
from itertools import product
import json


class QuadraticField:
    # u^2 = linear*u + constant; encode a+b*u as a+p*b.
    def __init__(self, p, linear, constant):
        assert all((x*x-linear*x-constant) % p for x in range(p))
        self.p, self.linear, self.constant = p, linear, constant
        self.elements = range(p*p)
        self.add = [[self.encode(x % p + y % p, x//p + y//p)
                     for y in self.elements] for x in self.elements]
        self.mul = [[self.encode((x % p)*(y % p) + constant*(x//p)*(y//p),
                                 (x % p)*(y//p)+(x//p)*(y % p)+linear*(x//p)*(y//p))
                     for y in self.elements] for x in self.elements]
        self.neg = [self.encode(-x % p, -(x//p)) for x in self.elements]

    def encode(self, a, b):
        return a % self.p + self.p*(b % self.p)

    def poly(self, coeffs, x):
        y = 0
        for c in reversed(coeffs):
            y = self.add[self.mul[y][x]][c % self.p]
        return y


def point_counts():
    f = QuadraticField(13, 0, 2)
    # The modular curve model and a separately tabulated isogenous Jacobian.
    models = {
        "modular_curve": ((-1, -1, 0, -1), (-2, 2, -3, 0, 0, -2)),
        "isogenous_jacobian": ((1, 1, 0, 1), (0, 0, 0, 0, 0, -1)),
    }
    out = {}
    for label, (h, rhs) in models.items():
        # Good reduction in odd characteristic: squarefree discriminant h^2+4*rhs.
        d = [sum(h[i]*h[j] for i, j in product(range(4), repeat=2) if i+j == k)
             + (4*rhs[k] if k < len(rhs) else 0) for k in range(7)]
        a, b = [v % 13 for v in d], [k*d[k] % 13 for k in range(1, 7)]
        while b:
            while a and a[-1] == 0:
                a.pop()
            while b and b[-1] == 0:
                b.pop()
            if not b:
                break
            while len(a) >= len(b):
                i = len(a)-len(b)
                c = a[-1]*pow(b[-1], -1, 13) % 13
                for j, bj in enumerate(b):
                    a[i+j] = (a[i+j]-c*bj) % 13
                while a and a[-1] == 0:
                    a.pop()
            a, b = b, a
        assert len(a) == 1
        counts = []
        for field in (range(13), f.elements):
            total = 2  # v^2+h_3*v=0 at u=1/x=0 has two distinct roots.
            for x in field:
                hx, rx = f.poly(h, x), f.poly(rhs, x)
                total += sum(f.add[f.mul[y][y]][f.mul[hx][y]] == rx for y in field)
            counts.append(total)
        assert counts == [8, 204]
        out[label] = counts
    s1, s2 = 14-8, 170-204
    # Sum of squares of the two GL2 traces is s2 + 4p.
    assert 2*(s2+4*13)-s1*s1 == 0
    assert s1//2 == 3
    return out


def rank(rows):
    pivots = {}
    for row in rows:
        row = list(row)
        for i in range(len(row)):
            v = row[i]
            if v == 0:
                continue
            if i not in pivots:
                pivots[i] = [(v*w) % 3 for w in row]
                break
            row = [(x-v*y) % 3 for x, y in zip(row, pivots[i])]
    return len(pivots)


def residual_checks():
    f = QuadraticField(3, 1, 1)
    add, mul, neg = f.add, f.mul, f.neg
    I, A = (1, 0, 0, 1), (0, 2, 1, 0)

    def mm(x, y):
        a, b, c, d = x
        e, g, h, k = y
        return (add[mul[a][e]][mul[b][h]], add[mul[a][g]][mul[b][k]],
                add[mul[c][e]][mul[d][h]], add[mul[c][g]][mul[d][k]])

    def det(x):
        return add[mul[x[0]][x[3]]][neg[mul[x[1]][x[2]]]]

    def tr(x):
        return add[x[0]][x[3]]

    def inv(x):  # determinant one
        return (x[3], neg[x[1]], neg[x[2]], x[0])

    def power(x, n):
        y = I
        for _ in range(n):
            y = mm(y, x)
        return y

    def closure(b):
        seen, todo = {I}, [I]
        while todo:
            x = todo.pop()
            for y in (mm(x, A), mm(x, b)):
                if y not in seen:
                    seen.add(y)
                    todo.append(y)
        return frozenset(seen)

    def proj_identity(x):
        return x[1] == x[2] == 0 and x[0] == x[3]

    def coords(x):
        return [z for v in x for z in (v % 3, v//3)]

    def difference_rows(group):
        for h in sorted(group):
            ha = mm(mm(h, A), inv(h))
            delta = tuple(add[x][neg[y]] for x, y in zip(ha, A))
            assert tr(delta) == 0
            yield coords(delta)

    sl = [x for x in product(f.elements, repeat=4) if det(x) == 1]
    assert len(sl) == 720
    eligible = [b for b in sl if not proj_identity(b)
                and proj_identity(power(b, 3))
                and not proj_identity(mm(A, b))
                and proj_identity(power(mm(A, b), 5))]
    groups = {closure(b) for b in eligible}
    assert groups
    i = next(z for z in f.elements if mul[z][z] == 2)
    rank_pairs, numbers = set(), set()
    for group in groups:
        assert len(group) == 120
        diffs = list(difference_rows(group))
        u = tuple(add[x][2 if j in (0, 3) else 0] for j, x in enumerate(A))
        assert tr(u) == 1
        pair = (rank(diffs), rank(diffs + [coords(u)]))
        assert pair == (6, 7)
        rank_pairs.add(pair)
        full = group | {tuple(mul[i][z] for z in x) for x in group}
        rational = [x for x in full if tr(x) < 3]
        n = (len(full), len(rational), sum(map(proj_identity, rational)))
        assert n == (240, 102, 2)
        numbers.add(n)
    assert Fraction(100, 240)*3 + Fraction(2, 240)*9*Fraction(13, 12) == Fraction(213, 160)
    return {
        "determinant_one_matrices": len(sl),
        "eligible_second_generators": len(eligible),
        "distinct_binary_groups_with_fixed_involution": len(groups),
        "adjoint_and_full_ranks": sorted(rank_pairs),
        "group_order_rational_trace_and_scalar_rational_trace": sorted(numbers),
        "local_factor_3": "213/160",
    }


if __name__ == "__main__":
    print(json.dumps({"point_counts": point_counts(), "residual": residual_checks()}, indent=2))
