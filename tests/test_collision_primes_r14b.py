"""R14b regression tests (numpy required; PARI cross-check runs only if cypari2 is installed).

python3 -m unittest discover -s tests -p test_collision_primes_r14b.py -v
"""
import unittest
from fractions import Fraction
import numpy as np
from experiments.collision_primes_r14b import (eigen_uv, law_audit, preimage_factor, gl2_over_O_mod_l,
                                               binary_icosahedral_G3, vh_inert, density_constant)

N = 300000


class CollisionPrimesR14b(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p, cls.u, cls.v, cls.B = eigen_uv(N)
        cls.audit = law_audit(cls.p, cls.u, cls.v)

    def test_01_R14_point_count_rows_reproduced(self):
        # (p, trace, v) from R14's independent F_p / F_p^2 point counts; trace = 2u + v
        rows = {3: (0, 2), 5: (-2, 2), 13: (6, 0), 19: (-4, 0), 29: (-6, 0), 31: (0, 6), 43: (0, 0), 97: (22, 6)}
        d = dict(zip(self.p.tolist(), zip(self.u.tolist(), self.v.tolist())))
        for q, (t, va) in rows.items():
            u, v = d[q]
            self.assertEqual(2 * u + v, t); self.assertEqual(abs(v), va)

    def test_02_weight_one_form_is_eta_product(self):
        # B = eta(z) eta(23z) = q - q^2 - q^3 + q^6 + q^8 - q^13 - ...
        self.assertEqual(self.B[:14].tolist(), [0, 1, -1, -1, 0, 0, 1, 0, 1, 0, 0, 0, 0, -1])

    def test_03_first_collisions(self):
        col = self.p[self.v == 0][:8].tolist()
        self.assertEqual(col, [13, 19, 29, 43, 109, 223, 229, 271])

    def test_04_theorem_and_lemmas_zero_exceptions(self):
        a = self.audit
        for key in ("L1_v_even_violations", "L2_eisenstein_P11_violations", "T_mod11_violations",
                    "T_parity_violations", "R1_violations", "R2_violations", "icosahedral_mod3_violations"):
            self.assertEqual(a[key], 0, key)
        self.assertGreater(a["collisions"], 150)

    def test_05_a_p_zero_needs_p_10_mod_11(self):
        z = self.audit["a_p_zero_primes"]
        self.assertEqual(z, [43, 109, 1033])
        self.assertTrue(all(q % 11 == 10 for q in z))

    def test_06_generic_local_factor_matches_van_hirtum_limit(self):
        full3 = gl2_over_O_mod_l(3)
        f1, _ = preimage_factor(full3, 3, 1); f2, _ = preimage_factor(full3, 3, 2)
        self.assertEqual((f1, f2), (Fraction(87, 80), Fraction(263, 240)))
        self.assertLess(abs(float(f2) - vh_inert(3)), 2e-3)

    def test_07_two_adic_and_icosahedral_factors(self):
        G2 = gl2_over_O_mod_l(2, pred=lambda e: e[1] == 0)          # mod-2 image inside GL2(F2)
        self.assertEqual(preimage_factor(G2, 2, 1)[0], 2)
        self.assertEqual(preimage_factor(G2, 2, 2)[0], Fraction(7, 3))
        G3 = binary_icosahedral_G3()
        self.assertEqual(len(G3), 240)
        f1, c1 = preimage_factor(G3, 3, 1, cond_mod=3)
        self.assertEqual(f1, Fraction(51, 40)); self.assertEqual(c1, {1: Fraction(3, 5), 2: Fraction(1, 4)})
        self.assertEqual(preimage_factor(G3, 3, 2)[0], Fraction(53, 40))

    def test_08_empirical_mod3_frequencies_match_2A5(self):
        m = (self.p % 3 == 1)
        f00 = float(((self.u % 3 == 0) & (self.v % 3 == 0))[m].mean())
        self.assertLess(abs(f00 - 0.25), 0.01)
        self.assertLess(abs(float((self.v % 3 == 0).mean()) * 3 - 1.275), 0.02)

    def test_09_density_constants(self):
        F1, K1, c1 = density_constant("M1", 10 ** 5); F2, K2, c2 = density_constant("M2", 10 ** 5)
        self.assertAlmostEqual(F1, 3.1729, delta=2e-3); self.assertAlmostEqual(c2, 4.864, delta=5e-3)

    def test_10_optional_pari_crosscheck(self):
        try:
            import cypari2
        except ImportError:
            self.skipTest("cypari2 not installed")
        pari = cypari2.Pari(); pari.allocatemem(4 * 10 ** 8)
        co = pari('mfcoefs(mfeigenbasis(mfinit([23,2],0))[1],1000)')
        d = dict(zip(self.p.tolist(), zip(self.u.tolist(), self.v.tolist())))
        for q in [x for x in d if x <= 1000]:
            c = pari.lift(co[q]); self.assertEqual((int(c.polcoef(0)), int(c.polcoef(1))), d[q])


if __name__ == "__main__":
    unittest.main()
