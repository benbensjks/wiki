"""Fast state-contract tests for the formal 34-state model."""
from __future__ import annotations

import copy
import unittest

import numpy as np

from model import BIT_NAMES, ZENG
from model_twobit34 import (N_STATE_34, STATE_NAMES_34,
                            CarryExpressionParameters, TwoBit34Model,
                            nominal_extension)


class TwoBit34Tests(unittest.TestCase):
    def setUp(self):
        self.zeng = copy.deepcopy(ZENG)
        self.model = TwoBit34Model()

    def tearDown(self):
        self.assertEqual(ZENG, self.zeng)

    def test_state_contract(self):
        self.assertEqual(N_STATE_34, 34)
        self.assertEqual(len(STATE_NAMES_34), len(set(STATE_NAMES_34)))
        for i in range(2):
            self.assertEqual(tuple(n.removeprefix(f'b{i}_') for n in STATE_NAMES_34[6 + 11*i:17 + 11*i]),
                             BIT_NAMES)
        self.assertEqual(STATE_NAMES_34[28:],
                         ('A0', 'F0', 'M_A0', 'A0_u', 'M_F0', 'F0_u'))
        forbidden = ('E0', 'E1', ' U', 'b0_U', 'b1_U')
        self.assertFalse(any(any(x in name for x in forbidden) for name in STATE_NAMES_34))

    def test_nominal_profile(self):
        e = nominal_extension()
        self.assertEqual(e.uM_per_au, 6.0)
        self.assertEqual(e.maturation_half_life_min, 20.0)
        self.assertEqual(e.complex_on_au_inv_h, 0.1)
        self.assertEqual(e.complex_off_h, 1.0)
        self.assertFalse(e.add_growth)
        c = CarryExpressionParameters()
        self.assertEqual(c.mrna_half_life_min, 2.0)
        self.assertEqual(c.activator_maturation_half_life_min, 30.0)
        self.assertEqual(c.repressor_maturation_half_life_min, 30.0)

    def test_initial_state_and_rhs(self):
        y = self.model.initial_state()
        self.assertEqual(y.shape, (34,))
        self.assertTrue(np.isfinite(y).all())
        self.assertGreaterEqual(y.min(), 0.0)
        d = self.model.rhs(0.0, y)
        self.assertEqual(d.shape, (34,))
        self.assertTrue(np.isfinite(d).all())
        cold = self.model.initial_state(cold=True)
        self.assertTrue(np.all(cold[28:34] == 0.0))

    def test_upstream_is_unloaded_reference(self):
        y = self.model.initial_state()
        d34 = self.model.rhs(0.0, y)
        dbase = self.model.base.rhs_twobit(0.0, y[:30])
        np.testing.assert_allclose(d34[:6], dbase[:6], rtol=0, atol=0)

    def test_short_integration_is_finite(self):
        sol = self.model.simulate(hours=4.0, sample_min=2.0, max_step_min=2.0)
        self.assertEqual(sol.y.shape[0], 34)
        self.assertTrue(np.isfinite(sol.y).all())
        self.assertGreaterEqual(sol.y.min(), -1e-8)


if __name__ == '__main__':
    unittest.main(verbosity=2)
