"""State contract and no-feedback projection tests for ThreeBit51Model."""
from __future__ import annotations

import copy
import json
import unittest
from dataclasses import replace

import numpy as np

from model import BIT_NAMES, ROOT, ZENG
from model_threebit51 import (N_STATE_51, STATE_NAMES_51, ThreeBit51Model,
                              selected_threebit_extension)
from model_twobit34 import (STATE_NAMES_34, CarryExpressionParameters,
                            TwoBit34Model, nominal_extension)


def frozen_twobit_model():
    profile = json.loads((ROOT/'twobit34_results'/'selected_profile.json').read_text(encoding='utf-8'))
    extension = replace(nominal_extension(), **profile['extension'])
    carry = CarryExpressionParameters(**profile['carry_expression'])
    return TwoBit34Model(extension, carry)


class ThreeBit51Tests(unittest.TestCase):
    def setUp(self):
        self.zeng = copy.deepcopy(ZENG)
        self.model = ThreeBit51Model()

    def tearDown(self):
        self.assertEqual(ZENG, self.zeng)

    def test_state_contract(self):
        self.assertEqual(N_STATE_51, 51)
        self.assertEqual(len(STATE_NAMES_51), len(set(STATE_NAMES_51)))
        self.assertEqual(STATE_NAMES_51[:34], STATE_NAMES_34)
        for level in range(2):
            names = STATE_NAMES_51[6+11*level:17+11*level]
            self.assertEqual(tuple(n.removeprefix(f'b{level}_') for n in names), BIT_NAMES)
        self.assertEqual(tuple(n.removeprefix('b2_') for n in STATE_NAMES_51[34:45]), BIT_NAMES)
        self.assertEqual(STATE_NAMES_51[45:],('A1','F1','M_A1','A1_u','M_F1','F1_u'))
        self.assertFalse(any(name in ('E0','E1','U') for name in STATE_NAMES_51))

    def test_confirmed_clock_gate(self):
        e = selected_threebit_extension()
        self.assertEqual(e.clock_K_au, 0.4)
        self.assertEqual(e.clock_n, 3.0)
        self.assertFalse(e.add_growth)

    def test_initial_state_and_rhs(self):
        y = self.model.initial_state()
        self.assertEqual(y.shape, (51,))
        self.assertTrue(np.isfinite(y).all())
        self.assertGreaterEqual(y.min(), 0.0)
        d = self.model.rhs(0.0, y)
        self.assertEqual(d.shape, (51,))
        self.assertTrue(np.isfinite(d).all())
        cold = self.model.initial_state(cold=True)
        self.assertTrue(np.all(cold[34:51] == 0.0))

    def test_lower_rhs_projection_is_exact(self):
        two = frozen_twobit_model()
        y51 = self.model.initial_state()
        y34 = self.model.project_twobit34(y51)
        np.testing.assert_allclose(y34, two.initial_state(), rtol=0, atol=1e-14)
        d51 = self.model.rhs(0.0, y51)
        d34 = two.rhs(0.0, y34)
        np.testing.assert_allclose(self.model.project_twobit34_derivative(d51), d34,
                                   rtol=0, atol=1e-13)

    def test_short_trajectory_preserves_lower_subsystem(self):
        two = frozen_twobit_model()
        a = self.model.simulate(hours=4.0, sample_min=2.0, max_step_min=2.0)
        b = two.simulate(hours=4.0, sample_min=2.0, max_step_min=2.0)
        np.testing.assert_allclose(self.model.project_twobit34(a.y), b.y,
                                   rtol=1e-7, atol=1e-9)

    def test_no_feedback_jacobian_block(self):
        rng = np.random.default_rng(20260920)
        for _ in range(3):
            y = self.model.initial_state()
            y[34:] = np.maximum(y[34:] * rng.uniform(.5, 1.5, 17), 1e-6)
            for j in range(34, 51):
                eps = 1e-6 * max(1.0, abs(y[j]))
                yp, ym = y.copy(), y.copy(); yp[j] += eps; ym[j] -= eps
                column = (self.model.rhs(0.0, yp)[:34] - self.model.rhs(0.0, ym)[:34])/(2*eps)
                self.assertLessEqual(float(np.max(np.abs(column))), 1e-12)


if __name__ == '__main__':
    unittest.main(verbosity=2)
