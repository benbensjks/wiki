"""Structural, unit, reference-equivalence and bounded smoke checks."""
import json
import sys
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from hybrid_model import (WIKI,STATE_NAMES,TAIL_NAMES,DONOR_NAMES,INDEX,
                          ZmhDonor,HbyReceiver,HybridModel,integrate,time_grid)
from diagnostics import analyse


class HybridTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.donor=ZmhDonor()
        cls.receiver=HbyReceiver()
        cls.model=HybridModel(cls.donor,cls.receiver)
        sys.path.insert(0,str(WIKI/'final_reconstruction'))
        from model import Extension
        from model_threebit51 import ThreeBit51Model,ThreeBitCarryParameters
        profile=json.loads((WIKI/'final_reconstruction/twobit34_results/selected_profile.json').read_text(encoding='utf-8'))
        cls.reference=ThreeBit51Model(extension=Extension(**profile['extension']),
                                      carry=ThreeBitCarryParameters(),n_A1_gate=6.)
        cls.ref_ids=[cls.reference.state_names.index(n) for n in TAIL_NAMES]

    def test_layout(self):
        self.assertEqual(len(STATE_NAMES),27)
        self.assertEqual(len(set(STATE_NAMES)),27)
        self.assertEqual(INDEX['b2_S'],26)
        grid=time_grid(2000/60,1.)
        self.assertEqual(len(grid),2001)
        self.assertTrue(np.all(np.diff(grid)>0))
        self.assertEqual(grid[-1],2000/60)

    def test_initial_reference(self):
        expected=self.reference.initial_state()[self.ref_ids]
        np.testing.assert_allclose(self.receiver.initial_state(),expected,rtol=1e-13,atol=1e-13)

    def test_receiver_rhs_reference(self):
        rng=np.random.default_rng(20260927)
        for _ in range(32):
            y=self.reference.initial_state().copy()
            y[self.ref_ids]=rng.uniform(.01,4,len(self.ref_ids))
            y[44]=rng.uniform(0,1); y[27]=rng.uniform(0,1); y[8]=rng.uniform(0,3)
            actual=self.receiver.rhs(1.,y[self.ref_ids],1-y[27],y[8])
            expected=self.reference.rhs(1.,y)[self.ref_ids]
            np.testing.assert_allclose(actual,expected,rtol=2e-12,atol=2e-12)
            sig=self.reference.diagnostic_signals(y)
            ours=self.receiver.signals(y[self.ref_ids],y[8])
            self.assertAlmostEqual(ours['g1'],float(sig['g1'][0]),places=13)

    def test_no_feedback(self):
        y=self.model.initial_state()
        reference=self.donor.rhs(1.,y[:10])
        for scale in (0,1,10,1000):
            z=y.copy();z[10:]=z[10:]*scale+.1
            np.testing.assert_array_equal(self.model.rhs(1.,z)[:10],reference)

    def test_hour_conversion_by_integration(self):
        t=np.linspace(0,1,61)
        minute=solve_ivp(self.donor.rhs_min,(0,60),self.donor.initial_state(),
                         t_eval=60*t,rtol=1e-9,atol=1e-11,method='DOP853',max_step=.5)
        hour=integrate(self.donor.rhs,self.donor.initial_state(),t,rtol=1e-9,atol=1e-11,max_step_min=.5)
        self.assertTrue(minute.success)
        np.testing.assert_allclose(hour.y,minute.y,rtol=2e-7,atol=2e-7)

    def test_mass_balance(self):
        tail=self.receiver.initial_state();tail[8]=1.2;tail[14]=2.3;tail[15]=.4
        d=self.receiver.rhs(0,tail,1,.8)
        p,c=self.receiver.p,self.receiver.c
        self.assertAlmostEqual(d[8]+d[15],self.receiver.kmb*tail[7]-p['gamma_int'][2]*tail[8]-c.complex_decay*tail[15])
        self.assertAlmostEqual(d[14]+d[15],self.receiver.kmb*tail[13]-p['gamma_rdf']*tail[14]-c.complex_decay*tail[15])

    def test_input_matches_original_formula_and_is_bounded(self):
        for t in (0,50,1000,2000,3000):
            expected=float(np.clip(np.interp(t,self.donor.upstream.t,self.donor.upstream.y[4])/self.donor.c31_max,0,2))
            self.assertEqual(self.donor.input_min(t),expected)
        with self.assertRaises(ValueError): self.donor.input_min(3001)
        with self.assertRaises(ValueError): time_grid(51,1)

    def test_short_run_and_conservative_readout(self):
        ts=time_grid(.5,1)
        sol=integrate(self.model.rhs,self.model.initial_state(),ts)
        self.assertEqual(sol.y.shape,(27,31))
        report=analyse(ts,sol.y,sol.y[1])
        self.assertEqual(report['status'],'EXPLORATORY_ONLY_NOT_CERTIFIED')
        self.assertFalse(report['steady_mod8']['enough_reads'])
        self.assertIsNone(report['causal_certification'])


if __name__=='__main__': unittest.main(verbosity=2)
