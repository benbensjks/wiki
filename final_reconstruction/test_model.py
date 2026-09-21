"""Structural, biochemical, upstream-reference and integration checks."""
import json
import unittest
from dataclasses import replace

import numpy as np
import pandas as pd

from model import Model, Extension, ZENG, STATE_NAMES, N_STATE, REFERENCE_PATH, ROOT, bit_slice


class ReconstructionTests(unittest.TestCase):
    def setUp(self):
        self.m = Model()

    def test_state_contract(self):
        self.assertEqual(N_STATE,43)
        self.assertEqual(len(set(STATE_NAMES)),43)
        self.assertNotIn('E0',STATE_NAMES)
        self.assertNotIn('E1',STATE_NAMES)
        self.assertEqual(STATE_NAMES[39:],('A0','F0','A1','F1'))
        self.assertEqual(self.m.initial_state().shape,(43,))

    def test_source_table(self):
        self.assertEqual(ZENG['k_fwd'],7)
        self.assertEqual(ZENG['k_rev'],5)
        self.assertEqual(ZENG['gamma_int'],(2,1.4,2.2))
        self.assertEqual(ZENG['alpha_Int'],(18,38))

    def test_QSS_expression_and_complex_mapping(self):
        m = self.m
        for J in (.01, 3.5, 6, 18, 38):
            mr = m.transcript_source(J)/m.lm
            pu = m.e.translation_h*mr/m.lu
            self.assertAlmostEqual(m.kmat*pu,J)
            np.testing.assert_allclose(m.expression(mr,pu,J/(2+m.growth),J,2),0,atol=1e-12)
        self.assertAlmostEqual(m.K_complex/m.q,ZENG['K_D_comp'])

    def test_AND_gate_and_autorepression(self):
        m = self.m
        y = m.initial_state()
        y[39:43] = [2,0,2,0]
        y[8] = 0
        g0,g1,_ = m.carry_promoters(y)
        self.assertGreater(g0,0)
        self.assertEqual(g1,0)
        self.assertEqual(m.rhs(0,y)[bit_slice(2)][0],0)
        y[8] = 5
        self.assertGreater(m.rhs(0,y)[bit_slice(2)][0],0)
        y[41] = 0
        self.assertEqual(m.carry_promoters(y)[1],0)
        y[41] = 2
        actual_source = m.rhs(0,y)[41]+(ZENG['gamma_A'][1]+m.growth)*2
        self.assertLess(actual_source,ZENG['alpha_A'][1])

    def test_binding_mass_balance(self):
        m = self.m
        y = m.initial_state()
        for i in range(3):
            b = y[bit_slice(i)]
            b[1],b[2],b[7],b[8],b[9] = .2,.7,.3,.6,.1
        d = m.rhs(0,y)
        for i in range(3):
            b,db = y[bit_slice(i)],d[bit_slice(i)]
            lossC=(m.e.complex_decay_h+m.growth)*b[9]
            self.assertAlmostEqual(db[2]+db[9],m.kmat*b[1]-(ZENG['gamma_int'][i]+m.growth)*b[2]-lossC)
            self.assertAlmostEqual(db[8]+db[9],m.kmat*b[7]-(ZENG['gamma_rdf']+m.growth)*b[8]-lossC)

    def test_invariant_boundaries(self):
        m=self.m
        rng=np.random.default_rng(20260919)
        for _ in range(10):
            y=m.initial_state()+rng.random(N_STATE)
            for i in range(3): y[bit_slice(i)][10]=rng.random()
            for k in range(N_STATE):
                z=y.copy();z[k]=0
                self.assertGreaterEqual(m.rhs(0,z)[k],-1e-12)
            for i in range(3):
                z=y.copy();z[bit_slice(i)][10]=1
                self.assertLessEqual(m.rhs(0,z)[bit_slice(i)][10],1e-12)

    def test_upstream_independent_of_downstream_and_units(self):
        m=self.m;y=m.initial_state();z=y.copy()
        z[7:]=2
        np.testing.assert_allclose(m.rhs(0,y)[:7],m.rhs(0,z)[:7])
        scaled=Model(replace(m.e,uM_per_au=2))
        self.assertAlmostEqual(m.flux(y),scaled.flux(scaled.initial_state()))

    def test_upstream_reference_and_numerical_convergence(self):
        m=self.m
        sol=m.simulate(100)
        strict=m.simulate(100,strict=True)
        self.assertTrue(np.isfinite(sol.y).all())
        self.assertGreaterEqual(sol.y.min(),-1e-9)
        errors=np.max(np.abs(sol.y-strict.y),axis=1)
        np.testing.assert_allclose(sol.y,strict.y,rtol=1e-4,atol=1e-5)
        ref=pd.read_csv(REFERENCE_PATH)
        flux=np.array([m.flux(y) for y in sol.y.T])
        ref_flux=np.interp(sol.t*60,ref.time_min,ref.actual_C31_translation_flux_uM_h)
        max_flux_error=float(np.max(abs(flux-ref_flux)))
        self.assertLess(max_flux_error,.002)
        # Check every retained upstream state, plus shared C31 mRNA.
        ref_cols=['TetR_mRNA','TetR_total','CI_mRNA','CI','LacI_mRNA','LacI','C31_mRNA']
        reconstructed=np.vstack([sol.y[:6],sol.y[6]*m.copies_per_au])
        target=np.array([np.interp(sol.t*60,ref.time_min,ref[col]) for col in ref_cols])
        relative_error=np.max(abs(reconstructed-target),axis=1)/np.maximum(np.max(abs(target),axis=1),1e-12)
        self.assertLess(relative_error.max(),.002)
        out=ROOT/'validation';out.mkdir(exist_ok=True)
        (out/'numerical_validation.json').write_text(json.dumps(dict(
            compared_hours=100,upstream_flux_max_abs_error_uM_h=max_flux_error,
            upstream_state_relative_errors=dict(zip(ref_cols,map(float,relative_error))),
            state_max_abs_convergence_error=dict(zip(STATE_NAMES,map(float,errors))),
            minimum_state=float(sol.y.min()),finite=True),indent=2),encoding='utf-8')


if __name__=='__main__':
    unittest.main(verbosity=2)
