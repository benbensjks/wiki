"""Check Han mappings without changing the original donor or receiver."""
import unittest
import numpy as np
from run_han_comparison import HanInput,HanDonor


class HanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.han=HanInput(50.)
        cls.norm=HanDonor(cls.han,'normalized_concentration')
        cls.flux=HanDonor(cls.han,'absolute_translation_flux')

    def test_normalized_mapping(self):
        for t in (0,100,1500,2999):
            expected=6*np.clip(np.interp(t,self.han.t,self.han.sol.y[7])/self.han.c31_max,0,2)
            self.assertAlmostEqual(self.norm.source_min(t),expected,places=13)

    def test_absolute_flux_no_extra_gain(self):
        y=self.flux.initial_state();y[1]=.3
        for t in (100,1500,2500):
            j=np.interp(t,self.han.t,self.han.flux_copies_min)/(602.214076*5.75)
            self.assertAlmostEqual(self.flux.rhs_min(t,y)[1],j-2*y[1],places=13)

    def test_only_source_changes_and_time_units(self):
        y=self.norm.initial_state();y[1]=.5
        a=self.norm.rhs_min(800,y);b=self.flux.rhs_min(800,y)
        np.testing.assert_array_equal(np.delete(a,1),np.delete(b,1))
        np.testing.assert_array_equal(self.flux.rhs(800/60,y),60*self.flux.rhs_min(800,y))

    def test_no_extrapolation(self):
        with self.assertRaises(ValueError):self.norm.input_min(3001)


if __name__=='__main__':unittest.main(verbosity=2)
