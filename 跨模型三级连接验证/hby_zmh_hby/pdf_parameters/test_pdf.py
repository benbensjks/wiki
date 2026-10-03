import unittest
import numpy as np
from model_pdf import PDFModel,PDF_VALUES,RATE_KEYS,HZHModel

class Tests(unittest.TestCase):
    def test_table_and_units(self):
        m=PDFModel('square')
        self.assertEqual(m.audit()['matching_count'],21)
        for k,v in PDF_VALUES.items():self.assertAlmostEqual(m.z[k]*(60 if k in RATE_KEYS else 1),v,places=13)
        self.assertEqual(m.z['n_int1'],4)

    def test_only_middle_changed(self):
        old=HZHModel('square');new=PDFModel('square')
        y=old.initial_state();a=old.rhs(8,y);b=new.rhs(8,y)
        np.testing.assert_array_equal(a[:11],b[:11]);np.testing.assert_array_equal(a[17:],b[17:])
        self.assertGreater(np.max(np.abs(a[11:17]-b[11:17])),.1)

    def test_pre_equilibria(self):
        m=PDFModel('square');y=m.initial_state();d=m.middle_rhs(0,y[11:17])
        np.testing.assert_allclose(d[[0,1,4]],0,atol=1e-13)
        self.assertAlmostEqual(y[15],3.5/.6)

if __name__=='__main__':unittest.main(verbosity=2)
