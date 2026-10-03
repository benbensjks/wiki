import unittest,sys,json
import numpy as np
from model_hzh import HZHModel,NAMES,IDX,PARENT
from hybrid_model import WIKI,donor_namespace
from run_han_comparison import HanInput


class Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.han=HanInput(50.)
        cls.model=HZHModel('han',cls.han)
        cls.square=HZHModel('square')
        sys.path.insert(0,str(WIKI/'final_reconstruction'))
        from model import Model,Extension
        prof=json.loads((WIKI/'final_reconstruction/twobit34_results/selected_profile.json').read_text(encoding='utf-8'))
        cls.ref=Model(Extension(**prof['extension']))

    def test_han_bit0_reference_rhs(self):
        rng=np.random.default_rng(27)
        for j in (0,600,1200,1800):
            old=self.ref.initial_state();old[:6]=self.han.sol.y[:6,j]
            b=rng.uniform(.01,2,11);b[10]=rng.uniform(0,1);old[6:17]=b
            expected=self.ref.rhs(self.han.t[j]/60,old)[6:17]
            actual=self.model.bit0_rhs(self.han.t[j]/60,b)
            np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-12)

    def test_middle_exact_source_rhs(self):
        ns=donor_namespace();ns['get_u_in']=lambda t:0.
        rng=np.random.default_rng(20)
        for _ in range(20):
            s=rng.uniform(0,1);middle=rng.uniform(.1,2,6);middle[2]=rng.uniform(0,1)
            y=np.r_[1-s,.5,1.,.2,middle]
            expected=60*np.asarray(ns['counter_2bit_ode'](0,y,ns['p']))[4:10]
            np.testing.assert_array_equal(self.square.middle_rhs(s,middle),expected)

    def test_structure_and_initial(self):
        self.assertEqual(len(set(NAMES)),34)
        y=self.model.initial_state();self.assertEqual(y.shape,(34,))
        np.testing.assert_allclose(y[:11],self.ref.initial_state()[6:17],atol=1e-13)
        before=self.model.rhs(0,y)
        y[17:]=y[17:]*3+.1
        np.testing.assert_array_equal(before[:17],self.model.rhs(0,y)[:17])
        y[11:]=y[11:]*2
        np.testing.assert_array_equal(before[:11],self.model.rhs(0,y)[:11])

    def test_square_spec(self):
        self.assertEqual(self.square.input_source(4.9),0)
        self.assertEqual(self.square.input_source(5),6)
        self.assertEqual(self.square.input_source(6.8),0)
        self.assertEqual(self.square.input_source(15),6)


if __name__=='__main__':unittest.main(verbosity=2)
