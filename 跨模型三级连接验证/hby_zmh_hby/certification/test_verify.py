import unittest
import numpy as np
from verify_hzh import associate,read_verdict,margins,segments


class VerifierTests(unittest.TestCase):
    def test_no_vacuous_causal_pass(self):
        x=associate([],[],[],100)
        self.assertFalse(x['passed']);self.assertFalse(x['alternating_directions'])
        self.assertFalse(associate([dict(start_h=1.,peak_h=1.5)],[],[],100)['passed'])

    def test_unassigned_and_order_fail(self):
        rev=[dict(start_h=k*10.,peak_h=k*10+.5) for k in range(1,6)]
        gates=[dict(start_h=k*10-.5,peak_h=k*10+1.) for k in range(1,6)]
        flips=[dict(time_h=k*10+2.,direction='up' if k%2 else 'down') for k in range(1,6)]
        x=associate(rev,gates,flips,100)
        self.assertTrue(x['one_to_one']);self.assertFalse(x['causal_order']);self.assertFalse(x['passed'])
        gates.append(dict(start_h=82.,peak_h=83.))
        self.assertFalse(associate(rev,gates,flips,100)['one_to_one'])

    def test_readout_requires_full_cycles(self):
        rr=[dict(value=i%8,covered=True,commitment=[1,1,1]) for i in range(12)]
        self.assertFalse(read_verdict(rr,8)['passed'])
        rr=[dict(value=i%8,covered=True,commitment=[1,1,1]) for i in range(24)]
        self.assertTrue(read_verdict(rr,8)['passed'])
        rr[-1]['value']=None
        self.assertFalse(read_verdict(rr,8)['passed'])

    def test_pulse_duration_and_dose(self):
        t=np.linspace(0,1,61);s=np.where((t>.2)&(t<.4),.2,0.)
        self.assertEqual(len(segments(t,s,.1,min_duration_h=.3)),0)
        self.assertEqual(len(segments(t,s,.1,min_duration_h=.1,min_dose=.1)),0)
        self.assertEqual(len(segments(t,s,.1,min_duration_h=.1,min_dose=.02)),1)

    def test_margins_without_events(self):
        m=margins([dict(start_h=1,end_h=2)],[])
        self.assertIsNone(m['min_setup_h']);self.assertIsNone(m['min_hold_h'])


if __name__=='__main__':unittest.main(verbosity=2)
