"""PDF-authoritative ZMH middle block; preserve the two HBY end modules."""
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from model_hzh import HZHModel,NAMES,IDX,PARENT
from hybrid_model import WIKI,sha256,require_hash

PDF_PATH=WIKI/'其他小组成员任务/前馈环脉冲进位级联模型参数分析 (1).pdf'
PDF_SHA='4EF770A748C77C977887B62EDDC1209CA87C070BDAD17C7EFB5A94D84A26473D'
# Values transcribed from PDF pages 1-2, verified visually.
# Internal middle_rhs is per-minute then multiplied by 60. Rate entries below
# are per-HOUR and must therefore be divided by 60 on assignment.
PDF_VALUES={
    'k_fwd':7., 'k_rev':5., 'K_D_comp':1.2, 'Kinh':.1,
    'K_D_int1':1., 'alpha_rep1':3.5, 'gamma_rep1':.6,
    'alpha_rdf1':6., 'gamma_rdf1':.8, 'K_rep':.85, 'n':3.9,
    'alpha_A0':8., 'gamma_A0':1.9, 'K_A0':.5, 'n_A0':2.,
    'alpha_R0':3.8, 'gamma_R0':.6, 'K_R0':.6, 'n_R0':4.,
    'alpha_Int1':18., 'gamma_int1':1.4,
}
RATE_KEYS={'k_fwd','k_rev','alpha_rep1','gamma_rep1','alpha_rdf1','gamma_rdf1',
           'alpha_A0','gamma_A0','alpha_R0','gamma_R0','alpha_Int1','gamma_int1'}


class PDFModel(HZHModel):
    def __init__(self,mode,han=None):
        require_hash(PDF_PATH,PDF_SHA)
        super().__init__(mode,han)
        self.old_z=self.z.copy()
        self.z.update({k:v/60 if k in RATE_KEYS else v for k,v in PDF_VALUES.items()})
        # n_int1=4 is absent from the parameter PDF. Preserve the actual code's
        # unspecified equation exponent, rather than inventing a PDF value.
        assert self.z['n_int1']==4.

    def audit(self):
        rows=[]
        for k,v in PDF_VALUES.items():
            scale=60 if k in RATE_KEYS else 1
            old=self.old_z[k]*scale;new=self.z[k]*scale
            rows.append(dict(parameter=k,old_code_per_hour_or_dimensionless=old,
                             pdf=v,used_after_conversion=new,
                             changed=not np.isclose(old,v,rtol=0,atol=1e-12),
                             matches_pdf=bool(np.isclose(new,v,rtol=0,atol=1e-12)),
                             stored_internal_value=self.z[k],is_rate=k in RATE_KEYS))
        return dict(pdf_sha256=PDF_SHA,pages=[1,2],rows=rows,
            matching_count=sum(r['matches_pdf'] for r in rows),changed_count=sum(r['changed'] for r in rows),
            not_in_pdf=dict(n_int1=self.z['n_int1']),
            retained_hby_extension=dict(n_A1_gate=self.tail.c.n_A1_gate,F1_production_exponent=self.p['n_A'][1],
                clock_K=self.tail.c.clock_K_au,clock_n=self.tail.c.clock_n),
            scope='Replace all PDF-listed parameters used by ZMH carry0/bit1. HBY bit0/carry1/bit2 remain unchanged; this is not a full reconstruction of the PDF circuit.',
            dimensional_note='K_D_comp numeric value 1.2 retained; because the original reverse Hill argument is I*R, its threshold has product-concentration dimensions despite the PDF a.u. label.')
