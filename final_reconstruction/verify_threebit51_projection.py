"""Verify that adding bit2 leaves the frozen 34-state subsystem unchanged."""
from __future__ import annotations

import json
from dataclasses import asdict, replace

import numpy as np

from model import ROOT
from model_threebit51 import ThreeBit51Model
from model_twobit34 import CarryExpressionParameters, TwoBit34Model, nominal_extension
from verify_twobit_causal import analyse_solution


OUT = ROOT/'threebit51_results'


def frozen_two():
    p=json.loads((ROOT/'twobit34_results'/'selected_profile.json').read_text(encoding='utf-8'))
    return TwoBit34Model(replace(nominal_extension(),**p['extension']),
                         CarryExpressionParameters(**p['carry_expression']))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    three=ThreeBit51Model(); two=frozen_two()
    a=three.simulate(hours=300,sample_min=2,max_step_min=2)
    b=two.simulate(hours=300,sample_min=2,max_step_min=2)
    projection=three.project_twobit34(a.y)
    error=np.abs(projection-b.y)
    rhs_errors=[]
    for k in np.linspace(0,a.t.size-1,21,dtype=int):
        d51=three.rhs(float(a.t[k]),a.y[:,k])
        d34=two.rhs(float(a.t[k]),projection[:,k])
        rhs_errors.append(float(np.max(np.abs(three.project_twobit34_derivative(d51)-d34))))
    aa=analyse_solution(three.base,a,asdict(three.e),300)
    bb=analyse_solution(two.base,b,asdict(two.e),300)
    report=dict(hours=300,samples=int(a.t.size),
                projected_state_max_abs_error=float(error.max()),
                per_state_max_abs_error={name:float(error[i].max()) for i,name in enumerate(two.state_names)},
                projected_rhs_max_abs_error=max(rhs_errors),
                threebit_projection_sequence=aa['cold_start']['sequence'],
                frozen_twobit_sequence=bb['cold_start']['sequence'],
                threebit_projection_certified=aa['certified'],
                frozen_twobit_certified=bb['certified'])
    report['passed']=bool(report['projected_state_max_abs_error']<1e-5 and
                          report['projected_rhs_max_abs_error']<1e-12 and
                          report['threebit_projection_sequence']==report['frozen_twobit_sequence'] and
                          report['threebit_projection_certified'] and report['frozen_twobit_certified'])
    (OUT/'projection_equivalence.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    if not report['passed']: raise SystemExit('51->34 projection equivalence failed')


if __name__=='__main__':main()
