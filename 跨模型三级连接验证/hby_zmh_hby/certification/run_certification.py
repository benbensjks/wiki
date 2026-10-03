"""Certify the current-code-parameter HZH Han-upstream 300 h candidate."""
from __future__ import annotations

import json
import platform
import sys
import time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp

from verify_hzh import HERE, NAMES, IDX, analyse, dump
from hybrid_model import WIKI, sha256, source_hashes
from model_hzh import HZHModel
from run_han_comparison import HanInput, HAN_PATH

OLD = HERE.parent/'results'/'20260927_233646_924908'/'han'
OUT_ROOT = HERE/'results'


def integrate(model,hours,sample_min,max_step_min,rtol,atol):
    count=round(hours*60/sample_min)
    t=np.linspace(0,hours,count+1)
    if len(t)<2 or not np.all(np.diff(t)>0):raise ValueError('Bad sampling grid')
    sol=solve_ivp(model.rhs,(0,hours),model.initial_state(),t_eval=t,
                  method='DOP853',rtol=rtol,atol=atol,max_step=max_step_min/60)
    if not sol.success or not np.all(np.isfinite(sol.y)):raise RuntimeError(sol.message)
    return sol.t,sol.y


def main():
    if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
    out=OUT_ROOT/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    out.mkdir(parents=True,exist_ok=False)
    dump(out/'status.json',dict(status='RUNNING'))
    sources={str(p):sha256(p) for p in (HERE/'verify_hzh.py',HERE/'run_certification.py',
              HERE.parent/'model_hzh.py',HERE.parent/'test_hzh.py',HERE.parent.parent/'hybrid_model.py',
              HERE.parent.parent/'run_han_comparison.py',HAN_PATH,OLD/'trajectory.npz',OLD/'parameters.json')}
    sources.update(source_hashes())
    try:
        source=json.loads((OLD/'parameters.json').read_text(encoding='utf-8'))
        z=np.load(OLD/'trajectory.npz')
        t0,y0=z['time_h'],z['states']
        assert tuple(z['state_names'])==tuple(NAMES)
        han=HanInput(300)
        m=HZHModel('han',han)
        if source['Zmh_parameters_per_min']!=m.z:raise RuntimeError('ZMH active parameters changed')
        if source['Hby_config']!=asdict(m.tail.c):raise RuntimeError('HBY active parameters changed')
        if not np.array_equal(np.array(list(source['initial_state'].values())),m.initial_state()):
            raise RuntimeError('Initial state differs from frozen candidate')
        cases=[('baseline_saved',t0,y0,dict(rtol=2e-7,atol=2e-9,max_step_min=1.,sample_min=1.))]
        settings=[('tight',1.,.5,1e-9,1e-11),('ultra',1.,.25,1e-11,1e-13),
                  ('sampling_0p5min',.5,.25,1e-11,1e-13)]
        for name,sample,step,rtol,atol in settings:
            print('Integrating',name,flush=True)
            tic=time.perf_counter()
            t,y=integrate(m,300,sample,step,rtol,atol)
            cases.append((name,t,y,dict(sample_min=sample,max_step_min=step,rtol=rtol,
                                        atol=atol,runtime_s=time.perf_counter()-tic)))
        summary=[]
        reference=None
        for name,t,y,cfg in cases:
            v,s=analyse(t,y,m.z,m.tail)
            folder=out/name;folder.mkdir()
            dump(folder/'verdict.json',v)
            if name!='baseline_saved':
                np.savez_compressed(folder/'trajectory.npz',time_h=t,states=y,state_names=NAMES)
            gap=None
            if len(t)==len(t0) and np.allclose(t,t0,rtol=0,atol=1e-10):
                gap=float(np.max(np.abs(y-y0)))
            else:
                idx=np.searchsorted(t,t0)
                if np.max(abs(t[idx]-t0))>1e-9:raise RuntimeError('Grid alignment failed')
                gap=float(np.max(np.abs(y[:,idx]-y0)))
            item=dict(name=name,settings=cfg,certified_v1=v['certified_v1'],
                      cold=v['cold'],steady=v['steady'],
                      event_chains={k:{j:x[j] for j in ('reverse_events','gate_events',
                            'flip_events','one_to_one','causal_order','alternating_directions','passed')}
                            for k,x in v['events'].items()},
                      bit_margins=v['bit_margins'],minimum_timing_margin_h=v['global_min_timing_margin_h'],
                      clock_peak_count=v['clock_peak_count'],
                      max_abs_state_gap_vs_baseline=gap,
                      source_window_ends=[v['read_windows'][0]['trough_h'],v['read_windows'][-1]['trough_h']])
            summary.append(item)
            if name=='baseline_saved':reference=item
            else:
                if item['cold']['sequence']!=reference['cold']['sequence'] or item['steady']['sequence']!=reference['steady']['sequence']:
                    raise RuntimeError('Code sequence changed under accuracy or grid check')
                if item['certified_v1']!=reference['certified_v1']:
                    raise RuntimeError('Certification changed under accuracy or grid check')
                for k in reference['event_chains']:
                    if item['event_chains'][k]!=reference['event_chains'][k]:
                        raise RuntimeError('Causal event verdict changed under numerical refinement')
            print(name,'certified',v['certified_v1'],'reads',v['steady']['reads'],
                  'minimum margin',v['global_min_timing_margin_h'],flush=True)
        # Reread source hashes after the entire experiment.
        after={k:sha256(k) for k in sources if Path(k).is_absolute()}
        after.update(source_hashes())
        if after!=sources:raise RuntimeError('Source file changed during certification')
        dump(out/'summary.json',dict(status='PASSED' if all(x['certified_v1'] for x in summary) else 'FAILED',
             architecture='HBY bit0 + ZMH current-code A0/F0 and bit1 + HBY A1/F1 and bit2',
             input='Han v53d shared C31 mRNA and true translation flux',
             criterion='HZH_MOD8_CAUSAL_V1, not legacy 51-state certified',
             all_numerical_levels_same_code_and_causal_verdict=True,
             run_results=summary,source_hashes=sources,
             environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform())))
        dump(out/'status.json',dict(status='COMPLETED'))
        dump(out/'SHA256SUMS.json',{str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file()})
        print('OUTPUT',out,flush=True)
    except Exception as exc:
        dump(out/'status.json',dict(status='FAILED',error=repr(exc)))
        raise


if __name__=='__main__':main()
