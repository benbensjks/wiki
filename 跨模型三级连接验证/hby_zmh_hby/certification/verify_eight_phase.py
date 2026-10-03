"""Eight complete biochemical orbit states, each with its own Han clock phase."""
from __future__ import annotations

import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from verify_hzh import HERE,NAMES,analyse,dump
from hybrid_model import sha256,source_hashes
from model_hzh import HZHModel
from run_han_comparison import HanInput,HAN_PATH

BASE=HERE.parent/'results'/'20260927_233646_924908'/'han'
ROOT=HERE/'eight_phase_results'


def main():
    if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
    out=ROOT/datetime.now().strftime('%Y%m%d_%H%M%S_%f');out.mkdir(parents=True)
    dump(out/'status.json',dict(status='RUNNING'))
    inputs=[BASE/'trajectory.npz',BASE/'parameters.json',BASE/'readout.json',HERE/'verify_hzh.py',
            HERE/'verify_eight_phase.py',HERE.parent/'model_hzh.py',HAN_PATH]
    hashes={str(p):sha256(p) for p in inputs};hashes.update(source_hashes())
    try:
        z=np.load(BASE/'trajectory.npz');t0=z['time_h'];y0=z['states']
        assert tuple(z['state_names'])==tuple(NAMES)
        par=json.loads((BASE/'parameters.json').read_text(encoding='utf-8'))
        read=json.loads((BASE/'readout.json').read_text(encoding='utf-8'))['reads']
        eight=read[8:16]
        if len(eight)!=8 or len({r['value'] for r in eight})!=8:raise RuntimeError('No full orbit phase set')
        end=max(r['trough_h'] for r in eight)+300
        han=HanInput(end/60*60+1) # hours; do not reset clock phase
        model=HZHModel('han',han)
        if model.z!=par['Zmh_parameters_per_min']:raise RuntimeError('Parameter mismatch')
        results=[]
        for r in eight:
            target=float(r['trough_h'])
            i=int(np.argmin(abs(t0-target)))
            if abs(t0[i]-target)>1/120:raise RuntimeError('Original orbit sample not aligned')
            start=float(t0[i]);t=np.linspace(start,start+300,18001)
            started=time.perf_counter()
            sol=solve_ivp(model.rhs,(start,start+300),y0[:,i],t_eval=t,method='DOP853',
                          rtol=2e-7,atol=2e-9,max_step=1/60)
            if not sol.success or not np.isfinite(sol.y).all():raise RuntimeError(sol.message)
            verdict,_=analyse(sol.t,sol.y,model.z,model.tail)
            expected=(int(r['value'])+1)%8
            observed=verdict['cold']['sequence'][:1]
            entry=dict(initial_read=int(r['value']),start_h=start,expected_first_new_read=expected,
                first_new_read=observed,first_read_matches=str(expected)==observed,
                certified_v1=verdict['certified_v1'],steady=verdict['steady'],
                events={k:{kk:x[kk] for kk in ('reverse_events','gate_events','flip_events',
                       'one_to_one','causal_order','alternating_directions','passed')}
                        for k,x in verdict['events'].items()},
                min_timing_h=verdict['global_min_timing_margin_h'],runtime_s=time.perf_counter()-started)
            if not entry['certified_v1'] or not entry['first_read_matches']:
                entry['attention']='Inspect full verdict; do not drop phase as startup'
            folder=out/f"phase_{int(r['value'])}";folder.mkdir()
            np.savez_compressed(folder/'trajectory.npz',time_h=sol.t,states=sol.y,state_names=NAMES)
            dump(folder/'verdict.json',verdict)
            results.append(entry)
            print('Phase',int(r['value']),'first',observed,'expected',expected,
                  'cert',verdict['certified_v1'],'runtime',entry['runtime_s'],flush=True)
        after={k:sha256(k) for k in hashes if Path(k).is_absolute()};after.update(source_hashes())
        if hashes!=after:raise RuntimeError('Input source changed during phases')
        dump(out/'summary.json',dict(status='PASSED' if all(r['certified_v1'] and r['first_read_matches'] for r in results) else 'PARTIAL',
             scope='Eight orbit-derived complete 34-state initial conditions at original absolute Han clock phase; new 300 h each',
             phase_count=len(results),passed=sum(r['certified_v1'] for r in results),
             results=results,source_hashes=hashes))
        dump(out/'status.json',dict(status='COMPLETED'))
        dump(out/'SHA256SUMS.json',{str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file()})
        print('OUTPUT',out,flush=True)
    except Exception as e:
        dump(out/'status.json',dict(status='FAILED',error=repr(e)))
        raise


if __name__=='__main__':main()
