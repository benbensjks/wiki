"""Han-upstream replacement: fixed-interface primary test and absolute-flux diagnostic.

Does not modify any existing model/parameter. Primary case retains ZMH Int0
production law. Absolute-flux case explicitly changes that source, not its loss.
"""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
import sys
import time
from dataclasses import asdict
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from hybrid_model import (ROOT,WIKI,STATE_NAMES,DONOR_NAMES,INDEX,ZmhDonor,donor_namespace,
                          HybridModel,HbyReceiver,integrate,source_hashes,sha256,require_hash)
from diagnostics import analyse,state_health

HAN_PATH=WIKI/'其他小组成员任务/Week4_振荡器-C31建模/code/Flux_Driven_Translation_Burden_Model.py'
HAN_SHA='456326A418C4C7EE23C75FA2DD4A86070F37CAF820B5F863BBAA44FE9DCE51BD'
MODES=('normalized_concentration','absolute_translation_flux')


def dump(path,obj):
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


class HanInput:
    def __init__(self,hours):
        require_hash(HAN_PATH,HAN_SHA)
        spec=importlib.util.spec_from_file_location('hybrid_han_reference',HAN_PATH)
        module=importlib.util.module_from_spec(spec)
        sys.modules[spec.name]=module
        spec.loader.exec_module(module)
        self.module=module
        self.p=module.make_reference_parameters()
        assert self.p.peak_load_fraction == 0.0
        duration=max(3000.,60*hours)
        self.t=np.r_[np.arange(0,duration,1.),duration]
        self.sol=solve_ivp(lambda t,y:module.rhs(t,y,self.p,1.),(0,duration),
                           module.initial_state(self.p),t_eval=self.t,
                           method='LSODA',rtol=1e-9,atol=1e-11,max_step=1.)
        if not self.sol.success or not np.isfinite(self.sol.y).all():
            raise RuntimeError(self.sol.message)
        # Normalize once on the original source's physical reference horizon;
        # extending the plotted experiment does NOT change this denominator.
        self.normalization_mask=(self.t>=1000)&(self.t<=3000)
        self.c31_max=float(self.sol.y[7,self.normalization_mask].max())
        self.flux_copies_min=self.p.c31_translation_per_mrna_per_min*self.sol.y[6]
        self.copies_per_uM=602.214076 # 1 fL, same explicit convention as frozen model

    def interp(self,t,signal):
        if t < -1e-9 or t > self.t[-1]+1e-9:
            raise ValueError('Han input extrapolation forbidden')
        return float(np.interp(t,self.t,signal))


class HanDonor(ZmhDonor):
    def __init__(self,han,mode):
        if mode not in MODES: raise ValueError(mode)
        self.ns=donor_namespace()
        self.p=self.ns['p'].copy()
        self.han,self.mode=han,mode
        self.receiver_scale=5.75
        self.ns['get_u_in']=self.input_min

    def source_min(self,t):
        if self.mode=='normalized_concentration':
            return self.p['k_int']*float(np.clip(self.han.interp(t,self.han.sol.y[7])/self.han.c31_max,0,2))
        return self.han.interp(t,self.han.flux_copies_min)/(self.han.copies_per_uM*self.receiver_scale)

    def input_min(self,t):
        # The existing donor RHS multiplies by k_int. Division here replaces
        # exactly its source for the absolute-flux diagnostic (no extra x6).
        return self.source_min(t)/self.p['k_int']


def run_case(out,han,mode,hours,sample_min):
    donor=HanDonor(han,mode)
    receiver=HbyReceiver()
    model=HybridModel(donor,receiver)
    t=np.r_[np.arange(0,hours,sample_min/60),hours]
    t=np.unique(t)
    print(f'Running {mode}, {hours:g} h...',flush=True)
    sol=integrate(model.rhs,model.initial_state(),t,max_step_min=.5)
    dsol=integrate(donor.rhs,donor.initial_state(),t,max_step_min=.5)
    signals=[]
    for j,tt in enumerate(t):
        s=receiver.signals(sol.y[10:,j],sol.y[1,j])
        s.update(Int0_source_au_per_h=60*donor.source_min(60*tt),
                 C31_concentration_copies=han.interp(60*tt,han.sol.y[7]),
                 C31_flux_uM_per_h=60*han.interp(60*tt,han.flux_copies_min)/han.copies_per_uM,
                 g0=donor.gate0(sol.y[:10,j]))
        signals.append(s)
    target=out/mode;target.mkdir()
    np.savez_compressed(target/'trajectory.npz',time_h=t,states=sol.y,donor=dsol.y,state_names=STATE_NAMES)
    with (target/'trajectory.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(['time_h',*STATE_NAMES,*signals[0]])
        for j,tt in enumerate(t): w.writerow([tt,*sol.y[:,j],*signals[j].values()])
    full=analyse(t,sol.y,sol.y[1])
    mask=t<=2000/60
    short=analyse(t[mask],sol.y[:,mask],sol.y[1,mask])
    extrema={n:dict(minimum=float(sol.y[i].min()),maximum=float(sol.y[i].max())) for i,n in enumerate(STATE_NAMES)}
    summary=dict(mode=mode,hours=hours,certified=None,
                 low_state_max_abs_gap=float(np.max(np.abs(sol.y[:10]-dsol.y))),
                 health=state_health(sol.y,STATE_NAMES,('pb0','pb1','b2_S')),
                 extrema=extrema,source_peak_au_per_h=max(s['Int0_source_au_per_h'] for s in signals),
                 source_trough_au_per_h=min(s['Int0_source_au_per_h'] for s in signals),
                 g1_peak=max(s['g1'] for s in signals),u2_peak_au_per_h=max(s['u2_target_au_per_h'] for s in signals),
                 matched_2000min={k:short[k] for k in ('cold_mod4','cold_mod8')},
                 extended={k:full[k] for k in ('cold_mod4','cold_mod8','steady_mod4','steady_mod8')},
                 event_counts={k:len(v) for k,v in full['crossings'].items()})
    dump(target/'readout.json',full);dump(target/'readout_2000min.json',short)
    dump(target/'summary.json',summary)
    dump(target/'parameters.json',dict(donor_parameters_per_min=donor.p,receiver=asdict(receiver.c),
              receiver_parameters=receiver.p,state_names=STATE_NAMES,initial_state=model.initial_state().tolist(),
              mapping=('u=C31/C31_refmax; dint0/dmin=6*u-2*int0' if mode==MODES[0] else
                       'dint0/dmin=J_C31_copies_min/(602.214076*5.75)-2*int0; no extra k_int'),
              assumptions=['clock_scale=1 uncalibrated','no downstream feedback to upstream','1 fL cell for flux conversion'],
              solver=dict(method='DOP853',rtol=2e-7,atol=2e-9,max_step_min=.5)))
    fig,axes=plt.subplots(4,1,figsize=(12,10),sharex=True,constrained_layout=True)
    axes[0].plot(t,sol.y[1],color='black');axes[0].set_ylabel('Int0 (donor a.u.)')
    for s,name in ((1-sol.y[0],'S0'),(1-sol.y[6],'S1'),(sol.y[26],'S2')):
        axes[1].plot(t,s,label=name)
    axes[1].axhspan(0,.3,alpha=.08,color='green');axes[1].axhspan(.7,1,alpha=.08,color='blue')
    axes[1].legend();axes[1].set_ylabel('LR fraction (raw)')
    rr=full['reads'];axes[2].plot([r['trough_h'] for r in rr],
        [np.nan if r['value'] is None else r['value'] for r in rr],'o-',drawstyle='steps-post')
    for r in rr:
        if r['value'] is None: axes[2].text(r['trough_h'],3.5,'x',color='red')
    axes[2].set_yticks(range(8));axes[2].set_ylabel('Exploratory readout')
    axes[3].plot(t,[s['u2_target_au_per_h'] for s in signals],label='Int2 target source')
    axes[3].legend();axes[3].set_ylabel('Receiver a.u./h');axes[3].set_xlabel('Time (h)')
    fig.suptitle(f'Han upstream / {mode}\nNew cross-model experiment, NOT certified')
    fig.savefig(target/'overview.png',dpi=150);fig.savefig(target/'overview.pdf');plt.close(fig)
    print(json.dumps({k:summary[k] for k in ('mode','source_peak_au_per_h','g1_peak','matched_2000min','extended','event_counts')},indent=2),flush=True)
    return summary


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--hours',type=float,default=120.)
    ap.add_argument('--sample-min',type=float,default=1.)
    args=ap.parse_args()
    if args.hours<2000/60 or args.sample_min<=0: raise ValueError('Need at least original 2000min and positive sampling')
    if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
    out=ROOT/'han_results'/datetime.now().strftime('%Y%m%d_%H%M%S_%f');out.mkdir(parents=True)
    paths=[HAN_PATH,ROOT/'hybrid_model.py',ROOT/'diagnostics.py',ROOT/'run_han_comparison.py']
    hashes={str(p):sha256(p) for p in paths};hashes.update(source_hashes())
    dump(out/'status.json',dict(status='RUNNING'))
    try:
        han=HanInput(args.hours)
        np.savez_compressed(out/'han_upstream.npz',time_min=han.t,states=han.sol.y,flux_copies_min=han.flux_copies_min)
        dump(out/'han_parameters.json',dict(parameters=asdict(han.p),source_sha256=HAN_SHA,
            solver=dict(method='LSODA',rtol=1e-9,atol=1e-11,max_step_min=1),
            normalization=dict(reference_start_min=1000,reference_end_min=3000,C31_max_copies=han.c31_max),
            state_names=['m_TetR','TetR','m_CI','CI','m_LacI','LacI','m_C31','C31'],
            initial_state=list(han.module.initial_state(han.p)),load_fraction=0.0))
        results=[run_case(out,han,mode,args.hours,args.sample_min) for mode in MODES]
        after={str(p):sha256(p) for p in paths};after.update(source_hashes())
        if hashes!=after:raise RuntimeError('Inputs changed during run')
        dump(out/'comparison.json',dict(cases=results,source_hashes=hashes,
             note='Primary test retains normalized-concentration interface. Absolute flux changes the source mapping; not an upstream-only comparison.'))
        dump(out/'status.json',dict(status='COMPLETED',certified=None))
        dump(out/'SHA256SUMS.json',{str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file()})
        print('OUTPUT:',out,flush=True)
    except Exception as e:
        dump(out/'status.json',dict(status='FAILED',error=repr(e)));raise


if __name__=='__main__':main()
