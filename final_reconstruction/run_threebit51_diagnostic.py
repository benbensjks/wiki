"""Untuned 600 h diagnostic run for the formal 51-state scaffold."""
from __future__ import annotations

import json
from dataclasses import asdict
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from model import ROOT
from model_threebit51 import STATE_NAMES_51, ThreeBit51Model
from verify_threebit51 import analyse_threebit


OUT=ROOT/'threebit51_results'/'diagnostic'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    model=ThreeBit51Model(); raw=OUT/'trajectory_raw.npz'
    if raw.exists():
        saved=np.load(raw); sol=SimpleNamespace(t=saved['t'],y=saved['y'])
    else:
        sol=model.simulate(hours=600,sample_min=2,max_step_min=2)
        np.savez_compressed(raw,t=sol.t,y=sol.y)
    analysis=analyse_threebit(model,sol,600);sig=model.diagnostic_signals(sol.y)
    flux=np.asarray([model.flux(sol.y[:,k]) for k in range(sol.y.shape[1])])
    data={'time_h':sol.t};data.update({n:sol.y[i] for i,n in enumerate(STATE_NAMES_51)})
    data['C31_flux']=flux;data.update(sig);pd.DataFrame(data).to_csv(OUT/'trajectories.csv',index=False)
    (OUT/'analysis.json').write_text(json.dumps(analysis,ensure_ascii=False,indent=2),encoding='utf-8')
    params=dict(extension=asdict(model.e),carry0=asdict(model.carry.carry0),
                carry1=asdict(model.carry.carry1),state_names=list(STATE_NAMES_51),
                note='carry1 parameters are initial uncalibrated assumptions')
    (OUT/'parameters.json').write_text(json.dumps(params,ensure_ascii=False,indent=2),encoding='utf-8')
    summary=dict(states=51,hours=600,samples=int(sol.t.size),cold_start=analysis['cold_start'],
                 steady_state=analysis['steady_state'],certified=analysis['certified'],
                 causal_verdict=analysis['causal_verdict'],gate_contrast=analysis['gate_contrast'],
                 signal_ranges=analysis['signal_ranges'],bit_margins=analysis['bit_margins'],
                 reverse_events=len(analysis['bit1_reverse_events']),
                 gate_events=len(analysis['carry1_gate_events']),bit2_crossings=len(analysis['bit2_crossings']))
    (OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')

    reads=analysis['read_windows'];rt=np.asarray([r['trough_h'] for r in reads]);
    rv=np.asarray([np.nan if r['value'] is None else r['value'] for r in reads])
    fig,ax=plt.subplots(4,1,figsize=(14,10),sharex=True,constrained_layout=True)
    ax[0].plot(sol.t,sol.y[16],label='S0');ax[0].plot(sol.t,sol.y[27],label='S1');ax[0].plot(sol.t,sol.y[44],label='S2')
    ax[0].axhline(.3,color='.6',ls='--');ax[0].axhline(.7,color='.6',ls='--');ax[0].set_ylabel('DNA states');ax[0].legend(ncol=3)
    ax[1].step(rt,rv,where='mid',color='black');ax[1].scatter(rt,rv,c=rv,cmap='viridis',vmin=0,vmax=7,s=20)
    ax[1].set_yticks(range(8));ax[1].set_ylabel('decoded 0..7')
    ax[2].plot(sol.t,sig['J_rev1'],label='bit1 reverse flux');ax[2].plot(sol.t,sig['g1'],label='g1');
    ax[2].plot(sol.t,sig['Int2_source']/max(sig['Int2_source'].max(),1e-12),label='Int2 source normalized')
    ax[2].set_ylabel('bit1→bit2 signals');ax[2].legend(ncol=3)
    cg=analysis['cycle_gate_metrics'];x=[reads[r['cycle']]['trough_h'] for r in cg]
    ax[3].plot(x,[r['g1_peak'] for r in cg],label='cycle peak');ax[3].plot(x,[r['g1_median'] for r in cg],label='cycle median')
    ax[3].plot(x,[r['g1_min'] for r in cg],label='cycle minimum');ax[3].set_ylabel('g1 by cycle');ax[3].set_xlabel('clock cycle / time (h)')
    ax[3].legend(ncol=3);ax[3].set_xlabel('time (h)');ax[0].set_title('Untuned 51-state three-bit diagnostic')
    fig.savefig(OUT/'threebit_diagnostic.png',dpi=220,bbox_inches='tight',facecolor='white')
    fig.savefig(OUT/'threebit_diagnostic.pdf',bbox_inches='tight',facecolor='white');plt.close(fig)
    print(json.dumps(summary,ensure_ascii=False,indent=2));print(f'wrote {OUT}')


if __name__=='__main__':main()
