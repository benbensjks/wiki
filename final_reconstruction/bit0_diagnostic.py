"""Bit0-only diagnosis and scan of extension assumptions.

Zeng table parameters are never varied. Only Extension fields listed in GRID
change. If a stable bit0 candidate exists, A0/F0 -> bit1 is tested afterwards.
"""
from __future__ import annotations

import argparse
import itertools
import json
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.signal import find_peaks, peak_widths

from model import Model, Extension, ZENG, ROOT, BIT0_STATE_NAMES, TWOBIT_STATE_NAMES

GRID = dict(
    uM_per_au=(0.1, 0.3, 1.0, 3.0, 10.0),
    maturation_half_life_min=(1.0, 2.5, 5.0, 10.0, 20.0),
    complex_on_au_inv_h=(0.1, 0.3, 1.0, 3.0),
    complex_off_h=(1.0, 3.0, 10.0, 30.0),
    add_growth=(False, True),
)


def source_peaks(t, flux):
    ids, _ = find_peaks(flux, prominence=max(1.0, .1*np.ptp(flux)), distance=int(5/(t[1]-t[0])))
    return ids


def cycle_metrics(model, sol):
    t = sol.t
    Y = sol.y.T
    flux = np.array([model.flux(model._expand_bit0(y)) for y in Y])
    peaks = source_peaks(t, flux)
    rows = []
    for n, (a, b) in enumerate(zip(peaks[:-1], peaks[1:]), 1):
        # Exclude first two source intervals from stable-alternation score.
        k = a + int(np.argmin(flux[a:b]))
        I, R, C, S = Y[:, 8], Y[:, 14], Y[:, 15], Y[:, 16]
        vf = ZENG['k_fwd']*(np.maximum(I,0)/ZENG['K_D_int'][0])**2
        vf = vf/(1+(np.maximum(I,0)/ZENG['K_D_int'][0])**2)*ZENG['K_inh']/(ZENG['K_inh']+np.maximum(R,0))
        vr = ZENG['k_rev']*(np.maximum(C,0)/model.K_complex)**2/(1+(np.maximum(C,0)/model.K_complex)**2)
        sl = slice(a,b+1)
        forward = vf[sl]*(1-S[sl]); reverse = vr[sl]*S[sl]
        active_f = forward > .1*max(float(forward.max()),1e-12)
        active_r = reverse > .1*max(float(reverse.max()),1e-12)
        rows.append(dict(cycle=n, left_h=float(t[a]), right_h=float(t[b]), sample_h=float(t[k]),
                         S_sample=float(S[k]), S_min=float(S[sl].min()), S_max=float(S[sl].max()),
                         delta_S=float(S[b]-S[a]), forward_integral=float(trapezoid(forward,t[sl])),
                         reverse_integral=float(trapezoid(reverse,t[sl])),
                         simultaneous_integral=float(trapezoid(np.minimum(forward,reverse),t[sl])),
                         overlap_h=float(np.sum(active_f & active_r)*(t[1]-t[0])),
                         I_max=float(I[sl].max()), R_max=float(R[sl].max()), C_max=float(C[sl].max())))
    return flux, peaks, pd.DataFrame(rows)


def score_cycles(cycles):
    stable = cycles.iloc[2:].copy()
    if len(stable) < 6:
        return dict(alternation_accuracy=0.0, confident_fraction=0.0, dynamic_range=0.0,
                    stable_alternation=False, codes='')
    values = stable.S_sample.to_numpy()
    codes = np.where(values >= .8, 1, np.where(values <= .2, 0, -1))
    confident = codes >= 0
    transitions = (codes[1:] == 1-codes[:-1]) & confident[1:] & confident[:-1]
    alt = float(transitions.mean()) if len(transitions) else 0.0
    conf = float(confident.mean())
    dynamic = float(np.median(stable.S_max-stable.S_min))
    late = codes[-5:]
    late_ok = bool(len(late) == 5 and np.all(late >= 0)
                   and np.all(late[1:] == 1-late[:-1]))
    return dict(alternation_accuracy=alt, confident_fraction=conf, dynamic_range=dynamic,
                stable_alternation=late_ok,
                codes=''.join('x' if q<0 else str(q) for q in codes),
                late_codes=''.join('x' if q<0 else str(q) for q in late))


def evaluate(cfg, hours=100.0):
    try:
        m = Model(Extension(**cfg))
        sol = m.simulate_bit0(hours=hours, sample_min=2, max_step_min=4)
        flux, peaks, cycles = cycle_metrics(m, sol)
        score = score_cycles(cycles)
        stable = cycles.iloc[2:]
        return dict(**cfg, **score, solver=True, minimum=float(sol.y.min()),
                    S_min=float(sol.y[16].min()), S_max=float(sol.y[16].max()),
                    median_forward=float(stable.forward_integral.median()),
                    median_reverse=float(stable.reverse_integral.median()),
                    median_simultaneous=float(stable.simultaneous_integral.median()),
                    median_overlap_h=float(stable.overlap_h.median()), error='')
    except Exception as exc:
        return dict(**cfg, solver=False, stable_alternation=False, error=repr(exc))


def diagnose_one(cfg, out, stem, hours=100.0):
    m = Model(Extension(**cfg)); sol = m.simulate_bit0(hours=hours)
    flux, peaks, cycles = cycle_metrics(m, sol); cycles.to_csv(out/f'{stem}_cycles.csv',index=False)
    df = pd.DataFrame(sol.y.T,columns=BIT0_STATE_NAMES); df.insert(0,'time_h',sol.t); df['flux_uM_h']=flux
    I,R,C,S=df.b0_I.to_numpy(),df.b0_R.to_numpy(),df.b0_C.to_numpy(),df.b0_S.to_numpy()
    vf=ZENG['k_fwd']*(np.maximum(I,0)/ZENG['K_D_int'][0])**2/(1+(np.maximum(I,0)/ZENG['K_D_int'][0])**2)*ZENG['K_inh']/(ZENG['K_inh']+np.maximum(R,0))
    vr=ZENG['k_rev']*(np.maximum(C,0)/m.K_complex)**2/(1+(np.maximum(C,0)/m.K_complex)**2)
    df['forward_flux']=vf*(1-S);df['reverse_flux']=vr*S;df.to_csv(out/f'{stem}_trajectory.csv',index=False)
    fig,ax=plt.subplots(5,1,figsize=(13,13),sharex=True,constrained_layout=True)
    ax[0].plot(sol.t,flux,c='black',label='real C31 translation flux')
    ax[1].plot(sol.t,I,label='free Int');ax[1].plot(sol.t,R,label='free RDF');ax[1].plot(sol.t,C,label='Int-RDF complex')
    ax[2].plot(sol.t,df.forward_flux,label='forward');ax[2].plot(sol.t,df.reverse_flux,label='reverse')
    ax[3].plot(sol.t,S,label='LR fraction',c='#0072B2');ax[3].axhspan(.2,.8,color='grey',alpha=.15,label='ambiguous')
    ax[4].plot(cycles.cycle,cycles.forward_integral,'o-',label='forward integral')
    ax[4].plot(cycles.cycle,cycles.reverse_integral,'o-',label='reverse integral')
    ax[4].plot(cycles.cycle,cycles.simultaneous_integral,'o-',label='simultaneous minimum integral')
    for a in ax: a.grid(alpha=.2);a.legend(loc='upper right')
    ax[-1].set_xlabel('time (h), except last panel: cycle');fig.savefig(out/f'{stem}_diagnostic.png',dpi=160);plt.close(fig)
    return m,sol,cycles,score_cycles(cycles)


def test_twobit(cfg, out, hours=120.0):
    m=Model(Extension(**cfg)); sol=m.simulate_twobit(hours=hours)
    df=pd.DataFrame(sol.y.T,columns=TWOBIT_STATE_NAMES);df.insert(0,'time_h',sol.t)
    expanded=[m._expand_twobit(y) for y in sol.y.T]
    df['upstream_flux_uM_h']=[m.flux(y) for y in expanded]
    df['carry_promoter_01']=[m.carry_promoters(y)[0] for y in expanded]
    df.to_csv(out/'twobit_candidate_trajectory.csv',index=False)
    S0,S1=df.b0_S.to_numpy(),df.b1_S.to_numpy(); carry=df.carry_promoter_01.to_numpy()
    # Falling crossings define LR->PB. Match carry peaks inside the next source period.
    fall=np.flatnonzero((S0[:-1]>=.8)&(S0[1:]<.8))+1
    p,_=find_peaks(carry,prominence=max(.02*np.ptp(carry),1e-6),distance=60)
    period=10.59; counts=[]
    for k in fall:
        counts.append(int(np.sum((sol.t[p]>=sol.t[k])&(sol.t[p]<sol.t[k]+period))))
    result=dict(falling_edges_h=sol.t[fall].tolist(),carry_peak_h=sol.t[p].tolist(),
                carry_peaks_per_falling_edge=counts,
                exactly_one_per_falling_edge=bool(len(fall)>=2 and counts and all(n==1 for n in counts)),
                bit1_S_range=[float(S1.min()),float(S1.max())])
    (out/'twobit_candidate_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    fig,ax=plt.subplots(4,1,figsize=(13,10),sharex=True,constrained_layout=True)
    ax[0].plot(sol.t,df.upstream_flux_uM_h,c='black');ax[1].plot(sol.t,S0,label='S0');ax[1].plot(sol.t,S1,label='S1')
    ax[2].plot(sol.t,df.A0,label='A0');ax[2].plot(sol.t,df.F0,label='R0 carry repressor')
    ax[3].plot(sol.t,carry,label='carry promoter');ax[3].plot(sol.t,df.b1_I,label='mature Int1')
    for a in ax:a.legend();a.grid(alpha=.2)
    ax[-1].set_xlabel('Time (h)');fig.savefig(out/'twobit_candidate.png',dpi=160);plt.close(fig)
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=6);ap.add_argument('--hours',type=float,default=100)
    ap.add_argument('--output',type=Path,default=ROOT/'bit0_results');args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    base=asdict(Extension()); scan_keys=list(GRID)
    configs=[]
    for vals in itertools.product(*(GRID[k] for k in scan_keys)):
        c=base.copy();c.update(dict(zip(scan_keys,vals)));configs.append(c)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows=list(pool.map(evaluate,configs,itertools.repeat(args.hours),chunksize=4))
    scan=pd.DataFrame(rows)
    scan['rank_score']=scan.alternation_accuracy.fillna(0)*100+scan.confident_fraction.fillna(0)*20+scan.dynamic_range.fillna(0)
    scan=scan.sort_values(['stable_alternation','rank_score'],ascending=False);scan.to_csv(out/'extension_parameter_scan.csv',index=False)
    stable=scan[scan.stable_alternation==True]
    baseline_cfg=base.copy();diagnose_one(baseline_cfg,out,'baseline',args.hours)
    top_cfg={k:scan.iloc[0][k].item() if hasattr(scan.iloc[0][k],'item') else scan.iloc[0][k] for k in base}
    _,_,_,top_score=diagnose_one(top_cfg,out,'top_candidate',args.hours)
    two=None
    if len(stable):
        chosen={k:stable.iloc[0][k].item() if hasattr(stable.iloc[0][k],'item') else stable.iloc[0][k] for k in base}
        two=test_twobit(chosen,out,max(args.hours,120))
    summary=dict(grid_size=len(scan),solved=int(scan.solver.sum()),stable_bit0_candidates=int(len(stable)),
                 varied_parameters=scan_keys,zeng_parameters_varied=False,
                 top_candidate=top_cfg,top_candidate_score=top_score,
                 twobit_test=two,
                 interpretation=('Stable bit0 found; two-bit carry test executed.' if len(stable) else
                                 'No stable bit0 in the declared extension grid; carry test correctly skipped.'))
    (out/'scan_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
