"""Descriptive F/R paired far-off leakage for the certified HZH bit2.

The denominator is the expected direction's gate-on dose summed across a
reverse and a forward carry. No threshold from the old 51-state verifier is used.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from verify_hzh import HERE, IDX, signals, dump
from hybrid_model import sha256


def windows(t,g,threshold=.01,min_duration_h=.05):
    mask=g>=threshold
    changes=np.diff(np.r_[False,mask,False].astype(np.int8))
    first=np.flatnonzero(changes==1);last=np.flatnonzero(changes==-1)-1
    return [(int(a),int(b)) for a,b in zip(first,last) if t[b]-t[a]>=min_duration_h]


def integral(x,t,a,b):
    return float(np.trapezoid(x[a:b+1],t[a:b+1])) if b>a else 0.0


def calc(t,y,par,tail_fraction):
    from model_hzh import HZHModel
    from run_han_comparison import HanInput
    # This object provides only the frozen receiver equations and constants.
    # No upstream simulation is necessary to analyze saved 34-state trajectories.
    from hybrid_model import HbyReceiver,HbyConfig
    receiver=HbyReceiver(HbyConfig(**par['Hby_config']))
    sig=signals(y,par['Zmh_parameters_per_min'],receiver)
    g,jr,jf=sig['g1'],sig['J_rev2'],sig['J_fwd2']
    i2=y[IDX['b2_I']];s2=sig['S2']
    wins=windows(t,g)
    rows=[]
    for k,(a,b) in enumerate(wins):
        peak=float(np.max(i2[a:b+1]));threshold=tail_fraction*peak
        nxt=wins[k+1][0] if k+1<len(wins) else len(t)-1
        tail=b
        while tail+1<nxt and i2[tail+1]>threshold:tail+=1
        far_h=float(t[nxt]-t[tail]) if nxt>tail else 0.
        good=bool(k+1<len(wins) and far_h>=.1 and tail<nxt)
        rows.append(dict(k=k,type='R' if s2[a]>=.5 else 'F',S2_at_open=float(s2[a]),
            start_h=float(t[a]),end_h=float(t[b]),
            gate_on_rev=integral(jr,t,a,b),gate_on_fwd=integral(jf,t,a,b),
            tail_rev=integral(jr,t,b,tail),tail_fwd=integral(jf,t,b,tail),
            far_off_rev=integral(jr,t,tail,nxt) if good else None,
            far_off_fwd=integral(jf,t,tail,nxt) if good else None,
            far_off_h=far_h,far_off_evaluable=good))
    pairs=[];anomalies=[];idx=0
    while idx+1<len(rows):
        a,b=rows[idx:idx+2]
        if a['type']!=b['type'] and a['far_off_evaluable'] and b['far_off_evaluable']:
            R,F=(a,b) if a['type']=='R' else (b,a)
            denom=R['gate_on_rev']+F['gate_on_fwd']
            pairs.append(dict(R=R['k'],F=F['k'],order=a['type']+b['type'],
                expected_dose=denom,wrong_direction_far_off_dose=F['far_off_rev']+R['far_off_fwd'],
                L_symmetric=(F['far_off_rev']+R['far_off_fwd'])/denom if denom>0 else None))
            idx+=2
        else:
            anomalies.append(dict(k1=a['k'],k2=b['k'],reason='same_type' if a['type']==b['type'] else 'not_evaluable'))
            idx+=1
    values=[x['L_symmetric'] for x in pairs if x['L_symmetric'] is not None]
    return dict(tail_fraction=tail_fraction,n_windows=len(wins),n_R=sum(r['type']=='R' for r in rows),
        n_F=sum(r['type']=='F' for r in rows),n_complete_pairs=len(pairs),
        n_evaluable_windows=sum(r['far_off_evaluable'] for r in rows),
        pattern=''.join(r['type'] for r in rows),pairs=pairs,rows=rows,anomalies=anomalies,
        L_symmetric_median=float(np.median(values)) if values else None,
        L_symmetric_min=float(min(values)) if values else None,
        L_symmetric_max=float(max(values)) if values else None,
        usage='descriptive; no predetermined acceptance threshold')


def main():
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('certification_run');args=ap.parse_args()
    source=Path(args.certification_run)
    out=source/'paired_leakage_supplement'
    out.mkdir(exist_ok=False)
    old=HERE.parent/'results'/'20260927_233646_924908'/'han'
    par=json.loads((old/'parameters.json').read_text(encoding='utf-8'))
    result={}
    sources={}
    for name in ('baseline_saved','tight','ultra','sampling_0p5min'):
        path=(old/'trajectory.npz') if name=='baseline_saved' else (source/name/'trajectory.npz')
        sources[str(path)]=sha256(path)
        z=np.load(path)
        t,y=z['time_h'],z['states']
        result[name]=dict(sample_min=float(np.median(np.diff(t))*60),
            by_cut={f'{int(100*f)}pct':calc(t,y,par,f) for f in (.01,.05,.10)})
    dump(out/'paired_leakage.json',dict(results=result,source_sha256=sources,
         definition='one R plus one F, far-off wrong-direction integral divided by expected on-window dose; tail excluded'))
    dump(out/'SHA256SUMS.json',{p.name:sha256(p) for p in out.iterdir() if p.is_file()})
    for name,entry in result.items():
        print(name,[(cut,v['n_complete_pairs'],v['pattern'],v['L_symmetric_median'])
                    for cut,v in entry['by_cut'].items()])


if __name__=='__main__':main()
