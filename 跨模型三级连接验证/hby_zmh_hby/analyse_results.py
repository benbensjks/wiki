"""Post-run diagnostics; preserve original readout and distinguish post-hoc windows."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from hybrid_model import sha256

ROOT=Path(__file__).resolve().parent

def dump(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')

def associations(events_low,events_high):
    triggers=[e['time_h'] for e in events_low if e['direction']=='down']
    rows=[]
    for j,t in enumerate(triggers):
        stop=triggers[j+1] if j+1<len(triggers) else float('inf')
        hits=[e for e in events_high if t<=e['time_h']<stop]
        rows.append(dict(trigger_h=t,count=len(hits),delays_h=[e['time_h']-t for e in hits]))
    assigned=sum(r['count'] for r in rows)
    return dict(rows=rows,one_per_observed_rollover=bool(rows and all(r['count']==1 for r in rows)),
                unassigned_high_events=len(events_high)-assigned,
                alternating_high_directions=all(a['direction']!=b['direction'] for a,b in zip(events_high[:-1],events_high[1:])),
                scope='0.5 DNA crossings only; not a gate-window causal certificate; last interval is right-censored')

def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('run');a=p.parse_args();run=Path(a.run)
    out=run/'follow_up_analysis';out.mkdir(exist_ok=False)
    result={}
    for mode in ('square','han'):
        r=json.loads((run/mode/'readout.json').read_text());ev=r['crossings']
        result[mode]=dict(level01=associations(ev['bit0'],ev['bit1']),level12=associations(ev['bit1'],ev['bit2']),
            min_original_commitment=min(min(w['commitment']) for w in r['reads']),
            unlabelled_by_bit=[sum(w['labels'][i] is None for w in r['reads']) for i in range(3)])
    data=np.load(run/'square/trajectory.npz');t=data['time_h'];y=data['states'];ss=[y[10],1-y[13],y[33]]
    # Named alternative protocol: input OFF-interval midpoint, not a phase sweep.
    pars=json.loads((run/'square/parameters.json').read_text())['square'];period=pars['period_h']
    rows=[]
    for start in np.arange(pars['start_h'],t[-1],period):
        centre=start+(pars['width_h']+period)/2
        left,right=centre-.1*period,centre+.1*period
        if right>t[-1]:continue
        mask=(t>=left)&(t<=right);labels=[];commit=[]
        for s in ss:
            lo=float(np.mean(s[mask]<=.3));hi=float(np.mean(s[mask]>=.7))
            labels.append(0 if lo>=.8 else 1 if hi>=.8 else None);commit.append(max(lo,hi))
        value=sum(2**i*v for i,v in enumerate(labels)) if all(v is not None for v in labels) else None
        rows.append(dict(centre_h=float(centre),start_h=float(left),end_h=float(right),labels=labels,commitment=commit,value=value))
    vals=[r['value'] for r in rows]
    result['square_alternative_readout']=dict(status='POST_HOC_DIAGNOSTIC_NOT_CERTIFICATION',
        rule='OFF interval midpoint; width=20% input period; no phase scan; original trough readout retained',
        reads=rows,sequence=''.join('x' if v is None else str(v) for v in vals),
        increments=all(v is not None for v in vals) and all((b-a)%8==1 for a,b in zip(vals[:-1],vals[1:])),
        min_commitment=min(min(r['commitment']) for r in rows))
    fig,ax=plt.subplots(3,1,figsize=(13,8),sharex=True,constrained_layout=True)
    for s,name in zip(ss,['S0 HBY','S1 ZMH','S2 HBY']):ax[0].plot(t,s,label=name)
    ax[0].legend();ax[0].set_ylabel('LR fraction')
    old=json.loads((run/'square/readout.json').read_text())['reads']
    ax[1].plot([r['trough_h'] for r in old],np.full(len(old),3.5),'rx');ax[1].set_ylabel('Original trough rule')
    ax[2].plot([r['centre_h'] for r in rows],vals,'o-',drawstyle='steps-post');ax[2].set_yticks(range(8));ax[2].set_ylabel('OFF-midpoint rule');ax[2].set_xlabel('Time (h)')
    fig.suptitle('Square wave: original readout fails; post-hoc OFF-midpoint diagnostic shown separately')
    fig.savefig(out/'square_readout_diagnostic.png',dpi=150);plt.close(fig)
    dump(out/'diagnostics.json',result)
    dump(out/'SHA256SUMS.json',{p.name:sha256(p) for p in out.iterdir() if p.is_file()})
    print(json.dumps({k:(v if k!='square_alternative_readout' else {kk:vv for kk,vv in v.items() if kk!='reads'}) for k,v in result.items()},indent=2))

if __name__=='__main__':main()
