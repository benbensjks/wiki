"""Read-only diagnostic of the existing 300 h ZMH-bit1 trajectories."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from model_pdf import HERE, NAMES
from hybrid_model import sha256

OLD = HERE.parent / 'results' / '20260927_233646_924908'
PDF = HERE / 'results' / '20260927_235006_013973'
DEST = HERE / 'diagnostic_bit1_20260928'


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def crossings(t, s):
    side = np.asarray(s) >= .5
    ids = np.flatnonzero(side[1:] != side[:-1])
    return [{'time_h': float(t[i] + (.5-s[i])*(t[i+1]-t[i])/(s[i+1]-s[i])),
             'direction': 'up' if side[i+1] else 'down'} for i in ids]


def load(root, mode):
    folder = root / mode
    z = np.load(folder / 'trajectory.npz')
    params = json.loads((folder / 'parameters.json').read_text(encoding='utf-8'))
    p = params['zmh_per_minute'] if 'zmh_per_minute' in params else params['Zmh_parameters_per_min']
    return z['time_h'], z['states'], p


def measures(root, mode, label):
    t,y,p=load(root,mode)
    assert y.shape == (34,len(t)) and tuple(NAMES)==tuple(np.load(root/mode/'trajectory.npz')['state_names'])
    s0,s1= y[10],1-y[13]
    a,f,pb,i,tr,r=y[11:17]
    act=a**p['n_A0']/(p['K_A0']**p['n_A0']+a**p['n_A0'])
    suppress=p['K_R0']**p['n_R0']/(p['K_R0']**p['n_R0']+f**p['n_R0'])
    gate=act*suppress
    u1=60*p['alpha_Int1']*gate
    positive=60*p['k_fwd']*pb*i**p['n_int1']/(p['K_D_int1']**p['n_int1']+i**p['n_int1'])*p['Kinh']/(p['Kinh']+r)
    ratio=i*r/p['K_D_comp']
    negative=60*p['k_rev']*(1-pb)*ratio**2/(1+ratio**2)
    # A transparent reconstruction check: d(1-pb1)/dt = J_fwd-J_rev.
    slope=np.gradient(s1,t)
    residual=float(np.max(np.abs(slope-(positive-negative))))
    events0=crossings(t,s0);events1=crossings(t,s1)
    starts=[event['time_h'] for event in events0 if event['direction']=='down']
    cycles=[]
    for j,start in enumerate(starts):
        if j+1>=len(starts):break # Last event is right-censored; never score as missing.
        stop=starts[j+1]
        mask=(t>=start)&(t<stop)
        if np.count_nonzero(mask)<2:continue
        inds=np.flatnonzero(mask); x=t[mask]
        hits=[e for e in events1 if start<=e['time_h']<stop]
        cycles.append(dict(label=label,mode=mode,cycle=j,start_h=start,end_h=stop,
            bit1_before=float(s1[inds[0]]),bit1_after=float(s1[inds[-1]]),
            bit1_min=float(s1[mask].min()),bit1_max=float(s1[mask].max()),
            bit1_midpoint_crossings=len(hits),crossing_directions=''.join('U' if e['direction']=='up' else 'D' for e in hits),
            A0_peak=float(a[mask].max()),F0_peak=float(f[mask].max()),
            gate_peak=float(gate[mask].max()),gate_median=float(np.median(gate[mask])),
            Int1_source_peak_per_h=float(u1[mask].max()),Int1_source_dose_au=float(np.trapezoid(u1[mask],x)),
            Int1_peak=float(i[mask].max()),Rep1_max=float(tr[mask].max()),RDF1_max=float(r[mask].max()),
            simultaneous_product_ratio_max=float(ratio[mask].max()),
            J_fwd1_integral=float(np.trapezoid(positive[mask],x)),
            J_rev1_integral=float(np.trapezoid(negative[mask],x))))
    after=t>=50
    summary=dict(label=label,mode=mode,n_bit0_down=len(starts),n_complete_intervals=len(cycles),
        n_bit1_midpoint_crossings=len(events1),
        interval_counts=dict(exactly_one=sum(c['bit1_midpoint_crossings']==1 for c in cycles),
                             zero=sum(c['bit1_midpoint_crossings']==0 for c in cycles),
                             multiple=sum(c['bit1_midpoint_crossings']>1 for c in cycles)),
        peak_gate=float(gate.max()),peak_Int1_source_per_h=float(u1.max()),peak_Int1=float(i.max()),
        peak_Rep1=float(tr.max()),peak_RDF1=float(r.max()),
        post50_min_S1=float(s1[after].min()),post50_max_S1=float(s1[after].max()),
        post50_min_PB1=float(pb[after].min()),post50_max_PB1=float(pb[after].max()),
        product_ratio_max=float(ratio.max()),
        total_J_fwd1=float(np.trapezoid(positive,t)),total_J_rev1=float(np.trapezoid(negative,t)),
        dna_balance_finite_diff_max_abs_residual=residual,
        params_per_min={k:p[k] for k in ('alpha_A0','gamma_A0','alpha_R0','gamma_R0','K_R0','n_R0',
                                        'alpha_Int1','gamma_int1','K_D_int1','n_int1','K_D_comp','alpha_rep1',
                                        'gamma_rep1','alpha_rdf1','gamma_rdf1','K_rep','n')})
    return summary,cycles,dict(t=t,S0=s0,S1=s1,A0=a,F0=f,gate=gate,u1=u1,I1=i,Rep1=tr,RDF1=r,Jf=positive,Jr=negative)


def main():
    if DEST.exists():raise FileExistsError(f'Preserving existing diagnostic: {DEST}')
    DEST.mkdir(parents=True)
    audit={'inputs':{}}
    rows=[];detail={}
    for root,name in ((OLD,'old_code_parameters'),(PDF,'pdf_table_parameters')):
        for mode in ('han','square'):
            for fn in ('trajectory.npz','parameters.json','readout.json'):
                p=root/mode/fn;audit['inputs'][str(p)]=sha256(p)
            summary,cycles,signals=measures(root,mode,name)
            rows.extend(cycles);detail[name+'_'+mode]=signals
            audit.setdefault('cases',[]).append(summary)
    dump(DEST/'summary.json',audit)
    with (DEST/'cycles.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    fig,axs=plt.subplots(4,2,figsize=(16,11),sharex='col',constrained_layout=True)
    for col,mode in enumerate(('han','square')):
        for key,color in [('old_code_parameters','tab:blue'),('pdf_table_parameters','tab:red')]:
            s=detail[key+'_'+mode]
            axs[0,col].plot(s['t'],s['S1'],color=color,label=key)
            axs[1,col].plot(s['t'],s['gate'],color=color)
            axs[2,col].plot(s['t'],s['I1'],color=color)
            axs[3,col].plot(s['t'],s['Jf']-s['Jr'],color=color)
        axs[0,col].set_title(mode);axs[0,col].legend(fontsize=8)
        axs[0,col].axhspan(0,.3,color='green',alpha=.07);axs[0,col].axhspan(.7,1,color='blue',alpha=.07)
        axs[0,col].set_ylabel('Bit1 LR fraction')
        axs[1,col].set_ylabel('g0')
        axs[2,col].set_ylabel('Int1 (a.u.)')
        axs[3,col].set_ylabel('J_fwd1-J_rev1 (/h)');axs[3,col].set_xlabel('Time (h)')
        axs[3,col].set_xlim(0,120)
    fig.suptitle('Existing trajectories: ZMH bit1 diagnosis by parameter version')
    fig.savefig(DEST/'diagnostic_first120h.png',dpi=150);fig.savefig(DEST/'diagnostic_first120h.pdf');plt.close(fig)
    dump(DEST/'SHA256SUMS.json',{p.name:sha256(p) for p in DEST.iterdir() if p.is_file()})
    for s in audit['cases']:
        print(s['label'],s['mode'],{k:s[k] for k in ('interval_counts','peak_gate','peak_Int1',
              'peak_RDF1','post50_min_S1','product_ratio_max','total_J_fwd1','total_J_rev1')})
    print('OUTPUT',DEST)


if __name__=='__main__':main()
