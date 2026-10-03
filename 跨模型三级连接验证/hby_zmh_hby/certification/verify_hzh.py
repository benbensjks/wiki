"""Versioned HBY-ZMH-HBY verifier; no legacy state indices or inherited verdict."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from model_hzh import NAMES,IDX
from hybrid_model import HbyReceiver,sha256

LOW=.3
HIGH=.7
OCC=.8
DROP=8
MIN_STEADY=16
PULSE_DURATION_H=.2
EVENT_TOL_H=2/60
JREV_THRESHOLD_H=.1
GATE_THRESHOLD=.05
GATE_MIN_DOSE_H=.02


def dump(path,obj):
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


def hill_array(x,k,n):
    z=np.maximum(np.asarray(x,dtype=float),0)/k
    return z**n/(1+z**n)


def crossings(t,s):
    z=np.asarray(s)-.5
    ids=np.flatnonzero(z[:-1]*z[1:]<0)
    return [dict(time_h=float(t[i]-z[i]*(t[i+1]-t[i])/(z[i+1]-z[i])),
                 direction='up' if z[i+1]>0 else 'down') for i in ids]


def segments(t,x,threshold,min_duration_h=PULSE_DURATION_H,min_dose=0):
    mask=np.asarray(x)>=threshold
    transitions=np.diff(np.r_[False,mask,False].astype(np.int8))
    starts=np.flatnonzero(transitions==1)
    stops=np.flatnonzero(transitions==-1)-1
    rows=[]
    for a,b in zip(starts,stops):
        if t[b]-t[a]<min_duration_h:continue
        dose=float(np.trapezoid(x[a:b+1],t[a:b+1]))
        if dose<min_dose:continue
        peak=a+int(np.argmax(x[a:b+1]))
        rows.append(dict(start_h=float(t[a]),end_h=float(t[b]),peak_h=float(t[peak]),
                         peak=float(x[peak]),dose=dose))
    return rows


def signals(y,zmh,receiver):
    p=receiver.p
    s0=y[IDX['b0_S']]; i0=y[IDX['b0_I']]; r0=y[IDX['b0_R']]; c0=y[IDX['b0_C']]
    a0=y[IDX['A0_zmh']]; f0=y[IDX['F0_zmh']]
    pb1=y[IDX['pb1_zmh']]; i1=y[IDX['I1_zmh']]; rdf1=y[IDX['RDF1_zmh']]
    a1=y[IDX['A1']]; f1=y[IDX['F1']]
    s2=y[IDX['b2_S']];i2=y[IDX['b2_I']];r2=y[IDX['b2_R']];c2=y[IDX['b2_C']]
    g0=hill_array(a0,zmh['K_A0'],zmh['n_A0'])*(1-hill_array(f0,zmh['K_R0'],zmh['n_R0']))
    clock=hill_array(receiver.c.clock_scale*i0,receiver.c.clock_K_au,receiver.c.clock_n)
    g1=hill_array(a1,p['K_A'][1],receiver.c.n_A1_gate)*(1-hill_array(f1,p['K_F'][1],p['n_F'][1]))*clock
    rev0=p['k_rev']*hill_array(c0,receiver.k_complex,2)*s0
    rev1=60*zmh['k_rev']*(1-pb1)*hill_array(i1*rdf1,zmh['K_D_comp'],2)
    rev2=p['k_rev']*s2*hill_array(c2,receiver.k_complex,2)
    fwd2=p['k_fwd']*(1-s2)*hill_array(i2,p['K_D_int'][2],2)*p['K_inh']/(p['K_inh']+np.maximum(r2,0))
    return dict(g0=g0,g1=g1,clock=clock,J_rev0=rev0,J_rev1=rev1,J_rev2=rev2,J_fwd2=fwd2,
                S0=s0,S1=1-pb1,S2=s2)


def read_windows(t,s):
    int0=s['int0']
    amp=np.ptp(int0)
    peaks=find_peaks(int0,prominence=.1*amp)[0] if amp>1e-12 else np.array([],dtype=int)
    reads=[]
    for cycle,(a,b) in enumerate(zip(peaks[:-1],peaks[1:])):
        trough=a+int(np.argmin(int0[a:b+1]))
        half=.1*(t[b]-t[a])
        start,end=t[trough]-half,t[trough]+half
        covered=bool(start>=t[a] and end<=t[b])
        mask=(t>=start)&(t<=end)
        labels=[];commit=[];extremes=[]
        for name in ('S0','S1','S2'):
            values=s[name][mask]
            low=float(np.mean(values<=LOW));high=float(np.mean(values>=HIGH))
            labels.append(0 if low>=OCC else 1 if high>=OCC else None)
            commit.append(max(low,high))
            extremes.append([float(values.min()),float(values.max())])
        value=sum((2**i)*label for i,label in enumerate(labels)) if covered and all(v is not None for v in labels) else None
        reads.append(dict(cycle=cycle,cycle_start_h=float(t[a]),cycle_end_h=float(t[b]),
                          trough_h=float(t[trough]),start_h=float(start),end_h=float(end),
                          covered=covered,samples=int(mask.sum()),labels=labels,
                          commitment=commit,extremes=extremes,value=value))
    return peaks,reads


def read_verdict(reads,drop):
    chosen=reads[drop:];vals=[r['value'] for r in chosen]
    enough=len(vals)>=MIN_STEADY if drop else len(vals)>=8
    complete=bool(enough and all(v is not None for v in vals))
    increment=bool(complete and all((b-a)%8==1 for a,b in zip(vals[:-1],vals[1:])))
    return dict(drop=drop,reads=len(vals),complete=complete,increments_mod8=increment,
                passed=increment,sequence=''.join('x' if v is None else str(v) for v in vals),
                minimum_commitment=min((min(r['commitment']) for r in chosen),default=None),
                boundary_clips=sum(not r['covered'] for r in chosen))


def margins(reads,cross):
    times=np.array([c['time_h'] for c in cross],dtype=float)
    setup=[];hold=[]
    for r in reads:
        before=times[times<=r['start_h']]
        after=times[times>=r['end_h']]
        if len(before):setup.append(float(r['start_h']-before[-1]))
        if len(after):hold.append(float(after[0]-r['end_h']))
    return dict(min_setup_h=min(setup) if setup else None,
                min_hold_h=min(hold) if hold else None,
                setup_samples=len(setup),hold_samples=len(hold))


def associate(reverse,gates,flips,hours):
    rows=[];used_g=set();used_f=set()
    for j,rev in enumerate(reverse):
        left=0 if j==0 else (reverse[j-1]['peak_h']+rev['peak_h'])/2
        right=hours if j+1==len(reverse) else (rev['peak_h']+reverse[j+1]['peak_h'])/2
        gs=[k for k,g in enumerate(gates) if left<=g['peak_h']<right]
        fs=[k for k,f in enumerate(flips) if left<=f['time_h']<right]
        if len(gs)==1:used_g.add(gs[0])
        if len(fs)==1:used_f.add(fs[0])
        gate=gates[gs[0]] if len(gs)==1 else None
        flip=flips[fs[0]] if len(fs)==1 else None
        rows.append(dict(reverse_id=j,reverse_start_h=rev['start_h'],reverse_peak_h=rev['peak_h'],
                         gate_count=len(gs),flip_count=len(fs),
                         gate_start_h=gate['start_h'] if gate else None,
                         flip_h=flip['time_h'] if flip else None,
                         gate_after_reverse=bool(gate and gate['start_h']>=rev['start_h']-EVENT_TOL_H),
                         flip_after_reverse=bool(flip and flip['time_h']>=rev['start_h']-EVENT_TOL_H),
                         flip_direction=flip['direction'] if flip else None))
    # Burn in two carry events, not two read windows; this is explicitly part of v1.
    late=rows[2:]
    unique=bool(len(late)>=3 and all(r['gate_count']==r['flip_count']==1 for r in late)
                and len(used_g)==len(gates) and len(used_f)==len(flips))
    order=bool(late and all(r['gate_after_reverse'] and r['flip_after_reverse'] for r in late))
    directions=[r['flip_direction'] for r in late]
    alt=bool(len(directions)>=3 and all(v is not None for v in directions)
             and all(a!=b for a,b in zip(directions[:-1],directions[1:])))
    return dict(reverse_events=len(reverse),gate_events=len(gates),flip_events=len(flips),
                associations=rows,late_associations=len(late),
                one_to_one=unique,causal_order=order,alternating_directions=alt,
                unassigned_gate_events=[i for i in range(len(gates)) if i not in used_g],
                unassigned_flip_events=[i for i in range(len(flips)) if i not in used_f],
                passed=bool(unique and order and alt))


def analyse(t,y,zmh,receiver):
    if y.shape!=(len(NAMES),len(t)):raise ValueError('Trajectory does not match HZH layout')
    sig=signals(y,zmh,receiver);sig['int0']=y[IDX['b0_I']]
    peaks,reads=read_windows(t,sig)
    cross={b:crossings(t,sig[b]) for b in ('S0','S1','S2')}
    events={}
    for stage,revname,gatename,bit in ((0,'J_rev0','g0','S1'),(1,'J_rev1','g1','S2')):
        rev=segments(t,sig[revname],JREV_THRESHOLD_H)
        gates=segments(t,sig[gatename],GATE_THRESHOLD,min_dose=GATE_MIN_DOSE_H)
        events[f'bit{stage}_to_bit{stage+1}']=associate(rev,gates,cross[bit],float(t[-1]))
    steady=read_verdict(reads,DROP)
    bit_margins={b:margins(reads,cross[b]) for b in cross}
    nonnull=[x for m in bit_margins.values() for x in (m['min_setup_h'],m['min_hold_h']) if x is not None]
    # Historical gate peak contrast is intentionally NOT an acceptance arm.
    verdict=dict(version='HZH_MOD8_CAUSAL_V1',hours=float(t[-1]),
       rule=dict(window='Int0 peak-to-peak trough, total width 20% period',
                 bands=[LOW,HIGH],occupancy=OCC,drop=DROP,min_steady_reads=MIN_STEADY,
                 reverse_flux_threshold_per_h=JREV_THRESHOLD_H,gate_threshold=GATE_THRESHOLD,
                 min_event_h=PULSE_DURATION_H,min_gate_dose_h=GATE_MIN_DOSE_H,
                 association='nearest reverse-event peak interval; causal order separately checked',
                 initial_carry_events_dropped=2,
                 acceptance='steady mod8 + both event chains one-to-one, ordered, direction alternating'),
       cold=read_verdict(reads,0),steady=steady,clock_peak_count=len(peaks),read_windows=reads,
       crossings=cross,events=events,bit_margins=bit_margins,
       global_min_timing_margin_h=min(nonnull) if nonnull else None,
       signal_ranges={k:[float(np.min(v)),float(np.max(v))] for k,v in sig.items()},
       certified_v1=bool(steady['passed'] and all(e['passed'] for e in events.values())),
       legacy_certified=None,
       scope='Deterministic 300 h HZH v1; no inheritance of 51-state gate-contrast predicate; no experimental reliability claim')
    return verdict,sig
