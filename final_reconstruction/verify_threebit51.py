"""Finite-window and causal verifier for the formal 51-state counter."""
from __future__ import annotations

import numpy as np

from verify_bit0_part2 import clock_cycles
from verify_twobit_causal import (_crossings, _segments, _setup_hold_margins,
                                  _window_label, EVENT_TOL_H, GATE_ON,
                                  JREV0_ON_PER_H, MIN_GATE_DOSE_H,
                                  READ_WINDOW_FRACTION)

MAX_OFF_ON_GATE_RATIO = 0.10


def read_windows_three(t, flux, states):
    peaks=clock_cycles(t,flux); rows=[]
    for cycle,(a,b) in enumerate(zip(peaks[:-1],peaks[1:])):
        trough=a+int(np.argmin(flux[a:b+1])); half=max(1,int(round(READ_WINDOW_FRACTION*(b-a)/2)))
        wa,wb=max(a,trough-half),min(b,trough+half)
        labels=[_window_label(s,wa,wb) for s in states]
        value=None if any(x['label']=='x' for x in labels) else sum((1<<i)*int(x['label']) for i,x in enumerate(labels))
        rows.append(dict(cycle=cycle,cycle_start_h=float(t[a]),cycle_end_h=float(t[b]),
                         trough_h=float(t[trough]),window_start_h=float(t[wa]),window_end_h=float(t[wb]),
                         bits=labels,value=value))
    return rows


def counter_verdict(reads,drop=0):
    chosen=reads[drop:]; vals=[r['value'] for r in chosen]
    valid=bool(len(vals)>=16 and all(v is not None for v in vals))
    increments=bool(valid and all((vals[i+1]-vals[i])%8==1 for i in range(len(vals)-1)))
    return dict(drop_reads=drop,reads=len(vals),valid_windows=valid,increments_mod8=increments,
                passed=bool(valid and increments),sequence=''.join('x' if v is None else str(v) for v in vals),
                minimum_commitment=(float(min(min(bit['commitment'] for bit in r['bits']) for r in chosen))
                                    if chosen else 0.0))


def analyse_threebit(model,sol,hours):
    t,y=sol.t,sol.y
    flux=np.asarray([model.flux(y[:,k]) for k in range(y.shape[1])])
    sig=model.diagnostic_signals(y); states=(y[16],y[27],y[44])
    reads=read_windows_three(t,flux,states)
    reverse=_segments(t,sig['J_rev1'],JREV0_ON_PER_H)
    gates=_segments(t,sig['g1'],GATE_ON,min_area=MIN_GATE_DOSE_H)
    s2cross=_crossings(t,y[44]); associations=[]; used_gates=set();used_cross=set()
    for i,event in enumerate(reverse):
        left=0.0 if i==0 else (reverse[i-1]['peak_h']+event['peak_h'])/2
        right=hours if i+1==len(reverse) else (event['peak_h']+reverse[i+1]['peak_h'])/2
        gs=[j for j,g in enumerate(gates) if left<=g['peak_h']<right]
        cs=[j for j,c in enumerate(s2cross) if left<=c['time_h']<right]
        gi=gs[0] if len(gs)==1 else None; ci=cs[0] if len(cs)==1 else None
        if gi is not None:used_gates.add(gi)
        if ci is not None:used_cross.add(ci)
        gate=gates[gi] if gi is not None else None; cross=s2cross[ci] if ci is not None else None
        associations.append(dict(reverse_event=i,gate_event=gi,bit2_crossing=ci,
                                 gate_count=len(gs),bit2_crossing_count=len(cs),
                                 gate_after_reverse_start=(gate['start_h']>=event['start_h']-EVENT_TOL_H if gate else None),
                                 bit2_after_reverse_start=(cross['time_h']>=event['start_h']-EVENT_TOL_H if cross else None),
                                 gate_delay_h=(gate['start_h']-event['start_h'] if gate else None),
                                 bit2_delay_h=(cross['time_h']-event['start_h'] if cross else None),
                                 bit2_direction=(cross['direction'] if cross else None)))
    late=associations[2:] if len(associations)>4 else associations
    one=bool(late and all(a['gate_count']==1 and a['bit2_crossing_count']==1 for a in late)
             and not(set(range(len(gates)))-used_gates) and not(set(range(len(s2cross)))-used_cross))
    order=bool(late and all(a['gate_after_reverse_start'] and a['bit2_after_reverse_start'] for a in late))
    dirs=[a['bit2_direction'] for a in late if a['bit2_direction']]
    alt=bool(len(dirs)==len(late) and all(dirs[i]!=dirs[i+1] for i in range(len(dirs)-1)))
    crossings=[_crossings(t,s) for s in states]
    cold=counter_verdict(reads,0);steady=counter_verdict(reads,8)
    cycle_gate=[]
    for r in reads:
        a=int(np.searchsorted(t,r['cycle_start_h'])); b=int(np.searchsorted(t,r['cycle_end_h']))
        g=sig['g1'][a:b+1]
        cycle_gate.append(dict(cycle=r['cycle'],bit1_label=r['bits'][1]['label'],
                               g1_min=float(g.min()),g1_median=float(np.median(g)),
                               g1_p90=float(np.percentile(g,90)),g1_peak=float(g.max()),
                               g1_dose=float(np.trapezoid(g,t[a:b+1])),
                               carry_cycle=bool(g.max()>=GATE_ON)))
    on=[r for r in cycle_gate if r['carry_cycle']]
    off=[r for r in cycle_gate if not r['carry_cycle']]
    on_peak=float(np.median([r['g1_peak'] for r in on])) if on else 0.0
    def ratio(field,mode='max'):
        if not off or on_peak<=0:return None
        vals=[r[field] for r in off]
        value=max(vals) if mode=='max' else float(np.median(vals))
        return float(value/on_peak)
    off_peak_ratio=ratio('g1_peak')
    contrast=dict(on_cycle_peak_median=on_peak,
                  off_cycle_min_max=(max(r['g1_min'] for r in off) if off else None),
                  off_cycle_median_max=(max(r['g1_median'] for r in off) if off else None),
                  off_cycle_peak_max=(max(r['g1_peak'] for r in off) if off else None),
                  off_min_to_on_peak=ratio('g1_min'),
                  off_median_to_on_peak=ratio('g1_median'),
                  off_peak_to_on_peak=off_peak_ratio,
                  threshold=MAX_OFF_ON_GATE_RATIO,
                  passed=bool(off_peak_ratio is not None and off_peak_ratio<=MAX_OFF_ON_GATE_RATIO))
    return dict(hours=hours,read_windows=reads,cold_start=cold,steady_state=steady,
                signal_ranges={k:[float(np.min(v)),float(np.max(v))] for k,v in sig.items()},
                bit_margins={f'bit{i}':_setup_hold_margins(reads,crossings[i]) for i in range(3)},
                bit1_reverse_events=reverse,carry1_gate_events=gates,bit2_crossings=s2cross,
                cycle_gate_metrics=cycle_gate,gate_contrast=contrast,
                associations=associations,causal_verdict=dict(
                    exactly_one_gate_and_flip_per_late_reverse=one,
                    causal_order_after_reverse_start=order,bit2_directions_alternate=alt,
                    unassigned_gate_events=sorted(set(range(len(gates)))-used_gates),
                    unassigned_bit2_crossings=sorted(set(range(len(s2cross)))-used_cross)),
                certified=bool(steady['passed'] and one and order and alt and contrast['passed']))
