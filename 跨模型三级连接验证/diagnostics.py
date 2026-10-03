"""Predeclared exploratory readout. Not the original certification predicate."""
import numpy as np
from scipy.signal import find_peaks
from hybrid_model import INDEX


def crossings(t, s):
    side = np.asarray(s) >= .5
    ids = np.flatnonzero(side[1:] != side[:-1])
    return [dict(time_h=float(t[i]+(.5-s[i])*(t[i+1]-t[i])/(s[i+1]-s[i])),
                 direction='up' if side[i+1] else 'down') for i in ids]


def state_health(y, names, fraction_names):
    return dict(finite=bool(np.isfinite(y).all()),
                minimum_by_state={n:float(y[i].min()) for i,n in enumerate(names)},
                maximum_by_state={n:float(y[i].max()) for i,n in enumerate(names)},
                fractions_within_tolerance=bool(all(y[names.index(n)].min()>=-1e-7 and
                                                    y[names.index(n)].max()<=1+1e-7 for n in fraction_names)),
                nonnegative_within_tolerance=bool(y.min()>=-1e-7))


def analyse(t,y,clock):
    states = [1-y[INDEX['pb0']],1-y[INDEX['pb1']],y[INDEX['b2_S']]]
    amplitude = float(np.ptp(clock))
    peaks = find_peaks(clock,prominence=.1*amplitude)[0] if amplitude>1e-10 else np.array([],dtype=int)
    reads=[]
    for a,b in zip(peaks[:-1],peaks[1:]):
        v = a+int(np.argmin(clock[a:b+1]))
        half=.1*(t[b]-t[a])
        left,right=t[v]-half,t[v]+half
        clipped=bool(left<t[a] or right>t[b])
        mask=(t>=left)&(t<=right)
        labels, commitments=[],[]
        for s in states:
            lo,hi=float(np.mean(s[mask]<=.3)),float(np.mean(s[mask]>=.7))
            labels.append(0 if lo>=.8 else 1 if hi>=.8 else None)
            commitments.append(max(lo,hi))
        usable=not clipped and all(x is not None for x in labels)
        value=int(sum((2**i)*x for i,x in enumerate(labels))) if usable else None
        reads.append(dict(trough_h=float(t[v]),window_start_h=float(left),window_end_h=float(right),
                          boundary_crossed=clipped,samples=int(mask.sum()),labels=labels,
                          commitment=commitments,value=value))
    def score(drop,modulus):
        rr=reads[drop:]
        vals=[r['value'] if modulus==8 else
              (r['labels'][0]+2*r['labels'][1] if not r['boundary_crossed'] and
               all(x is not None for x in r['labels'][:2]) else None) for r in rr]
        increments=bool(len(vals)>=2 and all(v is not None for v in vals) and
                        all((b-a)%modulus==1 for a,b in zip(vals[:-1],vals[1:])))
        return dict(drop=drop,reads=len(vals),sequence=''.join('x' if v is None else str(v) for v in vals),
                    observed_increments=increments,enough_reads=len(vals)>=2*modulus)
    return dict(status='EXPLORATORY_ONLY_NOT_CERTIFIED',
                rule='Int0 peak prominence=10% full range; trough between peaks; total window 20% period; no phase fitting',
                reads=reads,clock_peak_times_h=t[peaks].tolist(),
                cold_mod8=score(0,8),steady_mod8=score(8,8),
                cold_mod4=score(0,4),steady_mod4=score(4,4),
                crossings={f'bit{i}':crossings(t,s) for i,s in enumerate(states)},
                causal_certification=None,paired_leakage=None,
                limitation='Formal carry association and paired leakage are not implemented in this preliminary version.')
