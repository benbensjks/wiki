"""Stage-1 scan of new A1/F1 expression parameters in the 51-state model."""
from __future__ import annotations

import argparse,csv,itertools,json,time,traceback
from pathlib import Path
import numpy as np
import pandas as pd

from model import ROOT
from model_threebit51 import ThreeBit51Model,ThreeBitCarryParameters
from model_twobit34 import CarryExpressionParameters
from verify_threebit51 import analyse_threebit


OUT=ROOT/'threebit51_results'/'carry1_scan'
MRNA=(0.5,1.0,2.0,4.0,8.0)
MAT=(5.0,10.0,20.0,30.0,40.0,50.0,60.0)


def jobs():return [dict(carry1_mrna_min=a,carry1_maturation_min=b) for a,b in itertools.product(MRNA,MAT)]


def evaluate(job,hours=600.0):
    t0=time.perf_counter();row=dict(job)
    try:
        default=ThreeBitCarryParameters()
        c1=CarryExpressionParameters(mrna_half_life_min=job['carry1_mrna_min'],
             activator_maturation_half_life_min=job['carry1_maturation_min'],
             repressor_maturation_half_life_min=job['carry1_maturation_min'])
        model=ThreeBit51Model(carry=ThreeBitCarryParameters(carry0=default.carry0,carry1=c1))
        sol=model.simulate(hours=hours,sample_min=2,max_step_min=2)
        a=analyse_threebit(model,sol,hours);c=a['causal_verdict'];g=a['gate_contrast'];sig=a['signal_ranges']
        row.update(solver_success=True,error='',finite=bool(np.isfinite(sol.y).all()),
                   cold_passed=a['cold_start']['passed'],steady_passed=a['steady_state']['passed'],
                   cold_sequence=a['cold_start']['sequence'],steady_sequence=a['steady_state']['sequence'],
                   minimum_commitment=a['steady_state']['minimum_commitment'],certified=a['certified'],
                   one_to_one=c['exactly_one_gate_and_flip_per_late_reverse'],
                   causal_order=c['causal_order_after_reverse_start'],directions_alternate=c['bit2_directions_alternate'],
                   gate_contrast_passed=g['passed'],off_on_peak_ratio=g['off_peak_to_on_peak'],
                   on_g1_peak=g['on_cycle_peak_median'],g1_peak=sig['g1'][1],Int2_source_peak=sig['Int2_source'][1],
                   S2_min=sig['S2'][0],S2_max=sig['S2'][1],bit1_reverse_events=len(a['bit1_reverse_events']),
                   gate_events=len(a['carry1_gate_events']),bit2_crossings=len(a['bit2_crossings']),
                   bit2_setup_h=a['bit_margins']['bit2']['min_setup_h'],bit2_hold_h=a['bit_margins']['bit2']['min_hold_h'])
    except Exception as exc:  # noqa: BLE001
        row.update(solver_success=False,finite=False,certified=False,error=repr(exc),error_trace=traceback.format_exc(limit=3))
    row['runtime_s']=round(time.perf_counter()-t0,3);return row


def scan(args):
    OUT.mkdir(parents=True,exist_ok=True);allj=jobs();mine=allj[args.shard::args.nshards]
    tag=f'shard{args.shard:02d}of{args.nshards:02d}';path=OUT/f'carry1_{tag}.csv';rows=[]
    for i,j in enumerate(mine):
        r=evaluate(j,args.hours);rows.append(r);fields=sorted({k for x in rows for k in x})
        with path.open('w',newline='',encoding='utf-8') as fh:
            w=csv.DictWriter(fh,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
        print(f"[{tag}] {i+1}/{len(mine)} mRNA={j['carry1_mrna_min']} mat={j['carry1_maturation_min']} "
              f"cert={r['certified']} crossings={r.get('bit2_crossings','?')}",flush=True)
    (OUT/f'meta_{tag}.json').write_text(json.dumps(dict(shard=args.shard,nshards=args.nshards,hours=args.hours,
        points=len(mine),expected_total=len(allj))),encoding='utf-8')


def as_bool(s):return s.map(lambda x:str(x).lower() in ('true','1'))


def merge(args):
    frames=[]
    for i in range(args.nshards):
        p=OUT/f'carry1_shard{i:02d}of{args.nshards:02d}.csv'
        if not p.exists():raise SystemExit(f'missing {p}')
        frames.append(pd.read_csv(p))
    d=pd.concat(frames,ignore_index=True).sort_values(['carry1_mrna_min','carry1_maturation_min'])
    for c in ('solver_success','finite','cold_passed','steady_passed','certified','one_to_one','causal_order','directions_alternate','gate_contrast_passed'):
        d[c]=as_bool(d[c])
    d.to_csv(OUT/'carry1_all.csv',index=False)
    # If no point is certified, ranking still exposes the most complete bit2 response.
    d['rank_crossing_error']=(d.bit1_reverse_events-d.bit2_crossings).abs()
    best=d.sort_values(['certified','rank_crossing_error','minimum_commitment','off_on_peak_ratio'],
                       ascending=[False,True,False,True]).head(10)
    summary=dict(rows=len(d),expected=len(jobs()),solver_failures=int((~d.solver_success).sum()),
                 certified=int(d.certified.sum()),steady_passed=int(d.steady_passed.sum()),
                 gate_contrast_passed=int(d.gate_contrast_passed.sum()),best_points=best.to_dict(orient='records'))
    (OUT/'carry1_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='mode',required=True)
    s=sub.add_parser('scan');s.add_argument('--shard',type=int,required=True);s.add_argument('--nshards',type=int,default=7);s.add_argument('--hours',type=float,default=600)
    m=sub.add_parser('merge');m.add_argument('--nshards',type=int,default=7)
    a=ap.parse_args();scan(a) if a.mode=='scan' else merge(a)


if __name__=='__main__':main()
