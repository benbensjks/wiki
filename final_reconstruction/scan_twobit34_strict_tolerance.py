"""Strict-tolerance audit at the selected centre and pass/fail boundaries."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import time
import traceback
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from model import ROOT
from model_twobit34 import CarryExpressionParameters, TwoBit34Model, nominal_extension
from verify_twobit_causal import analyse_solution


PROFILE = ROOT/'twobit34_results'/'selected_profile.json'
OUT = ROOT/'twobit34_results'/'strict_tolerance'
POINTS = (
    dict(name='selected_centre', kind='centre', uM=5.75, mat=32.5, mrna=2.0),
    dict(name='pass_low_uM_edge', kind='pass_boundary', uM=5.5, mat=30.0, mrna=2.0),
    dict(name='pass_low_mat_edge', kind='pass_boundary', uM=6.0, mat=25.0, mrna=2.0),
    dict(name='pass_notch_lip', kind='pass_boundary', uM=6.25, mat=35.0, mrna=2.0),
    dict(name='pass_high_corner', kind='pass_boundary', uM=6.5, mat=25.0, mrna=4.0),
    dict(name='fail_low_uM', kind='fail_boundary', uM=5.5, mat=27.5, mrna=2.0),
    dict(name='fail_low_mat', kind='fail_boundary', uM=5.75, mat=25.0, mrna=2.0),
    dict(name='fail_notch_mid', kind='fail_boundary', uM=6.25, mat=32.5, mrna=2.0),
    dict(name='fail_notch_low', kind='fail_boundary', uM=6.25, mat=30.0, mrna=2.0),
)
SETTINGS = (
    dict(name='baseline', rtol=2e-7, atol=2e-9, max_step_min=2.0),
    dict(name='tight', rtol=1e-9, atol=1e-11, max_step_min=1.0),
    dict(name='ultra', rtol=1e-11, atol=1e-13, max_step_min=0.5),
)


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as fh:
        for chunk in iter(lambda:fh.read(1<<20),b''):h.update(chunk)
    return h.hexdigest().upper()


def jobs(): return [(i,p,s) for i,p in enumerate(POINTS) for s in SETTINGS]


def reason(a):
    if not a['steady_state']['valid_windows']: return 'readout_unlabelled_window'
    if not a['steady_state']['increments_mod4']: return 'readout_wrong_increment'
    c=a['causal_verdict']
    if not c['exactly_one_gate_and_flip_per_late_reverse']: return 'carry_not_one_to_one'
    if not c['causal_order_after_reverse_start']: return 'causal_order_failed'
    if not c['bit1_directions_alternate']: return 'bit1_direction_failed'
    return ''


def evaluate(point_id, point, setting, hours=300.0):
    t0=time.perf_counter(); row=dict(point_id=point_id,point_name=point['name'],kind=point['kind'],
                                     uM_per_au=point['uM'],carry_maturation_half_life_min=point['mat'],
                                     carry_mrna_half_life_min=point['mrna'],setting=setting['name'],
                                     rtol=setting['rtol'],atol=setting['atol'],max_step_min=setting['max_step_min'])
    try:
        e=replace(nominal_extension(),uM_per_au=point['uM'])
        c=CarryExpressionParameters(mrna_half_life_min=point['mrna'],
                                    activator_maturation_half_life_min=point['mat'],
                                    repressor_maturation_half_life_min=point['mat'])
        m=TwoBit34Model(e,c)
        sol=m.simulate(hours=hours,sample_min=2.0,rtol=setting['rtol'],atol=setting['atol'],
                       max_step_min=setting['max_step_min'])
        a=analyse_solution(m.base,sol,vars(e),hours); cv=a['causal_verdict']
        row.update(solver_success=True,error='',finite=bool(np.isfinite(sol.y).all()),
                   minimum_state=float(sol.y.min()),maximum_state=float(sol.y.max()),
                   cold_passed=bool(a['cold_start']['passed']),steady_passed=bool(a['steady_state']['passed']),
                   certified=bool(a['certified']),failure_reason=reason(a),
                   cold_sequence=a['cold_start']['sequence'],steady_sequence=a['steady_state']['sequence'],
                   minimum_commitment=float(a['steady_state']['minimum_commitment']),
                   one_to_one=bool(cv['exactly_one_gate_and_flip_per_late_reverse']),
                   causal_order=bool(cv['causal_order_after_reverse_start']),
                   directions_alternate=bool(cv['bit1_directions_alternate']),
                   bit0_min=float(a['signal_ranges']['S0'][0]),bit0_max=float(a['signal_ranges']['S0'][1]),
                   bit1_min=float(a['signal_ranges']['S1'][0]),bit1_max=float(a['signal_ranges']['S1'][1]),
                   Jrev0_peak=float(a['signal_ranges']['Jrev0'][1]),g0_peak=float(a['signal_ranges']['g0'][1]),
                   Int1_source_peak=float(a['signal_ranges']['Int1_source_peak']),
                   bit1_setup_h=a['margins']['bit1']['min_setup_h'],bit1_hold_h=a['margins']['bit1']['min_hold_h'])
    except Exception as exc:  # noqa: BLE001
        row.update(solver_success=False,finite=False,certified=False,failure_reason='integration_failed',
                   error=repr(exc),error_trace=traceback.format_exc(limit=3))
    row['runtime_s']=round(time.perf_counter()-t0,3); return row


def scan(args):
    OUT.mkdir(parents=True,exist_ok=True); all_jobs=jobs(); mine=all_jobs[args.shard::args.nshards]
    tag=f'shard{args.shard:02d}of{args.nshards:02d}'; path=OUT/f'strict_{tag}.csv'; rows=[]
    for i,(pid,p,s) in enumerate(mine):
        row=evaluate(pid,p,s,args.hours); rows.append(row); fields=sorted({k for r in rows for k in r})
        with path.open('w',newline='',encoding='utf-8') as fh:
            w=csv.DictWriter(fh,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
        print(f"[{tag}] {i+1}/{len(mine)} {p['name']} {s['name']} certified={row['certified']} "
              f"reason={row['failure_reason'] or 'pass'}",flush=True)
    (OUT/f'meta_{tag}.json').write_text(json.dumps(dict(shard=args.shard,nshards=args.nshards,
        hours=args.hours,points=len(mine),expected_total=len(all_jobs),profile_sha256=sha256(PROFILE),
        source_sha256=dict(model_py=sha256(ROOT/'model.py'),model34_py=sha256(ROOT/'model_twobit34.py'),
                           verifier_py=sha256(ROOT/'verify_twobit_causal.py'),scanner_py=sha256(Path(__file__)))),indent=2),encoding='utf-8')


def as_bool(s): return s.map(lambda x:str(x).strip().lower() in ('true','1'))


def merge(args):
    frames=[]
    for i in range(args.nshards):
        p=OUT/f'strict_shard{i:02d}of{args.nshards:02d}.csv'
        if not p.exists(): raise SystemExit(f'missing {p}')
        frames.append(pd.read_csv(p))
    d=pd.concat(frames,ignore_index=True).sort_values(['point_id','setting']).reset_index(drop=True)
    for c in ('solver_success','finite','cold_passed','steady_passed','certified','one_to_one','causal_order','directions_alternate'):
        d[c]=as_bool(d[c])
    d.to_csv(OUT/'strict_tolerance_all.csv',index=False)
    numeric=('minimum_state','maximum_state','minimum_commitment','bit0_min','bit0_max','bit1_min','bit1_max',
             'Jrev0_peak','g0_peak','Int1_source_peak','bit1_setup_h','bit1_hold_h')
    points=[]
    for pid,g in d.groupby('point_id'):
        g=g.set_index('setting'); ref=g.loc['ultra']; comparisons=[]
        for name in ('baseline','tight'):
            r=g.loc[name]; diffs={k:abs(float(r[k])-float(ref[k])) for k in numeric}
            comparisons.append(dict(setting=name,same_certified=bool(r.certified==ref.certified),
                                    same_reason=str(r.failure_reason)==str(ref.failure_reason),
                                    same_cold_sequence=str(r.cold_sequence)==str(ref.cold_sequence),
                                    same_steady_sequence=str(r.steady_sequence)==str(ref.steady_sequence),
                                    max_numeric_difference=max(diffs.values()),differences=diffs))
        converged=all(x['same_certified'] and x['same_reason'] and x['same_cold_sequence'] and
                      x['same_steady_sequence'] and x['max_numeric_difference']<=0.005 for x in comparisons)
        points.append(dict(point_id=int(pid),point_name=str(ref.point_name),kind=str(ref.kind),
                           uM=float(ref.uM_per_au),mat=float(ref.carry_maturation_half_life_min),
                           mrna=float(ref.carry_mrna_half_life_min),ultra_certified=bool(ref.certified),
                           ultra_reason='' if pd.isna(ref.failure_reason) else str(ref.failure_reason),
                           comparisons=comparisons,converged=bool(converged)))
    summary=dict(rows=len(d),expected=len(jobs()),solver_failures=int((~d.solver_success).sum()),
                 points=points,all_points_converged=all(p['converged'] for p in points))
    (OUT/'strict_tolerance_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    md=['# 中心与边界严格容差复核','',f"- 轨迹：{len(d)}/{len(jobs())}",
        f"- 积分失败：{summary['solver_failures']}",f"- 全部点结论收敛：**{summary['all_points_converged']}**",'',
        '| 点 | 类型 | uM/mat/mRNA | ultra结论 | 收敛 | baseline相对ultra最大数值差 |','|---|---|---|---|:---:|---:|']
    for p in points:
        md.append(f"| {p['point_name']} | {p['kind']} | {p['uM']:g}/{p['mat']:g}/{p['mrna']:g} | "
                  f"{'PASS' if p['ultra_certified'] else p['ultra_reason']} | {p['converged']} | "
                  f"{p['comparisons'][0]['max_numeric_difference']:.3g} |")
    (OUT/'strict_tolerance_summary.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    manifest={}
    for path in sorted(OUT.glob('*')):
        if path.is_file() and not path.name.startswith('SHA256'):
            manifest[str(path.relative_to(ROOT))]=sha256(path)
    for path in (PROFILE,ROOT/'model.py',ROOT/'model_twobit34.py',ROOT/'verify_twobit_causal.py',Path(__file__)):
        manifest[str(path.relative_to(ROOT))]=sha256(path)
    (OUT/'SHA256SUMS.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2,sort_keys=True),encoding='utf-8')
    (OUT/'SHA256SUMS.txt').write_text(''.join(f'{v}  {k}\n' for k,v in sorted(manifest.items())),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='mode',required=True)
    s=sub.add_parser('scan');s.add_argument('--shard',type=int,required=True);s.add_argument('--nshards',type=int,default=9);s.add_argument('--hours',type=float,default=300.0)
    m=sub.add_parser('merge');m.add_argument('--nshards',type=int,default=9)
    args=ap.parse_args();scan(args) if args.mode=='scan' else merge(args)


if __name__=='__main__':main()
