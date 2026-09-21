"""Sharded molecular-pool perturbation scan for the selected 34-state model.

Priority/order encoded in GROUP_ORDER:
S0/S1 -> RDF0/RDF1 -> A0/F0 -> Int1.

Initial conditions are the four phase-consistent 34-state branches already
validated at the selected robust centre.  The script supports independent
scan and merge modes and never changes ZENG or model source files.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import platform
import sys
import time
import traceback
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

from model import ROOT
from model_twobit34 import CarryExpressionParameters, TwoBit34Model, nominal_extension
from verify_twobit_causal import analyse_solution


PROFILE = ROOT / 'twobit34_results' / 'selected_profile.json'
INITIALS = ROOT / 'twobit34_results' / 'initial_states_selected' / 'four_initial_states.csv'
OUT = ROOT / 'twobit34_results' / 'perturbations'
GROUP_ORDER = ('S', 'RDF', 'AFFL', 'INT1')
S_ACTIONS = ('delta_-0.20', 'delta_-0.10', 'delta_+0.10', 'delta_+0.20', 'flip')
FACTORS = (0.5, 0.8, 1.2, 1.5, 2.0)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def load_configuration():
    profile = json.loads(PROFILE.read_text(encoding='utf-8'))
    e = profile['extension']; c = profile['carry_expression']
    extension = replace(nominal_extension(), **e)
    carry = CarryExpressionParameters(**c)
    table = pd.read_csv(INITIALS).set_index('state')
    if tuple(table.index) != tuple(TwoBit34Model.state_names):
        raise ValueError('initial-state CSV order does not match the formal 34-state contract')
    states = {value: table[f'initial_{value:02b}'].to_numpy(dtype=float) for value in range(4)}
    return profile, extension, carry, states


def jobs(group='all'):
    selected = GROUP_ORDER if group == 'all' else (group.upper(),)
    out = []
    for initial in range(4):
        if 'S' in selected:
            for target, action in itertools.product(('S0', 'S1'), S_ACTIONS):
                out.append(dict(group='S', target=target, action=action, initial_value=initial))
        if 'RDF' in selected:
            for target, factor in itertools.product(('RDF0', 'RDF1'), FACTORS):
                out.append(dict(group='RDF', target=target, action=f'factor_{factor:g}',
                                factor=factor, initial_value=initial))
        if 'AFFL' in selected:
            for target, factor in itertools.product(('A0_pool', 'F0_pool'), FACTORS):
                out.append(dict(group='AFFL', target=target, action=f'factor_{factor:g}',
                                factor=factor, initial_value=initial))
        if 'INT1' in selected:
            for factor in FACTORS:
                out.append(dict(group='INT1', target='Int1_pool', action=f'factor_{factor:g}',
                                factor=factor, initial_value=initial))
    return out


def perturb(y0, job):
    y = y0.copy(); target = job['target']; action = job['action']
    if target in ('S0', 'S1'):
        idx = 16 if target == 'S0' else 27
        before = float(y[idx])
        if action == 'flip':
            y[idx] = 1.0 - y[idx]
        else:
            y[idx] = np.clip(y[idx] + float(action.replace('delta_', '')), 0.0, 1.0)
        indices = [idx]
    else:
        factor = float(job['factor'])
        indices = {'RDF0': [12, 13, 14], 'RDF1': [23, 24, 25],
                   'A0_pool': [30, 31, 28], 'F0_pool': [32, 33, 29],
                   'Int1_pool': [17, 18, 19]}[target]
        before = float(np.sum(y[indices])); y[indices] *= factor
    after = float(np.sum(y[indices])) if len(indices) > 1 else float(y[indices[0]])
    return y, indices, before, after


def expected_sequence(initial_value, length):
    return ''.join(str((initial_value + i + 1) % 4) for i in range(length))


def phase_offset(observed, expected, drop=4):
    a, b = observed[drop:], expected[drop:]
    if not a or len(a) != len(b) or 'x' in a:
        return None
    offsets = {(int(x) - int(y)) % 4 for x, y in zip(a, b)}
    return offsets.pop() if len(offsets) == 1 else None


def evaluate(model, initial_states, job, hours=200.0):
    started = time.perf_counter(); row = dict(job)
    try:
        y0, indices, before, after = perturb(initial_states[job['initial_value']], job)
        sol = model.simulate(hours=hours, sample_min=2.0, max_step_min=2.0, initial_state=y0)
        analysis = analyse_solution(model.base, sol, vars(model.e), hours)
        observed = analysis['cold_start']['sequence']
        expected = expected_sequence(job['initial_value'], len(observed))
        offset = phase_offset(observed, expected)
        c = analysis['causal_verdict']
        if analysis['certified'] and offset == 0:
            outcome = 'recovered_same_phase'
        elif analysis['certified'] and offset is not None:
            outcome = f'stable_phase_shift_{offset}'
        elif analysis['steady_state']['passed']:
            outcome = 'readout_pass_causal_fail'
        else:
            outcome = 'lost_counting'
        row.update(solver_success=True, error='', finite=bool(np.isfinite(sol.y).all()),
                   perturbed_indices=','.join(map(str, indices)), pool_before=before, pool_after=after,
                   relative_pool_change=(after / before if before > 1e-15 else np.nan),
                   cold_passed=bool(analysis['cold_start']['passed']),
                   steady_passed=bool(analysis['steady_state']['passed']),
                   certified=bool(analysis['certified']), outcome=outcome,
                   phase_offset=offset if offset is not None else np.nan,
                   expected_sequence=expected, observed_sequence=observed,
                   minimum_commitment=float(analysis['steady_state']['minimum_commitment']),
                   one_to_one=bool(c['exactly_one_gate_and_flip_per_late_reverse']),
                   causal_order=bool(c['causal_order_after_reverse_start']),
                   directions_alternate=bool(c['bit1_directions_alternate']),
                   reverse_events=len(analysis['reverse_events']), gate_events=len(analysis['gate_events']),
                   bit1_setup_h=analysis['margins']['bit1']['min_setup_h'],
                   bit1_hold_h=analysis['margins']['bit1']['min_hold_h'])
    except Exception as exc:  # noqa: BLE001
        row.update(solver_success=False, finite=False, certified=False,
                   outcome='integration_failed', error=repr(exc),
                   error_trace=traceback.format_exc(limit=3))
    row['runtime_s'] = round(time.perf_counter() - started, 3)
    return row


def scan(args):
    OUT.mkdir(parents=True, exist_ok=True)
    profile, extension, carry, initial_states = load_configuration()
    model = TwoBit34Model(extension, carry)
    all_jobs = jobs(args.group); mine = all_jobs[args.shard::args.nshards]
    tag = f'{args.group.lower()}_shard{args.shard:02d}of{args.nshards:02d}'
    path = OUT / f'perturb_{tag}.csv'; rows = []
    for i, job in enumerate(mine):
        row = evaluate(model, initial_states, job, args.hours); rows.append(row)
        fields = sorted({k for r in rows for k in r})
        with path.open('w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore'); w.writeheader(); w.writerows(rows)
        print(f"[{tag}] {i+1}/{len(mine)} init={job['initial_value']:02b} "
              f"{job['target']} {job['action']} -> {row['outcome']}", flush=True)
    meta = dict(group=args.group, shard=args.shard, nshards=args.nshards, hours=args.hours,
                points=len(mine), expected_total=len(all_jobs), profile=profile['profile_name'],
                python=sys.version, platform=platform.platform(),
                versions=dict(numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
                sha256=dict(profile=sha256(PROFILE), initials=sha256(INITIALS),
                            model_py=sha256(ROOT/'model.py'), model34_py=sha256(ROOT/'model_twobit34.py'),
                            verifier_py=sha256(ROOT/'verify_twobit_causal.py'), scanner_py=sha256(Path(__file__)),
                            csv=sha256(path)))
    (OUT / f'meta_{tag}.json').write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'[{tag}] wrote {path}')


def _as_bool(s):
    return s.map(lambda x: str(x).strip().lower() in ('true', '1'))


def merge(args):
    all_jobs = jobs(args.group); frames = []
    for shard in range(args.nshards):
        tag = f'{args.group.lower()}_shard{shard:02d}of{args.nshards:02d}'
        path = OUT / f'perturb_{tag}.csv'
        if not path.exists(): raise SystemExit(f'missing {path}')
        frames.append(pd.read_csv(path))
    data = pd.concat(frames, ignore_index=True)
    for c in ('solver_success','finite','cold_passed','steady_passed','certified',
              'one_to_one','causal_order','directions_alternate'):
        data[c] = _as_bool(data[c])
    keys = ('group','target','action','initial_value')
    expected = {(j['group'],j['target'],j['action'],j['initial_value']) for j in all_jobs}
    got = {(r.group,r.target,r.action,int(r.initial_value)) for r in data.itertuples()}
    completeness = dict(rows=len(data), expected=len(expected), unique=len(got),
                        duplicates=len(data)-len(got), missing=sorted(expected-got), unexpected=sorted(got-expected))
    data = data.sort_values(['group','target','action','initial_value']).reset_index(drop=True)
    out_csv = OUT / f'perturbations_{args.group.lower()}_all.csv'; data.to_csv(out_csv,index=False)
    summary = dict(completeness=completeness, solver_failures=int((~data.solver_success).sum()),
                   outcomes={str(k):int(v) for k,v in data.outcome.value_counts().items()},
                   by_group={str(g):{str(k):int(v) for k,v in sub.outcome.value_counts().items()}
                             for g,sub in data.groupby('group')},
                   by_target={str(g):{str(k):int(v) for k,v in sub.outcome.value_counts().items()}
                              for g,sub in data.groupby('target')})
    (OUT / f'perturbations_{args.group.lower()}_summary.json').write_text(
        json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')

    targets=list(data.target.unique()); fig,axes=plt.subplots(len(targets),1,figsize=(12,max(3,2.6*len(targets))),
                                                              squeeze=False,constrained_layout=True)
    code={'recovered_same_phase':0,'stable_phase_shift_1':1,'stable_phase_shift_2':2,
          'stable_phase_shift_3':3,'readout_pass_causal_fail':4,'lost_counting':5,'integration_failed':6}
    for ax,target in zip(axes[:,0],targets):
        sub=data[data.target==target].copy(); actions=list(dict.fromkeys(sub.action))
        z=np.full((4,len(actions)),6.0)
        for r in sub.itertuples(): z[int(r.initial_value),actions.index(r.action)]=code.get(r.outcome,6)
        ax.imshow(z,aspect='auto',vmin=0,vmax=6,cmap='tab10')
        ax.set_yticks(range(4),('00','01','10','11')); ax.set_xticks(range(len(actions)),actions,rotation=25,ha='right')
        ax.set_ylabel('initial'); ax.set_title(target)
        for i in range(4):
            for j in range(len(actions)):
                outcome=sub[(sub.initial_value==i)&(sub.action==actions[j])].iloc[0].outcome
                ax.text(j,i,outcome.replace('recovered_','').replace('stable_',''),ha='center',va='center',fontsize=6,color='white')
    fig.suptitle(f'Molecular-pool perturbation outcomes: {args.group}',fontsize=14)
    fig.savefig(OUT/f'perturbations_{args.group.lower()}_heatmap.png',dpi=220,bbox_inches='tight',facecolor='white')
    fig.savefig(OUT/f'perturbations_{args.group.lower()}_heatmap.pdf',bbox_inches='tight',facecolor='white'); plt.close(fig)
    manifest = {}
    token = args.group.lower()
    for path in sorted(OUT.glob('*')):
        if path.is_file() and token in path.name and not path.name.startswith('SHA256'):
            manifest[str(path.relative_to(ROOT))] = sha256(path)
    for path in (PROFILE, INITIALS, ROOT/'model.py', ROOT/'model_twobit34.py',
                 ROOT/'verify_twobit_causal.py', Path(__file__)):
        manifest[str(path.relative_to(ROOT))] = sha256(path)
    (OUT/f'SHA256_{token}.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2,sort_keys=True),encoding='utf-8')
    (OUT/f'SHA256_{token}.txt').write_text(''.join(f'{v}  {k}\n' for k,v in sorted(manifest.items())),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2)); print(f'wrote {out_csv}')


def main():
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='mode',required=True)
    s=sub.add_parser('scan'); s.add_argument('--group',choices=('S','RDF','AFFL','INT1','all'),default='all')
    s.add_argument('--shard',type=int,required=True); s.add_argument('--nshards',type=int,default=20)
    s.add_argument('--hours',type=float,default=200.0)
    m=sub.add_parser('merge'); m.add_argument('--group',choices=('S','RDF','AFFL','INT1','all'),default='all')
    m.add_argument('--nshards',type=int,default=20)
    args=ap.parse_args(); scan(args) if args.mode=='scan' else merge(args)


if __name__=='__main__': main()
