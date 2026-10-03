"""Stage 2: bracket the OAT boundaries found in stage 1, and test whether
sharpening the carry gate can compensate a mis-placed concentration scale.

Imports the stage-1 harness (scan_oat.py) without modifying it, so the stage-1
provenance hashes stay valid. Adds two things stage 1 did not record:
  * signal_ranges (clock gate, g1, Int0, u2 target) from the verdict, which
    shows directly whether the second-stage AND gate is saturated;
  * a paired 2-D grid n_A1_gate x receiver_uM_per_au, which is the only design
    that can answer "can a steeper gate rescue a wrong a.u. scale?" -- a
    single-parameter scan cannot, because it never varies two levers together.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import scan_oat as S                                              # noqa: E402
from hybrid_model import HbyConfig                                # noqa: E402

FROZEN = S.FROZEN
ZENG = S.ZENG


def evaluate(prototype, label, cfg, zeng_patch=None, ref_traj=None, extra=None):
    """Stage-1 evaluate plus signal_ranges."""
    row = dict(label=label)
    if extra:
        row.update(extra)
    tic = time.perf_counter()
    try:
        m = S.variant(prototype, cfg, zeng_patch)
        row['work_point'] = S.work_point_readback(m)
        t, y = S.integrate(m)
        from verify_hzh import analyse
        v, sig = analyse(t, y, m.z, m.tail)
        bm = v['bit_margins']
        row.update(
            certified=bool(v['certified_v1']),
            steady_reads=int(v['steady']['reads']),
            increments_mod8=bool(v['steady']['increments_mod8']),
            sequence=str(v['steady']['sequence']),
            minimum_commitment=(None if v['steady']['minimum_commitment'] is None
                                else float(v['steady']['minimum_commitment'])),
            boundary_clips=int(v['steady']['boundary_clips']),
            margin_h=(None if v['global_min_timing_margin_h'] is None
                      else float(v['global_min_timing_margin_h'])),
            clock_peak_count=int(v['clock_peak_count']),
            margin_S1_h=S.min_of(bm.get('S1')),
            margin_S2_h=S.min_of(bm.get('S2')),
            bit_margins={b: {k: (None if x is None else float(x))
                             for k, x in mm.items()} for b, mm in bm.items()},
            signal_ranges={k: [float(x) for x in vv]
                           for k, vv in v['signal_ranges'].items()},
            stage0=dict(reverse=int(v['events']['bit0_to_bit1']['reverse_events']),
                        gate=int(v['events']['bit0_to_bit1']['gate_events']),
                        flip=int(v['events']['bit0_to_bit1']['flip_events']),
                        one_to_one=bool(v['events']['bit0_to_bit1']['one_to_one']),
                        order=bool(v['events']['bit0_to_bit1']['causal_order']),
                        alt=bool(v['events']['bit0_to_bit1']['alternating_directions'])),
            stage1=dict(reverse=int(v['events']['bit1_to_bit2']['reverse_events']),
                        gate=int(v['events']['bit1_to_bit2']['gate_events']),
                        flip=int(v['events']['bit1_to_bit2']['flip_events']),
                        one_to_one=bool(v['events']['bit1_to_bit2']['one_to_one']),
                        order=bool(v['events']['bit1_to_bit2']['causal_order']),
                        alt=bool(v['events']['bit1_to_bit2']['alternating_directions'])),
            error=None)
        if ref_traj is not None:
            row['max_abs_gap_vs_baseline'] = float(np.max(np.abs(y - ref_traj)))
        else:
            row['max_abs_gap_vs_baseline'] = 0.0
        row['_traj'] = y
    except Exception as exc:
        row.update(certified=False, error=repr(exc), steady_reads=0, sequence='',
                   increments_mod8=False, minimum_commitment=None,
                   boundary_clips=0, margin_h=None, clock_peak_count=0,
                   margin_S1_h=None, margin_S2_h=None, bit_margins=None,
                   signal_ranges=None, stage0=None, stage1=None,
                   max_abs_gap_vs_baseline=None)
    row['runtime_s'] = time.perf_counter() - tic
    return row


def cfg_with(**kw):
    return HbyConfig(**{**asdict(FROZEN), **kw})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', default='both', choices=['1', '2', 'both'])
    args = ap.parse_args()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    out = ROOT / 'results' / ('refine_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True, exist_ok=False)
    S.dump(out / 'status.json', dict(status='RUNNING'))

    print('building Han upstream once', flush=True)
    han = S.HanInput(S.HOURS)
    assert han.p.peak_load_fraction == 0.0
    prototype = S.HZHModel('han', han)

    rows = []
    base = evaluate(prototype, 'BASELINE', FROZEN)
    checks = dict(certified=base['certified'] is True,
                  reads=base['steady_reads'] == S.REF['steady_reads'],
                  sequence=base['sequence'] == S.REF['sequence'],
                  margin=abs(base['margin_h'] - S.REF['margin_h']) < 1e-6,
                  n_A1_gate_readback=base['work_point']['n_A1_gate_effective'] == 6.0)
    print('baseline guard:', json.dumps(checks), flush=True)
    if not all(checks.values()):
        S.dump(out / 'status.json', dict(status='FAILED', error='baseline guard',
                                         checks=checks))
        raise RuntimeError(f'Baseline guard failed: {checks}')
    rows.append(base)
    ref = base['_traj']

    # ---------------- Part 1: bracket the stage-1 boundaries --------------
    if args.part in ('1', 'both'):
        spec = []
        for v in (2.0, 3.0):
            spec.append((f'n_A1_gate={v:g}', 'B', 'config', 'n_A1_gate', v, None))
        for f in (0.9, 1.1, 1.5, 1.75):
            spec.append((f'receiver_uM_per_au=x{f:g}', 'A', 'config',
                         'receiver_uM_per_au', FROZEN.receiver_uM_per_au * f, None))
        for f in (1.5, 1.75):
            spec.append((f'K_A[1]=x{f:g}', 'B', 'zeng', ('K_A', 1),
                         ZENG['K_A'][1] * f, None))
        for f in (1.5, 1.75):
            spec.append((f'K_F[1]=x{f:g}', 'B', 'zeng', ('K_F', 1),
                         ZENG['K_F'][1] * f, None))
        for f in (0.25, 4.0):
            spec.append((f'clock_K_au=x{f:g}', 'A', 'config', 'clock_K_au',
                         FROZEN.clock_K_au * f, None))
        print(f'Part 1 rows: {len(spec)}', flush=True)
        for i, (label, tier, kind, target, value, _) in enumerate(spec, 1):
            if kind == 'config':
                cfg, patch = cfg_with(**{target: float(value)}), None
            else:
                cfg, patch = FROZEN, {target: float(value)}
            r = evaluate(prototype, label, cfg, patch, ref_traj=ref,
                         extra=dict(tier=tier, kind=kind, target=str(target),
                                    value=float(value)))
            rows.append(r)
            sr = r.get('signal_ranges') or {}
            cl = sr.get('clock') or [None, None]
            print(f'[P1 {i}/{len(spec)}] {label:28s} cert={str(r["certified"]):5s} '
                  f'reads={r["steady_reads"]:3d} S2m={r.get("margin_S2_h")} '
                  f'clock=[{cl[0]:.3g},{cl[1]:.3g}] gap={r["max_abs_gap_vs_baseline"]:.2e} '
                  f'{r["runtime_s"]:.0f}s', flush=True)

    # ---------------- Part 2: paired 2-D grid ----------------------------
    if args.part in ('2', 'both'):
        gates = [4.0, 6.0, 8.0]
        scales = [0.8, 1.0, 1.25, 2.0]
        spec = [(g, f) for g in gates for f in scales
                if not (g == 6.0 and f == 1.0)]
        print(f'Part 2 rows: {len(spec)} (paired grid, {len(gates)}x{len(scales)} '
              f'minus the frozen cell)', flush=True)
        for i, (g, f) in enumerate(spec, 1):
            label = f'GRID n_A1_gate={g:g} uM_per_au=x{f:g}'
            cfg = cfg_with(n_A1_gate=g,
                           receiver_uM_per_au=FROZEN.receiver_uM_per_au * f)
            r = evaluate(prototype, label, cfg, None, ref_traj=ref,
                         extra=dict(tier='GRID', kind='grid',
                                    n_A1_gate=g, uM_fold=f,
                                    target='n_A1_gate x receiver_uM_per_au',
                                    value=f))
            rows.append(r)
            print(f'[P2 {i}/{len(spec)}] n={g:g} uM=x{f:g}  '
                  f'cert={str(r["certified"]):5s} reads={r["steady_reads"]:3d} '
                  f'S2m={r.get("margin_S2_h")} seq={r["sequence"][:24]} '
                  f'{r["runtime_s"]:.0f}s', flush=True)

    # ---------------- persist -------------------------------------------
    slim = [{k: v for k, v in r.items() if not k.startswith('_')} for r in rows]
    S.dump(out / 'refine_all.json', dict(
        hours=S.HOURS, sample_min=S.SAMPLE_MIN, max_step_min=S.MAX_STEP_MIN,
        rtol=S.RTOL, atol=S.ATOL, frozen_config=asdict(FROZEN),
        reference=S.REF, baseline_guard=checks, rows=slim))

    cols = ['label', 'tier', 'kind', 'target', 'value', 'n_A1_gate', 'uM_fold',
            'certified', 'steady_reads', 'increments_mod8', 'sequence',
            'minimum_commitment', 'boundary_clips', 'margin_h',
            'margin_S1_h', 'margin_S2_h', 'clock_peak_count',
            'clock_min', 'clock_max', 'g1_min', 'g1_max',
            'int0_min', 'int0_max', 'g0_min', 'g0_max',
            's0_rev', 's0_gate', 's0_flip', 's1_rev', 's1_gate', 's1_flip',
            'max_abs_gap_vs_baseline', 'runtime_s', 'error']
    with (out / 'refine_all.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            s0 = r.get('stage0') or {}
            s1 = r.get('stage1') or {}
            sr = r.get('signal_ranges') or {}

            def rg(name, i):
                v = sr.get(name)
                return None if not v else v[i]
            w.writerow([r['label'], r.get('tier', ''), r.get('kind', ''),
                        r.get('target', ''), r.get('value', ''), r.get('n_A1_gate'),
                        r.get('uM_fold'), r['certified'], r['steady_reads'],
                        r['increments_mod8'], r['sequence'], r['minimum_commitment'],
                        r['boundary_clips'], r['margin_h'], r.get('margin_S1_h'),
                        r.get('margin_S2_h'), r['clock_peak_count'],
                        rg('clock', 0), rg('clock', 1), rg('g1', 0), rg('g1', 1),
                        rg('int0', 0), rg('int0', 1), rg('g0', 0), rg('g0', 1),
                        s0.get('reverse'), s0.get('gate'), s0.get('flip'),
                        s1.get('reverse'), s1.get('gate'), s1.get('flip'),
                        r.get('max_abs_gap_vs_baseline'),
                        round(r['runtime_s'], 3), r['error']])

    S.dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows)))
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
