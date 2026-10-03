"""Round 1b: fill the one cell round 1 left open.

Round 1 ran the clock probes on the AUTOREGULATION ON arm only, so the working
point it found is "A1 negative autoregulation ON + clock_K = 0.20". The source
file 前馈三级级联.py actually runs with A1 autoregulation OFF (its da1 line is
commented out) and clock_K = 0.40. The cell "autoregulation OFF + clock_K fixed"
was therefore never measured, and that is exactly the cell that answers:

    "does the source file's own tail configuration work once the clock gate is
     matched to the real Int0 scale?"

This script adds ONLY that cell (plus the K = 0.07 counterpart for symmetry).
Everything else is reused unchanged from diagnose_round1.py; no existing file is
modified and no existing result is overwritten.
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT.parent.parent))

import diagnose_round1 as D                                       # noqa: E402
import verify_hzh as V                                            # noqa: E402

# The missing cells. Round 1's own entries are re-listed for reference only.
NEW = [
    ('off', 0.20),   # <-- the cell the original configuration needs
    ('off', 0.07),
]
REFERENCE_FROM_ROUND1 = [
    ('on', 0.20), ('on', 0.07), ('on', 0.40), ('off', 0.40),
]


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = ROOT / 'results' / ('round1b_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir(parents=True, exist_ok=False)
    D.dump(out / 'status.json', dict(status='running'))

    print('build shared %g h Han input' % D.HOURS, flush=True)
    han = D.HanInput(D.HOURS)
    prefix = D.HZHModel('han', han)

    # prefix guard, so a bad prefix cannot silently invalidate the new cells
    tb, yb = D.integrate(prefix)
    bv, _ = V.analyse(tb, yb, prefix.z, prefix.tail)
    guard = dict(certified=bv['certified_v1'],
                 reads=bv['steady']['reads'] == 19,
                 sequence=bv['steady']['sequence'] == '1234567012345670123',
                 margin=abs(bv['global_min_timing_margin_h']
                            - 1.1784609018563117) < 1e-6)
    print('HZH prefix guard:', json.dumps(guard), flush=True)
    if not all(guard.values()):
        D.dump(out / 'status.json', dict(status='FAILED', guard=guard))
        raise RuntimeError(guard)

    rows = []
    for arm, clock_K in NEW:
        label = f'clk{int(round(clock_K*100)):03d}|{arm}'
        model = D.HZZModel(prefix, arm == 'on')
        model.clock_K = clock_K
        t, y = D.integrate(model)
        sig = model.signals(y)
        sig['_p'] = model.p
        if y.min() < -1e-7:
            raise RuntimeError(f'{label}: nonphysical state')
        res = D.evaluate(t, sig)
        row = dict(label=label, arm=arm, clock_K=clock_K, overrides={})
        row.update(D.summarise(t, y, sig, res, clock_K, {}))
        row['prefix_gap'] = float(np.max(np.abs(y[:17] - yb[:17])))
        row['per_bit'] = res['per_bit']
        rows.append(row)
        np.savez_compressed(out / f'traj_clk{int(round(clock_K*100)):03d}_{arm}.npz',
                            time_h=t, states=y, state_names=np.array(D.NAMES))
        D.dump(out / f'verdict_clk{int(round(clock_K*100)):03d}_{arm}.json', res)
        m = row['mechanism']
        print(f"[{label:12s}] cnt={int(row['counting_passed'])} "
              f"ev={int(row['event_causality_passed'])} "
              f"S2flips={row['s2_flips']} seq={row['steady']:20s} "
              f"g1pk={row['g1_peak']:.4f} "
              f"dur={[round(x,3) for x in (row['gate_duration_range'] or [])]} "
              f"dose={[round(x,4) for x in (row['gate_dose_range'] or [])]} "
              f"meetBoth={row['gate_segments_meeting_both']} "
              f"RDF2pk={m['RDF2_peak']:.4f} vr/vf={m['reverse_over_forward']:.2f}",
              flush=True)

    # re-read round 1 for the side-by-side table
    r1 = sorted((ROOT / 'results').glob('2026*'))
    r1_rows = []
    if r1:
        with (r1[-1] / 'round1_all.csv').open(encoding='utf-8') as fh:
            r1_rows = list(csv.DictReader(fh))

    def f(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return None
    table = []
    for r in r1_rows:
        table.append(dict(source='round1', label=r['label'], arm=r['arm'],
                          clock_K=f(r['clock_K']), counting=r['counting_passed'],
                          events=r['event_causality_passed'],
                          s2_flips=f(r['s2_flips']), sequence=r['steady'],
                          gate_dose_max=f(r['gate_dose_max']),
                          gate_duration_max=f(r['gate_dur_max']),
                          meets_both=f(r['gates_meet_both']),
                          I2_peak=f(r['I2_peak']), RDF2_peak=f(r['RDF2_peak']),
                          reverse_hill=f(r['reverse_hill']),
                          vr_over_vf=f(r['vr_over_vf'])))
    for r in rows:
        m = r['mechanism']
        dur = r['gate_duration_range'] or [None, None]
        dose = r['gate_dose_range'] or [None, None]
        table.append(dict(source='round1b', label=r['label'], arm=r['arm'],
                          clock_K=r['clock_K'],
                          counting=str(r['counting_passed']),
                          events=str(r['event_causality_passed']),
                          s2_flips=r['s2_flips'], sequence=r['steady'],
                          gate_dose_max=dose[1], gate_duration_max=dur[1],
                          meets_both=r['gate_segments_meeting_both'],
                          I2_peak=m['I2_peak'], RDF2_peak=m['RDF2_peak'],
                          reverse_hill=m['reverse_hill_value'],
                          vr_over_vf=m['reverse_over_forward']))

    D.dump(out / 'round1b_summary.json', dict(
        hours=D.HOURS, solver=D.SOLVER, prefix_guard=guard,
        new_cells=[dict(arm=a, clock_K=k) for a, k in NEW],
        reference_from_round1=[dict(arm=a, clock_K=k) for a, k in REFERENCE_FROM_ROUND1],
        question=('does the source file 前馈三级级联.py own tail configuration '
                  '(A1 autoregulation OFF, clock_K 0.40) work once the clock gate '
                  'is matched to the real HBY Int0 scale?'),
        rows=rows, side_by_side=table))

    cols = ['source', 'label', 'arm', 'clock_K', 'counting', 'events',
            's2_flips', 'sequence', 'gate_dose_max', 'gate_duration_max',
            'meets_both', 'I2_peak', 'RDF2_peak', 'reverse_hill', 'vr_over_vf']
    with (out / 'round1b_side_by_side.csv').open('w', newline='',
                                                 encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in table:
            w.writerow({k: r.get(k) for k in cols})

    D.dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows)))
    D.dump(out / 'SHA256SUMS.json',
           {str(p.relative_to(out)): D.sha256(p)
            for p in sorted(out.rglob('*')) if p.is_file()})
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
