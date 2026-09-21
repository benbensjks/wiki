"""Independent verification of the 34-state perturbation and strict-tolerance runs.

Reads the raw shard CSVs and meta JSONs (not only the merged files), recomputes
completeness, cross-shard environment and source-hash uniformity, recovery
limits, pool sensitivity, flip phase behaviour, dangerous readout-only points,
and strict-tolerance convergence.  Re-verifies every SHA256 manifest.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
PERT = ROOT / 'twobit34_results' / 'perturbations'
STRICT = ROOT / 'twobit34_results' / 'strict_tolerance'

GROUPS = dict(S=10, RDF=10, AFFL=10, INT1=5)
PER_SHARD = dict(S=4, RDF=4, AFFL=4, INT1=4)
SOURCE_FILES = ('model.py', 'model_twobit34.py', 'verify_twobit_causal.py',
                'twobit34_results/selected_profile.json',
                'twobit34_results/initial_states_selected/four_initial_states.csv')
SCANNERS = dict(perturbations='scan_twobit34_perturbations.py',
                strict='scan_twobit34_strict_tolerance.py')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def as_bool(s):
    return s.map(lambda x: str(x).strip().lower() in ('true', '1'))


def load_pert_group(name, nshards):
    frames = []
    for s in range(nshards):
        p = PERT / f'perturb_{name.lower()}_shard{s:02d}of{nshards:02d}.csv'
        if not p.exists():
            raise SystemExit(f'missing {p}')
        frames.append(pd.read_csv(p))
    d = pd.concat(frames, ignore_index=True)
    for c in ('solver_success', 'finite', 'cold_passed', 'steady_passed', 'certified',
              'one_to_one', 'causal_order', 'directions_alternate'):
        d[c] = as_bool(d[c])
    return d


def meta_audit(files):
    metas = [json.loads(p.read_text(encoding='utf-8')) for p in files]
    def hashmap(m):
        return m.get('sha256') or m.get('source_sha256') or {}
    keys = sorted(set().union(*[set(hashmap(m)) for m in metas]))
    return dict(shards=len(metas),
                python=sorted({m['python'] for m in metas}) if 'python' in metas[0] else None,
                platform=sorted({m['platform'] for m in metas}) if 'platform' in metas[0] else None,
                versions=({k: sorted({m['versions'][k] for m in metas})
                           for k in metas[0]['versions']} if 'versions' in metas[0] else None),
                hours=sorted({m['hours'] for m in metas}),
                points_per_shard=sorted({m['points'] for m in metas}),
                uniform=all(len({hashmap(m)[k] for m in metas}) == 1 for k in keys),
                hashes={k: sorted({hashmap(m)[k] for m in metas}) for k in keys})


def verify_manifest(path: Path):
    manifest = json.loads(path.read_text(encoding='utf-8'))
    bad, missing = [], []
    for name, want in manifest.items():
        p = ROOT / name
        if not p.exists():
            missing.append(name)
        elif sha256(p) != want:
            bad.append(name)
    return dict(entries=len(manifest), mismatches=bad, missing=missing, ok=not bad and not missing)


def main():
    report = {}
    pert_rows, pert_meta = {}, {}
    for g, n in GROUPS.items():
        d = load_pert_group(g, n)
        pert_rows[g] = d
        metas = sorted(PERT.glob(f'meta_{g.lower()}_shard*of*.json'))
        pert_meta[g] = meta_audit(metas)
        report[f'{g}_integrity'] = dict(
            rows=len(d),
            per_shard=[len(pd.read_csv(PERT / f'perturb_{g.lower()}_shard{s:02d}of{n:02d}.csv'))
                       for s in range(n)],
            solver_failures=int((~d.solver_success).sum()),
            nonfinite=int((~d.finite).sum()),
            cold_passed=int(d.cold_passed.sum()), steady_passed=int(d.steady_passed.sum()),
            certified=int(d.certified.sum()),
            outcomes={str(k): int(v) for k, v in d.outcome.value_counts().items()})

    allp = pd.concat(pert_rows.values(), ignore_index=True)
    report['overall_outcomes'] = {str(k): int(v) for k, v in allp.outcome.value_counts().items()}
    report['by_group_outcomes'] = {g: {str(k): int(v) for k, v in d.outcome.value_counts().items()}
                                   for g, d in pert_rows.items()}
    report['by_target_outcomes'] = {str(t): {str(k): int(v) for k, v in sub.outcome.value_counts().items()}
                                    for t, sub in allp.groupby('target')}
    report['dangerous_readout_pass_causal_fail'] = int((allp.outcome == 'readout_pass_causal_fail').sum())

    # per target x initial-state outcome matrices
    report['matrix'] = {}
    for t, sub in allp.groupby('target'):
        piv = sub.pivot_table(index='initial_value', columns='action', values='outcome', aggfunc='first')
        report['matrix'][str(t)] = {str(i): {str(c): str(v) for c, v in row.items()}
                                    for i, row in piv.iterrows()}

    # S deltas: largest recoverable magnitude per target
    s = pert_rows['S']
    deltas = {}
    for t in ('S0', 'S1'):
        sub = s[s.target == t]
        for a in sorted(sub.action.unique()):
            ok = int((sub[sub.action == a].outcome == 'recovered_same_phase').sum())
            deltas[f'{t}:{a}'] = dict(recovered=ok, total=int((sub.action == a).sum()))
    report['S_actions'] = deltas
    mag = {}
    for t in ('S0', 'S1'):
        sub = s[s.target == t]
        for a in ('delta_-0.20', 'delta_-0.10', 'delta_+0.10', 'delta_+0.20'):
            frac = float((sub[sub.action == a].outcome == 'recovered_same_phase').mean())
            mag.setdefault(t, {})[a] = frac
    report['S_recovery_fraction_by_magnitude'] = mag

    # Clip confound: the S deltas are applied to a state that is already near 0 or 1
    # at some phases, so np.clip turns a "delta" into a no-op.  Recovery must be
    # scored by the EFFECTIVE change, not by the requested one.
    clipped = []
    for r in s.itertuples():
        eff = float(r.pool_after) - float(r.pool_before)
        clipped.append(dict(target=str(r.target), action=str(r.action), initial=int(r.initial_value),
                            before=float(r.pool_before), after=float(r.pool_after),
                            effective_delta=eff, clipped=abs(eff) < 1e-12,
                            outcome=str(r.outcome)))
    eff_buckets = {}
    for row in clipped:
        key = ('direction=up' if row['effective_delta'] > 1e-12 else
               'direction=down' if row['effective_delta'] < -1e-12 else 'clipped_noop')
        b = eff_buckets.setdefault(key, dict(total=0, recovered=0, shifted=0, lost=0,
                                             magnitudes=[]))
        b['total'] += 1
        b['magnitudes'].append(abs(row['effective_delta']))
        if row['outcome'] == 'recovered_same_phase':
            b['recovered'] += 1
        elif row['outcome'].startswith('stable_phase_shift'):
            b['shifted'] += 1
        else:
            b['lost'] += 1
    for k, b in eff_buckets.items():
        b['effective_magnitudes'] = sorted({round(m, 4) for m in b['magnitudes'] if m > 1e-12})
        b.pop('magnitudes')
    report['S_clip_confound'] = dict(
        rows=clipped, by_effective_direction=eff_buckets,
        note=('Requested delta and effective delta differ whenever the DNA state is already '
              'near 0 or 1 at that phase, because the perturbation is clipped into [0,1]. '
              'Downward deltas applied at the low phase are no-ops and carry no tolerance '
              'information; only the effective change should be quoted as a limit.'))

    # F0 vulnerability conditioned on the pre-perturbation pool size
    f0 = pert_rows['AFFL']
    f0 = f0[f0.target == 'F0_pool']
    f0_rows = [dict(factor=float(r.factor), initial=int(r.initial_value),
                    pool_before=float(r.pool_before), pool_after=float(r.pool_after),
                    outcome=str(r.outcome)) for r in f0.itertuples()]
    lo = [r for r in f0_rows if r['pool_before'] < 5]
    hi = [r for r in f0_rows if r['pool_before'] >= 5]
    report['F0_vulnerability'] = dict(
        rows=sorted(f0_rows, key=lambda r: (r['factor'], r['initial'])),
        low_pool_phase=dict(n=len(lo), lost=sum(1 for r in lo if r['outcome'] == 'lost_counting'),
                            factors=sorted({r['factor'] for r in lo if r['outcome'] == 'lost_counting'})),
        high_pool_phase=dict(n=len(hi), lost=sum(1 for r in hi if r['outcome'] == 'lost_counting')),
        note=('F0 halving only breaks counting at the phases where the F0 pool is already low '
              '(~2.7 a.u. vs ~9.0 a.u.); the same factor at the high-F0 phases is harmless.'))

    # flip behaviour: expected permanent phase shift
    flips = s[s.action == 'flip']
    report['flip'] = {
        'rows': int(len(flips)),
        'outcomes': {str(k): int(v) for k, v in flips.outcome.value_counts().items()},
        'phase_offsets': {str(k): int(v) for k, v in flips.phase_offset.value_counts(dropna=False).items()},
        'detail': [dict(target=str(r.target), initial=int(r.initial_value), outcome=str(r.outcome),
                        phase_offset=(None if pd.isna(r.phase_offset) else int(r.phase_offset)),
                        pool_before=float(r.pool_before), pool_after=float(r.pool_after))
                   for r in flips.itertuples()]}

    # factor pools: largest upward / smallest downward factor that still recovers
    factor_table = {}
    for g in ('RDF', 'AFFL', 'INT1'):
        d = pert_rows[g]
        for t, sub in d.groupby('target'):
            rows = []
            for f in sorted(sub.factor.dropna().unique()):
                ss = sub[sub.factor == f]
                rows.append(dict(factor=float(f),
                                 recovered=int((ss.outcome == 'recovered_same_phase').sum()),
                                 total=int(len(ss)),
                                 lost=int((ss.outcome == 'lost_counting').sum())))
            factor_table[str(t)] = rows
    report['factor_recovery'] = factor_table

    # pool sensitivity ranking
    rank = []
    for g, d in pert_rows.items():
        lost = int((d.outcome == 'lost_counting').sum())
        shifted = int(d.outcome.str.startswith('stable_phase_shift').sum())
        rank.append(dict(group=g, points=len(d), lost_counting=lost, phase_shifted=shifted,
                         recovered=int((d.outcome == 'recovered_same_phase').sum()),
                         non_recovered=lost + shifted))
    rank.sort(key=lambda r: (-r['non_recovered'], r['group']))
    report['pool_sensitivity_ranking'] = rank

    # strict tolerance
    sd = pd.read_csv(STRICT / 'strict_tolerance_all.csv')
    for c in ('solver_success', 'finite', 'cold_passed', 'steady_passed', 'certified',
              'one_to_one', 'causal_order', 'directions_alternate'):
        sd[c] = as_bool(sd[c])
    strict_metas = sorted(STRICT.glob('meta_shard*of*.json'))
    report['strict_integrity'] = dict(
        rows=len(sd), expected=27, per_shard=[len(pd.read_csv(STRICT / f'strict_shard{s:02d}of09.csv'))
                                              for s in range(9)],
        solver_failures=int((~sd.solver_success).sum()), nonfinite=int((~sd.finite).sum()),
        settings=sorted(sd.setting.unique().tolist()))
    report['strict_meta'] = meta_audit(strict_metas) if strict_metas else None
    numeric = ('minimum_state', 'maximum_state', 'minimum_commitment', 'bit0_min', 'bit0_max',
               'bit1_min', 'bit1_max', 'Jrev0_peak', 'g0_peak', 'Int1_source_peak',
               'bit1_setup_h', 'bit1_hold_h')
    pts = []
    for pid, g in sd.groupby('point_id'):
        g = g.set_index('setting')
        ref = g.loc['ultra']
        cmp_rows = []
        for name in ('baseline', 'tight'):
            r = g.loc[name]
            diffs = {k: abs(float(r[k]) - float(ref[k])) for k in numeric}
            cmp_rows.append(dict(setting=name,
                                 same_certified=bool(r.certified == ref.certified),
                                 same_reason=str(r.failure_reason) == str(ref.failure_reason),
                                 same_cold_sequence=str(r.cold_sequence) == str(ref.cold_sequence),
                                 same_steady_sequence=str(r.steady_sequence) == str(ref.steady_sequence),
                                 max_numeric_difference=max(diffs.values()),
                                 worst_metric=max(diffs, key=diffs.get)))
        pts.append(dict(point_id=int(pid), name=str(ref.point_name), kind=str(ref.kind),
                        coords=[float(ref.uM_per_au), float(ref.carry_maturation_half_life_min),
                                float(ref.carry_mrna_half_life_min)],
                        ultra_certified=bool(ref.certified),
                        ultra_reason='' if pd.isna(ref.failure_reason) else str(ref.failure_reason),
                        cold_sequence=str(ref.cold_sequence), steady_sequence=str(ref.steady_sequence),
                        baseline_certified=bool(g.loc['baseline'].certified),
                        comparisons=cmp_rows,
                        verdicts_agree=all(c['same_certified'] and c['same_reason'] and
                                           c['same_cold_sequence'] and c['same_steady_sequence']
                                           for c in cmp_rows)))
    report['strict_points'] = pts
    report['strict_all_verdicts_agree'] = all(p['verdicts_agree'] for p in pts)
    report['strict_worst_numeric_difference'] = max(c['max_numeric_difference']
                                                    for p in pts for c in p['comparisons'])
    report['strict_criterion_changes'] = [p for p in pts if not p['verdicts_agree']]

    # manifests
    manifests = {}
    for g in GROUPS:
        p = PERT / f'SHA256_{g.lower()}.json'
        if p.exists():
            manifests[f'perturbations_{g}'] = verify_manifest(p)
    p = STRICT / 'SHA256SUMS.json'
    if p.exists():
        manifests['strict_tolerance'] = verify_manifest(p)
    report['manifests'] = manifests
    report['manifests_all_ok'] = all(m['ok'] for m in manifests.values())

    # environment / source hashes across every meta
    report['environment_uniform_across_all'] = all(m['uniform'] for m in pert_meta.values()) and \
        (report['strict_meta']['uniform'] if report['strict_meta'] else True)
    report['source_hashes_now'] = {f: sha256(ROOT / f) for f in SOURCE_FILES}
    report['scanner_hashes_now'] = {k: sha256(ROOT / v) for k, v in SCANNERS.items()}

    out = ROOT / 'twobit34_results' / 'experiment_verification.json'
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')

    print('== perturbation integrity ==')
    for g in GROUPS:
        r = report[f'{g}_integrity']
        print(f"  {g:5s} rows={r['rows']} solver_fail={r['solver_failures']} nonfinite={r['nonfinite']} "
              f"certified={r['certified']} outcomes={r['outcomes']}")
    print('== pool sensitivity ranking (most sensitive first) ==')
    for r in report['pool_sensitivity_ranking']:
        print(f"  {r['group']:5s} non_recovered={r['non_recovered']}/{r['points']} "
              f"(lost={r['lost_counting']}, shifted={r['phase_shifted']}, recovered={r['recovered']})")
    print('== dangerous readout-pass / causal-fail:',
          report['dangerous_readout_pass_causal_fail'])
    print('== S flip ==')
    print('  outcomes:', report['flip']['outcomes'], 'offsets:', report['flip']['phase_offsets'])
    print('== S clip confound (by EFFECTIVE change) ==')
    for k, b in report['S_clip_confound']['by_effective_direction'].items():
        print(f"  {k:16s} n={b['total']} recovered={b['recovered']} shifted={b['shifted']} "
              f"lost={b['lost']} |eff|={b['effective_magnitudes']}")
    print('== F0 vulnerability ==')
    f0r = report['F0_vulnerability']
    print(f"  low-pool phase: {f0r['low_pool_phase']}")
    print(f"  high-pool phase: {f0r['high_pool_phase']}")
    print('== strict ==')
    print(f"  rows={report['strict_integrity']['rows']} solver_fail={report['strict_integrity']['solver_failures']} "
          f"verdicts_agree={report['strict_all_verdicts_agree']} "
          f"worst_numeric_diff={report['strict_worst_numeric_difference']:.3g}")
    print('== manifests ==')
    for k, v in report['manifests'].items():
        print(f"  {k:22s} entries={v['entries']} mismatches={len(v['mismatches'])} missing={len(v['missing'])} ok={v['ok']}")
    print(f"wrote {out}")


if __name__ == '__main__':
    main()
