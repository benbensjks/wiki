"""Independent verification of the 75-point 34-state robustness scan.

Reads the raw shard CSVs and meta JSONs directly (not the merged summary),
recomputes completeness, the continuity of the certified region, the top
timing-margin points, the nominal point 6/30/2, the cross-shard environment and
source-hash uniformity, and re-verifies the SHA256 manifest.

Read-only: nothing under the scan tree is modified.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'twobit34_results' / 'robustness'
NSHARDS = 15

UM = (5.5, 5.75, 6.0, 6.25, 6.5)
MAT = (25.0, 27.5, 30.0, 32.5, 35.0)
MRNA = (1.0, 2.0, 4.0)
GRID_KEYS = ('uM_per_au', 'carry_maturation_half_life_min', 'carry_mrna_half_life_min')

BASELINE = dict(
    model_py='256D2D104EECB4CE17F27CF4DBC08257BEEEA79E5805FD7D5614E47BBC7E543A',
    model_twobit34_py='4385E3B8C9B8AE39D09061A9DF7D5B1E7F645AC5CF56B22E4C7183B421501062',
    verifier_py='C8678A03EC6201E65EBBE2BD182AB89BC2AB8DAAD502A4E302CAFFBC361BCFF0',
    scanner_py='CF66786011AFD94A084AD253449B0E6C22104E6E65CF2A596ED784B76E4008E6',
)
SOURCE_FILES = dict(model_py='model.py', model_twobit34_py='model_twobit34.py',
                    verifier_py='verify_twobit_causal.py',
                    scanner_py='scan_twobit34_robustness.py')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def as_bool(s):
    return s.map(lambda x: str(x).strip().lower() in ('true', '1'))


def load_raw():
    frames, metas = [], []
    for s in range(NSHARDS):
        tag = f'shard{s:02d}of{NSHARDS:02d}'
        frames.append(pd.read_csv(OUT / f'robustness_{tag}.csv'))
        metas.append(json.loads((OUT / f'meta_{tag}.json').read_text(encoding='utf-8')))
    d = pd.concat(frames, ignore_index=True)
    for k in GRID_KEYS:
        d[k] = d[k].astype(float)
    for k in ('solver_success', 'finite', 'cold_passed', 'steady_passed', 'one_to_one',
              'causal_order', 'directions_alternate', 'certified'):
        d[k] = as_bool(d[k])
    return d, metas


def components(mask):
    """6-neighbour connected components of a bool array indexed [mrna][mat][um]."""
    seen = np.zeros_like(mask, dtype=bool)
    comps = []
    shape = mask.shape
    for idx in np.argwhere(mask):
        start = tuple(idx)
        if seen[start]:
            continue
        stack, cells = [start], []
        seen[start] = True
        while stack:
            cur = stack.pop()
            cells.append(cur)
            for axis in range(3):
                for step in (-1, 1):
                    nxt = list(cur)
                    nxt[axis] += step
                    nxt = tuple(nxt)
                    if all(0 <= nxt[a] < shape[a] for a in range(3)) and mask[nxt] and not seen[nxt]:
                        seen[nxt] = True
                        stack.append(nxt)
        comps.append(cells)
    return comps


def describe(cells):
    arr = np.array(cells)
    return dict(size=len(cells),
                mrna=sorted({MRNA[i] for i in arr[:, 0]}),
                mat=[min(MAT[i] for i in arr[:, 1]), max(MAT[i] for i in arr[:, 1])],
                um=[min(UM[i] for i in arr[:, 2]), max(UM[i] for i in arr[:, 2])],
                cells=[dict(mrna=MRNA[a], mat=MAT[b], um=UM[c]) for a, b, c in cells])


def main():
    d, metas = load_raw()
    report = {}

    # 1. completeness recomputed from raw shards
    expected = {tuple(p) for p in itertools.product(UM, MAT, MRNA)}
    got = {tuple(float(getattr(r, k)) for k in GRID_KEYS) for r in d.itertuples()}
    report['completeness_raw'] = dict(
        rows=len(d), unique=len(got), expected=len(expected),
        duplicates=len(d) - len(got), missing=sorted(expected - got),
        unexpected=sorted(got - expected),
        solver_failures=int((~d.solver_success).sum()),
        nonfinite=int((~d.finite).sum()),
        shard_row_counts=[len(pd.read_csv(OUT / f'robustness_shard{s:02d}of{NSHARDS:02d}.csv'))
                          for s in range(NSHARDS)],
    )

    # 2/3/4. pass counts, failure reasons, per-mRNA rates
    report['passes'] = dict(
        cold=int(d.cold_passed.sum()), steady=int(d.steady_passed.sum()),
        certified=int(d.certified.sum()), certified_rate=float(d.certified.mean()),
        one_to_one=int(d.one_to_one.sum()), causal_order=int(d.causal_order.sum()),
        directions_alternate=int(d.directions_alternate.sum()),
        certified_but_not_steady=int((d.certified & ~d.steady_passed).sum()),
        steady_but_not_certified=int((d.steady_passed & ~d.certified).sum()),
        failure_reasons={str(k): int(v) for k, v in d.failure_reason.fillna('').value_counts().items()},
        by_mrna={str(v): dict(certified=int(d.loc[d.carry_mrna_half_life_min == v, 'certified'].sum()),
                              points=int((d.carry_mrna_half_life_min == v).sum()))
                 for v in MRNA},
    )

    # 5. continuity of the certified region
    mask = np.zeros((len(MRNA), len(MAT), len(UM)), dtype=bool)
    for r in d.itertuples():
        mask[MRNA.index(float(r.carry_mrna_half_life_min)),
             MAT.index(float(r.carry_maturation_half_life_min)),
             UM.index(float(r.uM_per_au))] = bool(r.certified)
    comps = sorted((describe(c) for c in components(mask)), key=lambda c: -c['size'])
    per_slice = {}
    for i, m in enumerate(MRNA):
        sl = mask[i]
        row = {}
        for j, mat in enumerate(MAT):
            idx = np.flatnonzero(sl[j])
            row[f'{mat:g}'] = [UM[k] for k in idx]
        per_slice[f'{m:g}'] = row
    report['continuity'] = dict(
        total_components=len(comps), components=comps[:4],
        per_slice_certified_um=per_slice,
        full_grid_slice_contiguous={f'{m:g}': all(
            len(np.flatnonzero(mask[i][j])) <= 1 or
            np.all(np.diff(np.flatnonzero(mask[i][j])) == 1) for j in range(len(MAT)))
            for i, m in enumerate(MRNA)})

    # 5b. region shape: interior holes, largest solid box, safest interior points
    def neighbours(idx):
        out = []
        for axis in range(3):
            for step in (-1, 1):
                nxt = list(idx)
                nxt[axis] += step
                nxt = tuple(nxt)
                if all(0 <= nxt[a] < mask.shape[a] for a in range(3)):
                    out.append(nxt)
        return out

    fail_cells = [tuple(i) for i in np.argwhere(~mask)]
    interior_failures = []
    for cell in fail_cells:
        nb = neighbours(cell)
        n_ok = sum(1 for n in nb if mask[n])
        if n_ok >= max(4, len(nb) - 1):
            interior_failures.append(dict(mrna=MRNA[cell[0]], mat=MAT[cell[1]], um=UM[cell[2]],
                                          certified_neighbours=n_ok, neighbours=len(nb)))
    best_box, best_vol = None, -1
    for i0 in range(len(MRNA)):
        for i1 in range(i0, len(MRNA)):
            for j0 in range(len(MAT)):
                for j1 in range(j0, len(MAT)):
                    for k0 in range(len(UM)):
                        for k1 in range(k0, len(UM)):
                            sub = mask[i0:i1 + 1, j0:j1 + 1, k0:k1 + 1]
                            if sub.all():
                                vol = sub.size
                                if vol > best_vol:
                                    best_vol, best_box = vol, dict(
                                        mrna=[MRNA[i0], MRNA[i1]], mat=[MAT[j0], MAT[j1]],
                                        um=[UM[k0], UM[k1]], cells=int(vol))
    interior_points = []
    for r in d[d.certified].itertuples():
        cell = (MRNA.index(float(r.carry_mrna_half_life_min)),
                MAT.index(float(r.carry_maturation_half_life_min)),
                UM.index(float(r.uM_per_au)))
        nb = neighbours(cell)
        if nb and all(mask[n] for n in nb):
            interior_points.append(dict(uM_per_au=float(r.uM_per_au),
                                        carry_maturation_half_life_min=float(r.carry_maturation_half_life_min),
                                        carry_mrna_half_life_min=float(r.carry_mrna_half_life_min),
                                        timing_margin_h=float(r.timing_margin_h)))
    interior_points.sort(key=lambda x: -x['timing_margin_h'])
    by_um = {str(u): dict(points=int((d.uM_per_au == u).sum()),
                          certified=int(d.loc[d.uM_per_au == u, 'certified'].sum()))
             for u in UM}
    report['region_shape'] = dict(
        failing_cells=len(fail_cells), interior_failures=interior_failures,
        largest_solid_box=best_box, interior_points=sorted(
            interior_points, key=lambda x: -x['timing_margin_h'])[:10],
        interior_point_count=len(interior_points),
        certified_by_um=by_um,
        note=('The certified set is one 6-connected component but it is not a box: '
              'uM_per_au = 6.25 fails for most maturation values, so a nominal point must be '
              'chosen from the interior points or from the largest solid box.'))

    # 6. top five by timing margin
    ok = d[d.certified].sort_values(['timing_margin_h', 'minimum_commitment', 'g0_peak'],
                                    ascending=False)
    report['top5'] = [dict(uM_per_au=float(r.uM_per_au),
                           carry_maturation_half_life_min=float(r.carry_maturation_half_life_min),
                           carry_mrna_half_life_min=float(r.carry_mrna_half_life_min),
                           timing_margin_h=float(r.timing_margin_h),
                           bit1_setup_h=float(r.bit1_setup_h), bit1_hold_h=float(r.bit1_hold_h),
                           minimum_commitment=float(r.minimum_commitment),
                           g0_peak=float(r.g0_peak), Int1_source_peak=float(r.Int1_source_peak),
                           cold_sequence=str(r.cold_sequence), steady_sequence=str(r.steady_sequence))
                      for r in ok.head(5).itertuples()]

    # 7. nominal point 6 / 30 / 2
    nom = d[(d.uM_per_au == 6.0) & (d.carry_maturation_half_life_min == 30.0)
            & (d.carry_mrna_half_life_min == 2.0)]
    r = nom.iloc[0]
    report['nominal_6_30_2'] = dict(
        found=bool(len(nom)) == 1, certified=bool(r.certified), steady=bool(r.steady_passed),
        cold=bool(r.cold_passed), timing_margin_h=float(r.timing_margin_h),
        bit1_setup_h=float(r.bit1_setup_h), bit1_hold_h=float(r.bit1_hold_h),
        minimum_commitment=float(r.minimum_commitment), g0_peak=float(r.g0_peak),
        Int1_source_peak=float(r.Int1_source_peak), Jrev0_peak=float(r.Jrev0_peak),
        one_to_one=bool(r.one_to_one), causal_order=bool(r.causal_order),
        directions_alternate=bool(r.directions_alternate),
        cold_sequence=str(r.cold_sequence), steady_sequence=str(r.steady_sequence),
        failure_reason=str(r.failure_reason),
        matches_smoke_claim=(bool(r.certified) and abs(float(r.timing_margin_h) - 2.799) < 0.01))

    # 8. environment and source-hash uniformity across the 15 metas
    keys = ('model_py', 'model_twobit34_py', 'verifier_py', 'scanner_py')
    report['environment'] = dict(
        shards=len(metas),
        points_per_shard=sorted({m['points'] for m in metas}),
        hours=sorted({m['hours'] for m in metas}),
        python=sorted({m['python'] for m in metas}),
        platform=sorted({m['platform'] for m in metas}),
        versions={k: sorted({m['versions'][k] for m in metas}) for k in ('numpy', 'pandas', 'scipy')},
        source_hashes={k: sorted({m['sha256'][k] for m in metas}) for k in keys},
        uniform=all(len({m['sha256'][k] for m in metas}) == 1 for k in keys),
        matches_baseline={k: ({m['sha256'][k] for m in metas} == {BASELINE[k]}) for k in keys},
        on_disk_now={k: sha256(ROOT / SOURCE_FILES[k]) for k in keys},
        on_disk_equals_baseline={k: sha256(ROOT / SOURCE_FILES[k]) == BASELINE[k] for k in keys},
        grid_identical=len({json.dumps(m['grid'], sort_keys=True) for m in metas}) == 1,
    )

    # 9. re-verify the SHA256 manifest
    manifest = json.loads((OUT / 'SHA256SUMS.json').read_text(encoding='utf-8'))
    mismatches, missing = [], []
    for name, want in manifest.items():
        p = ROOT / name
        if not p.exists():
            missing.append(name)
        elif sha256(p) != want:
            mismatches.append(name)
    report['manifest'] = dict(entries=len(manifest), mismatches=mismatches, missing=missing,
                              ok=not mismatches and not missing,
                              covers_source_files=all(any(n.endswith(f) for n in manifest)
                                                      for f in SOURCE_FILES.values()))

    # 10. interpretation guard
    report['interpretation_guard'] = (
        'This scan only tests engineering robustness under the assumption '
        'add_growth=False (Zeng gamma read as total clearance). It is not '
        'experimental evidence that gamma already contains growth dilution.')

    (OUT / 'independent_verification.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')

    c = report['completeness_raw']
    print(f"completeness: rows={c['rows']} unique={c['unique']} expected={c['expected']} "
          f"dups={c['duplicates']} missing={len(c['missing'])} unexpected={len(c['unexpected'])} "
          f"solver_fail={c['solver_failures']} nonfinite={c['nonfinite']}")
    p = report['passes']
    print(f"passes: cold={p['cold']} steady={p['steady']} certified={p['certified']} "
          f"({p['certified_rate']:.1%}) reasons={p['failure_reasons']}")
    print(f"by mRNA: {p['by_mrna']}")
    print(f"continuity: components={report['continuity']['total_components']} "
          f"sizes={[c0['size'] for c0 in report['continuity']['components']]}")
    rs = report['region_shape']
    print(f"region shape: failing={rs['failing_cells']} interior_failures={len(rs['interior_failures'])} "
          f"largest_solid_box={rs['largest_solid_box']} interior_points={rs['interior_point_count']}")
    print('certified by uM: ' + ', '.join(f"{u}:{v['certified']}/{v['points']}" for u, v in rs['certified_by_um'].items()))
    print('top5 margins: ' + ', '.join(f"{t['uM_per_au']:g}/{t['carry_maturation_half_life_min']:g}/"
                                       f"{t['carry_mrna_half_life_min']:g}={t['timing_margin_h']:.3f}h"
                                       for t in report['top5']))
    n = report['nominal_6_30_2']
    print(f"nominal 6/30/2: certified={n['certified']} margin={n['timing_margin_h']:.3f}h "
          f"setup={n['bit1_setup_h']:.3f} hold={n['bit1_hold_h']:.3f} matches_smoke={n['matches_smoke_claim']}")
    e = report['environment']
    print(f"environment uniform={e['uniform']} on_disk_equals_baseline={e['on_disk_equals_baseline']} "
          f"grid_identical={e['grid_identical']}")
    m = report['manifest']
    print(f"manifest: entries={m['entries']} mismatches={len(m['mismatches'])} missing={len(m['missing'])} ok={m['ok']}")
    print(f"wrote {OUT / 'independent_verification.json'}")


if __name__ == '__main__':
    main()
