"""How wide is the surviving bit0 window? Sequential region map (no process pool:
this sandbox blocks the multiprocessing pipe)."""
from __future__ import annotations

import itertools
import json

import numpy as np
import pandas as pd

from model import ROOT
from verify_bit0_part2 import OUT, CAND, _region_one  # noqa: F401
from dataclasses import asdict
from model import Extension

GRID = dict(uM_per_au=(6.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 18.0),
            maturation_half_life_min=(12.0, 16.0, 18.0, 20.0, 22.0, 26.0, 32.0))


def main(hours=120.0):
    OUT.mkdir(parents=True, exist_ok=True)
    base = asdict(Extension()); base.update(CAND)
    rows = []
    for u, mt in itertools.product(GRID['uM_per_au'], GRID['maturation_half_life_min']):
        c = base.copy(); c['uM_per_au'] = u; c['maturation_half_life_min'] = mt
        r = _region_one(c, hours)
        rows.append(r)
        print(f"uM_per_au={u:5.1f} t_mat={mt:5.1f} -> {r['codes']:>10s} commit={r['median_commitment']:.3f} "
              f"stored={r['stored_bit']}", flush=True)
        pd.DataFrame(rows).to_csv(OUT / 'bit0_region_map.csv', index=False)
    d = pd.DataFrame(rows)
    ok = d[d.stored_bit]
    pivot = d.pivot_table(index='uM_per_au', columns='maturation_half_life_min', values='stored_bit')
    summary = dict(hours=hours, points=len(d), stored_bit_points=int(len(ok)),
                   stored_region=(None if not len(ok) else dict(
                       uM_per_au=[float(ok.uM_per_au.min()), float(ok.uM_per_au.max())],
                       maturation_half_life_min=[float(ok.maturation_half_life_min.min()),
                                                  float(ok.maturation_half_life_min.max())])),
                   table={str(k): {str(kk): int(vv) for kk, vv in v.items()} for k, v in pivot.to_dict().items()})
    (OUT / 'bit0_region_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
