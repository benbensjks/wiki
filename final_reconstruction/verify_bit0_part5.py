"""Part 5: is there any toggle at all when growth dilution is switched back on?

Nothing else is changed: only add_growth and the a.u.->uM conversion are moved.
"""
from __future__ import annotations

import json
from dataclasses import asdict

import pandas as pd

from model import Extension, ROOT
from verify_bit0_part2 import OUT, CAND
from verify_bit0_part4 import run_one


def main(hours=120.0):
    OUT.mkdir(parents=True, exist_ok=True)
    base = asdict(Extension()); base.update(CAND)
    rows = []
    for u in (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0):
        for growth in (True, False):
            c = base.copy(); c['uM_per_au'] = u; c['add_growth'] = growth
            r = run_one(c, hours=hours)
            rows.append(r)
            print(f"uM_per_au={u:6.1f} add_growth={str(growth):5s} commit={r['codes_commitment']} "
                  f"({r['median_commitment']:.2f}) dwell={r['dwell_low']:.3f}/{r['dwell_high']:.3f} "
                  f"instant={r['codes_instant']}", flush=True)
            pd.DataFrame(rows).to_csv(OUT / 'bit0_growth_conversion_map.csv', index=False)
    d = pd.DataFrame(rows)
    ok = d[d.stored_commitment]
    return dict(hours=hours, points=len(d),
                stored_by_commitment=int(len(ok)),
                stored_with_growth=int(ok.add_growth.sum()),
                stored_without_growth=int((~ok.add_growth).sum()),
                stored_points=[dict(uM_per_au=float(r.uM_per_au), add_growth=bool(r.add_growth),
                                    codes=r.codes_commitment, dwell=[r.dwell_low, r.dwell_high])
                               for r in ok.itertuples()],
                frame=rows)


if __name__ == '__main__':
    print(json.dumps(main(), ensure_ascii=False, indent=2))
