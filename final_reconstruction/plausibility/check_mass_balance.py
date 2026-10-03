"""Check 8 - mass-balance audit for the explicit Int-RDF complex.

Our reconstruction differs from Zeng's reduced model in one structural way that
has never been audited: we carry an explicit complex C and let it *consume*
free Int and free RDF, whereas in her model I*R appears only inside the reverse
Hill and consumes nothing.  Two consequences must be quantified:

  1. implementation correctness - summing the implemented derivatives must give
         d(I + C)/dt = kmat*I_u - gamma_int*I - delta_C*C
         d(R + C)/dt = kmat*R_u - gamma_rdf*R - delta_C*C
     i.e. the binding / unbinding terms must cancel exactly between I (or R) and
     C.  A sign error or a missing term shows up here immediately.
  2. magnitude - what fraction of the Int / RDF inventory is sequestered in C.

Why the identity is checked at the RHS level
--------------------------------------------
A finite-difference residual from a 2-minute sampled trajectory is dominated by
O(h^2) truncation error: for the Int pool that is already ~1 %, which would
swamp the signal.  The identity is therefore verified by calling the model's own
right-hand side and combining its entries, which is exact; the finite-difference
residual is still reported, but only as a numerical diagnostic with a loose,
honest tolerance.

Pre-registered criteria
-----------------------
FLAG (implementation error) if the RHS-level residual exceeds 1e-9 of the
largest term in that balance.
FLAG (numerics) if the finite-difference residual exceeds 10 %; intermediate
values are reported as truncation noise.
The sequestration fraction is reported without pass/fail: it is a declaration
that has to go into the report.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT, build_threebit, frozen_carry0,  # noqa: E402
                                 restore, write_manifest)
from model import ZENG  # noqa: E402
from working_point import working_point_block  # noqa: E402

RHS_TOL = 1e-9
FD_TOL = 0.10
BIT_OFFSETS = {'bit0': (6, 0), 'bit1': (17, 1), 'bit2': (34, 2)}
N_PROBE_POINTS = 400


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hours', type=float, default=300.0)
    ap.add_argument('--sample-min', type=float, default=2.0)
    ap.add_argument('--max-step-min', type=float, default=2.0)
    args = ap.parse_args()

    model, patch = build_threebit(carry1=frozen_carry0())
    try:
        sol = model.simulate(hours=args.hours, sample_min=args.sample_min,
                             max_step_min=args.max_step_min)
        t, y = sol.t, sol.y
        kmat = model.base.kmat
        growth = model.base.growth
        delta_C = model.e.complex_decay_h
        stride = max(1, t.size // N_PROBE_POINTS)
        idxs = list(range(0, t.size, stride))

        rows = []
        for name, (off, index) in BIT_OFFSETS.items():
            gamma_I = ZENG['gamma_int'][index]
            gamma_R = ZENG['gamma_rdf']
            worst_rhs_I = worst_rhs_R = 0.0
            scale_I = scale_R = 0.0
            for k in idxs:
                d = model.rhs(float(t[k]), y[:, k])
                I_u, I, C = y[off + 1, k], y[off + 2, k], y[off + 9, k]
                R_u, R = y[off + 7, k], y[off + 8, k]
                ana_I = kmat * I_u - (gamma_I + growth) * I - (delta_C + growth) * C
                ana_R = kmat * R_u - (gamma_R + growth) * R - (delta_C + growth) * C
                worst_rhs_I = max(worst_rhs_I, abs((d[off + 2] + d[off + 9]) - ana_I))
                worst_rhs_R = max(worst_rhs_R, abs((d[off + 8] + d[off + 9]) - ana_R))
                scale_I = max(scale_I, abs(kmat * I_u), abs((gamma_I + growth) * I),
                              abs((delta_C + growth) * C))
                scale_R = max(scale_R, abs(kmat * R_u), abs((gamma_R + growth) * R),
                              abs((delta_C + growth) * C))
            for label, inv, source, loss, gamma, worst, scale in (
                    ('Int', y[off + 2] + y[off + 9], kmat * y[off + 1], y[off + 2], gamma_I, worst_rhs_I, scale_I),
                    ('RDF', y[off + 8] + y[off + 9], kmat * y[off + 7], y[off + 8], gamma_R, worst_rhs_R, scale_R)):
                fd = float(np.max(np.abs(np.gradient(inv, t) -
                                         (source - (gamma + growth) * loss -
                                          (delta_C + growth) * y[off + 9]))))
                Cmax = float(np.max(y[off + 9]))
                freemax = float(np.max(loss))
                rows.append(dict(bit=name, species=label,
                                 rhs_residual=worst, rhs_scale=scale,
                                 rhs_relative=worst / max(scale, 1e-30),
                                 rhs_flagged=bool(worst / max(scale, 1e-30) > RHS_TOL),
                                 fd_residual=fd, fd_relative=fd / max(scale, 1e-30),
                                 fd_flagged=bool(fd / max(scale, 1e-30) > FD_TOL),
                                 free_max=freemax, complex_max=Cmax,
                                 sequestered_fraction_of_free_max=Cmax / max(freemax, 1e-30)))
    finally:
        restore(patch)

    df = pd.DataFrame(rows)
    verdict = dict(
        implementation_ok=bool(not df.rhs_flagged.any()),
        worst_rhs_relative=float(df.rhs_relative.max()),
        rhs_tolerance=RHS_TOL,
        finite_difference_dominated_by_truncation=bool(not df.fd_flagged.any()),
        worst_fd_relative=float(df.fd_relative.max()),
        fd_note=('at %g min sampling the central-difference truncation error is O(h^2) '
                 'and roughly 1 %% for the Int pool; only the RHS-level identity is '
                 'conclusive about the implementation' % args.sample_min),
        max_sequestered_fraction=float(df.sequestered_fraction_of_free_max.max()),
        sequestered_by_pool={f'{r.bit}.{r.species}': round(float(r.sequestered_fraction_of_free_max), 4)
                             for r in df.itertuples()},
        declaration=("the explicit complex is a structural addition relative to Zeng's "
                     'reduced model; it sequesters up to %.1f %% of the peak free Int pool'
                     % (100 * float(df.sequestered_fraction_of_free_max.max()))),
        working_point=working_point_block(model))
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / 'mass_balance.csv', index=False, encoding='utf-8')
    (OUT / 'mass_balance.json').write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    pd.set_option('display.width', 240)
    print(df[['bit', 'species', 'rhs_relative', 'rhs_flagged', 'fd_relative', 'fd_flagged',
              'free_max', 'complex_max', 'sequestered_fraction_of_free_max']].to_string(index=False))
    print()
    print(json.dumps(verdict, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
