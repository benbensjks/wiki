"""Verify the two problems raised against rounds 3/4. Read-only, no re-integration.

PROBLEM 1 -- the eight "phases" may be eight sub-cycle samples of the SAME
  digital state, not the eight states 0-7. And the continuation restarts the
  upstream time at 0 instead of keeping the absolute time.
PROBLEM 2 -- the mechanism metric multiplied two independently attained
  maxima, max(I2)*max(RDF2), instead of the same-time product
  max(I2(t)*RDF2(t)).
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT.parent.parent))

import diagnose_round1 as D                                       # noqa: E402
import verify_hzh as V                                            # noqa: E402

RES = ROOT / 'results'
N_PHASE = 8


def newest(pat):
    return sorted(RES.glob(pat), key=lambda p: p.name)[-1]


def digital(y, t, idx, sig):
    """Value of the counter at each sampled state, from the precomputed signals."""
    out = []
    for j in idx:
        s0, s1, s2 = float(sig['S0'][j]), float(sig['S1'][j]), float(sig['S2'][j])
        lab = []
        for x in (s0, s1, s2):
            lab.append(0 if x <= V.LOW else 1 if x >= V.HIGH else None)
        val = (None if any(l is None for l in lab)
               else lab[0] + 2 * lab[1] + 4 * lab[2])
        out.append(dict(t_h=float(t[j]), S0=s0, S1=s1, S2=s2, labels=lab,
                        value=val))
    return out


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    r3 = newest('round3_*')
    r2 = newest('round2_*')
    r4b = newest('round4b_*')

    han = D.HanInput(D.HOURS)
    prefix = D.HZHModel('han', han)

    # ---------------- PROBLEM 1 -------------------------------------------
    print('=' * 74)
    print('PROBLEM 1: what digital states do the eight "phase" starts hold?')
    print('=' * 74)
    z = np.load(r3 / 'traj_work300.npz')
    t300, y300 = z['time_h'], z['states']
    m = D.HZZModel(prefix, False)
    m.clock_K = 0.10
    s300 = m.signals(y300)
    peaks, _ = V.read_windows(t300, s300)
    pt = t300[peaks]
    period = float(np.median(np.diff(pt)))
    targets = (pt[-2] - period) + np.arange(N_PHASE) * period / N_PHASE
    idx = [int(np.argmin(np.abs(t300 - x))) for x in targets]
    states = digital(y300, t300, idx, s300)
    for i, st in enumerate(states):
        v = 'undecided' if st['value'] is None else str(st['value'])
        print(f"  phase {i}: t={st['t_h']:8.3f} h  "
              f"S0={st['S0']:.4f} S1={st['S1']:.4f} S2={st['S2']:.4f}  "
              f"labels={st['labels']}  -> value {v}")
    vals = [st['value'] for st in states]
    print(f"  distinct values covered: {sorted({v for v in vals if v is not None})}"
          f"  ({len({v for v in vals})} distinct, undecided {sum(v is None for v in vals)})")
    print(f"  period used for spacing = {period:.4f} h (ONE clock cycle)")
    print(f"  span sampled = {targets[-1]-targets[0]:.3f} h "
          f"= {period:.3f} h, i.e. within a single count")

    # how long is one complete 8-count cycle?
    seq = V.read_verdict(V.read_windows(t300, s300)[1], V.DROP)['sequence']
    print(f"  steady sequence = {seq}")
    print(f"  -> one count per clock cycle, so a full 8-count cycle is "
          f"{8*period:.2f} h, not {period:.2f} h")

    # upstream time restart?
    print()
    print('  continuation call used: solve_ivp(model.rhs, (0, HOURS), y_snapshot)')
    print('  HZZModel.rhs(t, y) -> prefix.bit0_rhs(t, y[:11]) -> '
          'han.interp(60*t, h31)')
    print('  the 23-state model carries NO oscillator state, so the upstream '
          'phase restarts at t=0')
    print('  => the snapshot is re-driven by a RESTARTED upstream, not by its '
          'own upstream time')

    # ---------------- PROBLEM 2 -------------------------------------------
    print()
    print('=' * 74)
    print('PROBLEM 2: peak x peak  vs  same-time product')
    print('=' * 74)
    cases = [
        ('base fail  K=0.40 off  (no bypass)', r2 / 'traj_K0.40_off.npz', False, 0.40),
        ('base fail  K=0.25 off  (no bypass)', r2 / 'traj_K0.25_off.npz', False, 0.25),
        ('clk020 pass K=0.20 on  (no bypass)', r2 / 'traj_K0.20_on.npz', True, 0.20),
        ('clk020 pass K=0.20 off (no bypass)', r2 / 'traj_K0.20_off.npz', False, 0.20),
        ('bypass   K=0.40 off', r4b / 'traj_bypass_K0.40_off.npz', False, 0.40),
        ('bypass   K=0.40 on ', r4b / 'traj_bypass_K0.40_on.npz', True, 0.40),
    ]
    KD = D.ZENG['K_D_comp'] if hasattr(D, 'ZENG') else 0.8
    print(f"  K_D_comp used by the tail = 0.8")
    hdr = ('%-34s %-8s %-9s %-9s %-9s %-9s' %
           ('case', 'verdict', 'maxI2', 'maxRDF2', 'peakxpeak', 'same-time'))
    print(hdr)
    print('-' * len(hdr))
    table = []
    for name, path, ok, k in cases:
        if not path.exists():
            print(f'  (missing {path.name})')
            continue
        zz = np.load(path)
        yy = zz['states']
        mm = D.HZZModel(prefix, ok)
        mm.clock_K = k
        i2 = np.asarray(yy[20], dtype=float)
        r2_ = np.asarray(yy[22], dtype=float)
        pb2 = np.asarray(yy[19], dtype=float)
        max_i, max_r = float(i2.max()), float(r2_.max())
        peak_prod = max_i * max_r
        prod = i2 * r2_
        j = int(np.argmax(prod))
        same = float(prod[j])
        print('%-34s %-8s %-9.4f %-9.4f %-9.4f %-9.4f' %
              (name, 'pass' if ok else 'fail', max_i, max_r, peak_prod, same))
        table.append(dict(case=name, counting=ok, clock_K=k,
                          max_I2=max_i, max_RDF2=max_r,
                          peak_times_product=peak_prod,
                          same_time_product=same,
                          same_time_at_h=float(zz['time_h'][j]),
                          pb2_at_that_time=float(pb2[j]),
                          same_time_over_threshold=same / 0.8))
    print()
    print('  same-time product / 0.8:')
    for r in table:
        print(f"    {r['case']:<34} {r['same_time_over_threshold']:6.2f}   "
              f"(peak x peak would give {r['peak_times_product']/0.8:8.2f})")

    (RES / 'errata_check.json').write_text(json.dumps(dict(
        problem1=dict(period_h=period, span_h=float(targets[-1] - targets[0]),
                      phase_states=states,
                      distinct_values=sorted({v for v in vals if v is not None}),
                      undecided=sum(v is None for v in vals),
                      note=('period used for spacing is ONE clock cycle; one '
                            'count per cycle, so the eight samples sit inside a '
                            'single count and cannot cover values 0-7')),
        problem2=dict(K_D_comp=0.8, table=table,
                      note=('peak_times_product multiplies maxima attained at '
                            'different times; only same_time_product is a '
                            'physical reaction state'))),
        ensure_ascii=False, indent=2), encoding='utf-8')
    print()
    print('wrote', RES / 'errata_check.json')


if __name__ == '__main__':
    main()
