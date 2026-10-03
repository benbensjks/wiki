"""Cheap pre-flight for expr_corner.py: no integration, only construction.

(1) SplitCarryReceiver vs HbyReceiver equivalence at equal maturation rates.
(2) build() readback for all 25 point x config combinations, with an explicit
    check that the four previously-mislabelled rows really carry the high
    thresholds now.
"""
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import scan_oat as S                                              # noqa: E402
import expr_corner as EC                                          # noqa: E402
from hybrid_model import HbyConfig                                # noqa: E402
from dataclasses import asdict                                    # noqa: E402
from hybrid_model import HbyReceiver                              # noqa: E402

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

print('=== (1) SplitCarryReceiver equivalence at equal rates ===')
eq = EC.self_test_receiver_equivalence()
print(json.dumps(eq, indent=2))
assert eq['equivalent'], 'EQUIVALENCE FAILED'

print()
print('=== (1b) SplitCarryReceiver differs when rates really differ ===')
cfg = HbyConfig(**asdict(EC.FROZEN))
par = HbyReceiver(cfg)
chi = EC.SplitCarryReceiver(cfg, cfg.carry_maturation_min * 0.5,
                            cfg.carry_maturation_min * 2.0)
y = par.initial_state()
d = par.rhs(0.0, y, pb1=1.0, int0=0.5) - chi.rhs(0.0, y, pb1=1.0, int0=0.5)
print('max|rhs_parent - rhs_split| at a1=x0.5,f1=x2 =', float(abs(d).max()))
assert abs(d).max() > 1e-6, 'split receiver has no effect -> wiring broken'

print()
print('=== (2) build() readback over all 25 combinations ===')
han = S.HanInput(S.HOURS)
prototype = S.HZHModel('han', han)
bad = 0
for p in EC.POINTS:
    for c in EC.CONFIGS:
        m, meta = EC.build(prototype, p, c)
        rb = meta['readback']
        want_KA = EC.ZENG['K_A'][1] * EC.POINTS[p][1]
        want_KF = EC.ZENG['K_F'][1] * EC.POINTS[p][2]
        ok = (abs(rb['K_A1'] - want_KA) < 1e-12
              and abs(rb['K_F1'] - want_KF) < 1e-12
              and rb['n_tuple_keys'] == 0)
        if not ok:
            bad += 1
        flag = '' if ok else '   <-- MISMATCH'
        print(f"  {p:9s} {c:14s} split={int(meta['split_receiver'])} "
              f"K_A1={rb['K_A1']:.4f} (want {want_KA:.4f}) "
              f"K_F1={rb['K_F1']:.4f} (want {want_KF:.4f}) "
              f"kmc_A/kmc_F={rb['kmc_A']:.4f}/{rb['kmc_F']:.4f} "
              f"tuplekeys={rb['n_tuple_keys']}{flag}")

print()
print('=== (2b) the four previously-broken rows, explicitly ===')
for p in ('K_A_high', 'K_F_high'):
    for c in ('A1fast_F1slow', 'A1slow_F1fast'):
        m, meta = EC.build(prototype, p, c)
        rb = meta['readback']
        print(f"  {p:9s} {c:14s} K_A1={rb['K_A1']:.4f} K_F1={rb['K_F1']:.4f} "
              f"uM={rb['uM_per_au']:.4f} n={rb['n_A1_gate']:.1f} "
              f"carry_mat={rb['carry_maturation_min']:.3f} "
              f"kmc_A={rb['kmc_A']:.4f} kmc_F={rb['kmc_F']:.4f}")

print()
print('MISMATCHES:', bad)
assert bad == 0
print('ALL PRE-FLIGHT CHECKS PASSED')
