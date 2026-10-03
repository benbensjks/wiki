"""M-21: clock gate (Int0 AND) KEPT vs REMOVED, on the 34-state hybrid.

Answers Zeng's request for a comparison figure of the hybrid three-level
cascade with the Int0 gate removed vs retained.

Implemented per `时钟工作区扫描/复核_M21消融规格.md` (every finding independently
re-verified against the source before writing this file):

D1 (high)  Judgment must NOT go through `verify_hzh.analyse`: its `:164` calls
           `signals(y, zmh, receiver)` which RECOMPUTES clock/g1 from the
           receiver config, so a model-side bypass is invisible to it and the
           event arm would be judged on the un-bypassed gate.  Path used here:
           `common.sig_of` -> explicitly OVERWRITE clock/g1/_u2 -> `common.verdict`.
D2 (med)   The HZH certified clock is 0.30/2 (`hybrid_model.py:149`,
           `主结果_完整参数与来源.md:114`), NOT 0.10 -- 0.10 is the HZZ route's
           K_WORK.  (One combined clock_K for both arms: the gate is what the
           ablation removes, not a swept knob.)
D3 (med)   The originally proposed third arm is degenerate: once the clock is
           removed, clock_K no longer enters the dynamics.  Replaced by an
           explicit degeneracy control asserting max|dy| == 0.

Schema note: `common.verdict` returns `steady`/`events`/`certified_v1`/`cold`/
`read_windows`/`bit_margins` -- it does NOT return `counting_passed` /
`event_causality_passed` / `per_bit` (those live in `diagnose_round3.evaluate`),
so they are derived here from the same predicates.

Guards (all must pass or the ablation does not run):
  G1 certification reproduction on the KEPT arm, against `common.FROZEN_HZH`
  G2 `common.parity_guard` -- proves the injected-sig path equals `analyse`
     field-by-field when NOT bypassed
  G3 removed arm: the published clock signal is identically 1.0
  G4 degeneracy control: max|dy| == 0 between clock_K = 0.30 and 0.40
  G5 the clock-free gate reconstructed over the trajectory equals the model's
     own per-point `signals()['g1']` (guards the STATE-INDEX TRAP: on the full
     34-state array A1/F1 are at 17/18, not 0/1)
  G6 the removed arm's gate is non-trivial (peak > 0.1), i.e. a gate really
     does exist without the clock

G5/G6 were added after a real failure: the first 300 h run rebuilt the gate from
`y[:2]` of the full trajectory, reporting gate_peak = 0.0000 and a spurious
"no gate formed" event verdict, while the DYNAMICS were bypassed correctly.

Read-only with respect to every existing artifact: writes only into its own
results/<stamp>/ directory.  No global ZENG patching (H2); runs serial.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
SCAN = HERE.parent / '时钟工作区扫描'
# `common` lives in the scan directory and does not put its OWN directory on
# sys.path (it only adds the model dirs), so add it here before importing.
for _d in (str(SCAN), str(HERE)):
    if _d not in sys.path:
        sys.path.insert(0, _d)

import common as C                                                # noqa: E402
from hybrid_model import HbyReceiver, hill                        # noqa: E402
from model_hzh import HZHModel                                    # noqa: E402

# M21_HOURS lets a short smoke run exercise the whole path; M21_SMOKE=1 keeps
# going when the certification guards cannot pass at that duration.  Both
# default OFF, so a plain invocation is always the real 300 h run.
HOURS = float(os.environ.get('M21_HOURS', '300'))
SMOKE = os.environ.get('M21_SMOKE', '0') == '1'
K_CERT = 0.30          # D2: HZH certified clock threshold
K_ALT = 0.40           # degeneracy control only
DROP = int(C.VH.DROP)


def dump(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, ensure_ascii=False, default=str),
                    encoding='utf-8')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


class BypassReceiver(HbyReceiver):
    """Clock-gate REMOVED: only the third factor of g1 is replaced by 1.

    All three published keys are kept consistent -- `g1` (event arm),
    `u2_target_au_per_h` (DYNAMICS: `HbyReceiver.rhs` feeds this into the bit2
    I chain) and `clock_gate` (reporting).  HZH has NO `gate` key, so the
    round-4 key-name bug does not apply verbatim; the equivalent hazard here is
    leaving `u2_target_au_per_h` un-bypassed, which would bypass nothing at all.
    """

    def signals(self, y, int0):
        s = super().signals(y, int0)
        p, c = self.p, self.c
        a, f = y[:2]
        gate = (hill(a, p['K_A'][1], c.n_A1_gate)
                * (1.0 - hill(f, p['K_F'][1], p['n_F'][1])))
        s['clock_gate'] = 1.0
        s['g1'] = gate
        s['u2_target_au_per_h'] = p['alpha_Int'][1] * gate
        return s


def gate_no_clock_array(y, m):
    """The clock-free gate over a WHOLE 34-state trajectory (array version of
    `verify_hzh.py:69`).

    STATE INDEX TRAP: `HbyReceiver.signals(y, int0)` is called by the integrator
    with the 17-state TAIL (`model_hzh.py:90` passes `y[17:]`), so inside it
    `y[:2]` means [A1, F1].  Here `y` is the full 34-state trajectory, where A1
    and F1 sit at IDX['A1']=17 / IDX['F1']=18.  Writing `y[:2]` here silently
    rebuilds the gate from `b0_M_I`/`b0_I_u`, giving a near-zero gate that the
    event arm then reports as "no gate formed" -- a SPURIOUS failure.
    `gate_reconstruction_parity()` below asserts this cannot happen again.
    """
    p, c = m.tail.p, m.tail.c
    a, f = y[C.VH.IDX['A1']], y[C.VH.IDX['F1']]
    return (C.VH.hill_array(a, p['K_A'][1], c.n_A1_gate)
            * (1.0 - C.VH.hill_array(f, p['K_F'][1], p['n_F'][1])))


def gate_reconstruction_parity(y, m, step=97):
    """G5: the array reconstruction must equal the model's own per-point gate.

    Evaluates `BypassReceiver.signals` at sampled time points with the 17-state
    tail and compares to `gate_no_clock_array` on the full trajectory.  This is
    the guard for the state-index trap described above -- it is what turns that
    mistake from a silent wrong verdict into a hard failure.
    """
    worst, worst_i = 0.0, -1
    tail = m.tail
    for i in range(0, y.shape[1], step):
        ref = float(tail.signals(y[17:, i], y[C.VH.IDX['b0_I'], i])['g1'])
        mine = float(gate_no_clock_array(y[:, i:i + 1], m)[0])
        if abs(ref - mine) > worst:
            worst, worst_i = abs(ref - mine), i
    return worst, worst_i


def judge(t, y, m, bypassed: bool):
    """D1 path: sig_of, then explicit overwrite, then common.verdict."""
    sig = C.sig_of(m, y, 'hzh')
    if bypassed:
        g = gate_no_clock_array(y, m)
        sig['clock'] = np.ones_like(np.asarray(sig['clock'], dtype=float))
        sig['g1'] = g
        # sig_of derived _u2 from the OLD g1 -- re-derive after overwriting
        sig['_u2'] = g * float(m.tail.p['alpha_Int'][1])
    return sig, C.verdict(t, sig)


def summarise(t, y, m, sig, v, tag, clock_K, removed):
    s0 = v['events']['bit0_to_bit1']
    s1 = v['events']['bit1_to_bit2']
    segs = C.VH.segments(t, sig['g1'], C.VH.GATE_THRESHOLD,
                         min_dose=C.VH.GATE_MIN_DOSE_H)
    dose = [s['dose'] for s in segs] or [float('nan')]
    dur = [(s['end_h'] - s['start_h']) for s in segs] or [float('nan')]
    sel = v['read_windows'][DROP:]
    commit = {f'bit{i}': (min(r['commitment'][i] for r in sel) if sel else None)
              for i in (0, 1, 2)}
    unlab = {f'bit{i}': sum(r['labels'][i] is None for r in sel) for i in (0, 1, 2)}
    return dict(
        tag=tag, arm=('removed' if removed else 'kept'), clock_K=clock_K,
        hours=float(t[-1]), n_samples=int(len(t)),
        certified=bool(v['certified_v1']),
        counting=bool(v['steady']['passed']),
        events=bool(all(e['passed'] for e in v['events'].values())),
        steady_reads=int(v['steady']['reads']),
        sequence=str(v['steady']['sequence']),
        cold_reads=int(v['cold']['reads']),
        s2_crossings=len(v['crossings']['S2']),
        stage0=[s0['reverse_events'], s0['gate_events'], s0['flip_events']],
        stage1=[s1['reverse_events'], s1['gate_events'], s1['flip_events']],
        gate_segments=int(len(segs)),
        gate_dose_min=float(min(dose)), gate_dose_max=float(max(dose)),
        gate_dur_min=float(min(dur)), gate_dur_max=float(max(dur)),
        gate_peak=float(np.max(sig['g1'])),
        clock_signal_min=float(np.min(sig['clock'])),
        clock_signal_max=float(np.max(sig['clock'])),
        S0_min=float(sig['S0'].min()), S0_max=float(sig['S0'].max()),
        S1_min=float(sig['S1'].min()), S1_max=float(sig['S1'].max()),
        S2_min=float(sig['S2'].min()), S2_max=float(sig['S2'].max()),
        commitment=commit, unlabelled=unlab,
        margin_h=v['global_min_timing_margin_h'],
    )


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    out = HERE / 'results' / (('m21_SMOKE_' if SMOKE else 'm21_')
                              + datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir(parents=True, exist_ok=False)
    dump(out / 'status.json', dict(status='running', hours=HOURS, smoke=SMOKE,
                                   K_cert=K_CERT, K_alt=K_ALT))
    print('OUTPUT DIR', out, flush=True)
    if SMOKE:
        print('*** SMOKE RUN (%g h): numbers are NOT a result; the directory '
              'name is prefixed m21_SMOKE_ ***' % HOURS, flush=True)

    print('\n[0] shared Han upstream, %g h' % HOURS, flush=True)
    t0 = time.time()
    han = C.HanInput(HOURS)
    print('    HanInput ok  (%.1f s)' % (time.time() - t0), flush=True)

    rows, trajs = [], {}

    # ---------------- kept arm = the certification itself ----------------
    print('\n[1] KEPT arm: clock = H(b0_I; %.2f, 2.0), %g h' % (K_CERT, HOURS),
          flush=True)
    t0 = time.time()
    mk, _ = C.build('hzh', han)
    tk, yk = C.integrate(mk, hours=HOURS)
    sk, vk = judge(tk, yk, mk, bypassed=False)
    print('    integrated in %.1f s' % (time.time() - t0), flush=True)
    rk = summarise(tk, yk, mk, sk, vk, 'kept_K0.30', K_CERT, False)
    rows.append(rk)
    trajs['kept_K0.30'] = (tk, yk)
    print('    certified=%s counting=%s events=%s reads=%d seq=%s margin=%s'
          % (rk['certified'], rk['counting'], rk['events'], rk['steady_reads'],
             rk['sequence'], rk['margin_h']), flush=True)

    # ---------------- G1: certification reproduction ----------------
    F = C.FROZEN_HZH
    g1 = dict(certified=rk['certified'] is F['certified_v1'],
              reads=rk['steady_reads'] == 19,
              sequence=rk['sequence'] == F['steady'],
              stage0=tuple(rk['stage0']) == F['stage0'],
              stage1=tuple(rk['stage1']) == F['stage1'],
              margin=(rk['margin_h'] is not None
                      and abs(float(rk['margin_h']) - F['margin_h']) < 1e-6))
    print('    G1 certification reproduction:', json.dumps(g1), flush=True)

    # ---------------- G2: parity guard ----------------
    print('\n[2] G2 parity_guard (verdict == analyse when NOT bypassed, on the '
          '300 h trajectory)', flush=True)
    try:
        pg = C.parity_guard(han, traj=(tk, yk))
        g2 = True
        print('    parity OK:', json.dumps(pg, ensure_ascii=False), flush=True)
    except Exception as exc:                                     # noqa: BLE001
        g2 = False
        print('    PARITY FAILED:', exc, flush=True)

    guards = dict(G1_cert_reproduction=g1, G2_parity=bool(g2), G3=None, G4=None,
                  G5_gate_recon_parity=None, G6_gate_nontrivial=None)
    if not (all(g1.values()) and g2):
        if not SMOKE:
            dump(out / 'status.json', dict(status='FAILED_GUARD', guards=guards))
            print('\n*** GUARD FAILED -- ablation not run ***', flush=True)
            raise SystemExit(2)
        print('\n    (SMOKE: guard failure ignored, continuing to exercise the '
              'removed arm)', flush=True)

    # ---------------- removed arm ----------------
    print('\n[3] REMOVED arm: clock == 1, g1 = H(A1;1.2,6)*(1-H(F1;0.4,4))',
          flush=True)
    t0 = time.time()
    mr, _ = C.build('hzh', han)
    mr.tail = BypassReceiver(mr.tail.c)
    mr.copies_per_au = C.COPIES_PER_UM_PER_FL * mr.tail.c.receiver_uM_per_au
    tr, yr = C.integrate(mr, hours=HOURS)
    sr, vr = judge(tr, yr, mr, bypassed=True)
    print('    integrated in %.1f s' % (time.time() - t0), flush=True)
    rr = summarise(tr, yr, mr, sr, vr, 'removed_K0.30', K_CERT, True)
    rows.append(rr)
    trajs['removed_K0.30'] = (tr, yr)
    print('    certified=%s counting=%s events=%s reads=%d seq=%s'
          % (rr['certified'], rr['counting'], rr['events'], rr['steady_reads'],
             rr['sequence']), flush=True)
    print('    stage0=%s stage1=%s gate_segments=%d gate_peak=%.4f'
          % (rr['stage0'], rr['stage1'], rr['gate_segments'], rr['gate_peak']),
          flush=True)

    # ---------------- G3: the bypass is really published as 1 ----------------
    g3 = (rr['clock_signal_max'] == 1.0 and rr['clock_signal_min'] == 1.0)
    print('    G3 clock published as identically 1: %s (min=%.6f max=%.6f)'
          % (g3, rr['clock_signal_min'], rr['clock_signal_max']), flush=True)

    # ---------------- G5/G6: gate reconstruction sanity ----------------
    # Added after a real failure: the first 300 h run rebuilt the clock-free
    # gate from `y[:2]` on the FULL 34-state trajectory (states b0_M_I/b0_I_u)
    # instead of IDX['A1']/IDX['F1'], giving gate_peak = 0.0000 and a spurious
    # "no gate formed" event verdict.  G5 compares the reconstruction against
    # the model's own per-point signals; G6 asserts the gate is non-trivial.
    print('\n[3b] G5 gate-reconstruction parity (array vs model per-point)',
          flush=True)
    gw, gi = gate_reconstruction_parity(yr, mr)
    g5 = (gw == 0.0)
    print('    max|dGate| = %.3e at sample %d -> %s'
          % (gw, gi, 'OK' if g5 else 'MISMATCH'), flush=True)
    g6 = bool(rr['gate_peak'] > 0.1)
    print('    G6 removed-arm gate is non-trivial: gate_peak = %.4f -> %s'
          % (rr['gate_peak'], g6), flush=True)

    # ---------------- G4: degeneracy control ----------------
    print('\n[4] G4 degeneracy control: REMOVED arm at clock_K=%.2f vs %.2f'
          % (K_ALT, K_CERT), flush=True)
    t0 = time.time()
    ma, _ = C.build('hzh', han, clock_K=K_ALT)
    ma.tail = BypassReceiver(ma.tail.c)
    ta, ya = C.integrate(ma, hours=HOURS)
    dmax = float(np.max(np.abs(np.asarray(ya) - np.asarray(yr))))
    g4 = (dmax == 0.0)
    print('    integrated in %.1f s ; max|dy| = %.3e -> %s'
          % (time.time() - t0, dmax,
             'DEGENERATE as predicted' if g4 else 'NOT degenerate'), flush=True)
    trajs['removed_K0.40'] = (ta, ya)

    guards.update(G3=bool(g3), G4=bool(g4), G4_max_abs_dy=dmax,
                  G5_gate_recon_parity=bool(g5), G5_max_abs_dgate=gw,
                  G6_gate_nontrivial=bool(g6), G6_gate_peak=rr['gate_peak'])

    # G5/G6 guard the JUDGMENT of the removed arm: if the reconstructed gate is
    # wrong the event verdict is wrong, so do not publish the run.
    if not (g5 and g6) and not SMOKE:
        dump(out / 'status.json', dict(status='FAILED_GATE_SANITY', guards=guards))
        print('\n*** G5/G6 FAILED -- event verdict of the removed arm is NOT '
              'trustworthy; run not published ***', flush=True)
        raise SystemExit(3)

    # ---------------- save ----------------
    for name, (t, y) in trajs.items():
        np.savez_compressed(out / f'traj_{name}.npz', time_h=t, states=y)
    cols = ['tag', 'arm', 'clock_K', 'hours', 'certified', 'counting', 'events',
            'steady_reads', 'sequence', 'cold_reads', 's2_crossings', 'stage0',
            'stage1', 'gate_segments', 'gate_dose_min', 'gate_dose_max',
            'gate_dur_min', 'gate_dur_max', 'gate_peak', 'clock_signal_min',
            'clock_signal_max', 'S0_min', 'S0_max', 'S1_min', 'S1_max',
            'S2_min', 'S2_max', 'margin_h']
    with (out / 'ablation.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([('%.6g' % r[c]) if isinstance(r[c], float) else r[c]
                        for c in cols])
    dump(out / 'summary.json', dict(
        smoke=SMOKE,
        question=('on the 34-state hybrid HZH (HZH_MOD8_CAUSAL_V1), what changes '
                  'when the Int0 clock AND gate is removed from g1?'),
        note=('judgment via common.sig_of + explicit overwrite + common.verdict; '
              'NOT via verify_hzh.analyse (D1)'),
        K_cert=K_CERT, K_alt=K_ALT, hours=HOURS, guards=guards, cases=rows))
    dump(out / 'source_hashes.json', C.source_hashes())
    dump(out / 'SHA256SUMS.json',
         {str(p.relative_to(out)): sha256(p)
          for p in sorted(out.rglob('*')) if p.is_file()})
    dump(out / 'status.json', dict(status='COMPLETED', rows=len(rows),
                                   guards=guards))
    print('\nGUARDS:', json.dumps({k: v for k, v in guards.items()
                                   if k != 'G1_cert_reproduction'}),
          flush=True)
    print('OUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
