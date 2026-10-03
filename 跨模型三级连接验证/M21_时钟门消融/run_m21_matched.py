"""M-21 third arm: clock gate REMOVED, with the TOTAL Int2 source input matched.

Balancing rule chosen by the reviewer: option (c) -- integrate the Int2 source
u_2(t) = alpha_Int[1] * g1(t) over IDENTICAL time boundaries in both arms,
preferably one complete post-burn-in mod-8 cycle (8 read windows, which contain
exactly 2 bit1->bit2 carries plus the dwell periods between them), then adjust
ONLY the removed arm's alpha_Int[1] so the two total integrals agree.

Why not the other options:
  (a) matching only the event-threshold-identified gate windows leaves the
      outside-window input different;
  (b) matching peaks does not control the total input.
And crucially: the previously reported 0.2200 h / 0.3609 h are integrals INSIDE
the identified gate pulses, so their ratio must NOT be used as the balancing gain
-- the integral is recomputed here from the full saved g1(t) over whole cycles.

Fixed by the same rule as the 2-arm run:
  * judgment via common.sig_of + explicit overwrite of clock/g1/_u2 + common.verdict
    (never verify_hzh.analyse, which recomputes g1 WITH the clock -- D1);
  * clock working point 0.30/2 (HZH certified), not 0.10 (D2);
  * bit2's OWN setup/hold reported separately, since the global minimum timing
    margin is set by bit0 and therefore cannot discriminate these arms.

Reads the 2-arm trajectories for the matching window; runs exactly ONE new 300 h
integration (the matched arm).  Writes only into its own results/<stamp>/.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
SCAN = HERE.parent / '时钟工作区扫描'
for _d in (str(SCAN), str(HERE)):
    if _d not in sys.path:
        sys.path.insert(0, _d)

import common as C                                                # noqa: E402
from run_m21_clock_bypass import (BypassReceiver, gate_no_clock_array,  # noqa: E402
                                  judge, summarise, K_CERT, dump, sha256)

HOURS = 300.0
SRC_2ARM = 'm21_20261001_225706'
MATCH_WINDOWS = 8          # one mod-8 cycle = 8 read windows = 2 bit1 carries


def newest(root: Path, pat: str) -> Path:
    d = sorted([p for p in (root / 'results').glob(pat) if p.is_dir()])
    return d[-1]


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    src = HERE / 'results' / SRC_2ARM
    out = HERE / 'results' / ('m21matched_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir(parents=True, exist_ok=False)
    dump(out / 'status.json', dict(status='running', hours=HOURS, smoke=False))
    print('OUTPUT DIR', out, flush=True)
    print('2-arm source:', src.name, flush=True)

    print('\n[0] shared upstream + rebuild both 2-arm gates from saved states',
          flush=True)
    han = C.HanInput(HOURS)
    zk = np.load(src / 'traj_kept_K0.30.npz')
    zr = np.load(src / 'traj_removed_K0.30.npz')
    tk, yk = zk['time_h'], zk['states']
    tr, yr = zr['time_h'], zr['states']

    mk = C.HZHModel('han', han)
    mr = C.HZHModel('han', han)
    mr.tail = BypassReceiver(mr.tail.c)
    sk = C.sig_of(mk, yk, 'hzh')
    sr = C.sig_of(mr, yr, 'hzh')
    # D1: overwrite after sig_of, else g1 still carries the clock
    g_rem = gate_no_clock_array(yr, mr)
    sr['clock'] = np.ones_like(np.asarray(sr['clock'], dtype=float))
    sr['g1'] = g_rem
    sr['_u2'] = g_rem * float(mr.tail.p['alpha_Int'][1])
    g_kept = np.asarray(sk['g1'], dtype=float)

    _, reads = C.VH.read_windows(tk, sk)          # shared prefix -> shared windows
    DROP = int(C.VH.DROP)
    win_cycles = reads[DROP:DROP + MATCH_WINDOWS]
    t0h, t1h = win_cycles[0]['cycle_start_h'], win_cycles[-1]['cycle_end_h']
    print('    matching window = read windows %d..%d = %.3f .. %.3f h (%.3f h)'
          % (DROP, DROP + MATCH_WINDOWS - 1, t0h, t1h, t1h - t0h), flush=True)

    alpha0 = float(mk.tail.p['alpha_Int'][1])
    mech_k = C.mechanism(tk, sk, win_cycles)
    mech_r = C.mechanism(tk, sr, win_cycles)
    I_kept = sum(c['g1_dose'] for c in mech_k['per_cycle'])
    I_rem = sum(c['g1_dose'] for c in mech_r['per_cycle'])
    alpha_m = alpha0 * I_kept / I_rem
    print('    cycles = %d' % len(mech_k['per_cycle']), flush=True)
    print('    int g1 dt over window:  kept %.6f   removed %.6f   ratio %.6f'
          % (I_kept, I_rem, I_rem / I_kept), flush=True)
    print('    alpha_Int[1]:  %.6f -> matched %.6f  (x%.6f)'
          % (alpha0, alpha_m, alpha_m / alpha0), flush=True)
    seg_k = C.VH.segments(tk, g_kept, C.VH.GATE_THRESHOLD,
                          min_dose=C.VH.GATE_MIN_DOSE_H)
    seg_r = C.VH.segments(tk, g_rem, C.VH.GATE_THRESHOLD,
                          min_dose=C.VH.GATE_MIN_DOSE_H)
    print('    (in-pulse doses %.4f / %.4f h are NOT used as the gain -- they are '
          'integrals inside identified gate pulses only)'
          % (sum(s['dose'] for s in seg_k), sum(s['dose'] for s in seg_r)),
          flush=True)

    # ---------------- run the matched arm ----------------
    print('\n[1] MATCHED arm: clock == 1, alpha_Int[1] = %.6f' % alpha_m, flush=True)
    t_start = time.time()
    mm, _ = C.build('hzh', han, alpha_int2=alpha_m)
    recv = BypassReceiver(mm.tail.c)
    recv.p = dict(recv.p)                     # BypassReceiver re-reads the table
    seq = list(recv.p['alpha_Int'])
    seq[1] = float(alpha_m)
    recv.p['alpha_Int'] = tuple(seq)
    mm.tail = recv
    mm.copies_per_au = C.COPIES_PER_UM_PER_FL * recv.c.receiver_uM_per_au
    tm, ym = C.integrate(mm, hours=HOURS)
    print('    integrated in %.1f s' % (time.time() - t_start), flush=True)

    sm, vm = judge(tm, ym, mm, bypassed=True)
    rm = summarise(tm, ym, mm, sm, vm, 'matched_K0.30', K_CERT, True)
    print('    certified=%s counting=%s events=%s reads=%d seq=%s'
          % (rm['certified'], rm['counting'], rm['events'], rm['steady_reads'],
             rm['sequence']), flush=True)

    # ---------------- verify the matching really holds ----------------
    mech_m = C.mechanism(tm, sm, win_cycles)
    I_mat = sum(c['g1_dose'] for c in mech_m['per_cycle'])
    tot_kept = alpha0 * I_kept
    tot_rem = alpha0 * I_rem
    tot_mat = alpha_m * I_mat
    rel = abs(tot_mat - tot_kept) / abs(tot_kept)
    matched = rel < 1e-6
    print('\n[2] matching check over the same %d cycles:' % len(win_cycles), flush=True)
    print('    alpha*int g1 dt:  kept %.6f   removed %.6f   matched %.6f'
          % (tot_kept, tot_rem, tot_mat), flush=True)
    print('    relative gap matched-vs-kept = %.3e  -> %s'
          % (rel, 'MATCHED' if matched else 'NOT MATCHED'), flush=True)
    print('    (removed unmatched was %.2f %% higher than kept)'
          % (100 * (tot_rem / tot_kept - 1)), flush=True)

    # ---------------- bit2's OWN timing margins ----------------
    print('\n[3] timing margins -- global min is bit0-limited, so report bit2 too',
          flush=True)
    rows_margins = []
    # all three verdicts come from ONE judgment path so the comparison is fair
    vk = C.verdict(tk, sk)
    vr = C.verdict(tr, sr)
    for tag, v in (('kept', vk), ('removed', vr), ('matched', vm)):
        bm = v['bit_margins']
        rows_margins.append(dict(
            arm=tag,
            global_min_h=v['global_min_timing_margin_h'],
            bit0_setup=bm['S0']['min_setup_h'], bit0_hold=bm['S0']['min_hold_h'],
            bit1_setup=bm['S1']['min_setup_h'], bit1_hold=bm['S1']['min_hold_h'],
            bit2_setup=bm['S2']['min_setup_h'], bit2_hold=bm['S2']['min_hold_h']))
        print('    %-8s global %.4f | bit2 setup %s hold %s'
              % (tag, v['global_min_timing_margin_h'],
                 bm['S2']['min_setup_h'], bm['S2']['min_hold_h']), flush=True)

    # ---------------- save ----------------
    np.savez_compressed(out / 'traj_matched_K0.30.npz', time_h=tm, states=ym)
    cases = dict(kept=summarise(tk, yk, mk, sk, vk, 'kept_K0.30', K_CERT, False),
                 removed=summarise(tr, yr, mr, sr, vr, 'removed_K0.30', K_CERT, True),
                 matched=rm)
    dump(out / 'summary.json', dict(
        smoke=False,
        question=('does the clock gate act by TEMPORAL SCREENING or merely by '
                  'changing the total Int2 source input?'),
        balancing='option (c): equal total int u2 dt = int alpha_Int[1]*g1 dt over '
                  'identical boundaries (one post-burn-in mod-8 cycle = %d read '
                  'windows, 2 bit1->bit2 carries)' % MATCH_WINDOWS,
        match_window=dict(first_window=DROP, n_windows=MATCH_WINDOWS,
                          t_start_h=float(t0h), t_end_h=float(t1h)),
        alpha_Int1_default=alpha0, alpha_Int1_matched=float(alpha_m),
        integral_g1=dict(kept=I_kept, removed=I_rem, matched=I_mat),
        total_source=dict(kept=tot_kept, removed=tot_rem, matched_arm=tot_mat,
                          rel_gap_matched_vs_kept=float(rel),
                          matching_ok=bool(matched)),
        bit_margins=rows_margins,
        mechanism_per_cycle=dict(kept=mech_k['per_cycle'],
                                 removed=mech_r['per_cycle'],
                                 matched=mech_m['per_cycle']),
        cases=list(cases.values())))
    dump(out / 'source_hashes.json', C.source_hashes())
    dump(out / 'SHA256SUMS.json',
         {str(p.relative_to(out)): sha256(p)
          for p in sorted(out.rglob('*')) if p.is_file()})
    dump(out / 'status.json', dict(status='COMPLETED', matched=bool(matched),
                                  rel_gap=float(rel)))
    print('\nOUTPUT', out, flush=True)


if __name__ == '__main__':
    main()
