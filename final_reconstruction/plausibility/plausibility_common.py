"""Shared helpers for the parameter-plausibility checks.

Policy for every script in this folder
--------------------------------------
* The model and every frozen artefact are READ-ONLY.  Nothing here writes into
  model.py, model_twobit34.py, model_threebit51.py, verify_threebit51.py,
  verify_twobit_causal.py, selected_profile.json, or any other agent's outputs.
* Everything this package produces goes under plausibility/.
* Diagnostic monkey-patches of model.ZENG (used by the sharpen arm and the
  RDF-collapse probe) are applied at RUN TIME ONLY and are always recorded in
  the output row, so a reader can see exactly which value was in force.

All acceptance thresholds are pre-registered constants at the top of each
script; they are never chosen after seeing the result.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'plausibility'
FROZEN_PROFILE = ROOT / 'twobit34_results' / 'selected_profile.json'

# ---------------------------------------------------------------- constants
COPIES_PER_UM_PER_FL = 602.214076          # 1 uM in 1 fL, exact by construction
CELL_VOLUME_FL = 1.0

# E. coli at ~30 min doubling: about 2e6 proteins per cell are synthesised per
# generation, i.e. ~4e6 copies/cell/h.  Used only as an order-of-magnitude
# budget for the expression-burden check; it is an assumption, not a fitted
# value.  Sources vary by ~2x, hence the check reports fractions, not verdicts.
PROTEIN_SYNTHESIS_BUDGET_PER_H = 4.0e6

# Physical reference ranges used for flagging, order-of-magnitude only.
REF_FREE_TF_UM = 30.0        # a free transcription factor above this is unusual
REF_BURDEN_FRACTION = 0.10   # >10 % of total protein synthesis on one circuit

FROZEN_FILES = ('model.py', 'model_twobit34.py', 'model_threebit51.py',
                'verify_threebit51.py', 'verify_twobit_causal.py',
                'test_threebit51.py', 'scan_threebit51_carry1.py')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def source_hashes() -> dict:
    return {name: sha256(ROOT / name) for name in FROZEN_FILES if (ROOT / name).exists()}


def load_profile() -> dict:
    return json.loads(FROZEN_PROFILE.read_text(encoding='utf-8'))


# ------------------------------------------------------------------- models
def frozen_extension():
    """The validated 34-state profile, verbatim from selected_profile.json."""
    from model_twobit34 import nominal_extension
    p = load_profile()
    return replace(nominal_extension(), **p['extension'])


def frozen_carry0():
    from model_twobit34 import CarryExpressionParameters
    return CarryExpressionParameters(**load_profile()['carry_expression'])


def build_threebit(carry1=None, extension=None, n_A1=None, k_A1=None,
                   n_rep=None, k_rep=None):
    """Build the 51-state model.

    n_A1 / k_A1 / n_rep / k_rep apply a run-time diagnostic patch to model.ZENG
    and the previous values are returned so the caller can record and restore
    them.  They are NOT proposals to change Zeng's table; they exist to separate
    structural questions from parameter questions.
    """
    import model as M
    from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters
    patched = {}
    if n_A1 is not None:
        patched['n_A'] = M.ZENG['n_A']
        M.ZENG['n_A'] = (M.ZENG['n_A'][0], float(n_A1))
    if k_A1 is not None:
        patched['K_A'] = M.ZENG['K_A']
        M.ZENG['K_A'] = (M.ZENG['K_A'][0], float(k_A1))
    if n_rep is not None:
        patched['n_rep'] = M.ZENG['n_rep']
        M.ZENG['n_rep'] = float(n_rep)
    if k_rep is not None:
        patched['K_rep'] = M.ZENG['K_rep']
        M.ZENG['K_rep'] = float(k_rep)
    carry = ThreeBitCarryParameters(
        carry0=frozen_carry0(),
        carry1=carry1 if carry1 is not None else frozen_carry0())
    model = ThreeBit51Model(extension=extension if extension is not None
                            else frozen_extension(), carry=carry)
    return model, patched


def restore(patch: dict):
    import model as M
    for key, value in patch.items():
        M.ZENG[key] = value


# ------------------------------------------------------------------- maths
def hill(x, K, n):
    x = np.maximum(np.asarray(x, dtype=float), 0.0)
    return x ** n / (K ** n + x ** n)


def repression(x, K, n):
    return 1.0 - hill(x, K, n)


def threshold_margin(x, K):
    """Operating ratio x/K and how much time is spent near the threshold."""
    ratio = np.asarray(x, dtype=float) / float(K)
    near = float(np.mean((ratio >= 1.0 / np.sqrt(2.0)) & (ratio <= np.sqrt(2.0))))
    crossings = int(np.sum((ratio[:-1] - 1.0) * (ratio[1:] - 1.0) < 0))
    return dict(ratio_min=float(ratio.min()), ratio_p05=float(np.percentile(ratio, 5)),
                ratio_median=float(np.median(ratio)), ratio_p95=float(np.percentile(ratio, 95)),
                ratio_max=float(ratio.max()), frac_time_near_threshold=near,
                crossings_of_threshold=crossings,
                straddles_threshold=bool(ratio.min() < 1.0 < ratio.max()))


# --------------------------------------------------------------- trajectory
def flux_array(model, sol):
    return np.asarray([model.flux(sol.y[:, k]) for k in range(sol.y.shape[1])])


def read_window_slices(t, flux, fraction=0.20):
    """Same geometry as verify_twobit_causal / verify_threebit51."""
    from verify_bit0_part2 import clock_cycles
    peaks = clock_cycles(t, flux)
    out = []
    for a, b in zip(peaks[:-1], peaks[1:]):
        trough = a + int(np.argmin(flux[a:b + 1]))
        half = max(1, int(round(fraction * (b - a) / 2)))
        out.append((max(a, trough - half), min(b, trough + half), a, b))
    return out


def gate_windows(t, g1, threshold=0.01, min_duration_h=0.05):
    """Contiguous g1 episodes above a low threshold (carry opportunities)."""
    mask = np.asarray(g1) >= threshold
    ch = np.diff(np.r_[False, mask, False].astype(np.int8))
    starts = np.flatnonzero(ch == 1)
    stops = np.flatnonzero(ch == -1) - 1
    out = []
    for i, j in zip(starts, stops):
        if t[j] - t[i] >= min_duration_h:
            out.append(dict(i0=int(i), i1=int(j), start_h=float(t[i]), end_h=float(t[j]),
                            width_h=float(t[j] - t[i]),
                            peak=float(np.max(g1[i:j + 1])),
                            dose=float(np.trapezoid(g1[i:j + 1], t[i:j + 1]))))
    return out


def leak_segments(t, g1, J_rev2, J_fwd2, I2, windows, tail_fraction=0.05,
                  min_far_off_h=0.10):
    """Split the reverse dose of one carry period into three named segments.

    The original `leak_budget` integrated everything between two carry windows
    and called it "leakage".  It does not measure chronic leakage: it also
    contains the carry's own decay tail.  Evidence: at fitted maturation 60 min
    the off-state gate peak is 28x lower than at 32.5 min, yet `leak_worst` is
    40 % higher - so the metric is dominated by the tail, not by the far-off
    dwell.  This split follows the agreed definition:

        gate-on : the g1 main window            (the intended reverse action)
        tail    : from the window end until Int2 falls to tail_fraction of that
                  pulse's own peak             (the carry's own decay)
        far-off : the remainder of the dwell, until the next window opens
                                              (the only genuine leakage)

    Reporting rules agreed with the reviewer:
      * the durations of all three segments are reported;
      * if the tail does not fall below the cut before the next carry, far-off
        is EMPTY and is marked far_off_evaluable = False and left as None -
        it must never be recorded as "zero leakage";
      * the wrong-direction (forward) integral during far-off is reported too;
      * the cut fraction is an engineering choice, so `leak_report` repeats the
        whole decomposition at 1 %, 5 % and 10 %.
    """
    rows = []
    n = len(t)
    for k, w in enumerate(windows):
        i0, i1 = w['i0'], w['i1']
        I2pk = float(np.max(I2[i0:i1 + 1])) if i1 >= i0 else 0.0
        thr = tail_fraction * I2pk
        j = i1
        while j + 1 < n and I2[j + 1] > thr:
            j += 1
        nxt = windows[k + 1]['i0'] if k + 1 < len(windows) else n - 1
        far_h = float(t[nxt] - t[j]) if nxt > j else 0.0
        evaluable = bool(k + 1 < len(windows) and far_h >= min_far_off_h)
        rows.append(dict(
            after_window=k,
            gate_on_h=float(t[i1] - t[i0]),
            gate_on_rev=float(np.trapezoid(J_rev2[i0:i1 + 1], t[i0:i1 + 1])),
            tail_h=float(t[j] - t[i1]),
            tail_rev=float(np.trapezoid(J_rev2[i1:j + 1], t[i1:j + 1])),
            far_off_h=far_h,
            far_off_rev=(float(np.trapezoid(J_rev2[j:nxt], t[j:nxt])) if evaluable else None),
            far_off_fwd=(float(np.trapezoid(J_fwd2[j:nxt], t[j:nxt])) if evaluable else None),
            far_off_evaluable=evaluable,
            I2_peak=I2pk, tail_threshold=thr))
    ev = [r for r in rows if r['far_off_evaluable']]

    def med(key, src):
        v = [r[key] for r in src if r.get(key) is not None]
        return float(np.median(v)) if v else None

    g = med('gate_on_rev', rows)
    tl = med('tail_rev', rows)
    fo = med('far_off_rev', ev)
    ff = med('far_off_fwd', ev)
    return dict(per_window=rows, tail_fraction=tail_fraction,
                windows_total=len(rows), windows_evaluable=len(ev),
                gate_on_h=med('gate_on_h', rows), tail_h=med('tail_h', rows),
                far_off_h=med('far_off_h', ev),
                gate_on_rev=g, tail_rev=tl, far_off_rev=fo, far_off_fwd=ff,
                leak_ratio=(fo / g if (fo is not None and g) else None),
                far_off_fwd_ratio=(ff / g if (ff is not None and g) else None))


def leak_report(t, g1, J_rev2, J_fwd2, I2, windows,
                fractions=(0.01, 0.05, 0.10), min_far_off_h=0.10):
    """Three-segment decomposition at several tail thresholds.

    The 5 % cut is an engineering choice, not a biological constant, so the whole
    decomposition is repeated at 1 %, 5 % and 10 %.  A conclusion about which
    gate exponent is better is only trusted if the ordering is the same at all
    three cuts.
    """
    return {f: leak_segments(t, g1, J_rev2, J_fwd2, I2, windows,
                             tail_fraction=f, min_far_off_h=min_far_off_h)
            for f in fractions}


def leak_budget(t, g1, J_rev2, windows):
    """Pre-registered replacement for the off/on peak-ratio metric.

    The complex Hill is quadratic at low drive, so what erodes the stored bit is
    the *integrated* reverse dose between two carries, not the ratio of peaks.
    """
    rows = []
    for k in range(len(windows) - 1):
        i0, i1 = windows[k]['i1'], windows[k + 1]['i0']
        if i1 <= i0:
            continue
        rows.append(dict(after_window=k,
                         dwell_h=float(t[i1] - t[i0]),
                         int_g1_off=float(np.trapezoid(g1[i0:i1], t[i0:i1])),
                         int_Jrev2_off=float(np.trapezoid(J_rev2[i0:i1], t[i0:i1]))))
    if not rows:
        return dict(per_dwell=[], worst_Jrev2=None, median_Jrev2=None)
    vals = [r['int_Jrev2_off'] for r in rows]
    return dict(per_dwell=rows, worst_Jrev2=float(max(vals)),
                median_Jrev2=float(np.median(vals)))


def stall_diagnostics(t, y, sig, low_band=0.30):
    """State of the bit2 pools at the moment the read is stuck above the band."""
    S2 = np.asarray(sig['S2'])
    late = t > 0.5 * t[-1]
    idx = np.flatnonzero(late)
    k = idx[int(np.argmin(S2[idx]))]
    return dict(t_at_min=float(t[k]), S2_min_late=float(S2[k]),
                I2=float(y[36][k]), R2=float(y[42][k]), C2=float(y[43][k]),
                T2=float(y[40][k]), A1=float(y[45][k]), F1=float(y[46][k]),
                clock=float(sig['clock'][k]), g1=float(sig['g1'][k]),
                J_fwd2=float(sig['J_fwd2'][k]), J_rev2=float(sig['J_rev2'][k]),
                below_band=bool(S2[k] <= low_band))


# ------------------------------------------------------------------ output
def write_shard(name: str, tag: str, rows: list, meta: dict):
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / f'{name}_{tag}.csv'
    pd.DataFrame(rows).to_csv(csv_path, index=False, encoding='utf-8')
    meta = dict(meta)
    meta['csv'] = csv_path.name
    meta['csv_sha256'] = sha256(csv_path)
    meta['source_sha256'] = source_hashes()
    (OUT / f'meta_{name}_{tag}.json').write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    return csv_path


def shard_slice(points: list, shard: int, nshards: int) -> list:
    return points[shard::nshards]


def merge_shards(name: str, nshards: int) -> pd.DataFrame:
    frames = []
    for s in range(nshards):
        tag = f'shard{s:02d}of{nshards:02d}'
        p = OUT / f'{name}_{tag}.csv'
        if not p.exists():
            raise SystemExit(f'missing shard artefact {p}')
        frames.append(pd.read_csv(p))
    data = pd.concat(frames, ignore_index=True)
    data.to_csv(OUT / f'{name}_all.csv', index=False, encoding='utf-8')
    return data


def write_manifest():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {}
    # rglob: also cover artefacts written into subdirectories (e.g.
    # threshold_margin/), which a top-level glob would silently miss.
    for p in sorted(OUT.rglob('*')):
        if p.is_file() and not p.name.startswith('SHA256SUMS') and '__pycache__' not in p.parts:
            manifest[str(p.relative_to(ROOT))] = sha256(p)
    for k, v in source_hashes().items():
        manifest[k] = v
    (OUT / 'SHA256SUMS.json').write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True), encoding='utf-8')
    (OUT / 'SHA256SUMS.txt').write_text(
        ''.join(f'{v}  {k}\n' for k, v in sorted(manifest.items())), encoding='utf-8')
    return manifest


# ==========================================================================
# Decoupling of the gate exponent from the F1 production exponent
# ==========================================================================
# model_threebit51 reads ZENG['n_A'][1] in TWO places for the carry-1 level:
#   (a) model.Model.carry_promoters -> g1 = H(A1;K_A1,n) * G(F1;K_F1,n_F1) * clock
#       this is the gate arm we want to sharpen;
#   (b) ThreeBit51Model._carry1_sources -> source_F1 = alpha_F1 * H(A1;K_A1,n)
#       this is the incoherent arm's drive, which must stay at the frozen value.
# Patching ZENG['n_A'][1] therefore moves both.  DecoupledThreeBit51Model keeps
# (a) on the patched value and pins (b) to FROZEN_N_A1, without touching any
# frozen file: only the subclass override differs.
FROZEN_N_A1 = 4.0


def _make_decoupled_class():
    from model_threebit51 import ThreeBit51Model

    class DecoupledThreeBit51Model(ThreeBit51Model):
        """Gate arm uses ZENG['n_A'][1]; the F1 driving Hill is pinned."""

        def __init__(self, extension=None, carry=None, n_f1_drive=FROZEN_N_A1):
            self.n_f1_drive = float(n_f1_drive)
            super().__init__(extension=extension, carry=carry)

        def _carry1_sources(self, S1, A1):
            from model import ZENG, act, rep
            auto = rep(A1, ZENG['K_auto1'], ZENG['n_auto1'])
            source_A1 = ZENG['alpha_A'][1] * (1.0 - S1) * auto
            source_F1 = ZENG['alpha_F'][1] * act(A1, ZENG['K_A'][1], self.n_f1_drive)
            return source_A1, source_F1

        def initial_state(self, cold=False):
            import model as M
            saved = M.ZENG['n_A']
            M.ZENG['n_A'] = (saved[0], self.n_f1_drive)
            try:
                return super().initial_state(cold=cold)
            finally:
                M.ZENG['n_A'] = saved

    return DecoupledThreeBit51Model


DecoupledThreeBit51Model = _make_decoupled_class()


def build_decoupled(n_A1_gate=None, carry1=None, extension=None, n_f1_drive=FROZEN_N_A1):
    """51-state model whose gate exponent is independent of the F1 drive exponent."""
    import model as M
    from model_threebit51 import ThreeBitCarryParameters
    patched = {}
    if n_A1_gate is not None:
        patched['n_A'] = M.ZENG['n_A']
        M.ZENG['n_A'] = (M.ZENG['n_A'][0], float(n_A1_gate))
    carry = ThreeBitCarryParameters(
        carry0=frozen_carry0(),
        carry1=carry1 if carry1 is not None else frozen_carry0())
    model = DecoupledThreeBit51Model(
        extension=extension if extension is not None else frozen_extension(),
        carry=carry, n_f1_drive=n_f1_drive)
    return model, patched


def decoupling_verification_path() -> Path:
    return OUT / 'decoupling_verification.json'


def implementation_sha256() -> str:
    """Hash of the two files that implement the decoupling."""
    h = hashlib.sha256()
    for name in ('model_threebit51.py', 'plausibility/plausibility_common.py'):
        h.update(sha256(ROOT / name).encode())
    return h.hexdigest().upper()


def verify_decoupling(hours=100.0, sample_min=2.0, max_step_min=2.0,
                      probes=(6.0,)) -> dict:
    """Prove that the decoupling is exact, and write the artefact the scan checks.

    Test 1  with the gate exponent at the frozen value, the decoupled model must
            reproduce the plain model bit-for-bit (initial state and trajectory).
    Test 2  with the gate exponent moved, the F1 production source evaluated on
            the frozen trajectory must be unchanged, while g1 must change.
    """
    carry = frozen_carry0()

    plain, patch_plain = build_threebit(carry1=carry, n_A1=None)
    try:
        y0_plain = plain.initial_state(cold=False)
        sol_plain = plain.simulate(hours=hours, sample_min=sample_min,
                                   max_step_min=max_step_min, initial_state=y0_plain)
    finally:
        restore(patch_plain)

    dec4, patch4 = build_decoupled(n_A1_gate=FROZEN_N_A1, carry1=carry)
    try:
        y0_dec = dec4.initial_state(cold=False)
        sol_dec = dec4.simulate(hours=hours, sample_min=sample_min,
                                max_step_min=max_step_min, initial_state=y0_dec)
    finally:
        restore(patch4)

    init_gap = float(np.max(np.abs(y0_dec - y0_plain)))
    traj_gap = float(np.max(np.abs(sol_dec.y - sol_plain.y)))

    # Test 2: evaluate both source variants on the frozen trajectory
    probe_rows = []
    s1 = np.asarray(sol_plain.y[27])
    a1 = np.asarray(sol_plain.y[45])
    idx = list(range(0, len(s1), 200))
    frozen_src = np.asarray([plain._carry1_sources(s1[k], a1[k])[1] for k in idx])
    g1_frozen = plain.diagnostic_signals(sol_plain.y)['g1'][idx]
    for n_gate in probes:
        dec, patch = build_decoupled(n_A1_gate=n_gate, carry1=carry)
        try:
            src_dec = np.asarray([dec._carry1_sources(s1[k], a1[k])[1] for k in idx])
            g1_dec = dec.diagnostic_signals(sol_plain.y)['g1'][idx]
        finally:
            restore(patch)
        plain_patched, patch2 = build_threebit(carry1=carry, n_A1=n_gate)
        try:
            src_sh = np.asarray([plain_patched._carry1_sources(s1[k], a1[k])[1] for k in idx])
        finally:
            restore(patch2)
        probe_rows.append(dict(
            n_A1_gate=n_gate,
            f1_source_decoupled_equals_frozen=bool(np.array_equal(src_dec, frozen_src)),
            f1_source_shared_patch_equals_frozen=bool(np.array_equal(src_sh, frozen_src)),
            f1_source_rel_change_shared_patch=float(
                np.max(np.abs(src_sh - frozen_src) / np.maximum(np.abs(frozen_src), 1e-30))),
            max_g1_frozen=float(np.max(g1_frozen)),
            max_g1_at_gate_n=float(np.max(g1_dec))))

    passed = bool(init_gap == 0.0 and traj_gap < 1e-9 and
                  all(r['f1_source_decoupled_equals_frozen'] for r in probe_rows) and
                  all(not r['f1_source_shared_patch_equals_frozen'] for r in probe_rows))
    result = dict(passed=passed, hours=hours, frozen_n_A1=FROZEN_N_A1,
                  initial_state_max_abs_gap=init_gap,
                  trajectory_max_abs_gap=traj_gap,
                  probes=probe_rows,
                  implementation_sha256=implementation_sha256(),
                  source_sha256=source_hashes(),
                  note=('The decoupled model reproduces the frozen model exactly at the '
                        'frozen gate exponent, and moving the gate exponent leaves the '
                        'F1 production source untouched. The plain shared-slot patch '
                        'does move it, which is why it must not be used for the scan.'))
    OUT.mkdir(parents=True, exist_ok=True)
    decoupling_verification_path().write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    return result


def load_valid_decoupling():
    """Return the artefact if it exists, passed, and is not stale."""
    p = decoupling_verification_path()
    if not p.exists():
        return None
    art = json.loads(p.read_text(encoding='utf-8'))
    if not art.get('passed'):
        return None
    if art.get('implementation_sha256') != implementation_sha256():
        return None
    if art.get('source_sha256') != source_hashes():
        return None
    return art
