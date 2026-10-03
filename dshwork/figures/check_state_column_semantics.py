r"""Column -> model-state SEMANTICS checker for trajectory CSVs   (verification layer)

Why this file exists
--------------------
The existing figure self-checks proved that a plotted curve equals the CSV column it
was read from.  That is only COLUMN FIDELITY.  It cannot catch the case where the CSV
column itself holds the wrong model state -- which is exactly what happened to
`data/n4_matched_600h.csv`, whose columns were written as

    S0 = sol.y[44]   S1 = sol.y[45]   S2 = sol.y[46]

while the 51-state layout is S0 = y[16], S1 = y[27], S2 = y[44] (so the file actually
held S2, A1, F1 under the names S0, S1, S2).  Figure 4-3 Panel A therefore plotted bit
2's F1 protein in its n=4 panel and bit 2's real S2 in its n=6 panel: two different
physical quantities under one axis label, and a "S2 peaks at 2.25" note that was
really F1 -- a concentration in a.u., not a bounded fraction.

This tool checks the missing layer.  For each trajectory file it verifies, independently
of how the file was produced:

  1. the model's own `state_names` put the S states where the project says they are
     (b0_S = 16, b1_S = 27, b2_S = 44, A1 = 45, F1 = 46);
  2. every column that claims to be a DNA-conformation fraction (`*S`, or the derived
     `S0`/`S1`/`S2` in the wiki artefact) really lies in [0, 1] -- a concentration in
     a.u. blows straight through this, which is the cheapest possible tripwire;
  3. the claimed bit-2 fraction REPLAYS the recorded read-window bit-2 labels under
     `verify_twobit_causal._window_label` (band 0.30/0.70, occupancy 0.80).  This is the
     decisive check: it ties the plotted column to the same signal the certification
     read, so a wrong column cannot pass;
  4. when `A1`, `F1`, `g1`, `clock_gate` are all present, the gate identity
     g1 = act(A1, K_A[1], n_gate) * rep(F1, K_F[1], n_F[1]) * clock_gate
     is reproduced from the columns themselves.

Nothing here re-integrates the model: everything is read back off disk.

Run
---
    & 'D:\aconade\python.exe' .\check_state_column_semantics.py            # both files
    & 'D:\aconade\python.exe' .\check_state_column_semantics.py --csv X.csv --json X.json

Exit code 0 = every check passed; 1 = at least one failed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FR = HERE.parents[1] / 'final_reconstruction'
DATA = HERE / 'data'
N6DIR = FR / 'threebit51_results' / 'wiki_n6_20260924'

BAND_LOW, BAND_HIGH, MIN_OCCUPANCY = 0.30, 0.70, 0.80
EXPECTED_INDEX = {'b0_S': 16, 'b1_S': 27, 'b2_S': 44, 'A1': 45, 'F1': 46}
# column names seen in the wild that CLAIM to be a bit's LR fraction.  Whether the
# file lives up to the claim is what this tool tests.
FRACTION_ALIASES = {'S0': 'b0_S', 'S1': 'b1_S', 'S2': 'b2_S'}
FRACTION_NAMES = ('b0_S', 'b1_S', 'b2_S')


def window_label(sig: np.ndarray) -> str:
    """Verbatim re-implementation of verify_twobit_causal._window_label."""
    hi = float(np.mean(sig >= BAND_HIGH))
    lo = float(np.mean(sig <= BAND_LOW))
    return '1' if hi >= MIN_OCCUPANCY else ('0' if lo >= MIN_OCCUPANCY else 'x')


def act(x, K, n):
    x = np.maximum(np.asarray(x, dtype=float), 0.0)
    return x ** n / (K ** n + x ** n)


def rep(x, K, n):
    return 1.0 - act(x, K, n)


class Checker:
    def __init__(self, name: str):
        self.name = name
        self.fails: list[str] = []
        self.n_pass = 0

    def check(self, cond: bool, msg: str) -> bool:
        if cond:
            self.n_pass += 1
            print(f'  PASS  {msg}')
        else:
            self.fails.append(msg)
            print(f'  FAIL  {msg}')
        return bool(cond)


def state_index_map(n_A1_gate: float = 6.0) -> tuple[dict, object]:
    """Read the index map off the model itself, not off a hardcoded table.

    `n_A1_gate` matters because the gate identity in section 4 depends on it: pass the
    exponent the file under test was actually produced with.
    """
    sys.path.insert(0, str(FR))
    sys.path.insert(0, str(FR / 'plausibility'))
    from plausibility_common import build_threebit                     # noqa: E402
    model, _patched = build_threebit(n_A1_gate=n_A1_gate)
    idx = {n: i for i, n in enumerate(model.state_names)}
    return idx, model


def resolve_columns(cols: list[str]) -> dict[str, str]:
    """canonical state name -> the column in this file that claims to hold it.

    A canonical name always wins; an alias (`S2` -> `b2_S`) is only used when the
    canonical column is absent.
    """
    resolved: dict[str, str] = {}
    for name in FRACTION_NAMES + ('A1', 'F1'):
        if name in cols:
            resolved[name] = name
    for alias, canonical in FRACTION_ALIASES.items():
        if canonical not in resolved and alias in cols:
            resolved[canonical] = alias
    return resolved


def windows_from_json(json_path: Path) -> list[dict]:
    meta = json.loads(json_path.read_text(encoding='utf-8'))
    out = []
    for w in meta.get('read_windows') or []:
        raw = w.get('raw_labels')
        bit2 = raw[2] if raw else ('x' if w.get('bit2') is None else str(w['bit2']))
        out.append(dict(window_start_h=w['window_start_h'],
                        window_end_h=w['window_end_h'], bit2=str(bit2)))
    return out


def windows_from_csv(csv_path: Path) -> list[dict]:
    """The wiki artefact records the labels as `x` or a digit in a `bit2` column."""
    df = pd.read_csv(csv_path)
    return [dict(window_start_h=float(r.window_start_h),
                 window_end_h=float(r.window_end_h),
                 bit2=str(r.bit2)) for r in df.itertuples()]


def check_file(csv_path: Path, windows: list[dict], model, idx: dict,
               n_gate: float | None = None) -> int:
    c = Checker(csv_path.name)
    print(f'=== {csv_path.name}')
    df = pd.read_csv(csv_path)
    cols = list(df.columns)
    print(f'  {len(cols)} columns, {len(df)} rows')
    if len(cols) <= 12:
        print(f'  columns: {cols}')

    print("-- 1. the model's own state_names place the S states where the project says")
    for name, want in EXPECTED_INDEX.items():
        c.check(idx.get(name) == want,
                f'model.state_names[{name!r}] = {idx.get(name)} (expected {want})')

    resolved = resolve_columns(cols)
    print(f'  resolved claims: {resolved}')

    print('-- 2. every column claiming to be a bounded LR fraction lies in [0, 1]')
    seen = 0
    for canonical in FRACTION_NAMES:
        col = resolved.get(canonical)
        if col is None:
            continue
        seen += 1
        v = df[col].to_numpy(dtype=float)
        lo, hi = float(np.nanmin(v)), float(np.nanmax(v))
        c.check(lo >= -1e-12 and hi <= 1.0 + 1e-9,
                f'column {col!r} (claimed {canonical}) is within [0, 1]: '
                f'min {lo:.4f}, max {hi:.4f}')
        if hi > 1.0 + 1e-9:
            clue = [o for o in ('A1', 'F1', 'b2_I', 'b2_R')
                    if o in cols and np.allclose(df[o].to_numpy(dtype=float), v,
                                                 rtol=0, atol=1e-12)]
            print(f'        NOTE: {col!r} exceeds 1, so it is NOT an LR fraction.'
                  + (f' It is bit-identical to {clue} -- a species in a.u.' if clue
                     else ' It holds a species in a.u., not a fraction.'))
    c.check(seen == 3, f'all three bit fractions are resolvable (found {seen} of 3)')

    print('-- 2b. when both an alias and the canonical column exist they must agree')
    for alias, canonical in FRACTION_ALIASES.items():
        if alias in cols and canonical in cols:
            same = np.array_equal(df[alias].to_numpy(dtype=float),
                                  df[canonical].to_numpy(dtype=float))
            gap = float(np.abs(df[alias].to_numpy(dtype=float)
                               - df[canonical].to_numpy(dtype=float)).max())
            c.check(same, f'{alias!r} == {canonical!r} entry-wise (max|d| = {gap:.3g})')

    print('-- 3. the claimed b2 fraction replays the recorded read-window bit-2 labels')
    if not windows:
        print('  SKIP  no read windows were supplied for this file')
    else:
        col = resolved.get('b2_S')
        if col is None:
            c.check(False, 'read windows exist but no column claims to be b2_S')
        else:
            t = df['time_h'].to_numpy(dtype=float)
            sig = df[col].to_numpy(dtype=float)
            mism = []
            for w in windows:
                a = int(np.searchsorted(t, w['window_start_h']))
                b = int(np.searchsorted(t, w['window_end_h']))
                got = window_label(sig[a:b + 1])
                if got != w['bit2']:
                    mism.append((round(float(w['window_start_h']), 2), w['bit2'], got))
            c.check(not mism,
                    f'column {col!r} reproduces all {len(windows)} recorded bit-2 labels '
                    f'(mismatches: {len(mism)}; first few {mism[:4]})')

    print('-- 4. the gate identity holds for the columns as written')
    need = ('A1', 'F1', 'g1', 'clock_gate')
    if all(k in cols for k in need):
        from model import ZENG                                          # noqa: E402
        n_gate = float(model.n_A1_gate_effective if n_gate is None else n_gate)
        a1 = df['A1'].to_numpy(dtype=float)
        f1 = df['F1'].to_numpy(dtype=float)
        ck = df['clock_gate'].to_numpy(dtype=float)
        g1 = act(a1, ZENG['K_A'][1], n_gate) * rep(f1, ZENG['K_F'][1], ZENG['n_F'][1]) * ck
        gap = float(np.abs(g1 - df['g1'].to_numpy(dtype=float)).max())
        c.check(gap < 1e-12,
                f'g1 = act(A1, K_A[1], n_gate={n_gate:g}) * rep(F1, K_F[1], n_F[1]) '
                f'* clock_gate reproduces the g1 column (max|d| = {gap:.3g})')
        # the swapped reading must NOT reproduce it, otherwise the check is vacuous
        swapped = (act(f1, ZENG['K_A'][1], n_gate)
                   * rep(a1, ZENG['K_F'][1], ZENG['n_F'][1]) * ck)
        swap_gap = float(np.abs(swapped - df['g1'].to_numpy(dtype=float)).max())
        c.check(swap_gap > 1e-6,
                f'the swapped reading (F1 as A1, A1 as F1) does NOT reproduce g1 '
                f'(max|d| = {swap_gap:.3g}) -- so the identity above is discriminating')
    else:
        missing = [k for k in need if k not in cols]
        print(f'  SKIP  needs {need}; missing {missing}')

    print(f'  -> {c.n_pass} passed, {len(c.fails)} failed')
    return len(c.fails)


def gate_from_json(json_path: Path, default: float = 6.0) -> float:
    """The gate exponent the file was actually produced with, from its provenance."""
    try:
        meta = json.loads(json_path.read_text(encoding='utf-8'))
        return float(meta['gate_exponent']['realised_by_model'])
    except Exception:                                                   # noqa: BLE001
        return default


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:                                                   # noqa: BLE001
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', type=Path, default=None)
    ap.add_argument('--json', type=Path, default=None)
    ap.add_argument('--windows-csv', type=Path, default=None)
    ap.add_argument('--gate', type=float, default=None,
                    help='gate exponent the file was produced with; inferred from the '
                         'provenance JSON when omitted')
    args = ap.parse_args()

    if args.csv is not None:
        if args.json is not None and args.json.is_file():
            windows = windows_from_json(args.json)
            gate = args.gate if args.gate is not None else gate_from_json(args.json)
        elif args.windows_csv is not None and args.windows_csv.is_file():
            windows = windows_from_csv(args.windows_csv)
            gate = 6.0 if args.gate is None else args.gate
        else:
            windows = []
            gate = 6.0 if args.gate is None else args.gate
        targets = [(args.csv, windows, gate)]
    else:
        n4json = DATA / 'n4_matched_600h.json'
        targets = [
            (DATA / 'n4_matched_600h.csv',
             windows_from_json(n4json) if n4json.is_file() else [],
             gate_from_json(n4json) if n4json.is_file() else 4.0),
            (N6DIR / 'trajectories.csv',
             windows_from_csv(N6DIR / 'read_windows.csv')
             if (N6DIR / 'read_windows.csv').is_file() else [], 6.0),
        ]

    total = 0
    for csv_path, windows, gate in targets:
        if not csv_path.is_file():
            print(f'=== {csv_path.name}\n  SKIP  file not found: {csv_path}\n')
            continue
        idx, model = state_index_map(gate)
        print(f'model: {FR / "model_threebit51.py"}  (n_A1_gate = {gate:g})')
        print(f'state_names[b0_S/b1_S/b2_S/A1/F1] = '
              f'{[idx[k] for k in ("b0_S", "b1_S", "b2_S", "A1", "F1")]}')
        total += check_file(csv_path, windows, model, idx, n_gate=gate)
        print()
    print(f'TOTAL FAILURES: {total}')
    return 1 if total else 0


if __name__ == '__main__':
    sys.exit(main())
