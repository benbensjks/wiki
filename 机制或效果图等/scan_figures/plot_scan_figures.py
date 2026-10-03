"""Figures from the EXISTING scans: 2-bit single-parameter, 3-bit single-parameter,
and the 2-bit / 3-bit robustness analyses.

Nothing here simulates anything.  Every panel reads a completed scan artefact and
plots it; the input list with SHA256 and row counts is written next to the figures
as `figure_data_manifest.json`, so a reader can check exactly which run each panel
came from.

NON-NEGOTIABLE LABELLING RULES (the project's own discipline)
------------------------------------------------------------
1. A panel never mixes working points without saying so.  Several 51-state scans
   predate the `n_A1_gate` constructor parameter and were therefore run at the
   LEGACY shared-slot exponent (ZENG['n_A'][1] = 4.0), not at the frozen
   `n_A1_gate = 6.0`.  Those panels are titled "LEGACY" and state the value.
2. Leak is always the paired R/F super-period `L_symmetric` (H4), on the 2-minute
   output grid, and every panel says which tail cut it used.
3. Certified/uncertified is the frozen verifier's verdict, not a re-definition.
   Where a verdict is driven by the contrast arm alone, that is stated.
4. Figures use English labels (the project's existing figures do, and no CJK font
   is configured); the Chinese captions live in `README_图注与参数.md`.

    python plot_scan_figures.py            # writes all four figures
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                  # noqa: E402
import numpy as np                                               # noqa: E402
import pandas as pd                                              # noqa: E402

HERE = Path(__file__).resolve().parent
FR = HERE.parent                       # final_reconstruction/
OUT = HERE

BIT2 = FR / 'twobit34_results'
BIT3 = FR / 'threebit51_results'
PLAUS = FR / 'plausibility'

# Which artefacts are the figures allowed to read?  Explicit list; a missing file
# is a hard error, and anything listed as invalidated in
# plausibility/INVALIDATED_ARTIFACTS.json is refused outright (see _guard()).
SRC_BIT2_ROBUST = BIT2 / 'robustness' / 'robustness_all.csv'
SRC_BIT2_STRICT = BIT2 / 'strict_tolerance' / 'strict_tolerance_all.csv'
SRC_BIT2_PERTURB = {g: BIT2 / 'perturbations' / f'perturbations_{g}_all.csv'
                    for g in ('s', 'rdf', 'affl', 'int1')}
SRC_BIT3_CARRY1 = BIT3 / 'carry1_scan_dsh' / 'carry1_all.csv'
SRC_BIT3_SPLIT = BIT3 / 'carry1_split_scan' / 'split_all.csv'
SRC_BIT3_POSTFREEZE = BIT3 / 'postfreeze_peak_perturbation' / 'postfreeze_all.csv'
SRC_NA1_MAT = PLAUS / 'nA1_vs_maturation_all.csv'
SRC_PREAUDIT = PLAUS / 'preaudit_stochasticity_all.csv'
SRC_PREAUDIT_VERDICT = PLAUS / 'preaudit_stochasticity_verdict.json'
SRC_INVENTORY = PLAUS / 'preaudit_stochasticity_inventory.json'
SRC_PAIRING = PLAUS / 'carry_pairing_grid_verdict.json'
SRC_INVALIDATED = PLAUS / 'INVALIDATED_ARTIFACTS.json'

# the frozen working point, quoted in panel titles
WP = dict(n_A1_gate=6.0, mRNA_min=2.0, mat_min=32.5, uM_per_au=5.75)
C_OK, C_BAD = '#2b8a3e', '#c92a2a'
C_TXT = '0.25'

# The pool-perturbation outcome vocabulary, as ONE shared table: the figure and
# selfcheck_scan_figures.py both read it, so the check cannot drift from the plot.
# Note the semantics that a naive "recovered / not recovered" split would destroy:
# every `stable_phase_shift_*` row is certified=True - the counter still counts, the
# run is merely offset by N read windows.
PERTURB_CATS = (('recovered_same_phase', '#2b8a3e', 'same phase'),
                ('stable_phase_shift_1', '#74b816', 'phase +1'),
                ('stable_phase_shift_2', '#f59f00', 'phase +2'),
                ('stable_phase_shift_3', '#e8590c', 'phase +3'),
                ('lost_counting', '#c92a2a', 'lost counting'))

USED: list = []


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest().upper()


def _invalidated_basenames() -> set:
    doc = json.loads(SRC_INVALIDATED.read_text(encoding='utf-8'))
    names = set()
    for group in doc.get('invalidated_artifacts', []):
        for p in group.get('paths', []):
            names.add(Path(p).name)
    return names


INVALIDATED = _invalidated_basenames()


def _guard(path: Path):
    if not path.is_file():
        raise SystemExit(f'FATAL: required scan artefact is missing: {path}')
    if path.name in INVALIDATED:
        raise SystemExit(
            f'FATAL: refusing to plot {path.name}: it is listed as INVALID in '
            f'plausibility/INVALIDATED_ARTIFACTS.json (wrong working point).')


def load_csv(path: Path) -> pd.DataFrame:
    _guard(path)
    df = pd.read_csv(path)
    USED.append(dict(path=str(path.relative_to(FR)).replace('\\', '/'),
                     sha256=sha256(path), rows=int(len(df)), kind='csv'))
    return df


def load_json(path: Path):
    _guard(path)
    doc = json.loads(path.read_text(encoding='utf-8'))
    USED.append(dict(path=str(path.relative_to(FR)).replace('\\', '/'),
                     sha256=sha256(path), kind='json'))
    return doc


def verdict_markers(ax, x, certified, y, dy=0.0, size=170, annotate_legend=True):
    """Filled = certified, open = not certified (the frozen verifier's verdict)."""
    x = np.asarray(x, dtype=float)
    cert = np.asarray(certified, dtype=bool)
    y = np.asarray(y, dtype=float) + dy
    ok_lab = 'certified' if annotate_legend else '_nolegend_'
    bad_lab = 'not certified' if annotate_legend else '_nolegend_'
    ax.scatter(x[cert], y[cert], s=size, marker='o', facecolor=C_OK,
               edgecolor='k', linewidth=0.8, zorder=5, label=ok_lab)
    ax.scatter(x[~cert], y[~cert], s=size, marker='X', facecolor=C_BAD,
               edgecolor='k', linewidth=0.8, zorder=5, label=bad_lab)


# =====================================================================  FIG 1
def fig_2bit_single_param():
    d = load_csv(SRC_BIT2_ROBUST)
    centre = dict(mrna=2.0, mat=32.5, uM=5.75)
    slices = {
        'carry mRNA half-life (min)': ('carry_mrna_half_life_min', centre['mat'], centre['uM']),
        'carry maturation half-life (min)': ('carry_maturation_half_life_min',
                                             centre['mrna'], centre['uM']),
        'uM per a.u.': ('uM_per_au', centre['mrna'], centre['mat']),
    }
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for ax, (label, (col, a, b)) in zip(axes.ravel()[:3], slices.items()):
        if col == 'carry_mrna_half_life_min':
            sub = d[(d.carry_maturation_half_life_min == a) & (d.uM_per_au == b)]
        elif col == 'carry_maturation_half_life_min':
            sub = d[(d.carry_mrna_half_life_min == a) & (d.uM_per_au == b)]
        else:
            sub = d[(d.carry_mrna_half_life_min == a)
                    & (d.carry_maturation_half_life_min == b)]
        sub = sub.sort_values(col)
        if sub.empty:
            raise SystemExit(f'FATAL: empty centre slice for {col}')
        ax.plot(sub[col], sub.bit1_setup_h, 'o-', color='#1c7ed6', lw=1.6,
                ms=6, label='bit1 setup margin (h)')
        ax.plot(sub[col], sub.bit1_hold_h, 's-', color='#e8590c', lw=1.6,
                ms=6, label='bit1 hold margin (h)')
        for _, r in sub.iterrows():
            ax.annotate(('ok' if r.certified else 'FAIL'),
                        xy=(r[col], max(r.bit1_setup_h, r.bit1_hold_h)),
                        xytext=(0, 7), textcoords='offset points',
                        ha='center', fontsize=8,
                        color=(C_OK if r.certified else C_BAD))
        ax.set_xlabel(label)
        ax.set_ylabel('margin (h)')
        ax.set_title(f'2-bit single-parameter slice through the selected centre\n'
                     f'(other two fixed at mRNA={a:g} min, mat={b:g} min, '
                     f'uM={centre["uM"]:g})', fontsize=10)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc='best')
    # marginal pass rate per parameter value over the whole 75-point grid
    ax = axes.ravel()[3]
    groups = [('mRNA (min)', 'carry_mrna_half_life_min'),
              ('maturation (min)', 'carry_maturation_half_life_min'),
              ('uM per a.u.', 'uM_per_au')]
    cols = ['#1c7ed6', '#e8590c', '#7048e8']
    xt, ticks, ticklabels = 0.0, [], []
    for (name, col), c in zip(groups, cols):
        g = d.groupby(col).certified.agg(['sum', 'count']).sort_index()
        for val, row in g.iterrows():
            frac = float(row['sum']) / float(row['count'])
            ax.bar(xt, frac, width=0.72, color=c, edgecolor='k', linewidth=0.6)
            ax.text(xt, frac + 0.02, f'{int(row["sum"])}/{int(row["count"])}',
                    ha='center', fontsize=8)
            ticks.append(xt)
            ticklabels.append(f'{val:g}')
            xt += 1.0
        xt += 0.8
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels, fontsize=8)
    ax.set_ylim(0, 1.18)
    ax.set_ylabel('certified fraction')
    ax.set_title('marginal robustness over the whole 75-point grid\n'
                 '(3 mRNA x 5 maturation x 5 uM_per_au, 600 h each) - note the '
                 'NON-MONOTONE uM dependence', fontsize=10)
    ax.grid(alpha=0.25, axis='y')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=c, edgecolor='k', label=n)
                       for (n, _), c in zip(groups, cols)], fontsize=8,
              title='sliced axis', title_fontsize=8)
    fig.suptitle('2-bit (34-state) single-parameter scans — frozen two-bit profile, '
                 '600 h, finite read-window criterion', fontsize=12)
    fig.savefig(OUT / 'fig_2bit_single_param.png', dpi=190, bbox_inches='tight')
    plt.close(fig)


# =====================================================================  FIG 2
def fig_3bit_single_param():
    legacy = load_csv(SRC_BIT3_CARRY1)
    frozen = load_csv(SRC_NA1_MAT)
    pre = load_csv(SRC_PREAUDIT)
    fig, axes = plt.subplots(3, 2, figsize=(15, 12.5), constrained_layout=True)

    # (a) LEGACY carry1 mRNA slice
    ax = axes[0, 0]
    sub = legacy[legacy.carry1_maturation_min == 60.0].sort_values('carry1_mrna_min')
    ratio = sub.bit2_crossings / sub.cycles
    ax.plot(sub.carry1_mrna_min, ratio, 'o-', color='#c92a2a', lw=1.6, ms=6,
            label='bit2 crossings / clock cycles')
    ax.set_ylim(-0.03, 1.05)
    ax.set_xscale('log')
    ax.set_xlabel('carry1 mRNA half-life (min)')
    ax.set_ylabel('crossings / cycles')
    ax2 = ax.twinx()
    ax2.plot(sub.carry1_mrna_min, sub.g1_peak, 's--', color='#1c7ed6', lw=1.4,
             ms=5, label='gate peak g1')
    ax2.set_ylabel('gate peak g1 (a.u.)', color='#1c7ed6')
    ax.set_title('LEGACY 51-state scan — carry1 mRNA (maturation fixed 60 min)\n'
                 '0 of 35 points certified: the carry gate forms but bit2 never '
                 'counts', fontsize=10)
    ax.grid(alpha=0.25)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc='center right')

    # (b) LEGACY carry1 maturation slice + the split-scan verdict
    ax = axes[0, 1]
    sub = legacy[legacy.carry1_mrna_min == 4.0].sort_values('carry1_maturation_min')
    ax.plot(sub.carry1_maturation_min, sub.bit2_crossings / sub.cycles, 'o-',
            color='#c92a2a', lw=1.6, ms=6, label='bit2 crossings / clock cycles')
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel('carry1 maturation half-life (min)')
    ax.set_ylabel('crossings / cycles')
    ax2 = ax.twinx()
    # gate_dose_median is NaN exactly where gate_events == 0, i.e. where no gate
    # formed at all.  The line is therefore drawn WITH A BREAK: a missing value
    # means "not measurable here", and joining across it would invent a dose.
    ax2.plot(sub.carry1_maturation_min, sub.gate_dose_median, 's--',
             color='#e8590c', lw=1.4, ms=5, label='gate dose (a.u.*h)')
    n_nan = int(sub.gate_dose_median.isna().sum())
    if n_nan:
        first = sub[sub.gate_dose_median.isna()].carry1_maturation_min.max()
        ax2.annotate(f'break: no gate formed\nat maturation <= {first:g} min\n'
                     f'(dose not measurable, {n_nan} point(s))',
                     xy=(first, 0.01), xytext=(14, 26), textcoords='offset points',
                     fontsize=8, color='#e8590c',
                     arrowprops=dict(arrowstyle='->', color='#e8590c', lw=1.0))
    ax2.set_ylabel('gate dose (a.u.·h)', color='#e8590c')
    split = load_csv(SRC_BIT3_SPLIT)
    n_cls = split.classification.value_counts().to_dict()
    ax.text(0.03, 0.62, 'split A1/F1 maturation scan\n%d points, certified %d\n%s'
            % (len(split), int(split.certified.sum()),
               ', '.join(f'{k}: {v}' for k, v in n_cls.items())),
            transform=ax.transAxes, fontsize=8.5, va='top',
            bbox=dict(fc='#fff4e6', ec='0.6'))
    ax.set_title('LEGACY 51-state scan — carry1 maturation (mRNA fixed 4 min)\n'
                 'longer maturation moves the dose but no point counts', fontsize=10)
    ax.grid(alpha=0.25)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc='center right')

    # (c) FROZEN n_A1_gate slice
    ax = axes[1, 0]
    sub = frozen[frozen.a1_mat == 32.5].sort_values('n_A1_gate')
    ax.plot(sub.n_A1_gate, sub.crossings, 'o-', color='#2b8a3e', lw=1.6, ms=6,
            label='bit2 crossings (expect 14)')
    ax.plot(sub.n_A1_gate, sub.gates, '^--', color='0.5', lw=1.2, ms=5,
            label='carry gate events')
    verdict_markers(ax, sub.n_A1_gate, sub.certified,
                    np.full(len(sub), 4.0), size=120)
    ax.set_ylim(2, 15.5)
    ax.set_xlabel('n_A1_gate  (A1 Hill exponent of the carry gate)')
    ax.set_ylabel('events in 400 h')
    ax2 = ax.twinx()
    ax2.plot(sub.n_A1_gate, sub.leak_symmetric, 'd-.', color='#c92a2a', lw=1.4,
             ms=5, label='paired leak L_symmetric (5 %)')
    ax2.set_yscale('log')
    ax2.set_ylabel('L_symmetric (paired R/F, 5 % cut)', color='#c92a2a')
    ax.set_title('FROZEN 51-state working point — n_A1_gate (A1 maturation 32.5 min)\n'
                 'n=4 (legacy shared slot) fails outright; n>=5 counts', fontsize=10)
    ax.grid(alpha=0.25)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc='center right')

    # (d) FROZEN A1 maturation slice
    ax = axes[1, 1]
    sub = frozen[frozen.n_A1_gate == 6.0].sort_values('a1_mat')
    ax.plot(sub.a1_mat, sub.crossings, 'o-', color='#2b8a3e', lw=1.6, ms=6,
            label='bit2 crossings')
    verdict_markers(ax, sub.a1_mat, sub.certified, np.full(len(sub), 8.0), size=120)
    ax.set_ylim(6, 15.5)
    ax.set_xlabel('A1 maturation half-life (min)   [F1 arm frozen at 32.5]')
    ax.set_ylabel('bit2 crossings in 400 h')
    ax2 = ax.twinx()
    ax2.plot(sub.a1_mat, sub.gate_dose_median, 's--', color='#e8590c', lw=1.4,
             ms=5, label='gate dose (a.u.·h)')
    ax2.set_ylabel('gate dose (a.u.·h)', color='#e8590c')
    ax.set_title('FROZEN 51-state working point — A1 maturation at n_A1_gate = 6\n'
                 'all 5 values certify; dose changes by <5 %, leak by <=2.7 %',
                 fontsize=10)
    ax.grid(alpha=0.25)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc='center right')

    # (e) FROZEN static knobs: conc_scale and clock_K / K_A1
    for ax, knobs, title in (
            (axes[2, 0], ('conc_scale',),
             'FROZEN — global concentration scale f x uM_per_au\n'
             'certified only at f = 0.95 and 1.00'),
            (axes[2, 1], ('clock_K', 'K_A1'),
             'FROZEN — clock threshold and gate threshold K_A1\n'
             'clock_K certifies across the whole ±20 % grid')):
        for i, knob in enumerate(knobs):
            sub = pre[pre.knob == knob].sort_values('factor')
            style = dict(lw=1.6, ms=6)
            ax.plot(sub.factor, sub.off_on_gate_peak_ratio,
                    'o-' if i == 0 else 's-',
                    color='#1c7ed6' if i == 0 else '#7048e8', label=knob, **style)
            verdict_markers(ax, sub.factor, sub.certified,
                            np.full(len(sub), 0.268), size=90)
        ax.axhline(0.10, color=C_BAD, ls='--', lw=1.4)
        ax.annotate('contrast limit 0.10  (verify_threebit51.MAX_OFF_ON_GATE_RATIO)',
                    xy=(0.02, 0.10), xycoords=('axes fraction', 'data'),
                    xytext=(2, 4), textcoords='offset points', fontsize=8,
                    color=C_BAD)
        ax.set_ylim(0, 0.30)
        ax.set_xlabel('factor f (1.0 = frozen value)')
        ax.set_ylabel('off / on gate peak ratio')
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc='upper left')
    fig.suptitle('3-bit (51-state) single-parameter scans — LEGACY shared-slot '
                 'exponent (n_A[1] = 4.0) vs FROZEN n_A1_gate = 6.0, 600 h/400 h, '
                 'finite read-window criterion', fontsize=12)
    fig.savefig(OUT / 'fig_3bit_single_param.png', dpi=190, bbox_inches='tight')
    plt.close(fig)


# =====================================================================  FIG 3
def fig_2bit_robustness():
    d = load_csv(SRC_BIT2_ROBUST)
    strict = load_csv(SRC_BIT2_STRICT)
    perturb = {g: load_csv(p) for g, p in SRC_BIT2_PERTURB.items()}
    fig, axes = plt.subplots(2, 3, figsize=(15, 9.5), constrained_layout=True)

    for ax, mrna in zip(axes.ravel()[:3], (1.0, 2.0, 4.0)):
        sub = d[d.carry_mrna_half_life_min == mrna]
        ok = sub[sub.certified]
        bad = sub[~sub.certified]
        ax.scatter(ok.uM_per_au, ok.carry_maturation_half_life_min, s=150,
                   marker='o', facecolor=C_OK, edgecolor='k', linewidth=0.8,
                   label='certified')
        ax.scatter(bad.uM_per_au, bad.carry_maturation_half_life_min, s=150,
                   marker='X', facecolor=C_BAD, edgecolor='k', linewidth=0.8,
                   label='not certified')
        if mrna == 2.0:
            ax.scatter([WP['uM_per_au']], [WP['mat_min']], s=380, marker='*',
                       facecolor='#f59f00', edgecolor='k', linewidth=0.9,
                       zorder=6, label='selected centre')
        ax.set_xlabel('uM per a.u.')
        ax.set_ylabel('carry maturation half-life (min)')
        ax.set_title(f'carry mRNA = {mrna:g} min   '
                     f'({int(sub.certified.sum())}/{len(sub)} certified)',
                     fontsize=10)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc='lower right')

    ax = axes[1, 0]
    sub = d[d.carry_mrna_half_life_min == 2.0]
    piv = sub.pivot_table(index='carry_maturation_half_life_min',
                          columns='uM_per_au', values='timing_margin_h')
    im = ax.imshow(piv.to_numpy(), origin='lower', aspect='auto', cmap='viridis',
                   extent=[piv.columns.min(), piv.columns.max(),
                           piv.index.min(), piv.index.max()])
    for _, r in sub.iterrows():
        ax.scatter(r.uM_per_au, r.carry_maturation_half_life_min, s=90,
                   marker='o' if r.certified else 'X',
                   facecolor='none', edgecolor='w' if r.certified else C_BAD,
                   linewidth=1.4, zorder=5)
    ax.set_xlabel('uM per a.u.')
    ax.set_ylabel('carry maturation half-life (min)')
    ax.set_title('timing margin (h) at mRNA = 2 min\n'
                 'white circles = certified, red X = not certified', fontsize=10)
    fig.colorbar(im, ax=ax, label='timing margin (h)')

    ax = axes[1, 1]
    order = ['selected_centre', 'pass_low_uM_edge', 'pass_low_mat_edge',
             'pass_notch_lip', 'pass_high_corner', 'fail_low_uM', 'fail_low_mat',
             'fail_notch_mid', 'fail_notch_low']
    g = strict.groupby('point_name').certified.agg(['sum', 'count'])
    g = g.reindex([o for o in order if o in g.index])
    y = np.arange(len(g))
    ax.barh(y, g['sum'] / g['count'], color=[C_OK if s == c else '#f59f00'
                                             for s, c in zip(g['sum'], g['count'])],
            edgecolor='k', linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([f'{n}\n({int(s)}/{int(c)} settings)' for n, s, c
                        in zip(g.index, g['sum'], g['count'])], fontsize=8)
    ax.set_xlim(0, 1.15)
    ax.set_xlabel('fraction of solver settings that certify')
    ax.set_title('strict solver-tolerance re-check (rtol / atol / max_step)\n'
                 'no verdict flips: 4 pass points stay pass, 4 fail points stay fail',
                 fontsize=10)
    ax.grid(alpha=0.25, axis='x')

    ax = axes[1, 2]
    # The first version of this panel collapsed every perturbation to a single
    # "fraction recovered" bar.  The data does not support that: the outcome
    # vocabulary has FIVE values, and three of them are not failures - a
    # phase-shifted run is still certified and still counts.  Collapsing them hid
    # both the phase-shift results and the seven genuine `lost_counting` rows.
    CATS = PERTURB_CATS
    tally = {}
    for gname, df in perturb.items():
        xcol = 'factor' if 'factor' in df.columns else 'action'
        for xval, sub in df.groupby(xcol):
            label = f'{gname} {str(xval).replace("delta_", "")}'
            tally[label] = {c: int((sub.outcome == c).sum()) for c, _, _ in CATS}
    names = list(tally)
    bottoms = np.zeros(len(names))
    for cat, col, lab in CATS:
        vals = np.array([tally[n][cat] for n in names], dtype=float)
        ax.bar(np.arange(len(names)), vals, bottom=bottoms, color=col,
               edgecolor='k', linewidth=0.5, label=lab)
        bottoms += vals
    ax.set_xticks(np.arange(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_ylabel('initial states out of 4')
    ax.set_title('pool perturbation outcome, by category (4 initial states each)\n'
                 'green = same phase, orange = phase shifted but STILL CERTIFIED, '
                 'red = counting lost', fontsize=10)
    ax.legend(fontsize=7.5, ncol=2, loc='upper left', framealpha=0.95)
    ax.grid(alpha=0.25, axis='y')
    fig.suptitle('2-bit (34-state) robustness — 75-point joint parameter grid, '
                 'strict solver tolerance, and molecular-pool perturbations',
                 fontsize=12)
    fig.savefig(OUT / 'fig_2bit_robustness.png', dpi=190, bbox_inches='tight')
    plt.close(fig)


# =====================================================================  FIG 4
def fig_3bit_robustness():
    frozen = load_csv(SRC_NA1_MAT)
    pre = load_csv(SRC_PREAUDIT)
    verdict = load_json(SRC_PREAUDIT_VERDICT)
    inventory = load_json(SRC_INVENTORY)
    pairing = load_json(SRC_PAIRING)
    post = load_csv(SRC_BIT3_POSTFREEZE)
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 10), constrained_layout=True)

    # (a) certification map n x A1 maturation
    ax = axes[0, 0]
    for cert, mk, col, lab in ((True, 'o', C_OK, 'certified'),
                               (False, 'X', C_BAD, 'not certified')):
        sub = frozen[frozen.certified == cert]
        ax.scatter(sub.n_A1_gate, sub.a1_mat, s=180, marker=mk, facecolor=col,
                   edgecolor='k', linewidth=0.8, label=lab)
    ax.scatter([WP['n_A1_gate']], [WP['mat_min']], s=420, marker='*',
               facecolor='#f59f00', edgecolor='k', linewidth=0.9, zorder=6,
               label='frozen working point')
    ax.set_xlabel('n_A1_gate')
    ax.set_ylabel('A1 maturation half-life (min)')
    ax.set_title('FROZEN 25-point grid, 400 h — all 25 certified\n'
                 'the gate exponent is what decides, not the maturation', fontsize=10)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc='lower right')

    # (b) leakage knee from the 600 h pairing grid
    ax = axes[0, 1]
    by_n = pairing['by_n_A1_gate']
    ns = sorted(by_n, key=float)
    for cut, col in (('1pct', '#1c7ed6'), ('5pct', '#2b8a3e'), ('10pct', '#e8590c')):
        ys = [by_n[n][f'L_symmetric_{cut}_median'] for n in ns]
        ax.plot([float(n) for n in ns], ys, 'o-', color=col, lw=1.6, ms=6,
                label=f'{cut.replace("pct", " %")} tail cut')
    ax.set_yscale('log')
    ax.set_xlabel('n_A1_gate')
    ax.set_ylabel('L_symmetric (paired R/F, median)')
    ax.annotate('order 8<7<6<5<4\nholds at all three cuts',
                xy=(0.03, 0.60), xycoords='axes fraction', fontsize=8.5,
                bbox=dict(fc='#ebfbee', ec='0.6'))
    ax.annotate('5 -> 6 = 3.08x\n6 -> 7 = 1.22x\n7 -> 8 = 1.06x',
                xy=(0.55, 0.62), xycoords='axes fraction', fontsize=8.5,
                bbox=dict(fc='#fff9db', ec='0.6'))
    ax.set_title('FROZEN 20-point grid, 600 h — the leakage knee\n'
                 'per-n median over mRNA x maturation', fontsize=10)
    ax.grid(alpha=0.25, which='both')

    # (c) per-point ordering across mRNA/maturation
    ax = axes[0, 2]
    pp = pairing['per_point']
    styles = {('2', '32.5'): ('o-', '#2b8a3e'), ('2', '60'): ('s--', '#e8590c'),
              ('4', '32.5'): ('^-.', '#1c7ed6'), ('4', '60'): ('d:', '#7048e8')}
    for (mrna, mat), (ls, col) in styles.items():
        xs, ys = [], []
        for n in '45678':
            key = f'n={n}|mRNA={mrna}|mat={mat}'
            if key in pp:
                xs.append(float(n))
                ys.append(pp[key]['5pct_L_symmetric_median'])
        ax.plot(xs, ys, ls, color=col, lw=1.5, ms=6,
                label=f'mRNA={mrna} min, mat={mat} min')
    ax.set_yscale('log')
    ax.set_xlabel('n_A1_gate')
    ax.set_ylabel('L_symmetric (5 % cut)')
    ax.set_title('the ordering is not an artefact of one setting\n'
                 'all four mRNA x maturation combinations agree', fontsize=10)
    ax.grid(alpha=0.25, which='both')
    ax.legend(fontsize=8)

    # (d) static tolerance widths + contrast arm
    ax = axes[1, 0]
    widths = verdict['symmetric_certified_width']
    names = list(widths)
    vals = [widths[k] if widths[k] is not None else np.nan for k in names]
    bars = ax.bar(names, vals, color=['#c92a2a' if (v != v or v < 0.05) else
                                      '#2b8a3e' for v in vals],
                  edgecolor='k', linewidth=0.6)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, (0 if v != v else v) + 0.004,
                'n/a' if v != v else f'{v:.2f}', ha='center', fontsize=9)
    ax.axhline(0.05, color='0.4', ls='--', lw=1.2)
    ax.annotate('pre-registered SMALL threshold 0.05', xy=(0.02, 0.05),
                xycoords=('axes fraction', 'data'), xytext=(2, 4),
                textcoords='offset points', fontsize=8, color='0.35')
    ax.set_ylim(0, 0.24)
    ax.set_ylabel('symmetric certified width')
    ax.set_title('FROZEN static-tolerance scan, 21 points x 600 h\n'
                 'worst knob = global concentration scale (0.00): only '
                 'f=0.95 and 1.00 certify', fontsize=10)
    ax.grid(alpha=0.25, axis='y')

    # (e) which arm failed (the reason the above is not a counting failure)
    ax = axes[1, 1]
    arms = verdict['certification_arms']['failures_by_arm']
    ax.bar(list(arms), list(arms.values()), color='#f59f00', edgecolor='k',
           linewidth=0.6)
    for i, (k, v) in enumerate(arms.items()):
        ax.text(i, v + 0.06, str(v), ha='center', fontsize=9)
    ax.set_ylim(0, max(arms.values()) + 1)
    ax.set_ylabel('non-certified points failing this arm')
    ax.set_title('why those 7 points fail: predicate =\n'
                 'steady ∧ one ∧ order ∧ alt ∧ contrast', fontsize=10)
    ax.annotate('4 of the 7 fail the CONTRAST arm ALONE\n'
                '(off/on gate peak ratio > 0.10) while their\n'
                'digital readout is perfect (14/14),\n'
                'and 1 of those has BETTER leak than the frozen point',
                xy=(0.03, 0.96), xycoords='axes fraction', fontsize=8,
                va='top', bbox=dict(fc='#fff4e6', ec='0.6'))
    ax.grid(alpha=0.25, axis='y')

    # (f) copy-number inventory (why stochasticity is the open question)
    ax = axes[1, 2]
    reg = inventory['noise_regime_per_pool']
    labels, vals, cols, hatch, absent = [], [], [], [], 0
    for pool, entries in reg.items():
        for stat in ('gate_phase_minimum', 'whole_cycle_trough'):
            if stat not in entries:
                continue
            e = entries[stat]
            is_absent = e.get('implied_cv_pct') is None
            absent += int(is_absent)
            labels.append(f'{pool}\n{stat.replace("_", " ")}')
            vals.append(max(float(e['copies']), 1e-3))
            cols.append('#868e96' if is_absent else '#1c7ed6')
            hatch.append('//' if is_absent else None)
    bars = ax.barh(np.arange(len(vals)), vals, color=cols, edgecolor='k',
                   linewidth=0.5)
    for b, h in zip(bars, hatch):
        b.set_hatch(h)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(100, color=C_BAD, ls='--', lw=1.4)
    ax.annotate('100 copies = the registered\nsmall-number threshold',
                xy=(100, len(vals) - 0.6), xytext=(4, -14),
                textcoords='offset points', fontsize=8, color=C_BAD)
    ax.set_xscale('log')
    ax.set_xlim(1e-3, 1e5)
    ax.set_xlabel('copies per cell (log scale)')
    ax.set_title('FROZEN copy-number inventory (600 h, uM_per_au = 5.75)\n'
                 'the write is carried out by ~20 Int2 molecules', fontsize=10)
    ax.annotate('grey hatch = STRUCTURALLY ABSENT\n'
                '(below 1 copy: a zero, not a level);\n'
                'bars are floored at 1e-3 to stay on the log axis',
                xy=(0.30, 0.30), xycoords='axes fraction', fontsize=7.5,
                bbox=dict(fc='#f1f3f5', ec='0.6'))
    ax.grid(alpha=0.25, axis='x')

    pf = post[post.outcome != 'baseline']
    n_rec = int((pf.outcome == 'recovered_original').sum())
    fig.suptitle('3-bit (51-state) robustness — frozen point n_A1_gate = 6, '
                 '600 h pairing grid, 400 h maps, static tolerance, copy numbers '
                 f'| peak-phase pool perturbation: {n_rec}/{len(pf)} recover',
                 fontsize=11.5)
    fig.savefig(OUT / 'fig_3bit_robustness.png', dpi=190, bbox_inches='tight')
    plt.close(fig)


def png_dimensions(path: Path):
    with open(path, 'rb') as fh:
        head = fh.read(24)
    if head[:8] != b'\x89PNG\r\n\x1a\n':
        raise SystemExit(f'FATAL: {path.name} is not a PNG')
    import struct
    return struct.unpack('>II', head[16:24])


def main():
    which = sys.argv[1:] or ['2bit_sp', '3bit_sp', '2bit_rob', '3bit_rob']
    made = []
    if '2bit_sp' in which:
        fig_2bit_single_param()
        made.append('fig_2bit_single_param.png')
    if '3bit_sp' in which:
        fig_3bit_single_param()
        made.append('fig_3bit_single_param.png')
    if '2bit_rob' in which:
        fig_2bit_robustness()
        made.append('fig_2bit_robustness.png')
    if '3bit_rob' in which:
        fig_3bit_robustness()
        made.append('fig_3bit_robustness.png')
    outputs = {}
    for name in made:
        p = OUT / name
        w, h = png_dimensions(p)
        outputs[name] = dict(sha256=sha256(p), bytes=p.stat().st_size,
                             width=int(w), height=int(h))
    manifest = dict(
        purpose=('provenance for every scan figure: which completed artefact each '
                 'panel was drawn from, with SHA256 and row counts, plus the SHA256 '
                 'of every figure produced by this run'),
        figures=made, outputs=outputs, inputs=USED,
        working_points=dict(
            frozen_51state=WP,
            legacy_51state=dict(note=('these scans predate the n_A1_gate constructor '
                                      'parameter and therefore ran at the shared-slot '
                                      'ZENG["n_A"][1] = 4.0, NOT the frozen 6.0'),
                                files=['threebit51_results/carry1_scan_dsh/carry1_all.csv',
                                       'threebit51_results/carry1_split_scan/split_all.csv']),
            twobit_34state=dict(carry_mrna_half_life_min=2.0,
                                carry_maturation_half_life_min=32.5, uM_per_au=5.75)),
        metric_notes=('leak is always the paired R/F super-period L_symmetric (H4) on the '
                      '2-minute output grid, with the tail cut stated per panel; the '
                      'legacy 51-state scans report crossings/cycles and gate dose, not '
                      'L_symmetric'),
        verification=('figures were NOT visually inspected (the model has no image '
                      'input); they are verified numerically by '
                      'selfcheck_scan_figures.py, which re-renders every panel and '
                      'asserts the plotted values against the source artefacts'),
        not_plotted=[str(p) for p in sorted(INVALIDATED)][:40])
    (OUT / 'figure_data_manifest.json').write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    print('wrote:', ', '.join(made))
    print('inputs:', len(USED))
    print('wrote figure_data_manifest.json')


if __name__ == '__main__':
    main()
