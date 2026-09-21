"""Run the formal 34-state model and generate report-ready outputs/figures."""
from __future__ import annotations

import json
import hashlib
from dataclasses import asdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd

from model import ROOT, ZENG
from model_twobit34 import STATE_NAMES_34, TwoBit34Model
from verify_twobit_causal import analyse_solution


OUT = ROOT / 'twobit34_results'
FIG = OUT / 'figures'


def savefig(fig, name):
    fig.savefig(FIG / f'{name}.png', dpi=220, bbox_inches='tight', facecolor='white')
    fig.savefig(FIG / f'{name}.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)


def overview_figure(t, y, flux, sig, reads):
    fig, ax = plt.subplots(4, 1, figsize=(13, 10), sharex=True, constrained_layout=True)
    ax[0].plot(t, flux, color='#6A3D9A', lw=1.2)
    ax[0].set_ylabel('C31 flux\n(a.u./h)')
    ax[0].set_title('Formal 34-state two-bit counter: 300 h trajectory')

    ax[1].plot(t, y[16], label='bit0 S0', color='#1F78B4', lw=1.3)
    ax[1].plot(t, y[27], label='bit1 S1', color='#E31A1C', lw=1.3)
    ax[1].axhspan(0, .3, color='#33A02C', alpha=.08)
    ax[1].axhspan(.7, 1, color='#1F78B4', alpha=.08)
    ax[1].axhline(.3, color='0.6', ls='--', lw=.8)
    ax[1].axhline(.7, color='0.6', ls='--', lw=.8)
    ax[1].set_ylim(-.03, 1.03); ax[1].set_ylabel('DNA state S')
    ax[1].legend(ncol=2, loc='upper right')

    rt = np.asarray([r['trough_h'] for r in reads])
    rv = np.asarray([np.nan if r['value'] is None else r['value'] for r in reads])
    ax[2].step(rt, rv, where='mid', color='#111111', lw=1.6)
    ax[2].scatter(rt, rv, c=rv, cmap='viridis', vmin=0, vmax=3,
                  edgecolor='white', linewidth=.4, s=35, zorder=3)
    ax[2].set_yticks((0, 1, 2, 3)); ax[2].set_ylim(-.35, 3.35)
    ax[2].set_ylabel('decoded value')
    ax[2].text(.01, .9, ''.join(str(int(v)) for v in rv if np.isfinite(v)),
               transform=ax[2].transAxes, fontsize=9, family='monospace')

    ax[3].plot(t, sig['J_rev0'], label='bit0 reverse flux', color='#FF7F00', lw=1.1)
    ax[3].plot(t, sig['Int1_source'] / max(sig['Int1_source'].max(), 1e-12),
               label='Int1 source (normalized)', color='#33A02C', lw=1.1)
    ax[3].plot(t, sig['g0'] / max(sig['g0'].max(), 1e-12),
               label='carry gate g0 (normalized)', color='#B15928', lw=1.0, alpha=.8)
    ax[3].set_ylabel('causal signals'); ax[3].set_xlabel('time (h)')
    ax[3].legend(ncol=3, loc='upper right', fontsize=8)
    return fig


def causal_zoom_figure(t, y, sig, analysis):
    event = analysis['reverse_events'][5]
    lo, hi = event['start_h'] - 3.0, event['end_h'] + 6.0
    mask = (t >= lo) & (t <= hi)
    fig, ax = plt.subplots(4, 1, figsize=(11, 9), sharex=True, constrained_layout=True)
    ax[0].plot(t[mask], sig['J_fwd0'][mask], label='forward recombination', color='#1F78B4')
    ax[0].plot(t[mask], sig['J_rev0'][mask], label='reverse recombination', color='#E31A1C')
    ax[0].set_ylabel('DNA flux (/h)'); ax[0].legend(loc='upper right')

    ax[1].plot(t[mask], sig['g0'][mask], label='g0', color='#FF7F00')
    ax1r = ax[1].twinx()
    ax1r.plot(t[mask], sig['Int1_source'][mask], label='Int1 source', color='#33A02C')
    ax[1].set_ylabel('g0'); ax1r.set_ylabel('Int1 source (a.u./h)')

    ax[2].plot(t[mask], y[16, mask], label='S0', color='#1F78B4')
    ax[2].plot(t[mask], y[27, mask], label='S1', color='#E31A1C')
    ax[2].axhline(.5, color='0.5', ls='--', lw=.8)
    ax[2].set_ylabel('DNA state'); ax[2].legend(loc='upper right')

    ax[3].plot(t[mask], y[28, mask], label='A0', color='#6A3D9A')
    ax[3].plot(t[mask], y[29, mask], label='F0', color='#B15928')
    ax[3].plot(t[mask], y[31, mask], label='A0 immature', color='#CAB2D6', ls='--')
    ax[3].plot(t[mask], y[33, mask], label='F0 immature', color='#FDBF6F', ls='--')
    ax[3].set_ylabel('I1-FFL states'); ax[3].set_xlabel('time (h)')
    ax[3].legend(ncol=4, fontsize=8, loc='upper right')

    for a in ax:
        a.axvspan(event['start_h'], event['end_h'], color='#FB9A99', alpha=.12)
        a.axvline(event['peak_h'], color='#E31A1C', ls=':', lw=1)
    ax[0].set_title('Causal zoom: reverse recombination → carry gate → bit1 switch')
    return fig


def grouped_states_figure(t, y):
    fig, ax = plt.subplots(4, 1, figsize=(13, 11), sharex=True, constrained_layout=True)

    def normalized(indices):
        vals = y[indices]
        scale = np.maximum(np.max(np.abs(vals), axis=1, keepdims=True), 1e-12)
        return vals / scale

    for v, label in zip(normalized([1, 3, 5]), ('TetR total', 'CI', 'LacI')):
        ax[0].plot(t, v, label=label, lw=1)
    ax[0].set_ylabel('upstream\nnormalized'); ax[0].legend(ncol=3, fontsize=8)

    for v, label in zip(normalized([8, 11, 14, 15]), ('Int0', 'Rep0', 'RDF0', 'Int-RDF0')):
        ax[1].plot(t, v, label=label, lw=1)
    ax[1].plot(t, y[16], label='S0', color='black', lw=1.3)
    ax[1].set_ylabel('bit0\nnormalized'); ax[1].legend(ncol=5, fontsize=8)

    for v, label in zip(normalized([19, 22, 25, 26]), ('Int1', 'Rep1', 'RDF1', 'Int-RDF1')):
        ax[2].plot(t, v, label=label, lw=1)
    ax[2].plot(t, y[27], label='S1', color='black', lw=1.3)
    ax[2].set_ylabel('bit1\nnormalized'); ax[2].legend(ncol=5, fontsize=8)

    for idx, label, ls in ((28, 'A0', '-'), (29, 'F0', '-'), (30, 'M_A0', '--'),
                           (31, 'A0_u', '--'), (32, 'M_F0', ':'), (33, 'F0_u', ':')):
        v = y[idx] / max(np.max(np.abs(y[idx])), 1e-12)
        ax[3].plot(t, v, label=label, lw=1, ls=ls)
    ax[3].set_ylabel('carry module\nnormalized'); ax[3].set_xlabel('time (h)')
    ax[3].legend(ncol=6, fontsize=8)
    ax[0].set_title('All functional state groups (each molecular trace normalized to its own maximum)')
    return fig


def readout_figure(reads):
    bits = np.asarray([[int(r['bit1']['label']), int(r['bit0']['label'])] for r in reads])
    values = np.asarray([r['value'] for r in reads])
    fig, ax = plt.subplots(figsize=(13, 3.3), constrained_layout=True)
    ax.imshow(bits.T, aspect='auto', cmap='Blues', vmin=0, vmax=1,
              extent=(-.5, len(reads)-.5, 1.5, -.5))
    ax.set_yticks((0, 1), ('bit1', 'bit0')); ax.set_xlabel('clock read index')
    ax.set_title('Finite-window digital readout (all windows have commitment = 1.0)')
    for j, value in enumerate(values):
        ax.text(j, 1.82, str(value), ha='center', va='center', fontsize=8, clip_on=False)
    ax.text(-1.2, 1.82, 'value', ha='right', va='center', fontsize=9, clip_on=False)
    ax.set_xticks(range(len(reads)))
    ax.tick_params(axis='x', labelsize=7)
    return fig


def architecture_figure():
    fig, ax = plt.subplots(figsize=(14, 5.5), constrained_layout=True)
    ax.set_xlim(0, 14); ax.set_ylim(0, 6); ax.axis('off')

    boxes = [
        (0.4, 2.1, 2.2, 1.7, '#CAB2D6', 'Han oscillator\n6 states\ntrue C31 flux'),
        (3.1, 1.7, 2.4, 2.5, '#A6CEE3', 'bit0\n11 states\nM/Iu/I, M/Tu/T,\nM/Ru/R, C, S'),
        (6.1, 1.5, 2.5, 2.9, '#FDBF6F', 'A0/F0 I1-FFL\n6 states\nM_A0 → A0_u → A0\nM_F0 → F0_u → F0'),
        (9.3, 1.7, 2.4, 2.5, '#B2DF8A', 'bit1\n11 states\ncarry enters existing\nM_I / I_u / I'),
        (12.1, 2.1, 1.4, 1.7, '#FB9A99', 'mod-4\nreadout\nS1 S0'),
    ]
    for x, y, w, h, color, text in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.12',
                                    facecolor=color, edgecolor='0.25', lw=1.2))
        ax.text(x+w/2, y+h/2, text, ha='center', va='center', fontsize=10)
    for x0, x1, label in ((2.6, 3.1, 'C31'), (5.5, 6.1, 'PB=1-S0'),
                          (8.6, 9.3, 'Int1 source'), (11.7, 12.1, 'decode')):
        ax.add_patch(FancyArrowPatch((x0, 3), (x1, 3), arrowstyle='-|>',
                                     mutation_scale=14, lw=1.4, color='0.25'))
        ax.text((x0+x1)/2, 3.25, label, ha='center', va='bottom', fontsize=8)
    ax.text(7, 5.35, 'Formal two-bit biochemical counter: 34 continuous states',
            ha='center', fontsize=15, weight='bold')
    ax.text(7, .45,
            'ZENG unchanged • E0/E1 removed • no stateless U • downstream gamma treated as total clearance',
            ha='center', fontsize=10, color='0.25')
    return fig


def main():
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    model = TwoBit34Model()
    sol = model.simulate(hours=300.0, sample_min=2.0, max_step_min=2.0)
    t, y = sol.t, sol.y
    flux = np.asarray([model.flux(y[:, k]) for k in range(y.shape[1])])
    sig = model.diagnostic_signals(y)
    cfg = asdict(model.e)
    analysis = analyse_solution(model.base, sol, cfg, 300.0)
    reads = analysis['read_windows']

    data = {'time_h': t}
    data.update({name: y[i] for i, name in enumerate(STATE_NAMES_34)})
    data['C31_flux_au_h'] = flux
    data.update(sig)
    pd.DataFrame(data).to_csv(OUT / 'trajectories.csv', index=False)

    flat_reads = []
    for r in reads:
        flat_reads.append(dict(cycle=r['cycle'], trough_h=r['trough_h'],
                               window_start_h=r['window_start_h'], window_end_h=r['window_end_h'],
                               bit0=r['bit0']['label'], bit1=r['bit1']['label'], value=r['value'],
                               bit0_commitment=r['bit0']['commitment'],
                               bit1_commitment=r['bit1']['commitment'],
                               bit0_median=r['bit0']['median'], bit1_median=r['bit1']['median']))
    pd.DataFrame(flat_reads).to_csv(OUT / 'read_windows.csv', index=False)

    parameters = dict(extension=asdict(model.e), carry_expression=asdict(model.carry),
                      zeng=ZENG, state_names=list(STATE_NAMES_34),
                      interpretation=('Zeng gamma is treated as total downstream clearance; '
                                      'growth dilution is not added again downstream.'))
    (OUT / 'parameters.json').write_text(json.dumps(parameters, ensure_ascii=False, indent=2), encoding='utf-8')
    summary = dict(states=34, hours=300.0, samples=int(t.size),
                   cold_start=analysis['cold_start'], steady_state=analysis['steady_state'],
                   certified=analysis['certified'], causal_verdict=analysis['causal_verdict'],
                   signal_ranges=analysis['signal_ranges'], margins=analysis['margins'],
                   reverse_events=len(analysis['reverse_events']), gate_events=len(analysis['gate_events']))
    (OUT / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

    savefig(overview_figure(t, y, flux, sig, reads), '01_counter_overview')
    savefig(causal_zoom_figure(t, y, sig, analysis), '02_carry_causality_zoom')
    savefig(grouped_states_figure(t, y), '03_state_groups')
    savefig(readout_figure(reads), '04_digital_readout')
    savefig(architecture_figure(), '05_model_architecture')

    md = f"""# 正式34状态两比特模型结果

- 300 h冷启动序列：`{analysis['cold_start']['sequence']}`
- 稳态序列：`{analysis['steady_state']['sequence']}`
- 完整因果认证：**{analysis['certified']}**
- 反向重组事件 / carry事件：{len(analysis['reverse_events'])} / {len(analysis['gate_events'])}
- bit1最小建立/保持裕量：{analysis['margins']['bit1']['min_setup_h']:.2f} h / {analysis['margins']['bit1']['min_hold_h']:.2f} h

图：

1. `figures/01_counter_overview.png`：300 h总体计数；
2. `figures/02_carry_causality_zoom.png`：反向重组到bit1翻转的局部因果链；
3. `figures/03_state_groups.png`：34状态按功能分组；
4. `figures/04_digital_readout.png`：有限读窗bit图和模4值；
5. `figures/05_model_architecture.png`：报告用模型结构图。

参数与假设见 `parameters.json`；逐状态轨迹见 `trajectories.csv`；读窗见 `read_windows.csv`；
与临时验证变体的逐状态等价结果见 `equivalence.json`。
"""
    (OUT / '结果说明.md').write_text(md, encoding='utf-8')

    def sha256(path):
        h = hashlib.sha256()
        with open(path, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b''):
                h.update(chunk)
        return h.hexdigest().upper()

    manifest = {}
    for path in sorted(OUT.rglob('*')):
        if path.is_file() and not path.name.startswith('SHA256SUMS'):
            manifest[str(path.relative_to(ROOT))] = sha256(path)
    for name in ('model.py', 'model_twobit34.py', 'test_twobit34.py',
                 'verify_twobit34_equivalence.py', 'run_twobit34.py',
                 'verify_twobit_causal.py'):
        path = ROOT / name
        manifest[str(path.relative_to(ROOT))] = sha256(path)
    (OUT / 'SHA256SUMS.json').write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True), encoding='utf-8')
    (OUT / 'SHA256SUMS.txt').write_text(
        ''.join(f'{value}  {name}\n' for name, value in sorted(manifest.items())),
        encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f'wrote {OUT}')


if __name__ == '__main__':
    main()
