"""Reproduce the archived successful n=6 configuration and plot actual ODE data.

IMPORTANT: the archived successful scan uses frozen_extension() from the 34-state
profile (clock K=0.3,n=2), not selected_threebit_extension() (K=0.4,n=3).
This script makes that historical configuration explicit without modifying it.
"""
from __future__ import annotations
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from model import ROOT, Extension, UPSTREAM_PATH, ZENG
from model_twobit34 import CarryExpressionParameters
from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters
from verify_threebit51 import analyse_threebit

OUT = ROOT / 'threebit51_results' / 'wiki_n6_20260924'
PROFILE34 = ROOT / 'twobit34_results' / 'selected_profile.json'
PROFILE51 = ROOT / 'plausibility' / 'threebit51_selected_v1.json'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest().upper()

def save(fig, name):
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(OUT / f'{name}.{ext}', dpi=220, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def overview(sol, sig, flux, reads, left, right, name):
    t, y = sol.t, sol.y
    mask = (t >= left) & (t <= right)
    chosen = [r for r in reads if left <= r['trough_h'] <= right]
    fig, ax = plt.subplots(4, 1, figsize=(14, 10.5), sharex=True,
                           gridspec_kw={'height_ratios': [1, 1.55, 1.15, 1]},
                           layout='constrained')
    ax[0].plot(t[mask], flux[mask], color='#6A3D9A', lw=1.45)
    ax[0].set_ylabel('C31 flux\n(' + r'$\mu$M/h)')
    ax[0].set_title(f'51-state three-bit counter | {left:g}–{right:g} h | independent A1 gate exponent = 6',
                    loc='left', fontsize=15, pad=14)
    colors = ('#2479B5', '#D75B34', '#318A63')
    for idx, color, label in zip((16, 27, 44), colors, ('bit0 / S0', 'bit1 / S1', 'bit2 / S2')):
        ax[1].plot(t[mask], y[idx, mask], color=color, lw=1.6, label=label)
    ax[1].axhspan(0, .3, color='#DCE5E7', alpha=.5)
    ax[1].axhspan(.7, 1, color='#DCEBEC', alpha=.5)
    for h in (.3, .7): ax[1].axhline(h, color='0.55', ls='--', lw=.85)
    ax[1].set_ylim(-.045, 1.045)
    ax[1].set_ylabel('DNA LR fraction')
    ax[1].legend(loc='upper left', bbox_to_anchor=(0, 1.17), ncol=3, frameon=False, fontsize=10)
    rt = np.asarray([r['trough_h'] for r in chosen])
    vals = np.asarray([r['value'] for r in chosen], dtype=float)
    ax[2].step(rt, vals, where='mid', color='#33383A', lw=1.55, zorder=2)
    ax[2].scatter(rt, vals, c=vals, cmap='viridis', vmin=0, vmax=7,
                  edgecolor='white', linewidth=.6, s=42, zorder=3)
    if right-left < 150:
        for r, v in zip(chosen, vals):
            ax[1].axvspan(r['window_start_h'], r['window_end_h'], color='#304B53', alpha=.055)
            ax[2].annotate(str(int(v)), (r['trough_h'], v), xytext=(0, 9),
                           textcoords='offset points', ha='center', fontsize=10)
    ax[2].set_yticks(range(8)); ax[2].set_ylim(-.5, 7.7)
    ax[2].set_ylabel('Decoded value\n4 bit2 + 2 bit1 + bit0')
    ax[2].set_title('Finite-window reads; connecting steps are a visual guide', fontsize=10, loc='left')
    ax[3].plot(t[mask], sig['Int1_source'][mask], color=colors[1], lw=1.45,
                label='Int1 target production rate')
    ax[3].plot(t[mask], sig['Int2_source'][mask], color=colors[2], lw=1.5,
                label='Int2 target production rate')
    ax[3].set_ylabel('Carry drive\n(a.u./h)'); ax[3].set_xlabel('Time (h)')
    ax[3].legend(loc='upper left', bbox_to_anchor=(0, 1.19), ncol=2, frameon=False, fontsize=9)
    for i, a in enumerate(ax):
        a.text(-.075, 1.04, chr(65+i), transform=a.transAxes, weight='bold', fontsize=13)
        a.spines[['top', 'right']].set_visible(False)
        a.set_xlim(left, right)
    fig.supxlabel('Deterministic simulation | frozen scan configuration | readout is defined within windows',
                   fontsize=9, color='0.4')
    save(fig, name)

def digital(reads):
    values = [r['value'] for r in reads]
    bits = np.asarray([[int(r['bits'][i]['label']) for r in reads] for i in (2, 1, 0)])
    fig, ax = plt.subplots(figsize=(15, 3.7), layout='constrained')
    ax.imshow(bits, cmap='Blues', aspect='auto', vmin=0, vmax=1, interpolation='nearest')
    ax.set_yticks(range(3), ('bit2 (MSB)', 'bit1', 'bit0 (LSB)'))
    ax.set_xticks(range(0, len(reads), 2), range(0, len(reads), 2))
    ax.set_xlabel('Read index'); ax.set_title('56 valid read windows | 1 2 3 4 5 6 7 0 repeating', loc='left')
    for j, v in enumerate(values):
        ax.text(j, 2.92, str(v), ha='center', va='center', fontsize=8)
    ax.text(-1.2, 2.92, 'Value', ha='right', va='center', fontsize=9)
    ax.set_ylim(3.2, -.5)
    save(fig, '03_mod8_readout')

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p34 = json.loads(PROFILE34.read_text(encoding='utf-8'))
    p51 = json.loads(PROFILE51.read_text(encoding='utf-8'))
    f = p51['frozen_parameters']
    ext = Extension(**p34['extension'])
    c0 = CarryExpressionParameters(**p34['carry_expression'])
    c1 = CarryExpressionParameters(mrna_half_life_min=f['carry1_mrna_half_life_min'],
            activator_maturation_half_life_min=f['a1_f1_maturation_half_life_min'],
            repressor_maturation_half_life_min=f['a1_f1_maturation_half_life_min'])
    model = ThreeBit51Model(ext, ThreeBitCarryParameters(c0, c1), n_A1_gate=f['n_A1_gate'])
    sources = [ROOT/'model.py', ROOT/'model_twobit34.py', ROOT/'model_threebit51.py',
               ROOT/'verify_threebit51.py', ROOT/'verify_twobit_causal.py',
               PROFILE34, PROFILE51, UPSTREAM_PATH]
    hashes = {str(p): sha(p) for p in sources}
    cfg = dict(extension=asdict(ext), carry0=asdict(c0), carry1=asdict(c1),
               n_A1_gate=model.n_A1_gate_effective, F1_production_exponent=ZENG['n_A'][1],
               hours=600, sample_min=2, max_step_min=2, rtol=2e-7, atol=2e-9,
               initial_condition='PB DNA with pre-equilibrated regulatory pools (cold=False)',
               source_hashes=hashes)
    key = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
    cache = OUT/'trajectory.npz'
    metadata = OUT/'parameters.json'
    if cache.exists() and metadata.exists() and json.loads(metadata.read_text(encoding='utf-8'))['cache_key'] == key:
        from types import SimpleNamespace
        with np.load(cache) as z: sol = SimpleNamespace(t=z['t'], y=z['y'])
        print('Reusing hash-matched trajectory.', flush=True)
    else:
        print('Integrating frozen scan configuration: n_gate=6, clock K=0.3,n=2, 600 h.', flush=True)
        sol = model.simulate(hours=600, sample_min=2, max_step_min=2)
        if not np.isfinite(sol.y).all(): raise RuntimeError('Nonfinite state')
        np.savez_compressed(cache, t=sol.t, y=sol.y, names=np.array(model.state_names))
        metadata.write_text(json.dumps(dict(cache_key=key, **cfg), ensure_ascii=False, indent=2), encoding='utf-8')
    analysis = analyse_threebit(model, sol, cfg['hours'])
    expected = ''.join(str((i+1)%8) for i in range(56))
    assert analysis['cold_start']['sequence'] == expected, analysis['cold_start']
    assert analysis['cold_start']['passed'] and analysis['certified']
    assert hashes == {str(p):sha(p) for p in sources}, 'Source changed during run'
    (OUT/'verification.json').write_text(json.dumps(analysis, ensure_ascii=False, indent=2), encoding='utf-8')
    sig = model.diagnostic_signals(sol.y)
    flux = np.asarray([model.flux(y) for y in sol.y.T])
    data = {'time_h': sol.t, **dict(zip(model.state_names, sol.y)), 'C31_flux_uM_h': flux, **sig}
    pd.DataFrame(data).to_csv(OUT/'trajectories.csv', index=False)
    rows = []
    for r in analysis['read_windows']:
        rows.append(dict(time_h=r['trough_h'], window_start_h=r['window_start_h'],
                         window_end_h=r['window_end_h'], value=r['value'],
                         **{f'bit{i}':r['bits'][i]['label'] for i in range(3)}))
    pd.DataFrame(rows).to_csv(OUT/'read_windows.csv', index=False)
    plt.rcParams.update({'font.size':11, 'axes.labelsize':11, 'axes.titlesize':13,
                          'svg.fonttype':'none', 'pdf.fonttype':42})
    overview(sol,sig,flux,analysis['read_windows'],0,600,'01_threebit_overview_600h')
    overview(sol,sig,flux,analysis['read_windows'],90,190,'02_threebit_zoom_90_190h')
    digital(analysis['read_windows'])
    summary = dict(sequence=analysis['cold_start']['sequence'], certified=analysis['certified'],
                   reads=analysis['cold_start']['reads'], minimum_commitment=analysis['cold_start']['minimum_commitment'],
                   bit2_setup_hold=analysis['bit_margins']['bit2'],
                   actual_clock_gate=dict(K=ext.clock_K_au,n=ext.clock_n))
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    note = '''# 51状态成功模8计数图（2026-09-24）

这组图按已归档成功扫描的实际构造复算600 h，未改变模型、ZENG或冻结参数文件。

## 选择哪张图

- `01_threebit_overview_600h.png`：全程56个有效读窗，适合完整结果展示。
- `02_threebit_zoom_90_190h.png`：真实稳态片段，完整显示1→2→…→7→0，建议Wiki正文优先用这张。
- `03_mod8_readout.png`：三位数字条带，适合紧凑展示。

每张同时提供PDF与SVG；所有曲线来自真实ODE轨迹，整数阶梯仅连接读窗解码结果，不表示转变期间的连续数字输出。

## 本次确认的配置差异

成功扫描 `scan_gate_segments.build` / `check_read_commitment.build` 使用 `frozen_extension()`，从34状态selected_profile继承clock_K_au=0.3、clock_n=2。
51状态默认函数和曾同学参考代码中的0.4/3是另一套值。因此本次图使用实际已成功的0.3/2，不声称已复现0.4/3版本。

共同参数：uM_per_au=5.75；bit mRNA半衰期2 min、成熟半衰期20 min；carry0/1 mRNA半衰期2 min、A/F成熟半衰期32.5 min；独立n_A1_gate=6；F1产生指数4；add_growth=False。

## 可直接使用的图注

**中文：** 51状态三级计数器的确定性ODE轨迹。上游C31翻译通量驱动bit0，两级进位依次驱动bit1与bit2。600 h仿真产生56个有效读窗，按模8重复输出1–2–3–4–5–6–7–0。放大图取90–190 h真实片段；0.3/0.7虚线为读带边界，浅色竖带为周期20%的读窗。独立A1门指数为6；本图实际时钟门为K=0.3、n=2。结果为模型预测，尚未经湿实验验证。

**English:** Deterministic simulation of the 51-state three-bit counter. The C31 translation flux drives bit0, while two carry modules drive bit1 and bit2. The 600 h trajectory yields 56 valid finite-window reads following 1–2–3–4–5–6–7–0. The zoom shows the actual 90–190 h segment. The independent A1 gate exponent is 6; the clock gate used here has K=0.3 and n=2. Steps connect decoded reads only.

## 数据与复现

- `trajectory.npz`：全51状态、时间和状态名。
- `trajectories.csv`：全状态及物理信号。
- `read_windows.csv`：56个读窗、三个位和值。
- `parameters.json`：完整实际配置、求解器设置及源文件SHA256。
- `verification.json` / `summary.json`：实际验收结果。
- `make_wiki_threebit_n6_figures.py`（模型根目录）：生成脚本；仅复用参数与源文件哈希匹配的缓存。

C31 flux轴为µM/h；底部carry轴为下一位整合酶的目标产生速率a.u./h，经过现有转录与成熟链后才贡献成熟整合酶，不是游离Int浓度。
'''
    (OUT/'README_图注与参数.md').write_text(note, encoding='utf-8')
    manifest = {p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='SHA256SUMS.json'}
    manifest['generator']=sha(__file__)
    (OUT/'SHA256SUMS.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)

if __name__ == '__main__': main()
