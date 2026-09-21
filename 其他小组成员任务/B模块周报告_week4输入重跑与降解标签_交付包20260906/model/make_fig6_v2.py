"""make_fig6_v2.py — 选型裕度 fig6(P 判据, RBS0.45 单线).
主曲线: scan_joint_v2.json 的 RBS0.45 列(tag 1..28, 双初态 P 最小值);
标记: 下界精扫 edge_scan.json 2.4-2.8(含初态分裂), 上缘 25/26 易碎点.
色带: 黄=原生 ssrA 3-8, 绿=LAA+4/SsrA2X 8-12(估), 红=LAA-LAA 14-38,
      橙=悲观 Int×1/3 1-2.7, 暗红=AAV 0.7-1.7(淘汰).
"""
import json, subprocess, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE = os.path.dirname(os.path.abspath(__file__))
joint = json.load(open(os.path.join(HERE, 'scan_joint_v2.json')))
edge = json.load(open(os.path.join(HERE, 'edge_scan.json')))
rbs45 = [r for r in joint if r['rbs'] == 0.45]
tags = sorted({int(r['tag']) for r in rbs45})
vals = {t: min(r['P'] for r in rbs45 if int(r['tag']) == t) for t in tags}
edges = sorted({round(float(r['tag']), 1) for r in edge})
evals = {u: min(r['P'] for r in edge if round(float(r['tag']), 1) == u) for u in edges}

for fam in ('FandolSong-Regular.otf', 'FandolHei-Regular.otf'):
    p = subprocess.run(['kpsewhich', fam], capture_output=True, text=True).stdout.strip()
    if p:
        font_manager.fontManager.addfont(p)
plt.rcParams.update({'font.size': 10.5, 'font.family': 'sans-serif',
                     'font.sans-serif': ['FandolSong', 'FandolHei', 'DejaVu Sans'],
                     'axes.unicode_minus': False})

fig, ax = plt.subplots(figsize=(8.6, 4.8))
ax.axvspan(0.7, 1.7, color='#9b59b6', alpha=0.10, label='AAV 慢档（淘汰）')
ax.axvspan(1.0, 2.7, color='#e67e22', alpha=0.10, label='悲观：Int 慢 3×')
ax.axvspan(3, 8, color='#f1c40f', alpha=0.18, label='原生 ssrA / LAA+4（3--8）')
ax.axvspan(14, 38, color='#e74c3c', alpha=0.16, label='LAA-LAA（SI 换算）')
tt = np.array(tags)
ax.plot(tt, [vals[t] for t in tt], 'o-', color='#15396e', ms=4, lw=1.6, label='RBS0.45 P 判据')
ev = sorted(evals)
ax.plot(ev, [evals[u] for u in ev], 's', color='#c0392b', ms=5, mew=0,
        label='下界/上缘精扫')
ax.axvline(2.8, color='#8e44ad', ls=':', lw=1.4)
ax.text(2.35, 0.36, '下缘 k≈2.8', fontsize=8, color='#8e44ad', va='center')
ax.axvline(24.0, color='#c0392b', ls=':', lw=1.3)
ax.text(24.0, 0.42, '24', fontsize=8, color='#c0392b')
ax.axvline(27.0, color='#c0392b', ls=':', lw=1.3)
ax.text(27.3, 0.42, '27（失效带起点）', fontsize=8, color='#c0392b')
ax.axhline(0.5, color='#7f8c8d', ls='--', lw=0.8, label='判据一半（判果分界）')
ax.set_xlabel('Int 降解标签速率 $k_{tag,int}$ (h$^{-1}$)')
ax.set_ylabel('双初态最小计数 score（P 判据）')
ax.set_xlim(0, 30)
ax.set_ylim(-0.03, 1.06)
ax.grid(alpha=0.25, ls=':')
ax.legend(fontsize=8, ncol=2, loc='lower left')
fig.tight_layout()
out = os.path.join(os.path.join(HERE, '..', 'figures'), 'fig6_candidate_margin.pdf')
fig.savefig(out)
print('saved', out)
print('RBS0.45 P-min:', {t: round(v, 3) for t, v in sorted(vals.items())})
