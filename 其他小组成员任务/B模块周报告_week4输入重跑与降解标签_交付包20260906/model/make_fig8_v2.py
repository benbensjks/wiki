"""make_fig8_v2.py — 异质性蒙特卡洛 fig8(P/副判据, 12 档).
数据: robustness_mc_v2.json (30 实例/档, 修正相位采样).
柱状图: P 判据均值±SD + 失败实例数; 附 M 判据细线.
"""
import json, subprocess, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'robustness_mc_v2.json')))
st = d['stats']
tags = sorted(int(k) for k in st)
means = [st[str(t)]['P_mean'] for t in tags]
sds = [st[str(t)]['P_sd'] for t in tags]
fails = [st[str(t)]['P_fail'] for t in tags]
mmeans = [st[str(t)]['M_mean'] for t in tags]

for fam in ('FandolSong-Regular.otf', 'FandolHei-Regular.otf'):
    p = subprocess.run(['kpsewhich', fam], capture_output=True, text=True).stdout.strip()
    if p:
        font_manager.fontManager.addfont(p)
plt.rcParams.update({'font.size': 10.5, 'font.family': 'sans-serif',
                     'font.sans-serif': ['FandolSong', 'FandolHei', 'DejaVu Sans'],
                     'axes.unicode_minus': False})
fig, ax = plt.subplots(figsize=(8.2, 4.6))
colors = ['#e74c3c' if m < 0.9 else ('#f1c40f' if m < 0.98 else '#27ae60') for m in means]
ax.bar([str(t) for t in tags], means, yerr=sds, color=colors, capsize=5, width=0.55,
       label='P 判据')
ax.plot([str(t) for t in tags], mmeans, 'o--', color='#7f8c8d', ms=4, lw=1.0, label='M 判据')
for i, (m, f, sd) in enumerate(zip(means, fails, sds)):
    ax.text(i, m + sd + 0.02, f'{m:.2f}', ha='center', fontsize=9)
    ax.text(i, 0.02, f'稳 {30 - f}/30', ha='center', fontsize=8,
            color='white' if colors[i] == '#27ae60' else '#333')
ax.axhline(0.9, color='#15396e', ls='--', lw=1.2)
ax.text(0.2, 0.905, '阈值 0.9', fontsize=8, color='#15396e', va='bottom')
ax.set_ylim(0, 1.2)
ax.set_xlabel('Int 降解标签速率 $k_{tag,int}$ (h$^{-1}$)')
ax.set_ylabel('异质性下计数保真度（30 实例）')
ax.set_title('单细胞/异质性鲁棒性蒙特卡洛（P 判据）', pad=8)
ax.grid(axis='y', alpha=0.25)
ax.legend(fontsize=9)
fig.tight_layout()
out = os.path.join(FIGDIR := os.path.join(HERE, '..', 'figures'),
                   'fig8_stochastic_fidelity.pdf')
fig.savefig(out)
print('saved', out)
