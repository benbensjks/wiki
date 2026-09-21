"""make_fig7_v2.py — 联合工作区 fig7(P 判据重扫版).
数据: scan_joint_v2.json (tag 1..28 × RBS 0.30..0.60, 双初态 P 判据最小值).
画: 热图 + 0.95 等高线 + 候选速率带(黄=原生 ssrA 3-8, 红=高速 14-38),
    平台区 tag 6-18 用白框标出, 标注 RBS0.45 处上缘与下缘.
"""
import json, subprocess, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, '..', 'figures')
rows = json.load(open(os.path.join(HERE, 'scan_joint_v2.json')))
tags = sorted({int(r['tag']) for r in rows})
rbss = sorted({r['rbs'] for r in rows})
Z = np.full((len(rbss), len(tags)), np.nan)
cell = {}
for r in rows:
    k = (r['rbs'], int(r['tag']))
    if k not in cell or r['P'] < cell[k]:
        cell[k] = r['P']
for (rbs, tag), v in cell.items():
    Z[rbss.index(rbs), tags.index(tag)] = v

for fam in ('FandolSong-Regular.otf', 'FandolHei-Regular.otf'):
    p = subprocess.run(['kpsewhich', fam], capture_output=True, text=True).stdout.strip()
    if p:
        font_manager.fontManager.addfont(p)
plt.rcParams.update({'font.size': 10.5, 'font.family': 'sans-serif',
                     'font.sans-serif': ['FandolSong', 'FandolHei', 'DejaVu Sans'],
                     'axes.unicode_minus': False})

fig, ax = plt.subplots(figsize=(9.2, 4.9))
im = ax.pcolormesh(tags, rbss, Z, cmap='RdYlGn', vmin=0, vmax=1.0, shading='auto')
cs = ax.contour(tags, rbss, Z, levels=[0.95], colors='#15396e', linewidths=1.5)
ax.clabel(cs, inline=True, fontsize=8, fmt='%.2f')
ax.axvspan(3, 8, color='#f1c40f', alpha=0.16)
ax.axvspan(14, 38, color='#e74c3c', alpha=0.16)
# 全平台框: tag 7-18 全域通过
ax.add_patch(plt.Rectangle((7 - 0.5, rbss[0] - 0.01), 12, rbss[-1] - rbss[0] + 0.02,
                           fc='none', ec='#15396e', lw=1.6, ls='--'))
ax.text(13, rbss[-1] + 0.02, '全稳平台 tag 7–18', ha='center', va='bottom', fontsize=9,
        color='#15396e')
ax.text(2.6, 0.545, '下缘\nk≈2.8', fontsize=8, color='#8e44ad', ha='right', va='center')
ax.axvline(2.8, color='#8e44ad', ls=':', lw=1.4)
ax.set_xlabel('Int 降解标签速率 $k_{tag,int}$ (h$^{-1}$)')
ax.set_ylabel('C31 RBS scale')
ax.set_xticks(np.arange(1, 29, 2))
ax.set_yticks(rbss)
ax.grid(alpha=0.2, ls=':')
cb = plt.colorbar(im, ax=ax, pad=0.02)
cb.set_label('双初态最小计数 score（P 判据）')
fig.tight_layout()
out = os.path.join(FIGDIR, 'fig7_rbs_tag_region.pdf')
fig.savefig(out)
print('saved', out)

print('=== 可行 RBS 范围(score>=0.95) ===')
for t in tags:
    ok = [r for r in rbss if Z[rbss.index(r), tags.index(t)] >= 0.95]
    print(f'  tag={t:>3}: RBS {min(ok):.2f}--{max(ok):.2f}' if ok else f'  tag={t:>3}: 无')
