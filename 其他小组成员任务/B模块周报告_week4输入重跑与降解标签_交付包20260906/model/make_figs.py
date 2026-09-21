"""make_figs.py — 阶段3: 生成报告用图 → ../figures/*.png (dpi 150)
fig1 输入波形(通量+蛋白) | fig2 RBS×tag 矩阵 | fig3 krep×krdf(tag8/12)
fig4 工作点全程轨迹 | fig5 标签机制(间隙/峰 Int vs tag)
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import zhao_core as zc
import couple_oscillator as co

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
os.makedirs(FIG, exist_ok=True)
# ---- 中文字体(Fandol, TeX Live 自带) + PDF 矢量输出 ----
from matplotlib import font_manager
for f in ('FandolSong-Regular.otf', 'FandolHei-Regular.otf', 'FandolFang-Regular.otf'):
    p = os.popen('kpsewhich ' + f).read().strip()
    if p:
        font_manager.fontManager.addfont(p)
plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 11,
    'font.sans-serif': ['FandolSong', 'FandolHei', 'DejaVu Sans'],
    'font.family': 'sans-serif',
    'axes.unicode_minus': False,
})
K = {'format': 'pdf'}

# ---------- fig1 输入波形 ----------
t_h, flux, c31 = co._load()
fig, ax = plt.subplots(2, 1, figsize=(7.2, 4.6), sharex=True)
ax[0].plot(t_h, flux, lw=0.8, color='#1a6f9c')
ax[0].set_ylabel('C31 通量 (μM/h)'); ax[0].set_title('week4 定稿 A 模块输出（unloaded，RBS=1 基准）')
ax[1].plot(t_h, c31/602, lw=0.8, color='#b5651d')
ax[1].set_ylabel('C31 蛋白 (μM)'); ax[1].set_xlabel('time (h)')
pk, tr = co.pulse_times(mode='flux', tail_frac=0.85)
ax[0].plot(pk, np.interp(pk, t_h, flux), 'rv', ms=3)
ax[0].plot(tr, np.interp(tr, t_h, flux), 'g^', ms=3)
fig.tight_layout(); fig.savefig(f'{FIG}/fig1_input.pdf'); plt.close(fig)

def min_both(rows):
    d = {r['init']: r for r in rows}
    return min(d['PB']['score'], d['LR']['score']) if len(d) == 2 else 0.0

# ---------- fig2 RBS×tag 细矩阵 (scan_tag_w4) ----------
out = json.load(open('scan_tag_w4.json'))
rbs_list = [0.3, 0.45, 0.6]; tag_list = list(range(0, 25, 2))
Z = np.array([[min_both([r for r in out if abs(r['rbs']-rb)<1e-9 and abs(r['tag']-tg)<1e-9])
               for tg in tag_list] for rb in rbs_list])
fig, ax = plt.subplots(figsize=(7.6, 2.9))
im = ax.imshow(Z, aspect='auto', cmap='viridis', vmin=0, vmax=1,
               extent=[tag_list[0]-1, tag_list[-1]+1, len(rbs_list)-0.5, -0.5])
ax.set_yticks(range(len(rbs_list))); ax.set_yticklabels([f'{r:.2f}' for r in rbs_list])
ax.set_ylabel('C31 RBS scale'); ax.set_xlabel('Int 降解标签 $k_{tag,int}$ (h$^{-1}$)')
ax.set_title('B 模块双初态最小 toggle score（krep=15, krdf=200, K=18.6nM, n=3.4）')
for i in range(len(rbs_list)):
    for j in range(len(tag_list)):
        ax.text(tag_list[j], i, f'{Z[i,j]:.2f}', ha='center', va='center',
                color='w' if Z[i,j] < 0.55 else 'k', fontsize=7)
ax.axvline(0.5, color='red', lw=1.5, ls='--')
fig.colorbar(im, label='score'); fig.tight_layout()
fig.savefig(f'{FIG}/fig2_RBS_tag_matrix.pdf'); plt.close(fig)

# ---------- fig3 krep×krdf (轮1) ----------
out = json.load(open('scan_workpoint.json'))
krep_list = [6, 10, 15, 20, 30]; krdf_list = [50, 100, 200, 400]
fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.2), sharey=True)
for ax, tag in zip(axes, (8.0, 12.0)):
    Z = np.array([[min_both([r for r in out if r['rbs']==0.45 and abs(r['tag']-tag)<1e-9
                             and r['krep']==k and r['krdf']==kd])
                   for kd in krdf_list] for k in krep_list])
    im = ax.imshow(Z, aspect='auto', cmap='viridis', vmin=0, vmax=1,
                   extent=[0, len(krdf_list)-1, len(krep_list)-0.5, -0.5])
    ax.set_xticks(range(len(krdf_list))); ax.set_xticklabels(krdf_list)
    ax.set_yticks(range(len(krep_list))); ax.set_yticklabels(krep_list)
    ax.set_xlabel('$k_{rdf,tsl}$ (h$^{-1}$)'); ax.set_title(f'RBS0.45, tag={tag:g}')
    for i in range(len(krep_list)):
        for j in range(len(krdf_list)):
            ax.text(j, i, f'{Z[i,j]:.2f}', ha='center', va='center',
                    color='w' if Z[i,j] < 0.55 else 'k', fontsize=8)
axes[0].set_ylabel('$k_{BM3,tsl}$ (h$^{-1}$)')
fig.suptitle('轮1: krep × krdf（双初态最小 score，RBS=0.45）', y=1.0, fontsize=10)
fig.tight_layout(); fig.savefig(f'{FIG}/fig3_krep_krdf.pdf'); plt.close(fig)

# ---------- fig4 工作点轨迹 ----------
P = zc.default_params(); P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=15.0,
                                  krdf_tsl=200.0, k_tag_int=12.0)
C = zc._rate_constants(P); src = co.make_flux_source(scale=0.45)
y0 = zc.y0_PB(P, rep_mrna=P['k_tscr']*P['Dtot']/P['k_rna'],
              rep=15.0*P['k_tscr']*P['Dtot']/P['k_rna']/P['k_dil'])
t_eval = np.arange(0, 100, 0.02)
sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src), (0, 100), y0,
                method='LSODA', rtol=1e-7, atol=1e-12, max_step=0.1, t_eval=t_eval)
Y = sol.y.T
fig, axes = plt.subplots(4, 1, figsize=(7.6, 7.4), sharex=True)
axes[0].plot(sol.t, np.maximum(0, co.load_flux_spline(0.45)(sol.t)), lw=0.6, color='#1a6f9c')
axes[0].set_ylabel('Int 产生通量\nRBS0.45×flux (μM/h)')
axes[1].plot(sol.t, zc.LR_total(Y)/P['Dtot'], lw=0.9, color='#9c1a2f')
axes[1].axhline(0.5, color='gray', lw=0.5, ls=':')
axes[1].set_ylabel('LR 分数 (DNA 状态)'); axes[1].set_ylim(-0.05, 1.05)
# 面板3: Int 总 —— 独立轴(峰 ~0.23 μM, 相对 RDF/BM3R1 小两~三个数量级, 共享轴会贴 0)
axes[2].plot(sol.t, zc.int_total(Y), lw=0.7, color='#2f7d32')
axes[2].set_ylabel('Int 总 (μM)'); axes[2].set_ylim(0, 0.30)
# 面板4: RDF 总(左轴) 与 BM3R1(右轴) 尺度差 ~10 倍 → 双轴
ax4 = axes[3]
l1, = ax4.plot(sol.t, zc.rdf_total(Y), lw=0.7, color='#7d2f9c', label='RDF 总 (μM, 左轴)')
ax4.set_ylabel('RDF 总 (μM)', color='#7d2f9c'); ax4.tick_params(axis='y', colors='#7d2f9c')
ax4r = ax4.twinx()
l2, = ax4r.plot(sol.t, Y[:, 37], lw=0.7, color='#b5651d', label='BM3R1 (μM, 右轴)')
ax4r.set_ylabel('BM3R1 (μM)', color='#b5651d'); ax4r.tick_params(axis='y', colors='#b5651d')
ax4.legend(handles=[l1, l2], fontsize=7, loc='upper right')
ax4.set_xlabel('time (h)')
fig.suptitle('工作点: RBS0.45 / tag=12 / krep=15 / krdf=200（PB 初态）', y=0.995, fontsize=10)
fig.tight_layout(); fig.savefig(f'{FIG}/fig4_operating_point.pdf'); plt.close(fig)

# ---------- fig5 标签机制 ----------
tags = np.array([0, 4, 8, 12, 16, 20.0])
gap = np.array([0.05605, 0.00099, 0.00042, 0.00027, 0.0002, 0.00018])
peak = np.array([0.6090, 0.2234, 0.1472, 0.1140, 0.0953, 0.0829])
fig, ax = plt.subplots(figsize=(5.6, 3.4))
ax.semilogy(tags, gap, 'o-', color='#9c1a2f', label='间隙游离 Int (模拟中值)')
ax.semilogy(tags, peak, 's-', color='#1a6f9c', label='脉冲 Int 峰 (模拟)')
ax.axhline(0.1, color='gray', ls=':', lw=1)
ax.text(0.3, 0.11, '翻转阈值 ~0.1 μM', fontsize=7, color='gray')
ax.axvspan(4, 20, color='#2f7d32', alpha=0.12)
ax.text(7.5, 6e-3, 'RBS0.45 工作窗口', fontsize=7, color='#2f7d32')
ax.set_xlabel('$k_{tag,int}$ (h$^{-1}$)'); ax.set_ylabel('Int 浓度 (μM, log)')
ax.legend(fontsize=7); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(f'{FIG}/fig5_tag_mechanism.pdf'); plt.close(fig)
print('figures 已生成:', sorted(os.listdir(FIG)))
