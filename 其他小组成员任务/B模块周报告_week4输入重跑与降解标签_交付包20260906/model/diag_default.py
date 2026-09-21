"""diag_default.py — 阶段1.1: 默认/未调参组合接 week4 定稿波形,失效诊断.
情景(均 k_dil=0.8318 定稿):
  A. 原文 TetR 参数  (K=0.01, n=2, krep_tsl=0.3, krdf_tsl=4.0, tag=0, RBS=1.0)
  B. BM3R1 未调参    (K=0.0186, n=3.4, krep_tsl=0.3, krdf_tsl=4.0, tag=0, RBS=1.0)
  C. BM3R1 弱延迟     (同上但 krep_tsl=0.3, krdf_tsl=4.0, tag=0, RBS=0.45)
双初态 PB / LR,100h;报告每个通量谷(脉冲间隙)的 LR 分数序列与终态。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import zhao_core as zc
import couple_oscillator as co

T_END = 100.0
_, troughs = co.pulse_times(mode='flux')
troughs = troughs[(troughs > 30) & (troughs < T_END)]

def run(label, P_upd, rbs, y0):
    P = zc.default_params()
    P.update(P_upd)
    C = zc._rate_constants(P)
    src = co.make_flux_source(scale=rbs)
    t_eval = np.arange(0, T_END, 0.05)
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src),
                    (0, T_END), y0(P), method='LSODA',
                    rtol=1e-7, atol=1e-12, max_step=0.1, t_eval=t_eval)
    Y = sol.y.T
    LRf = zc.LR_total(Y) / P['Dtot']
    cons = np.abs(LRf + zc.PB_total(Y)/P['Dtot'] - 1.0).max()
    samp = np.interp(troughs, sol.t, LRf)
    seq = ' '.join(f'{v:.2f}' for v in samp)
    print(f'{label}: LR@谷={seq} | 终态={LRf[-1]:.2f} | 守恒{cons:.1e}')

def y0_PB(P):
    return zc.y0_PB(P, rep_mrna=P['k_tscr']*P['Dtot']/P['k_rna'],
                    rep=P['krep_tsl']*P['k_tscr']*P['Dtot']/P['k_rna']/P['k_dil'])
def y0_LR(P):
    y = np.zeros(38); y[0] = P['Dtot']
    # LR 稳态:rdf/rdf_mrna 取典型;rep 蛋白随稀释
    y[22] = 1.0; y[32] = 0.5
    y[36] = 0.001; y[37] = 0.001
    return y

print('情景A 原文TetR 默认 (RBS1.0)')
run('  PB初态', dict(krep_tsl=0.3, krdf_tsl=4.0, k_tag_int=0.0), 1.0, y0_PB)
run('  LR初态', dict(krep_tsl=0.3, krdf_tsl=4.0, k_tag_int=0.0), 1.0, y0_LR)
print('情景B BM3R1未调参 (RBS1.0)')
run('  PB初态', dict(K_rep=0.0186, n_rep=3.4, krep_tsl=0.3, krdf_tsl=4.0, k_tag_int=0.0), 1.0, y0_PB)
run('  LR初态', dict(K_rep=0.0186, n_rep=3.4, krep_tsl=0.3, krdf_tsl=4.0, k_tag_int=0.0), 1.0, y0_LR)
print('情景C BM3R1弱延迟 (RBS0.45)')
run('  PB初态', dict(K_rep=0.0186, n_rep=3.4, krep_tsl=0.3, krdf_tsl=4.0, k_tag_int=0.0), 0.45, y0_PB)
run('  LR初态', dict(K_rep=0.0186, n_rep=3.4, krep_tsl=0.3, krdf_tsl=4.0, k_tag_int=0.0), 0.45, y0_LR)
