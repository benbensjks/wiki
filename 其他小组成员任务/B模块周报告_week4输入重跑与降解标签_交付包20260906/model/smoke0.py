"""smoke0: 新输入(k_dil=0.8318 + week4 unloaded 通量)下 B 模块探针.
基准情景 = 韩亚轩 v53g 推荐 RBS 0.45 / tag 8;对照 tag 12 (用户旧工作点标签).
非正式扫描,只验证管线可跑 + DNA 守恒 + 翻转行为初观.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import zhao_core as zc
import couple_oscillator as co

T_END = 100.0
print('zhao_core k_dil =', zc.default_params()['k_dil'])
print('波形摘要:', co.waveform_summary())

def run(rbs, tag, y0_fn=zc.y0_PB):
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=15.0, krdf_tsl=200.0,
             k_tag_int=tag)
    C = zc._rate_constants(P)
    src = co.make_flux_source(scale=rbs)
    y0 = y0_fn(P, rep_mrna=P['k_tscr']*P['Dtot']/P['k_rna'],
               rep=15.0*P['k_tscr']*P['Dtot']/P['k_rna']/P['k_dil'])
    t_eval = np.arange(0, T_END, 0.05)
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src),
                    (0, T_END), y0, method='LSODA',
                    rtol=1e-7, atol=1e-12, max_step=0.1, t_eval=t_eval)
    assert sol.success, sol.message
    Y = sol.y.T
    LRf = zc.LR_total(Y) / P['Dtot']
    cons = np.abs(LRf + zc.PB_total(Y)/P['Dtot'] - 1.0)
    return sol.t, LRf, cons.max()

for tag in (8.0, 12.0):
    t, LRf, cons = run(0.45, tag)
    _, troughs = co.pulse_times(mode='flux')
    troughs = troughs[(troughs > 30) & (troughs < T_END)]
    samp = np.interp(troughs, t, LRf)
    alt = ' '.join(f'{v:.2f}' for v in samp)
    print(f'RBS0.45 tag={tag:g}: DNA守恒max|Δ|={cons:.2e} | LR@通量谷(30h后): {alt}')
