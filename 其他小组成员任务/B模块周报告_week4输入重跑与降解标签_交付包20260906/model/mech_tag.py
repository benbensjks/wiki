"""mech_tag.py — 阶段2.2: 降解标签机制定量(RBS0.45 工作点).
对 tag ∈ {0,4,8,12,16,20}: 输出
  - 间隙(通量谷)游离 Int 稳态 ≈ 谷通量/kd_int(解析) 与 模拟中值
  - 脉冲期 Int 峰(模拟)
  - 每周期 Int 产生剂量 = RBS·通量峰间积分
  - 总损失率 λ=k_dil+tag, 内禀 t½=ln2/tag, 总 t½=ln2/λ
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import zhao_core as zc
import couple_oscillator as co

RBS = 0.45
_, troughs = co.pulse_times(mode='flux', tail_frac=0.85)
troughs = troughs[(troughs>50)]
peaks, _ = co.pulse_times(mode='flux', tail_frac=0.85)
peaks = peaks[(peaks>50)]

def run(tag):
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=15.0, krdf_tsl=200.0,
             k_tag_int=tag)
    C = zc._rate_constants(P)
    src = co.make_flux_source(scale=RBS)
    y0 = zc.y0_PB(P, rep_mrna=P['k_tscr']*P['Dtot']/P['k_rna'],
                  rep=15.0*P['k_tscr']*P['Dtot']/P['k_rna']/P['k_dil'])
    t_eval = np.arange(0, 100, 0.02)
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src), (0, 100), y0,
                    method='LSODA', rtol=1e-7, atol=1e-12, max_step=0.1,
                    t_eval=t_eval)
    Y = sol.y.T
    itot = zc.int_total(Y)   # µM 总 Int 当量
    ifree = Y[:, 1]          # 游离 Int
    # 谷窗口 (谷±0.5h) 中值
    gap_int = []
    for tg in troughs:
        m = (sol.t > tg-0.5) & (sol.t < tg+0.5)
        gap_int.append(np.median(ifree[m]))
    pk_int = []
    for tp in peaks:
        m = (sol.t > tp-1.0) & (sol.t < tp+1.0)
        pk_int.append(np.max(ifree[m]))
    kd = P['k_dil'] + tag
    flux_trough = 0.036*RBS; flux_peak = 6.639*RBS
    return dict(tag=tag, kd=round(kd,2),
                t12_total_min=round(np.log(2)/kd*60,1),
                t12_tag_min=round(np.log(2)/max(tag,1e-9)*60,1),
                gap_ss_analytic=flux_trough/kd,
                gap_int_med=float(np.median(gap_int)),
                peak_int=float(np.max(pk_int)),
                peak_ss_analytic=flux_peak/kd)

print('tag | λ总 | 总t½ | 内禀t½ | 间隙Int 解析 | 间隙Int 模拟 | Int峰模拟 | 峰准稳态解析')
for tag in (0.0, 4.0, 8.0, 12.0, 16.0, 20.0):
    r = run(tag)
    print(f"{r['tag']:4.0f} | {r['kd']:5.2f} | {r['t12_total_min']:5.1f}min | "
          f"{r['t12_tag_min']:6.1f}min | {r['gap_ss_analytic']:.5f} | "
          f"{r['gap_int_med']:.5f} | {r['peak_int']:.4f} | {r['peak_ss_analytic']:.4f}")
