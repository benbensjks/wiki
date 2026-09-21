"""verify_kdil2_control.py — 对照:把 k_dil 临时改回 2.0 (20min 倍增,原文 Zhao2019 语境),
复跑 verify_tetr FigS11A 与 verify_bm3r1_square 同款方波。
目的:证明 0.8318 下的行为差异纯属 k_dil 效应(机制移植本身无 bug)。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
import zhao_core as zc

def run_sq(name, K, n, krep_tsl=0.3, k_dil=2.0):
    P = zc.default_params()
    P.update(K_rep=K, n_rep=n, krep_tsl=krep_tsl, k_dil=k_dil)
    C = zc._rate_constants(P)
    src = lambda t: P['k_int'] * zc.square_pulse(t, P)
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src), (0, 24*5),
                    zc.y0_PB_ss(P), method='LSODA', rtol=1e-6, atol=1e-11,
                    max_step=0.05, t_eval=np.arange(0, 24*5, 0.02))
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    samples = [LRf[np.argmin(np.abs(sol.t - (k*24 + 12)))] for k in range(1, 6)]
    cons = np.abs(LRf + zc.PB_total(sol.y.T)/P['Dtot'] - 1.0).max()
    print(f'{name}: LR@mid-gap c1-5 = ' + ' '.join(f'{v:.2f}' for v in samples)
          + f' | DNA守恒 {cons:.1e}')

print('--- k_dil=2.0 对照(应接近 SummerWeek3 原值 0.20/0.86 交替档) ---')
run_sq('TetR(paper) k=0.01 n=2 ', 0.01, 2.0)
run_sq('BM3R1      k=0.0186 n=3.4', 0.0186, 3.4)
print('--- 附加:工作点参数(krep=15)下 k_dil 2.0 vs 0.8318 ---')
for kd in (2.0, 0.8318):
    run_sq(f'BM3R1 krep=15 k_dil={kd}', 0.0186, 3.4, krep_tsl=15.0, k_dil=kd)
