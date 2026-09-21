"""exp1_upper_diag.py — 上界非单调诊断: RBS0.45, tag ∈ {24,26,28,30,32,36}.
对每个 tag 打印双判据、谷序、以及每个通量峰窗口内的翻转沿, 以判定
“超高通量降速后系统发生了什么”。假设失败应从“峰 Int 削到阈值以下 → 不再翻转”
单调延伸到全体; 若出现部分复原, 需查明是“隔脉冲翻转”(period doubling)
还是“翻转发生但采样相位偏移”。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import zhao_core as zc
import couple_oscillator as co
from twind import run_both, both_scores, trough_samples

def y0s(P):
    pb = zc.y0_PB(P, rep_mrna=P['k_tscr']*P['Dtot']/P['k_rna'],
                  rep=15.0*P['k_tscr']*P['Dtot']/P['k_rna']/P['k_dil'])
    lr = np.zeros(38); lr[0] = P['Dtot']; lr[22] = 1.0
    lr[32] = 0.5; lr[36] = 1e-4; lr[37] = 1e-4
    return {'PB': pb, 'LR': lr}

def mk(tag):
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=15.0, krdf_tsl=200.0, k_tag_int=tag)
    return P

import multiprocessing as mp

def one(job):
    tag, init = job
    P = mk(tag)
    src = co.make_flux_source(scale=0.45)
    t, LRf, cons = run_both(P, src, y0s(P)[init], dt=0.05)
    s = trough_samples(t, LRf)
    sc = both_scores(s)
    # 峰窗口内 LR 过 0.5 的沿(相对峰时刻), 区分"每次脉冲翻转"与"隔次翻转"
    _, peaks = co.pulse_times(mode='flux', tail_frac=0.85)
    peaks = peaks[(peaks > 25) & (peaks < t[-1])]
    edges = []
    LRfull = LRf[np.argsort(t)] if False else None
    ts = t
    for p in peaks:
        seg = (ts >= p - 1.2) & (ts <= p + 1.2)
        lr = LRf[seg]; tt = ts[seg]
        up = tt[np.where((lr[:-1] < 0.5) & (lr[1:] >= 0.5))[0]]
        dn = tt[np.where((lr[:-1] > 0.5) & (lr[1:] <= 0.5))[0]]
        edges.append((round(float(p), 2), [round(float(x - p), 2) for x in up],
                      [round(float(x - p), 2) for x in dn]))
    return dict(tag=tag, init=init, **sc, cons=cons, edges=edges)

if __name__ == '__main__':
    jobs = [(t, i) for t in (24.0, 26.0, 28.0, 30.0, 32.0, 36.0) for i in ('PB', 'LR')]
    with mp.Pool(12) as p:
        res = p.map(one, jobs)
    for r in sorted(res, key=lambda r: (r['tag'], r['init'])):
        print(f"tag={r['tag']:4.0f} {r['init']}: P={r['P']} M={r['M']} fid={r['fid']} "
              f"H={r['H']} L={r['L']} seq={r['seq']} cons={r['cons']:.1e}")
        print(f"    edges(pk, up@pk+, dn@pk+): {r['edges']}")
