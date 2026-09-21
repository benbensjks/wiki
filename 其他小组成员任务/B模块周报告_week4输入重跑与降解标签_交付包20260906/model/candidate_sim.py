"""candidate_sim.py — 选型数值验证: 文献候选标签 → 模型双初态 100h 审计.
把每个候选标签的可实现 k_tag 区间(含悲观打折)代入, 报告 score 与失效边界距离.
固定 B 侧工作点 krep15/krdf200/K18.6nM/n3.4, RBS0.45(主)+RBS0.3/0.6(参考).
"""
import sys, os, json
sys.path.insert(0, '/home/steven/workspace/2026/IGEM/Week4_RDF/model')
import numpy as np
from scipy.integrate import solve_ivp
import zhao_core as zc
import couple_oscillator as co

def one(job):
    rbs, tag, init = job
    P = zc.default_params(); P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=15.0,
                                      krdf_tsl=200.0, k_tag_int=float(tag))
    C = zc._rate_constants(P)
    src = co.make_flux_source(float(rbs))
    if init == 'PB':
        y0 = zc.y0_PB(P, rep_mrna=P['k_tscr']*P['Dtot']/P['k_rna'],
                      rep=15.0*P['k_tscr']*P['Dtot']/P['k_rna']/P['k_dil'])
    else:
        y0 = np.zeros(38); y0[0] = P['Dtot']; y0[22] = 1.0
        y0[32] = 0.5; y0[36] = 1e-4; y0[37] = 1e-4
    t_eval = np.arange(0, 100, 0.02)
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src), (0, 100), y0,
                    method='LSODA', rtol=1e-7, atol=1e-12, max_step=0.1, t_eval=t_eval)
    Y = sol.y.T
    LRf = zc.LR_total(Y)/P['Dtot']
    _, tr = co.pulse_times(mode='flux', tail_frac=0.85)
    tr = tr[tr > 25]
    s = np.interp(tr, sol.t, LRf)
    s = s[2:]
    # 谷序列应为交替 1,0,1,0.. 或 0,1,0,1; score = 与理想序列匹配率
    exp = np.array([1 if i % 2 == 0 else 0 for i in range(len(s))]) if s[0] > 0.5 \
          else np.array([0 if i % 2 == 0 else 1 for i in range(len(s))])
    match = np.mean((s > 0.5) == exp)
    return {'rbs': rbs, 'tag': tag, 'init': init, 'score': round(float(match), 3),
            'seq': [round(float(x), 2) for x in s]}

# 候选 → 测试 k_tag 集(文献预期区间 + 悲观打折)
CAND = {
    'A1 原生ssrA (Szydlo 3-8 h⁻¹ 区间)':      [2.5, 3, 4, 5, 6, 8],
    'A2 LAA-LAA (>原生, 假设快2-4×=6-24)':   [6, 8, 12, 16, 20, 24],
    'A3 +ClpX共表达 (再上探)':                [12, 16, 20],
    'C SENY变体 (≈原生±)':                    [4, 6, 8],
    '悲观: Int融合比GFP慢3× (A1打折)':        [1.0, 1.5, 2.0, 2.7],
    '对照: 慢档 AAV (0.7-1.7 h⁻¹)':            [0.7, 1.0, 1.7],
}
jobs = []
for cand, tags in CAND.items():
    for tag in tags:
        for init in ('PB', 'LR'):
            jobs.append((0.45, tag, init, cand))

from multiprocessing import Pool
def wrap(j):
    rbs, tag, init, cand = j
    r = one((rbs, tag, init))
    r['cand'] = cand
    return r

if __name__ == '__main__':
    with Pool(12) as p:
        res = p.map(wrap, jobs)
    # 汇总: 每候选×tag 取双初态 min score
    agg = {}
    for r in res:
        k = (r['cand'], r['tag'])
        agg.setdefault(k, []).append(r['score'])
    out = []
    for (cand, tag), sc in sorted(agg.items()):
        out.append({'cand': cand, 'tag': tag, 'min_score': round(min(sc), 3),
                    'PB': sc[0] if len(sc) > 0 else None})
    # 按候选分组打印
    cur = None
    for o in out:
        if o['cand'] != cur:
            print(f"\n== {o['cand']} =="); cur = o['cand']
        print(f"  k_tag={o['tag']:>5}: 双初态min={o['min_score']}")
    json.dump(out, open('/home/steven/workspace/2026/IGEM/Week4_RDF/model/candidate_sim.json', 'w'), indent=1)
