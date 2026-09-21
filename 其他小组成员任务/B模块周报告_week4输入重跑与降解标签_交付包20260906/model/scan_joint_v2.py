"""scan_joint_v2.py — 联合工作区(修订): tag 1..28 × RBS 0.30..0.60(步进0.05) × 双初态.
主判据 P(banded toggle score), 副判据 M(>0.5 匹配率) 同时记录, 另记 n_mid/seq 供边界描述.
输出 scan_joint_v2.json; 打印: 每 tag 的最小可行 RBS(P≥0.95), 平台形貌, 低速段行为.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import multiprocessing as mp
import zhao_core as zc
import couple_oscillator as co
from twind import run_both, both_scores, trough_samples, trough_pts

def y0s(P):
    pb = zc.y0_PB(P, rep_mrna=P['k_tscr']*P['Dtot']/P['k_rna'],
                  rep=15.0*P['k_tscr']*P['Dtot']/P['k_rna']/P['k_dil'])
    lr = np.zeros(38); lr[0] = P['Dtot']; lr[22] = 1.0
    lr[32] = 0.5; lr[36] = 1e-4; lr[37] = 1e-4
    return {'PB': pb, 'LR': lr}

def one(job):
    rbs, tag, init = job
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=15.0, krdf_tsl=200.0, k_tag_int=float(tag))
    src = co.make_flux_source(scale=float(rbs))
    t, LRf, cons = run_both(P, src, y0s(P)[init], t_pts=trough_pts())
    s = trough_samples(t, LRf)
    sc = both_scores(s)
    return dict(rbs=rbs, tag=float(tag), init=init, **sc, cons=cons)

if __name__ == '__main__':
    tags = list(range(1, 29))
    rbss = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60]
    jobs = [(r, t, i) for r in rbss for t in tags for i in ('PB', 'LR')]
    print(len(jobs), 'runs', flush=True)
    with mp.Pool(int(os.environ.get('NPROC', '12'))) as p:
        res = p.map(one, jobs)
    json.dump(res, open('scan_joint_v2.json', 'w'), ensure_ascii=False, indent=1)
    agg = {}
    for r in res:
        k = (r['tag'], r['rbs'])
        agg.setdefault(k, []).append(r)
    def minP(tag, rbs):
        v = [r['P'] for r in agg[(tag, rbs)]]
        return min(v)
    print('\nP-metric 双初态最小矩阵:  tag\\rbs ' + ' '.join(f'{x:.2f}' for x in rbss))
    for tag in tags:
        row = ' '.join(f'{minP(tag, r):5.2f}' for r in rbss)
        print(f'tag={tag:>3}: {row}')
    print('\n=== 平台/边界结论 ===')
    for tag in tags:
        ok = [r for r in rbss if minP(tag, r) >= 0.95]
        maj = [r for r in rbss if 0.9 <= minP(tag, r) < 0.95]
        if not ok and not maj:
            print(f'tag={tag:>3}: 全域 <0.9')
        else:
            print(f'tag={tag:>3}: P>=0.95 {ok}  0.9-0.95 {maj}')
