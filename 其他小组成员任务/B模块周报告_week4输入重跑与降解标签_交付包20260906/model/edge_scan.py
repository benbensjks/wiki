"""edge_scan.py — 上下界精扫: RBS0.45 × 双初态 × 指定 tag 列表, 双判据 + 中间态统计.
用法: python edge_scan.py --tags "2.0 2.2 ... 3.2"   或 --tags "24 25 26 27 28 29"
"""
import sys, os, json, argparse
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
    tag, init = job
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=15.0, krdf_tsl=200.0, k_tag_int=tag)
    src = co.make_flux_source(scale=0.45)
    t, LRf, cons = run_both(P, src, y0s(P)[init], t_pts=trough_pts())
    s = trough_samples(t, LRf)
    sc = both_scores(s)
    n_mid = int(np.sum((np.abs(s - 0.5) < 0.05) & (s > 0.4) & (s < 0.6)))
    return dict(tag=tag, init=init, **sc, n_mid=n_mid, cons=cons)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--tags', nargs='+', type=float, required=True)
    a = ap.parse_args()
    jobs = [(float(t), i) for t in a.tags for i in ('PB', 'LR')]
    with mp.Pool(int(os.environ.get('NPROC', '12'))) as p:
        res = p.map(one, jobs)
    for r in sorted(res, key=lambda r: (r['tag'], r['init'])):
        print(f"tag={r['tag']:5.2f} {r['init']}: P={r['P']} M={r['M']} fid={r['fid']} "
              f"H={r['H']} L={r['L']} n_mid={r['n_mid']} cons={r['cons']:.1e}")
    json.dump(res, open('edge_scan.json', 'w'), ensure_ascii=False, indent=1)
