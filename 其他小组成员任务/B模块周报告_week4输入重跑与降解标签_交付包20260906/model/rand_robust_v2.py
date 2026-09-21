"""rand_robust_v2.py — 蒙特卡洛(修订): tag {3,5,6,8,10,12,14,16,18,22,25,26}.
异质性=参数+输入扰动(与原版一致): RBS±0.04(N), K±15%, n±0.3, 相位均匀[0,1.04]周期,
幅度±8%; 30 实例/tag(双初态各半); 判据=双指标(P 主, M 副), 统计 min 序列的得分。
输出 robustness_mc_v2.json + 控制台表(mean±sd, 失败率 P<0.9 与 M<0.9)。
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import multiprocessing as mp
import zhao_core as zc, couple_oscillator as co
from twind import both_scores

def run_inst(job):
    tag, k = job
    rng = np.random.default_rng(hash((tag, k)) % (2**32))
    rbs = float(np.clip(0.45 + rng.normal(0, 0.04), 0.25, 0.65))
    K = 0.0186 * (1 + rng.normal(0, 0.15))
    n = 3.4 + rng.normal(0, 0.3)
    phase = rng.uniform(0, 1.04)
    amp = 1 + rng.normal(0, 0.08)
    init = 'PB' if rng.random() < 0.5 else 'LR'
    P = zc.default_params()
    P.update(K_rep=float(K), n_rep=float(n), krep_tsl=15.0,
             krdf_tsl=200.0, k_tag_int=float(tag))
    C = zc._rate_constants(P)
    t_h, flux, _ = co._load()
    f = rbs * amp * np.asarray(flux, dtype=float)
    sp = CubicSpline(t_h, f, bc_type='natural')
    _ph = phase * 10.59
    def src(t):
        return max(0.0, float(sp(np.asarray(t, dtype=float) - _ph)))
    if init == 'PB':
        y0 = zc.y0_PB(P, rep_mrna=P['k_tscr']*P['Dtot']/P['k_rna'],
                      rep=15.0*P['k_tscr']*P['Dtot']/P['k_rna']/P['k_dil'])
    else:
        y0 = np.zeros(38); y0[0] = P['Dtot']; y0[22] = 1.0
        y0[32] = 0.5; y0[36] = 1e-4; y0[37] = 1e-4
    _, tr0 = co.pulse_times(mode='flux', tail_frac=0.85)
    tr = tr0 + _ph
    tr = tr[(tr > 25) & (tr < 100.0)]
    t_pts = np.concatenate([tr, [25.0]])
    sol = solve_ivp(lambda t, y: zc.rhs(t, y, P, C, src), (0, 100), y0,
                    method='LSODA', rtol=1e-7, atol=1e-12, max_step=0.1,
                    t_eval=np.sort(np.unique(t_pts)))
    Y = sol.y.T
    LRf = zc.LR_total(Y) / P['Dtot']
    s = LRf[np.searchsorted(sol.t, tr)] if False else np.interp(tr, sol.t, LRf)[2:]
    if len(s) < 4:
        return {'tag': tag, 'init': init, 'P': None, 'M': None}
    sc = both_scores(s)
    return {'tag': tag, 'init': init, 'P': sc['P'], 'M': sc['M']}

if __name__ == '__main__':
    tags = [3, 5, 6, 8, 10, 12, 14, 16, 18, 22, 25, 26]
    n_inst = 30
    jobs = [(t, k) for t in tags for k in range(n_inst)]
    with mp.Pool(int(os.environ.get('NPROC', '12'))) as p:
        res = p.map(run_inst, jobs)
    out = {}
    print('tag | n | P: mean±sd  fail<0.9 | M: mean±sd  fail<0.9')
    per_inst = []
    for t in tags:
        vP = [r['P'] for r in res if r['tag'] == t and r['P'] is not None]
        vM = [r['M'] for r in res if r['tag'] == t and r['M'] is not None]
        mP, sP = float(np.mean(vP)), float(np.std(vP))
        mM, sM = float(np.mean(vM)), float(np.std(vM))
        fP = int(sum(1 for x in vP if x < 0.9))
        fM = int(sum(1 for x in vM if x < 0.9))
        out[t] = dict(n=len(vP), P_mean=round(mP, 3), P_sd=round(sP, 3), P_fail=fP,
                      M_mean=round(mM, 3), M_sd=round(sM, 3), M_fail=fM)
        per_inst += [{'tag': t, 'i': i, 'P': r['P'], 'M': r['M']}
                     for i, r in enumerate(res) if r['tag'] == t]
        print(f'{t:>4} | {len(vP):>2} | P {mP:.3f}±{sP:.3f} {fP}/{len(vP)} | '
              f'M {mM:.3f}±{sM:.3f} {fM}/{len(vM)}')
    json.dump({'stats': out, 'instances': per_inst},
              open('robustness_mc_v2.json', 'w'), indent=1)
