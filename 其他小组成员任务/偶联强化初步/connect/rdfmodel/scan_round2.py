"""scan_round2.py — §7.2 决定性补充实验（2026-07-30 二轮）的可复现重建。

Round A (36 combos): uncorrected waveform (bsub=0),
    krep {15,30,60} × krdf {100,200,400} × tag {8,12,16,20}  -> all score=0
Round B (32 combos): tag=0, uncorrected waveform,
    krep {15,30,60,120} × krdf {100,150,200,300,400,500,600,800} -> all score=0

Claim: with A-module v36 waveform taken strictly as-is (Int baseline 0.28 µM),
NO counter-side modification (expression strength + degradation tag) rescues
alternation. Results saved to scan_round2.npz, table to scan_round2.log.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from multiprocessing import Pool
import zhao_core as zc
from analysis_lib import simulate, make_source_baseline, toggle_score
from couple_oscillator import pulse_times

_, TROUGHS = pulse_times()

ROUND_A = [(krep, krdf, tag)
           for krep in [15.0, 30.0, 60.0]
           for krdf in [100.0, 200.0, 400.0]
           for tag in [8.0, 12.0, 16.0, 20.0]]
ROUND_B = [(krep, krdf, 0.0)
           for krep in [15.0, 30.0, 60.0, 120.0]
           for krdf in [100.0, 150.0, 200.0, 300.0, 400.0, 500.0, 600.0, 800.0]]


def one(job):
    krep, krdf, tag = job
    P = zc.default_params()
    P.update(K_rep=0.0186, n_rep=3.4, krep_tsl=krep, krdf_tsl=krdf,
             k_tag_int=tag)
    src = make_source_baseline(P['k_dil'], baseline_sub=0.0, scale=1.0)
    y0 = zc.y0_PB(P, rep_mrna=P['k_tscr'] * P['Dtot'] / P['k_rna'],
                  rep=krep * P['k_tscr'] * P['Dtot'] / P['k_rna'] / P['k_dil'])
    sol = simulate(P, src, y0)
    LRf = zc.LR_total(sol.y.T) / P['Dtot']
    score, H, L, fid, _ = toggle_score(sol.t, LRf, TROUGHS)
    return krep, krdf, tag, score, H, L, fid


if __name__ == '__main__':
    jobs = ROUND_A + ROUND_B
    with Pool(12) as p:
        out = p.map(one, jobs)

    lines = []
    header = f"{'round':>5} {'krep':>5} {'krdf':>6} {'tag':>4} {'score':>7} {'H':>6} {'L':>6} {'fid':>5}"
    lines.append(header)
    for i, (krep, krdf, tag, score, H, L, fid) in enumerate(out):
        rnd = 'A' if i < len(ROUND_A) else 'B'
        Hs = f"{H:.2f}" if H == H else "  - "
        Ls = f"{L:.2f}" if L == L else "  - "
        lines.append(f"{rnd:>5} {krep:5.0f} {krdf:6.0f} {tag:4.0f} "
                     f"{score:7.3f} {Hs:>6} {Ls:>6} {fid:5.2f}")
    table = '\n'.join(lines)
    print(table)
    with open('scan_round2.log', 'w') as f:
        f.write(table + '\n')
    nA = sum(1 for i, r in enumerate(out) if i < len(ROUND_A) and r[3] == 0.0)
    nB = sum(1 for i, r in enumerate(out) if i >= len(ROUND_A) and r[3] == 0.0)
    print(f"\nRound A: {nA}/{len(ROUND_A)} score=0   Round B: {nB}/{len(ROUND_B)} score=0")
    np.savez('scan_round2.npz',
             out=np.array([[k, r, t, s] for k, r, t, s, H, L, f in out]),
             n_round_a=len(ROUND_A))
    print('saved scan_round2.npz / scan_round2.log')
