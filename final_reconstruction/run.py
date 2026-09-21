"""Run reconstruction, save all states/fluxes and report measured performance."""
import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from model import Model, Extension, ZENG, STATE_NAMES, ROOT, UPSTREAM_PATH, bit_slice


def run(hours=150.0, cold=False, output=None, strict=False):
    out = Path(output) if output else ROOT/'results'
    out.mkdir(parents=True, exist_ok=True)
    m = Model()
    sol = m.simulate(hours, cold=cold, strict=strict)
    df = pd.DataFrame(sol.y.T, columns=STATE_NAMES)
    df.insert(0, 'time_h', sol.t)
    df['C31_translation_flux_uM_h'] = [m.flux(y) for y in sol.y.T]
    gates = np.array([m.carry_promoters(y) for y in sol.y.T])
    df['carry_promoter_01'] = gates[:, 0]
    df['carry_promoter_12'] = gates[:, 1]
    df['clock_AND_activity'] = gates[:, 2]
    for i in range(3):
        df[f'b{i}_mature_appearance_au_h'] = m.kmat*df[f'b{i}_I_u']
        # These are observations ONLY, never fed into the ODE.
        S = df[f'b{i}_S'].to_numpy()
        df[f'b{i}_decoded'] = np.where(S >= .7, 1, np.where(S <= .3, 0, -1))
    peak, _ = find_peaks(df.C31_translation_flux_uM_h, prominence=1.0, distance=5*60)
    # Sample at source-flux minimum between consecutive peaks: fixed criterion.
    rows = []
    for n, (left, right) in enumerate(zip(peak[:-1], peak[1:]), 1):
        k = left+int(np.argmin(df.C31_translation_flux_uM_h.iloc[left:right].to_numpy()))
        states = [float(df[f'b{i}_S'].iloc[k]) for i in range(3)]
        decoded = [int(df[f'b{i}_decoded'].iloc[k]) for i in range(3)]
        code = sum(q*2**i for i,q in enumerate(decoded)) if min(decoded) >= 0 else -1
        rows.append(dict(pulse=n, peak_h=float(sol.t[left]), sample_h=float(sol.t[k]),
                         S0=states[0], S1=states[1], S2=states[2], code=code,
                         expected_from_zero=n%8, confident=min(decoded)>=0))
    samples = pd.DataFrame(rows)
    observed = [r['code'] for r in rows]
    valid = np.array([r['confident'] for r in rows], bool)
    increments = [bool(b >= 0 and a >= 0 and b == (a+1)%8) for a,b in zip(observed[:-1], observed[1:])]
    peak_rows = []
    for name in ('C31_translation_flux_uM_h', 'b1_mature_appearance_au_h', 'b2_mature_appearance_au_h'):
        values = df[name].to_numpy()
        ids, _ = find_peaks(values, prominence=max(.05*np.ptp(values), 1e-6), distance=60)
        for k in ids:
            peak_rows.append(dict(signal=name, time_h=float(sol.t[k]), value=float(values[k])))
    summary = dict(n_states=len(STATE_NAMES), total_hours=hours, initial_condition='cold' if cold else 'PB regulators pre-equilibrated',
                   finite=bool(np.isfinite(sol.y).all()), minimum_state=float(sol.y.min()),
                   S_ranges={f'b{i}':[float(df[f'b{i}_S'].min()),float(df[f'b{i}_S'].max())] for i in range(3)},
                   upstream_period_h=float(np.median(np.diff(sol.t[peak][sol.t[peak]>25]))) if sum(sol.t[peak]>25)>1 else None,
                   observed_codes=observed, expected_codes=[r['expected_from_zero'] for r in rows],
                   confident_samples=int(valid.sum()), sample_count=len(rows),
                   sequence_accuracy_from_zero=float(np.mean([r['code']==r['expected_from_zero'] for r in rows])) if rows else None,
                   increment_accuracy=float(np.mean(increments)) if increments else None,
                   binary_counter_pass=bool(len(rows)>=8 and all(r['code']==r['expected_from_zero'] for r in rows)),
                   note='Numerical/structural validation is distinct from binary-counter performance. No parameter fitting performed.')
    df.to_csv(out/'trajectories.csv', index=False)
    samples.to_csv(out/'clock_samples.csv', index=False)
    pd.DataFrame(peak_rows).to_csv(out/'pulse_peaks.csv', index=False)
    (out/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    params = dict(upstream=asdict(m.a), zeng=ZENG, extension_assumptions=asdict(m.e),
                  derived=dict(mu_h=m.mu, mRNA_loss_h=m.lm, maturation_h=m.kmat,
                               complex_QSS_factor=m.q, K_complex_au=m.K_complex),
                  source=dict(path=str(UPSTREAM_PATH), sha256=hashlib.sha256(UPSTREAM_PATH.read_bytes()).hexdigest()),
                  units=dict(time='h', oscillator='copies/cell', downstream='a.u.', input_flux='uM/h'),
                  unused_source_parameter='k_int=6: square-input synthesis replaced by mechanistic upstream translation',
                  state_names=STATE_NAMES)
    (out/'parameters.json').write_text(json.dumps(params,ensure_ascii=False,indent=2),encoding='utf-8')
    plt.rcParams.update({'font.size':10})
    fig, ax = plt.subplots(5, 1, figsize=(12, 13), sharex=True, constrained_layout=True)
    t = sol.t
    ax[0].plot(t, df.C31_translation_flux_uM_h, c='black', label='Han upstream C31 translation')
    ax[0].set_ylabel('uM / h'); ax[0].legend(loc='upper right')
    for j,c in enumerate(('#0072B2','#D55E00','#009E73')):
        ax[1].plot(t, df[f'b{j}_S'], c=c, label=f'Bit {j}: LR fraction')
        ax[3].plot(t, df[f'b{j}_I'], c=c, label=f'Int {j} (free mature)')
    ax[1].set_ylim(-.05,1.05); ax[1].set_ylabel('LR fraction'); ax[1].legend(ncol=3)
    for n,style in [('A0','-'),('F0','--'),('A1','-'),('F1','--')]:
        ax[2].plot(t, df[n], style, label=n+' (carry R)' if n.startswith('F') else n)
    ax[2].set_ylabel('a.u.'); ax[2].legend(ncol=4)
    ax[3].set_ylabel('a.u.'); ax[3].legend(ncol=3)
    ax[4].plot(t, gates[:,0], label='A0 activation x F0 repression')
    ax[4].plot(t, gates[:,1], label='A1 activation x F1 repression x clock')
    ax[4].plot(t, gates[:,2], ':', label='clock gate')
    ax[4].set_ylabel('promoter activity'); ax[4].set_xlabel('Time (h)'); ax[4].legend(ncol=2)
    for a in ax: a.grid(alpha=.2)
    fig.suptitle('43-state reconstruction | '+('binary count PASS' if summary['binary_counter_pass'] else 'binary count NOT validated'))
    fig.savefig(out/'overview.png',dpi=160); plt.close(fig)
    fig, axes = plt.subplots(11,3,figsize=(14,22),sharex=True,constrained_layout=True)
    for i in range(3):
        for j,s in enumerate(STATE_NAMES[6:17]):
            name = s.replace('b0_',f'b{i}_')
            axes[j,i].plot(t,df[name],lw=.8)
            axes[j,i].set_ylabel(name)
            axes[j,i].grid(alpha=.2)
        axes[-1,i].set_xlabel('Time (h)')
    fig.savefig(out/'all_33_states.png',dpi=120); plt.close(fig)
    (out/'结果说明.md').write_text(
        '# 本次实际运行结果\n\n'
        f'- 总状态：43（6 个振荡器状态 + 33 个 bit 状态 + 4 个 I1-FFL 状态）。\n'
        f'- 上游稳态周期：{summary["upstream_period_h"]} h。\n'
        f'- 逐拍状态码：{observed}；-1 表示至少一位在 0.3–0.7 中间区。\n'
        f'- 从初态 000 起顺序正确率：{summary["sequence_accuracy_from_zero"]}。\n'
        f'- 相邻采样模 8 加一正确率：{summary["increment_accuracy"]}。\n'
        f'- 完整二进制计数通过：{summary["binary_counter_pass"]}。\n\n'
        '结构实现完成不代表功能自动成功；真实宽脉冲、显式复合物、转录成熟延迟和生长稀释改变了原粗粒化模型的动态。'
        '本次保留曾墨涵表内数值，没有为获得理想计数轨迹重新拟合。详见上一级 README 的参数映射与限制。\n',encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--hours',type=float,default=150)
    parser.add_argument('--cold',action='store_true')
    parser.add_argument('--strict',action='store_true')
    parser.add_argument('--output',type=Path,default=ROOT/'results')
    args = parser.parse_args()
    run(args.hours,args.cold,args.output,args.strict)
