"""Run matched Han-driven HZZ A1-autoregulation OFF/ON and an HZH control."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import csv
import json
import hashlib
from dataclasses import asdict
from datetime import datetime
import platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from model_hzz import ROOT, PARENT, NAMES, IDX, HZZModel, HZHModel, TAIL_SOURCE, structural_checks
sys.path.insert(0, str(PARENT / 'hby_zmh_hby' / 'certification'))
from run_han_comparison import HanInput, HAN_PATH
from hybrid_model import sha256, source_hashes
import verify_hzh as V

HOURS = 300.0
SOLVER = dict(method='DOP853', rtol=2e-7, atol=2e-9, max_step=1/60)


def dump(path, x):
    path.write_text(json.dumps(x, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def integrate(model):
    t = np.linspace(0, HOURS, int(HOURS*60)+1)
    sol = solve_ivp(model.rhs, (0, HOURS), model.initial_state(), t_eval=t, **SOLVER)
    if not sol.success or not np.isfinite(sol.y).all():
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def evaluate(t, sig):
    peaks, reads = V.read_windows(t, sig)
    crosses = {b: V.crossings(t, sig[b]) for b in ('S0','S1','S2')}
    events = {}
    for j in (0,1):
        rev = V.segments(t, sig[f'J_rev{j}'], V.JREV_THRESHOLD_H)
        gate = V.segments(t, sig[f'g{j}'], V.GATE_THRESHOLD, min_dose=V.GATE_MIN_DOSE_H)
        events[f'bit{j}_to_bit{j+1}'] = V.associate(rev, gate, crosses[f'S{j+1}'], HOURS)
    steady = V.read_verdict(reads, V.DROP)
    event_pass = all(x['passed'] for x in events.values())
    # Include unlabelled windows, which must not disappear from per-bit metrics.
    perbit = {}
    for i, b in enumerate(('S0','S1','S2')):
        selected = reads[V.DROP:]
        perbit[b] = dict(minimum_commitment=min(r['commitment'][i] for r in selected),
                         unlabelled=sum(r['labels'][i] is None for r in selected),
                         minimum=float(sig[b].min()), maximum=float(sig[b].max()),
                         timing=V.margins(selected, crosses[b]))
    return dict(version='HZZ_EXPLORATORY_V1', certified=None,
                counting_passed=steady['passed'], event_causality_passed=event_pass,
                diagnostic_passed=bool(steady['passed'] and event_pass),
                cold=V.read_verdict(reads,0), steady=steady, per_bit=perbit,
                crossings=crosses, events=events, read_windows=reads, clock_peaks=len(peaks),
                rule=dict(clock='own bit0 Int0 trough',window_total_fraction=.2,
                          bands=[V.LOW,V.HIGH],occupancy=V.OCC,drop=V.DROP,
                          min_steady_reads=V.MIN_STEADY,reverse_threshold_per_h=V.JREV_THRESHOLD_H,
                          gate_threshold=V.GATE_THRESHOLD,min_event_h=V.PULSE_DURATION_H,
                          min_gate_dose_h=V.GATE_MIN_DOSE_H,drop_carry_events=2),
                scope='New 23-state diagnostic: same scalar read/event rules, no inherited certification; no threshold tuning')


def draw(t, y, sig, result, target, title):
    fig, ax = plt.subplots(5,1,figsize=(14,12),sharex=True,layout='constrained')
    ax[0].plot(t,sig['int0'],color='#7550a1'); ax[0].set_ylabel('Int0\n(a.u.)')
    for b,c in zip(('S0','S1','S2'),('#2685ad','#d97726','#43875a')):
        ax[1].plot(t,sig[b],label=b,color=c,lw=1.4)
    for a in (ax[1],):
        a.axhspan(0,.3,color='#43875a',alpha=.08);a.axhspan(.7,1,color='#2685ad',alpha=.08)
        a.axhline(.3,ls='--',color='gray',lw=.7);a.axhline(.7,ls='--',color='gray',lw=.7)
        a.set_ylim(-.04,1.04);a.set_ylabel('DNA LR fraction');a.legend(loc='upper right',ncol=3)
    reads=result['read_windows']; rt=np.array([r['trough_h'] for r in reads])
    values=np.array([np.nan if r['value'] is None else r['value'] for r in reads])
    ax[2].step(rt,values,where='post',color='#33383a');ax[2].scatter(rt,values,c=values,cmap='viridis',vmin=0,vmax=7,s=25)
    bad=np.isnan(values);ax[2].scatter(rt[bad],np.full(bad.sum(),-.55),marker='x',color='crimson',label='unlabelled')
    ax[2].set_yticks(range(8));ax[2].set_ylim(-.8,7.6);ax[2].set_ylabel('Decoded value')
    ax[2].text(.01,.96,'Steady: '+result['steady']['sequence'],transform=ax[2].transAxes,va='top',fontsize=10)
    ax[3].plot(t,sig['clock'],color='gray',alpha=.7,label='clock factor')
    ax[3].plot(t,sig['g1'],color='#d97726',label='g1');ax[3].set_ylabel('Gate activity');ax[3].legend(loc='upper right',ncol=2)
    ax[4].plot(t,y[IDX['I2_zmh']],color='#2685ad',label='mature Int2 (a.u.)');ax[4].set_ylabel('Int2 (a.u.)')
    a2=ax[4].twinx();a2.plot(t,sig['int2_source'],color='#d97726',alpha=.65,label='Int2 source (a.u./h)');a2.set_ylabel('Source (a.u./h)')
    lines=ax[4].get_lines()+a2.get_lines();ax[4].legend(lines,[l.get_label() for l in lines],loc='upper right',ncol=2)
    ax[-1].set_xlabel('Time (h)')
    for a in ax:a.grid(alpha=.15);a.set_xlim(0,HOURS)
    fig.suptitle(title+'\nHan input; 23-state HBY-ZMH-ZMH; diagnostic, not inherited certification',fontsize=13)
    for ext in ('png','svg','pdf'):fig.savefig(target/f'overview.{ext}',dpi=180)
    plt.close(fig)


def comparison(t,cases,out):
    fig, axes=plt.subplots(4,2,figsize=(15,10),sharex=True,sharey='row',layout='constrained')
    for col,(tag,(y,sig,res)) in enumerate(cases.items()):
        axes[0,col].set_title('A1 negative autoregulation '+tag.upper())
        for b,c in zip(('S0','S1','S2'),('#2685ad','#d97726','#43875a')):
            axes[0,col].plot(t,sig[b],label=b,color=c,lw=1.4)
        axes[0,col].axhspan(0,.3,color='#43875a',alpha=.08);axes[0,col].axhspan(.7,1,color='#2685ad',alpha=.08)
        axes[0,col].legend(ncol=3,fontsize=8);axes[0,col].set_ylim(-.03,1.03)
        axes[1,col].plot(t,y[17],label='A1');axes[1,col].plot(t,y[18],label='F1');axes[1,col].legend()
        axes[2,col].plot(t,sig['g1'],label='g1',color='#d97726');axes[2,col].plot(t,sig['clock'],label='clock',color='gray',alpha=.6);axes[2,col].legend()
        axes[3,col].plot(t,sig['J_fwd2'],label='forward',color='#2685ad');axes[3,col].plot(t,sig['J_rev2'],label='reverse',color='#d97726');axes[3,col].legend()
        for row in range(4):axes[row,col].set_xlim(100,220);axes[row,col].grid(alpha=.15)
        axes[3,col].set_xlabel('Time (h)')
    for a,label in zip(axes[:,0],('LR fraction','A1 / F1 (a.u.)','Gate activity','Bit2 flux (/h)')):a.set_ylabel(label)
    fig.suptitle('Matched tail comparison: only A1 autoregulation and its consistent initial state differ')
    for ext in ('png','svg','pdf'):fig.savefig(out/f'comparison.{ext}',dpi=180)
    plt.close(fig)


def main():
    if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
    out=ROOT/'results'/datetime.now().strftime('%Y%m%d_%H%M%S');out.mkdir(parents=True)
    paths=[ROOT/'model_hzz.py',Path(__file__),TAIL_SOURCE,PARENT/'hybrid_model.py',
           PARENT/'hby_zmh_hby/model_hzh.py',PARENT/'hby_zmh_hby/certification/verify_hzh.py',
           PARENT/'run_han_comparison.py',HAN_PATH]
    before={str(p):sha256(p) for p in paths};before.update(source_hashes())
    dump(out/'status.json',dict(status='running'))
    print('Build shared 300 h Han input',flush=True)
    han=HanInput(HOURS);prefix=HZHModel('han',han)
    checks=structural_checks(prefix);dump(out/'structural_checks.json',checks)
    print('Structural checks:',checks,flush=True)
    print('Run HZH reference',flush=True)
    t,base=integrate(prefix);bv,bs=V.analyse(t,base,prefix.z,prefix.tail)
    guard=dict(certified=bv['certified_v1'],reads=bv['steady']['reads']==19,
               sequence=bv['steady']['sequence']=='1234567012345670123',
               margin=abs(bv['global_min_timing_margin_h']-1.1784609018563117)<1e-6)
    dump(out/'baseline_guard.json',dict(checks=guard,verdict=bv))
    if not all(guard.values()):raise RuntimeError(guard)
    print('Baseline guard passed',flush=True)
    cases={};summaries={}
    for tag,enabled in [('off',False),('on',True)]:
        print('Run autoregulation '+tag,flush=True)
        target=out/tag;target.mkdir()
        model=HZZModel(prefix,enabled);t,y=integrate(model);sig=model.signals(y)
        if y.min() < -1e-7 or any(sig[b].min() < -1e-7 or sig[b].max()>1+1e-7 for b in ('S0','S1','S2')):
            raise RuntimeError('Nonphysical state trajectory')
        res=evaluate(t,sig)
        res['prefix_trajectory_gap_vs_HZH']=float(np.max(np.abs(y[:17]-base[:17])))
        res['clock_distribution']=dict(p05=float(np.quantile(sig['clock'],.05)),median=float(np.median(sig['clock'])),
                                       maximum=float(sig['clock'].max()),duty_gt_05=float(np.mean(sig['clock']>.5)),
                                       duty_gt_01=float(np.mean(sig['clock']>.1)))
        res['state_extrema']={n:[float(y[i].min()),float(y[i].max())] for i,n in enumerate(NAMES)}
        dump(target/'verdict.json',res)
        dump(target/'parameters.json',dict(autoregulation=enabled,tail_hour_parameters=model.p,
             middle_minute_parameters=prefix.z,middle_time_conversion=60,hby_bit0_table=prefix.p,
             hby_config=asdict(prefix.tail.c),clock_K=.4,clock_n=3.,clock_scale=1.,
             clock_source='own HBY b0_I; no amplitude retuning',state_names=NAMES,
             initial_state=model.initial_state().tolist(),solver=SOLVER,sample_min=1,
             tail_source_sha256=before[str(TAIL_SOURCE)],initialization='A1/F1 equilibrium at PB1=1; no-feedback analytical, feedback brentq'))
        np.savez_compressed(target/'trajectory.npz',time_h=t,states=y,state_names=np.array(NAMES),**sig)
        with (target/'trajectory.csv').open('w',encoding='utf-8',newline='') as f:
            w=csv.writer(f);w.writerow(['time_h',*NAMES,*sig])
            for k,tt in enumerate(t):w.writerow([tt,*y[:,k],*[v[k] for v in sig.values()]])
        draw(t,y,sig,res,target,'A1 negative autoregulation '+tag.upper())
        cases[tag]=(y,sig,res)
        summaries[tag]=dict(counting=res['counting_passed'],events=res['event_causality_passed'],
             diagnostic_passed=res['diagnostic_passed'],cold=res['cold']['sequence'],steady=res['steady']['sequence'],
             per_bit=res['per_bit'],prefix_gap=res['prefix_trajectory_gap_vs_HZH'],
             flips={b:len(c) for b,c in res['crossings'].items()},
             gate_events={k:[v['reverse_events'],v['gate_events'],v['flip_events']] for k,v in res['events'].items()},
             A1_initial=float(model.initial_state()[17]),F1_initial=float(model.initial_state()[18]),
             S2_range=[float(sig['S2'].min()),float(sig['S2'].max())])
        print(tag,json.dumps(summaries[tag],ensure_ascii=False),flush=True)
    comparison(t,cases,out)
    after={str(p):sha256(p) for p in paths};after.update(source_hashes())
    if before!=after:raise RuntimeError('Input source changed during simulation')
    dump(out/'summary.json',dict(hours=HOURS,cases=summaries,source_hashes=before,source_unchanged=True,
         python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()))
    text=['# HBY–ZMH–ZMH：A1负自馈对照结果','',
          '韩上游300 h；23状态。前17状态复用原HBY–ZMH–HBY的低两位。末级来自《前馈三级级联.py》，属于探索参数版本。',
          '','| A1负自馈 | 稳态码串（drop8） | 计数 | 事件臂 | bit2翻转数 |', '|---|---|---|---|---|']
    for tag,r in summaries.items():text.append(f"| {tag} | `{r['steady']}` | {r['counting']} | {r['events']} | {r['flips']['S2']} |")
    text += ['', '## 实现与解释', '',
      '- 两组只改变A1产生项是否乘负自馈因子，并分别设置方程一致的A1/F1初值；其余参数完全相同。',
      '- 保留新ZMH末级clock_K=0.4、clock_n=3；时钟读取本次HBY bit0的成熟游离Int0。没有调幅或修改门指数求成功。',
      '- A1/F1为直接蛋白方程，bit2为PB、Int、Rep、RDF四态，无显式复合物。当前ZMH中间位按分钟参数乘60转换为小时；新末级本来以小时积分。',
      '- 新末级保留k_fwd/k_rev=7/7、K_D_comp=0.8、K_rep=0.45、alpha_Int2=28、gamma_int2=2、K_D_int2=1；不能称为现行二级参数的单位换算。',
      '- 原三级脚本关闭A1负自馈，但以含负自馈迭代式生成初值；本实验off采用A1=6.4，on采用brentq解约1.2301，并据此平衡F1。',
      '- 沿用原验证器的通用读窗/事件函数，按本23态命名信号重新评估；certified字段为null，diagnostic_passed为本轮组合判据，不继承旧模型认证。',
      '- x表示未满足0.3/0.7、80%占用率读窗要求，不绘成0，也不以预设数字代替。事件臂使用既有阈值，失败可有阈值敏感性。',
      '- 本轮是单初值300 h探索。未做八初态、长期扫描、随机性或实验可靠性验证。',
      '', '## 文件', '', '每组overview.png/svg/pdf、trajectory.csv/npz、verdict.json、parameters.json；根目录comparison.png/svg/pdf为100–220 h相同时间窗对比。',
      'structural_checks.json记录原ZMH末级方程一致性、负自馈开关定位、低17态无反馈和初值平衡；baseline_guard.json记录已认证HZH基线复现。']
    (out/'结果说明.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    dump(out/'status.json',dict(status='completed',runs=3))
    files={str(p.relative_to(out)):sha256(p) for p in sorted(out.rglob('*')) if p.is_file()}
    dump(out/'SHA256SUMS.json',files)
    print('OUTPUT',out,flush=True)


if __name__=='__main__':main()
