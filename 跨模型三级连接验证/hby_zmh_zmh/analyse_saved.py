"""Postprocess saved trajectories: no integration, preserve raw data."""
import sys
sys.dont_write_bytecode=True
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model_hzz import ROOT, IDX, NAMES, PARENT
sys.path.insert(0,str(PARENT/'hby_zmh_hby/certification'))
import verify_hzh as V
from hybrid_model import sha256


def main():
    out=Path(sys.argv[1]) if len(sys.argv)>1 else sorted((ROOT/'results').iterdir())[-1]
    records={};data={}
    for tag in ('off','on'):
        a=np.load(out/tag/'trajectory.npz',allow_pickle=False)
        v=json.loads((out/tag/'verdict.json').read_text(encoding='utf-8'))
        t=a['time_h'];y=a['states'];late=t>=200
        assert list(a['state_names'])==list(NAMES)
        assert np.array_equal(a['S2'],1-y[IDX['pb2_zmh']])
        assert np.array_equal(a['int2_source'],28*a['g1'])
        raw=V.segments(t,a['g1'],V.GATE_THRESHOLD,min_duration_h=0,min_dose=0)
        rows=[dict(**g,duration_h=g['end_h']-g['start_h'],
                   meets_duration=g['end_h']-g['start_h']>=V.PULSE_DURATION_H,
                   meets_dose=g['dose']>=V.GATE_MIN_DOSE_H) for g in raw]
        records[tag]=dict(g1_peak=float(a['g1'].max()),Int2_source_peak_au_per_h=float(a['int2_source'].max()),
             mature_Int2_peak_au=float(y[IDX['I2_zmh']].max()),
             S2_last100h=[float(a['S2'][late].min()),float(a['S2'][late].max())],
             final_S2=float(a['S2'][-1]),raw_gate_segments=rows,
             qualifying_gates=v['events']['bit1_to_bit2']['gate_events'],
             J_fwd2_total=float(np.trapezoid(a['J_fwd2'],t)),J_rev2_total=float(np.trapezoid(a['J_rev2'],t)),
             counted_flips=len(v['crossings']['S2']),steady=v['steady'])
        records[tag]['DNA_balance_residual']=abs(records[tag]['J_fwd2_total']-records[tag]['J_rev2_total']
                                                -float(a['S2'][-1]-a['S2'][0]))
        assert records[tag]['DNA_balance_residual']<1e-4
        data[tag]=(a,v)
    prefixgap=float(np.max(np.abs(data['off'][0]['states'][:17]-data['on'][0]['states'][:17])))
    records['prefix_gap_off_on']=prefixgap
    assert prefixgap<1e-6
    (out/'diagnostics.json').write_text(json.dumps(records,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    events=[e['time_h'] for e in data['off'][1]['crossings']['S1'] if e['direction']=='down' and 200<e['time_h']<280]
    centre=events[-1]
    fig,axes=plt.subplots(5,2,figsize=(14,11),sharex=True,sharey='row',layout='constrained')
    for col,tag in enumerate(('off','on')):
        a,v=data[tag];t=a['time_h'];y=a['states'];x=t-centre
        axes[0,col].set_title('A1 autoregulation '+tag.upper())
        for key,color in [('activation','#2685ad'),('repression','#43875a'),('clock','#888888')]:
            axes[0,col].plot(x,a[key],label=key,color=color)
        axes[0,col].legend(fontsize=8);axes[0,col].set_ylim(-.02,1.02)
        axes[1,col].plot(x,a['g1'],color='#d97726',label='g1')
        axes[1,col].axhline(.05,ls='--',color='gray',label='event threshold 0.05');axes[1,col].legend(fontsize=8)
        axes[2,col].plot(x,y[IDX['I2_zmh']],color='#2685ad')
        axes[3,col].plot(x,a['J_fwd2'],color='#2685ad',label='forward')
        axes[3,col].plot(x,a['J_rev2'],color='#d97726',label='reverse');axes[3,col].legend(fontsize=8)
        axes[4,col].plot(x,a['S2'],color='#43875a');axes[4,col].axhline(.7,ls='--',color='gray')
        axes[4,col].set_xlabel('Time relative to bit1 reset (h)')
        for ax in axes[:,col]:ax.set_xlim(-1,4);ax.axvline(0,color='gray',ls=':',lw=.8);ax.grid(alpha=.15)
    for ax,label in zip(axes[:,0],('Gate factors','g1','Mature Int2 (a.u.)','Bit2 flux (/h)','S2 (LR fraction)')):ax.set_ylabel(label)
    fig.suptitle(f'One matched bit1 reset at t={centre:.3f} h; narrow gate pulses and incomplete bit2 writing')
    for ext in ('png','svg','pdf'):fig.savefig(out/f'pulse_zoom.{ext}',dpi=180)
    plt.close(fig)
    lines=['# 初步诊断：两种末级均未实现模8','',
           '这次运行保留了曾同学新三级文件的探索性末级参数与时钟门0.4/3。负自馈仅作用于A1产生项。',
           '', '| 指标 | 无负自馈 | 有负自馈 |', '|---|---|---|']
    for title,key in [('g1峰值','g1_peak'),('Int2产生源峰值（a.u./h）','Int2_source_peak_au_per_h'),
                      ('实际成熟Int2峰值（a.u.）','mature_Int2_peak_au'),('bit2穿越0.5次数','counted_flips'),
                      ('J_fwd2全程积分','J_fwd2_total'),('J_rev2全程积分','J_rev2_total')]:
        lines.append(f"| {title} | {records['off'][key]:.6g} | {records['on'][key]:.6g} |")
    for tag,label in [('off','无负自馈'),('on','有负自馈')]:
        r=records[tag];segments=r['raw_gate_segments']
        lines += ['',f'## {label}', '',f"- 稳态码串：`{r['steady']['sequence']}`。",
            f"- 200–300 h的S2范围为{r['S2_last100h'][0]:.5f}–{r['S2_last100h'][1]:.5f}；没有回到0带。",
            f"- g1超过0.05的原始片段有{len(segments)}个；时长范围{min(s['duration_h'] for s in segments):.4f}–{max(s['duration_h'] for s in segments):.4f} h，片段剂量范围{min(s['dose'] for s in segments):.5f}–{max(s['dose'] for s in segments):.5f} h。",
            '- 既有事件规则另要求持续至少0.2 h且片段积分至少0.02 h。没有片段同时满足，所以识别门事件为0；不能据此说没有进位信号。']
    lines += ['', '## 可作出的结论', '',
      f'- 前17状态与原已认证模型的轨迹一致；两组之间前17状态最大差{prefixgap:.3g}。低两位正常，末级替换后没有实现反复高低交替。',
      '- 开启A1负自馈改变了门波形，并产生部分反向写入，但当前参数下仍不足以让bit2完整回零。不能推断负自馈在其他参数或时钟条件下无效。',
      '- 这次同时替换了末级表达层、DNA重组响应、参数表及clock门0.3/2→0.4/3；它是末级整模块兼容性试验，不能将失败唯一归因于其中某个参数。',
      '- 本轮没有调参救活。后续若诊断脉冲不足，应同时检查持续时间、正反向剂量、RDF池及门和时钟的相对时序，不能只凭峰值下结论。',
      '- 负自馈开关两组的模型初值各自平衡；此结果不等同于逐字运行原三级文件（原文件初值与关闭负自馈方程不一致）。',
      '', '图：`off/overview.png`、`on/overview.png`为完整300 h；`comparison.png`为同时间窗对照；`pulse_zoom.png`为单次进位放大。']
    (out/'诊断补充.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (out/'analysis_metadata.json').write_text(json.dumps(dict(script=str(Path(__file__)),sha256=sha256(Path(__file__)),
                no_new_integration=True,
                parameter_note='parameters.json preserves the full source p dictionary; only the listed tail keys are used in the new six-state module',
                active_tail_parameter_keys=['alpha_A1','gamma_A1','K_auto1','n_auto1','alpha_R1','gamma_R1',
                    'K_A1','n_A1','K_R1','n_R1','alpha_Int2','gamma_int2','K_D_int2','k_fwd','k_rev',
                    'Kinh','K_D_comp','alpha_rep','gamma_rep','alpha_rdf','gamma_rdf','K_rep','n'],
                hardcoded_source_clock=dict(K=.4,n=3,scale=1)),indent=2),encoding='utf-8')
    files={str(p.relative_to(out)):sha256(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.json'}
    (out/'SHA256SUMS.json').write_text(json.dumps(files,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:{kk:vv for kk,vv in r.items() if kk!='raw_gate_segments'} for k,r in records.items() if isinstance(r,dict)},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
