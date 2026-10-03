"""Render the certified saved HZH trajectory for the Wiki, without integration."""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'hby_zmh_hby/results/20260927_233646_924908/han/trajectory.npz'
VERDICT=ROOT/'hby_zmh_hby/certification/results/20260928_141723_738780/baseline_saved/verdict.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest().upper()

def main():
    z=np.load(SOURCE);t,y=z['time_h'],z['states'];idx={str(n):i for i,n in enumerate(z['state_names'])}
    v=json.loads(VERDICT.read_text(encoding='utf-8'))
    assert v['certified_v1'] and len(v['read_windows'])==27
    state=[y[idx['b0_S']],1-y[idx['pb1_zmh']],y[idx['b2_S']]]
    assert all(s.min()>=-1e-7 and s.max()<=1+1e-7 for s in state)
    reads=v['read_windows'];rt=np.array([r['trough_h'] for r in reads]);rv=np.array([r['value'] for r in reads])
    assert ''.join(map(str,rv))=='123456701234567012345670123'
    plt.rcParams.update({'svg.fonttype':'none','font.size':11})
    colors=['#2878b5','#dc8b22','#329b78']
    fig,ax=plt.subplots(4,1,figsize=(12,10),sharex=True,constrained_layout=True,
                        gridspec_kw={'height_ratios':[.8,1.3,1.0,.9]})
    ax[0].plot(t,30*y[idx['b0_M_I']],color='#6853a3',lw=1.3)
    ax[0].set_ylabel('Int0 translation\n(a.u./h)')
    for s,c,n in zip(state,colors,['Bit 0: biochemical','Bit 1: reduced','Bit 2: biochemical']):ax[1].plot(t,s,color=c,lw=1.3,label=n)
    ax[1].axhspan(0,.3,color='#329b78',alpha=.06);ax[1].axhspan(.7,1,color='#2878b5',alpha=.06)
    ax[1].axhline(.3,color='gray',ls='--',lw=.7);ax[1].axhline(.7,color='gray',ls='--',lw=.7)
    ax[1].legend(loc='upper center',ncol=3,fontsize=9);ax[1].set_ylabel('LR fraction');ax[1].set_ylim(-.04,1.13)
    ax[2].step(rt,rv,where='post',color='#303030',lw=1.4)
    ax[2].scatter(rt,rv,c=rv,cmap='viridis',s=27,zorder=3)
    ax[2].set_yticks(range(8));ax[2].set_ylabel('Decoded value')
    for n,c in zip(['b0_I','I1_zmh','b2_I'],colors):ax[3].plot(t,y[idx[n]],color=c,lw=1.2,label={'b0_I':'Int0','I1_zmh':'Int1','b2_I':'Int2'}[n])
    ax[3].legend(loc='upper right',ncol=3);ax[3].set_ylabel('Integrase (a.u.)');ax[3].set_xlabel('Time (h)')
    for a in ax:a.grid(alpha=.12);a.set_xlim(0,300)
    fig.suptitle('Cross-model three-bit counting under repressilator input\n34 downstream states; 27 finite-window reads over 300 h')
    for ext in ('png','svg'):fig.savefig(HERE/f'hybrid_threebit_300h.{ext}',dpi=180)
    plt.close(fig)
    meta=dict(source=str(SOURCE.relative_to(ROOT)),source_sha256=sha(SOURCE),
              verdict=str(VERDICT.relative_to(ROOT)),verdict_sha256=sha(VERDICT),
              state_index=idx,read_count=27,sequence=''.join(map(str,rv)),
              translation_signal='30/h * shared b0_M_I; actual immature Int0 production term',
              note='Saved states and certified read windows only; no integration, smoothing or resampling.',
              outputs={p.name:sha(p) for p in (HERE/'hybrid_threebit_300h.png',HERE/'hybrid_threebit_300h.svg')})
    (HERE/'hybrid_threebit_300h.provenance.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Rendered certified 27-read trajectory')

if __name__=='__main__':main()
