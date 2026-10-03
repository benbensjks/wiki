"""Wiki figure from saved original-RK45 donor output; no new ODE integration."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'results/20260927_222727_239025/trajectories.npz'
PARAMS=ROOT/'results/20260927_222727_239025/parameters.json'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def main():
    z=np.load(SOURCE)
    t=z['original_time_min']/60
    y=z['original_donor']
    if y.shape!=(10,len(t)) or not np.isfinite(y).all():raise ValueError('Bad donor trajectory')
    s0,s1=1-y[0],1-y[6]
    if min(s0.min(),s1.min())<-1e-6 or max(s0.max(),s1.max())>1+1e-6:raise ValueError('DNA range failed')
    plt.rcParams.update({'svg.fonttype':'none','font.size':11})
    fig,axs=plt.subplots(3,1,figsize=(11,8),sharex=True,constrained_layout=True)
    axs[0].plot(t,s0,label='Bit 0: LR fraction',color='#2878b5')
    axs[0].plot(t,s1,label='Bit 1: LR fraction',color='#dc8b22')
    axs[0].set_ylabel('DNA state');axs[0].legend(loc='upper left',ncol=2)
    axs[0].set_ylim(-.03,1.07)
    axs[1].plot(t,y[4],label='A0: activator',color='#329b78')
    axs[1].plot(t,y[5],label='F0: repressor (R0 in source)',color='#b44764')
    axs[1].set_ylabel('Regulators (a.u.)');axs[1].legend(loc='upper right')
    axs[2].plot(t,y[1],label='Int0',color='#333333',linestyle='--')
    axs[2].plot(t,y[7],label='Int1: carry to bit 1',color='#6853a3')
    axs[2].set_ylabel('Integrase (a.u.)');axs[2].set_xlabel('Time (h)');axs[2].legend(loc='upper left')
    for ax in axs:ax.grid(alpha=.16);ax.set_xlim(0,t[-1])
    fig.suptitle('ZMH two-bit model: current Python parameters\nNormalized C31-protein input; no explicit maturation layer')
    for ext in ('png','svg'):fig.savefig(HERE/f'zmh_twobit_current.{ext}',dpi=180)
    plt.close(fig)
    meta=dict(source=str(SOURCE.relative_to(ROOT)),source_sha256=sha(SOURCE),
              parameters=str(PARAMS.relative_to(ROOT)),parameters_sha256=sha(PARAMS),
              array='original_donor',time='original_time_min / 60',
              rows=len(t),hours=float(t[-1]),state_indices=dict(pb0=0,pb1=6,A0=4,F0_source_r0=5,Int0=1,Int1=7),
              note='Read-only rendering of original default RK45 donor trajectory already saved by run_initial.py; no reintegration.',
              outputs={p.name:sha(p) for p in (HERE/'zmh_twobit_current.png',HERE/'zmh_twobit_current.svg')})
    (HERE/'zmh_twobit_current.provenance.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Rendered',len(t),'saved samples')


if __name__=='__main__':main()
