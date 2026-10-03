"""Run two untuned HBY-ZMH-HBY cases; save evidence, not inherited certification."""
import argparse,csv,json,sys,time
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.integrate import solve_ivp
from model_hzh import HZHModel,NAMES,IDX,PARENT
from hybrid_model import source_hashes,sha256
from run_han_comparison import HanInput,HAN_PATH
from diagnostics import analyse,state_health

ROOT=Path(__file__).resolve().parent
def dump(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')

def integrate_case(model,t):
    if model.mode=='han':
        sol=solve_ivp(model.rhs,(0,t[-1]),model.initial_state(),t_eval=t,method='DOP853',rtol=2e-7,atol=2e-9,max_step=1/60)
        if not sol.success:raise RuntimeError(sol.message)
        return sol.y
    # Split exactly at square-wave jumps; use a constant on/off source in each piece.
    q=model.square
    starts=np.arange(q['start_h'],t[-1],q['period_h'])
    edges=np.unique(np.r_[0,t[-1],starts,starts+q['width_h']]);edges=edges[(edges>=0)&(edges<=t[-1])]
    saved=model.input_source
    y0=model.initial_state();y=np.empty((34,len(t)))
    try:
        for a,b in zip(edges[:-1],edges[1:]):
            value=saved((a+b)/2)
            model.input_source=lambda _t,v=value:v
            sol=solve_ivp(model.rhs,(a,b),y0,method='DOP853',rtol=2e-7,atol=2e-9,max_step=1/60,dense_output=True)
            if not sol.success:raise RuntimeError(sol.message)
            mask=(t>=a)&(t<=b);y[:,mask]=sol.sol(t[mask]);y0=sol.y[:,-1]
    finally:model.input_source=saved
    return y


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--hours',type=float,default=300.);args=ap.parse_args()
    if not np.isfinite(args.hours) or args.hours<=0:raise ValueError('positive duration')
    if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
    out=ROOT/'results'/datetime.now().strftime('%Y%m%d_%H%M%S_%f');out.mkdir(parents=True)
    hashes=source_hashes();hashes.update({str(p):sha256(p) for p in [HAN_PATH,ROOT/'model_hzh.py',ROOT/'run_hzh.py',PARENT/'hybrid_model.py',PARENT/'diagnostics.py',PARENT/'run_han_comparison.py']})
    dump(out/'status.json',dict(status='RUNNING'))
    results=[]
    try:
        print('Preparing Han upstream...',flush=True)
        han=HanInput(args.hours)
        np.savez_compressed(out/'han_upstream.npz',time_min=han.t,states=han.sol.y)
        dump(out/'han_parameters.json',dict(parameters=asdict(han.p),source_sha256=sha256(HAN_PATH),solver=dict(method='LSODA',rtol=1e-9,atol=1e-11,max_step_min=1.)))
        t=np.r_[np.arange(0,args.hours,1/60),args.hours];t=np.unique(t)
        for mode in ('square','han'):
            print('Running',mode,args.hours,'h...',flush=True)
            tic=time.perf_counter();model=HZHModel(mode,han if mode=='han' else None)
            y=integrate_case(model,t)
            if not np.isfinite(y).all():raise RuntimeError('Nonfinite trajectory')
            case=out/mode;case.mkdir()
            source=np.array([model.input_source(tt) for tt in t])
            u2=np.array([model.tail.signals(y[17:,j],y[2,j])['u2_target_au_per_h'] for j in range(len(t))])
            np.savez_compressed(case/'trajectory.npz',time_h=t,states=y,state_names=NAMES,input_source=source,u2=u2)
            with (case/'trajectory.csv').open('w',newline='',encoding='utf-8') as f:
                w=csv.writer(f);w.writerow(['time_h',*NAMES,'input_source_au_per_h','u2_au_per_h'])
                for j,tt in enumerate(t):w.writerow([tt,*y[:,j],source[j],u2[j]])
            read=analyse(t,model.diagnostic_projection(y),y[2])
            ss=[y[10],1-y[13],y[33]]
            peaks=find_peaks(y[2],prominence=.1*np.ptp(y[2]))[0]
            summary=dict(mode=mode,hours=args.hours,runtime_s=time.perf_counter()-tic,
                certified=None,readout_status=read['status'],
                health=state_health(y,NAMES,('b0_S','pb1_zmh','b2_S')),
                cold_mod8=read['cold_mod8'],steady_mod8=read['steady_mod8'],
                cold_mod4=read['cold_mod4'],steady_mod4=read['steady_mod4'],
                events={k:len(v) for k,v in read['crossings'].items()},
                post_50h_state_ranges=[dict(min=float(s[t>=min(50,args.hours)].min()),max=float(s[t>=min(50,args.hours)].max())) for s in ss],
                Int0_peak=float(y[2].max()),Int1_peak=float(y[14].max()),Int2_peak=float(y[25].max()),
                input_source_peak=float(source.max()),u2_peak=float(u2.max()),
                median_Int0_period_h=float(np.median(np.diff(t[peaks]))) if len(peaks)>1 else None)
            dump(case/'summary.json',summary);dump(case/'readout.json',read)
            dump(case/'parameters.json',dict(mode=mode,square=model.square if mode=='square' else None,
                Hby_config=asdict(model.tail.c),Hby_parameters=model.p,Zmh_parameters_per_min=model.z,
                initial_state=dict(zip(NAMES,model.initial_state().tolist())),
                units='hours and module-specific a.u.; PB is dimensionless; clock_scale=1 uncalibrated',
                bit0_input=('target production 6 a.u./h through mRNA/maturation' if mode=='square' else
                            'shared C31 mRNA from Han promoter; true translation into immature Int0; no normalized concentration or x6'),
                solver=dict(method='DOP853',rtol=2e-7,atol=2e-9,max_step_min=1,square_discontinuities_segmented=True)))
            fig,ax=plt.subplots(4,1,figsize=(13,10),sharex=True,constrained_layout=True)
            ax[0].plot(t,source,label='Input translation/target source');ax[0].set_ylabel('a.u./h');ax[0].legend()
            for s,n in zip(ss,['S0 HBY','S1 ZMH','S2 HBY']):ax[1].plot(t,s,label=n)
            ax[1].axhspan(0,.3,alpha=.07,color='green');ax[1].axhspan(.7,1,alpha=.07,color='blue');ax[1].legend();ax[1].set_ylabel('LR fraction')
            rr=read['reads'];ax[2].plot([r['trough_h'] for r in rr],[np.nan if r['value'] is None else r['value'] for r in rr],'o-',drawstyle='steps-post')
            for r in rr:
                if r['value'] is None:ax[2].text(r['trough_h'],3.5,'x',color='red')
            ax[2].set_yticks(range(8));ax[2].set_ylabel('Window value')
            for idx,n in ((2,'Int0'),(14,'Int1'),(25,'Int2')):ax[3].plot(t,y[idx],label=n)
            ax[3].legend();ax[3].set_xlabel('Time (h)');ax[3].set_ylabel('a.u.')
            fig.suptitle(f'HBY bit0 -> ZMH bit1 -> HBY bit2 / {mode}\nExploratory - no inherited certification')
            fig.savefig(case/'overview.png',dpi=150);fig.savefig(case/'overview.pdf');plt.close(fig)
            results.append(summary)
            print(json.dumps({k:summary[k] for k in ['mode','cold_mod8','steady_mod8','cold_mod4','steady_mod4','events','Int0_peak','Int1_peak','Int2_peak']},indent=2),flush=True)
        after=source_hashes();after.update({k:sha256(k) for k in hashes if Path(k).is_absolute()})
        if hashes!=after:raise RuntimeError('Source changed during run')
        dump(out/'comparison.json',dict(source_hashes=hashes,cases=results))
        dump(out/'status.json',dict(status='COMPLETED',certified=None))
        dump(out/'SHA256SUMS.json',{str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file()})
        print('OUTPUT:',out,flush=True)
    except Exception as exc:
        dump(out/'status.json',dict(status='FAILED',error=repr(exc)));raise


if __name__=='__main__':main()
