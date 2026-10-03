import json,csv,sys
from pathlib import Path
from datetime import datetime
from dataclasses import asdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model_pdf import PDFModel,HERE,PDF_PATH,PDF_VALUES,RATE_KEYS,NAMES,IDX,PARENT
from model_hzh import HZHModel
from run_hzh import integrate_case
from run_han_comparison import HanInput,HAN_PATH
from hybrid_model import sha256,source_hashes
from diagnostics import analyse,state_health

def dump(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')

def square_off_reads(t,states,params):
    out=[]
    for start in np.arange(params['start_h'],t[-1],params['period_h']):
        centre=start+(params['width_h']+params['period_h'])/2
        left=centre-.1*params['period_h'];right=centre+.1*params['period_h']
        if right>t[-1]:continue
        mask=(t>=left)&(t<=right);labels=[];commit=[]
        for s in states:
            lo=float(np.mean(s[mask]<=.3));hi=float(np.mean(s[mask]>=.7))
            labels.append(0 if lo>=.8 else 1 if hi>=.8 else None);commit.append(max(lo,hi))
        val=sum(2**i*x for i,x in enumerate(labels)) if all(x is not None for x in labels) else None
        out.append(dict(centre_h=float(centre),labels=labels,value=val,commitment=commit))
    v=[r['value'] for r in out]
    return dict(rule='OFF midpoint; total width20% period; chosen before this rerun from previous diagnostic',
        rows=out,sequence=''.join('x' if x is None else str(x) for x in v),
        observed_increments=bool(len(v)>1 and all(x is not None for x in v) and all((b-a)%8==1 for a,b in zip(v[:-1],v[1:]))))

def main():
    if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
    out=HERE/'results'/datetime.now().strftime('%Y%m%d_%H%M%S_%f');out.mkdir(parents=True)
    paths=[PDF_PATH,HAN_PATH,HERE/'model_pdf.py',HERE/'run_pdf.py',HERE.parent/'model_hzh.py',HERE.parent/'run_hzh.py',PARENT/'hybrid_model.py',PARENT/'diagnostics.py',PARENT/'run_han_comparison.py']
    before=source_hashes();before.update({str(p):sha256(p) for p in paths})
    dump(out/'status.json',dict(status='RUNNING'))
    try:
        han=HanInput(300);t=np.linspace(0,300,18001)
        np.savez_compressed(out/'han_upstream.npz',time_min=han.t,states=han.sol.y)
        dump(out/'han_parameters.json',dict(parameters=asdict(han.p),solver=dict(method='LSODA',rtol=1e-9,atol=1e-11,max_step_min=1)))
        results=[]
        for mode in ('square','han'):
            print('Running PDF middle parameters:',mode,flush=True)
            model=PDFModel(mode,han if mode=='han' else None)
            audit=model.audit();assert audit['matching_count']==len(PDF_VALUES)
            if mode=='square':dump(out/'parameter_audit.json',audit)
            y=integrate_case(model,t)
            if not np.isfinite(y).all():raise RuntimeError('Nonfinite')
            folder=out/mode;folder.mkdir()
            np.savez_compressed(folder/'trajectory.npz',time_h=t,states=y,state_names=NAMES)
            with (folder/'trajectory.csv').open('w',encoding='utf-8',newline='') as f:
                w=csv.writer(f);w.writerow(['time_h',*NAMES]);w.writerows(zip(t,*y))
            r=analyse(t,model.diagnostic_projection(y),y[2]);dump(folder/'readout.json',r)
            ss=[y[10],1-y[13],y[33]]
            extra=square_off_reads(t,ss,model.square) if mode=='square' else None
            if extra is not None:dump(folder/'square_off_readout.json',extra)
            summary=dict(mode=mode,hours=300,certified=None,parameter_audit=audit,
                cold_mod8=r['cold_mod8'],steady_mod8=r['steady_mod8'],
                cold_mod4=r['cold_mod4'],steady_mod4=r['steady_mod4'],
                events={k:len(v) for k,v in r['crossings'].items()},
                health=state_health(y,NAMES,('b0_S','pb1_zmh','b2_S')),
                post50_ranges=[dict(min=float(s[t>=50].min()),max=float(s[t>=50].max())) for s in ss],
                square_off_sequence=extra['sequence'] if extra else None,
                square_off_increments=extra['observed_increments'] if extra else None)
            dump(folder/'summary.json',summary)
            dump(folder/'parameters.json',dict(zmh_per_minute=model.z,hby=model.p,hby_config=asdict(model.tail.c),
                 initial_state=dict(zip(NAMES,model.initial_state().tolist())),square=model.square if mode=='square' else None,
                 solver=dict(method='DOP853',rtol=2e-7,atol=2e-9,max_step_min=1,square_edges_segmented=True)))
            fig,ax=plt.subplots(3,1,figsize=(13,9),sharex=True,constrained_layout=True)
            for s,n in zip(ss,['S0 HBY','S1 ZMH PDF','S2 HBY']):ax[0].plot(t,s,label=n)
            ax[0].axhspan(0,.3,color='green',alpha=.08);ax[0].axhspan(.7,1,color='blue',alpha=.08);ax[0].legend();ax[0].set_ylabel('LR fraction')
            reads=extra['rows'] if extra else r['reads'];key='centre_h' if extra else 'trough_h'
            ax[1].plot([q[key] for q in reads],[np.nan if q['value'] is None else q['value'] for q in reads],'o-',drawstyle='steps-post')
            for q in reads:
                if q['value'] is None:ax[1].text(q[key],3.5,'x',color='red')
            ax[1].set_yticks(range(8));ax[1].set_ylabel('OFF-midpoint read' if extra else 'Trough read')
            for idx,n in [(2,'Int0'),(14,'Int1'),(25,'Int2')]:ax[2].plot(t,y[idx],label=n)
            ax[2].legend();ax[2].set_ylabel('a.u.');ax[2].set_xlabel('Time (h)')
            fig.suptitle(f'HBY-ZMH-HBY / PDF-listed middle parameters / {mode}\nn_int1=4 retained (not in table); no inherited certification')
            fig.savefig(folder/'overview.png',dpi=150);fig.savefig(folder/'overview.pdf');plt.close(fig)
            results.append(summary)
            print(json.dumps({k:summary[k] for k in ['mode','cold_mod8','steady_mod8','events','post50_ranges','square_off_sequence','square_off_increments']},indent=2),flush=True)
        after=source_hashes();after.update({str(p):sha256(p) for p in paths})
        if before!=after:raise RuntimeError('Source changed')
        dump(out/'comparison.json',dict(source_hashes=before,cases=results))
        dump(out/'status.json',dict(status='COMPLETED',certified=None))
        dump(out/'SHA256SUMS.json',{str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file()})
        print('OUTPUT',out,flush=True)
    except Exception as e:
        dump(out/'status.json',dict(status='FAILED',error=repr(e)));raise

if __name__=='__main__':main()
