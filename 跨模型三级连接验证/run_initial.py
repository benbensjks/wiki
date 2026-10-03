"""Run a reproducible preliminary donor/hybrid comparison, never overwrite a run."""
import argparse
import csv
import json
import platform
import sys
import time
from dataclasses import asdict
from datetime import datetime

import numpy as np
import scipy

from hybrid_model import (ROOT,DONOR_NAMES,STATE_NAMES,INDEX,ZmhDonor,HbyConfig,HbyReceiver,
                          HybridModel,integrate,time_grid,source_hashes,sha256)
from diagnostics import analyse,state_health


def write_json(path,data):
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hours',type=float,default=2000/60)
    parser.add_argument('--sample-min',type=float,default=1.)
    parser.add_argument('--max-step-min',type=float,default=1.)
    parser.add_argument('--upstream-mode',choices=['original','tight'],default='original')
    parser.add_argument('--n-gate',type=float,default=6.)
    parser.add_argument('--clock-scale',type=float,default=1.)
    parser.add_argument('--no-plot',action='store_true')
    args=parser.parse_args()
    if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
    ts=time_grid(args.hours,args.sample_min)
    before=source_hashes()
    out=ROOT/'results'/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    out.mkdir(parents=True,exist_ok=False)
    write_json(out/'status.json',{'status':'RUNNING','arguments':vars(args)})
    start=time.perf_counter()
    try:
        print('Preparing original 0..3000 minute upstream...',flush=True)
        donor=ZmhDonor(args.upstream_mode)
        receiver=HbyReceiver(HbyConfig(n_A1_gate=args.n_gate,clock_scale=args.clock_scale))
        model=HybridModel(donor,receiver)
        print('Integrating independent donor and hybrid...',flush=True)
        dsol=integrate(donor.rhs,donor.initial_state(),ts,max_step_min=args.max_step_min)
        hsol=integrate(model.rhs,model.initial_state(),ts,max_step_min=args.max_step_min)
        # Original counter solver/grid is independently retained for comparison.
        orig_t=np.linspace(0,args.hours*60,2000)
        original=scipy.integrate.solve_ivp(donor.rhs_min,(0,args.hours*60),donor.initial_state(),
                                          t_eval=orig_t,method='RK45')
        if not original.success or not np.isfinite(original.y).all():
            raise RuntimeError(f'Original donor solver failed: {original.message}')
        np.savez_compressed(out/'trajectories.npz',time_h=ts,hybrid=hsol.y,donor=dsol.y,
                            state_names=np.array(STATE_NAMES),original_time_min=original.t,
                            original_donor=original.y,upstream_time_min=donor.upstream.t,
                            upstream=donor.upstream.y)
        signals={n:[] for n in ('u_in','g0','u1_target_au_per_h','g1','clock_gate','u2_target_au_per_h','J_fwd2_per_h','J_rev2_per_h')}
        for k,t in enumerate(ts):
            yy=hsol.y[:,k]
            sig=receiver.signals(yy[10:],yy[1])
            sig.update(u_in=donor.input_min(60*t),g0=donor.gate0(yy[:10]))
            sig['u1_target_au_per_h']=60*donor.p['alpha_Int1']*sig['g0']
            for n in signals: signals[n].append(float(sig[n]))
        with (out/'trajectories.csv').open('w',encoding='utf-8',newline='') as f:
            w=csv.writer(f);w.writerow(['time_h',*STATE_NAMES,'S0','S1','S2',*signals])
            for k,t in enumerate(ts):
                yy=hsol.y[:,k]
                w.writerow([t,*yy,1-yy[0],1-yy[6],yy[INDEX['b2_S']],*[signals[n][k] for n in signals]])
        verdict=analyse(ts,hsol.y,hsol.y[1])
        gap=float(np.max(np.abs(hsol.y[:10]-dsol.y)))
        write_json(out/'readout.json',verdict)
        summary=dict(status='PRELIMINARY_COMPLETED_NOT_CERTIFIED',hours=args.hours,
                     runtime_s=time.perf_counter()-start,prefix_trajectory_max_abs_gap=gap,
                     donor_health=state_health(dsol.y,DONOR_NAMES,('pb0','pb1')),
                     hybrid_health=state_health(hsol.y,STATE_NAMES,('pb0','pb1','b2_S')),
                     donor_and_hybrid_low_states_same_within_1e_4=gap<1e-4,
                     cold_mod8=verdict['cold_mod8'],steady_mod8=verdict['steady_mod8'],
                     cold_mod4=verdict['cold_mod4'],steady_mod4=verdict['steady_mod4'],
                     crossing_counts={k:len(v) for k,v in verdict['crossings'].items()})
        write_json(out/'parameters.json',dict(arguments=vars(args),donor=donor.metadata(),
                   receiver=asdict(receiver.c),receiver_table=receiver.p,initial_state=dict(zip(STATE_NAMES,model.initial_state().tolist())),
                   state_names=STATE_NAMES,source_hashes=before,
                   input_assumption='clock_scale=1 assumes compatible donor Int0 and receiver clock-threshold a.u.; uncalibrated',
                   solver=dict(method='DOP853',rtol=2e-7,atol=2e-9,max_step_min=args.max_step_min),
                   original_counter_solver=dict(method='RK45',rtol=1e-3,atol=1e-6),
                   environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform())))
        if source_hashes()!=before: raise RuntimeError('Source hash changed during run')
        write_json(out/'summary.json',summary)
        if not args.no_plot:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig,axes=plt.subplots(4,1,figsize=(12,10),sharex=True,constrained_layout=True)
            axes[0].plot(ts,hsol.y[1],color='black',label='ZMH Int0 (a.u.)');axes[0].legend()
            for s,label in zip((1-hsol.y[0],1-hsol.y[6],hsol.y[INDEX['b2_S']]),('S0 donor','S1 donor','S2 receiver')):
                axes[1].plot(ts,s,label=label)
            axes[1].axhspan(0,.3,color='green',alpha=.08);axes[1].axhspan(.7,1,color='blue',alpha=.08)
            axes[1].legend();axes[1].set_ylabel('LR fraction (raw)')
            rr=verdict['reads'];xs=[r['trough_h'] for r in rr];ys=[np.nan if r['value'] is None else r['value'] for r in rr]
            axes[2].plot(xs,ys,'o-',drawstyle='steps-post');axes[2].set_yticks(range(8));axes[2].set_ylabel('Window readout')
            for r in rr:
                if r['value'] is None: axes[2].text(r['trough_h'],3.5,'x',color='red')
            axes[3].plot(ts,signals['u1_target_au_per_h'],label='Int1 source (donor a.u./h)')
            axes[3].plot(ts,signals['u2_target_au_per_h'],label='Int2 source (receiver a.u./h)')
            axes[3].legend();axes[3].set_xlabel('Time (h)')
            fig.suptitle('Preliminary cross-model connection - not certified')
            fig.savefig(out/'overview.png',dpi=160);fig.savefig(out/'overview.pdf');plt.close(fig)
        write_json(out/'status.json',{'status':'COMPLETED','certified':None})
        paths=list(out.iterdir())+[ROOT/'hybrid_model.py',ROOT/'diagnostics.py',ROOT/'run_initial.py']
        write_json(out/'SHA256SUMS.json',{str(p.relative_to(ROOT)):sha256(p) for p in paths if p.is_file()})
        print(json.dumps({'output':str(out),**summary},ensure_ascii=False,indent=2),flush=True)
    except Exception as exc:
        write_json(out/'status.json',{'status':'FAILED','error':repr(exc)})
        raise


if __name__=='__main__': main()
