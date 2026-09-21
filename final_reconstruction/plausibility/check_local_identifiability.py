"""Local identifiability of the frozen three-bit working point.

What question this answers
--------------------------
NOT "which parameters matter most" (that is sensitivity ranking).  It answers
"which parameters can be DISTINGUISHED FROM EACH OTHER, and by which
observables", i.e. it looks for directions in parameter space that the data
cannot separate.

Two observation models are kept strictly apart, because pooling them is how an
analysis overstates experimental identifiability:

  theoretical (full state)  every quantity the simulation can reach:
      S0 S1 S2, I0 R0 C0, I1 R1 C1, I2 R2 C2, A1 F1, g1, J_rev1, J_rev2
  experimental (realistic)  only what an experiment could plausibly report:
      S0 S1 S2 (DNA-configuration reporters), I0 (C31/Int reporter), A1 (output),
      sampled every 30 min instead of every 2 min

The experimental set deliberately EXCLUDES the RDF pools, the complexes, the
immature/mRNA states, the internal gate g1 and all fluxes.  A parameter that is
identifiable only in the theoretical model is NOT experimentally identifiable.

Method
------
Central-difference sensitivity of the observable vector to each LOG parameter:

    S[i,j] = (y_i(p_j(1+h)) - y_i(p_j(1-h))) / (scale_i * (ln(1+h) - ln(1-h)))

with scale_i = max_t |y_i| on the baseline trajectory, so every observable
contributes in relative units.  Then

    S = U Sigma V^T
    singular spectrum, effective rank (#{sigma > rank_tol * sigma_max}),
    condition number sigma_max/sigma_min,
    F = S^T S, correlation matrix from pinv(F),
    weakest identifiable direction = right singular vector of sigma_min,
    leave-one-observable-out effect on sigma_min,
    finite-difference convergence over h in {1 %, 2 %, 5 %}.

Parameters: the 9 named in the task spec.  `--with-translation` adds
translation_h as a 10th, because the predicted uM_per_au <-> translation-strength
confusion cannot be tested inside the 9-parameter set.

Assumptions this result is conditional on (recorded with the verdict)
    perturbation scales; observable sets; sampling cadence; an IDEALIZED
    NOISELESS observation model; the assumption that states are absolutely
    quantified in model a.u. (scale_i uses the model's own magnitudes); the
    initial condition (the frozen default, not fitted).  These are assumptions of
    the EXPERIMENT, not of the model, and identifiability only holds within them.

    python check_local_identifiability.py run   --shard S --nshards N
    python check_local_identifiability.py merge --nshards N
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (OUT as PLAUS_OUT, frozen_carry0,  # noqa: E402
                                 frozen_extension, sha256, source_hashes,
                                 write_manifest)
from model_threebit51 import ThreeBit51Model, ThreeBitCarryParameters  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = PLAUS_OUT / 'local_identifiability'
JOBS = OUT / 'jobs'
MODEL_FILE = ROOT / 'model_threebit51.py'

HOURS = 350.0
ANALYSIS_START_H = 200.0
SAMPLE_MIN = 2.0
MAX_STEP_MIN = 2.0
STEPS = (0.01, 0.02, 0.05)
RANK_TOL = 1e-3
EXPERIMENTAL_SAMPLE_MIN = 30.0

PARAMS_CORE = (
    ('uM_per_au', 'extension', 'uM_per_au'),
    ('n_A1_gate', 'ctor', 'n_A1_gate'),
    ('carry1_mrna_half_life_min', 'carry1', 'mrna_half_life_min'),
    ('A1_maturation_half_life_min', 'carry1', 'activator_maturation_half_life_min'),
    ('F1_maturation_half_life_min', 'carry1', 'repressor_maturation_half_life_min'),
    ('complex_on', 'extension', 'complex_on_au_inv_h'),
    ('complex_off', 'extension', 'complex_off_h'),
    ('complex_decay', 'extension', 'complex_decay_h'),
    ('K_D_comp_scale', 'zeng', 'K_D_comp'),
)
PARAMS_EXTRA = (('translation_h', 'extension', 'translation_h'),)

NOMINAL = dict(uM_per_au=5.75, n_A1_gate=6.0, carry1_mrna_half_life_min=2.0,
               A1_maturation_half_life_min=32.5, F1_maturation_half_life_min=32.5,
               complex_on=0.1, complex_off=1.0, complex_decay=1.0, K_D_comp_scale=1.0,
               translation_h=30.0)

THEORETICAL_OBS = (('S0', 16), ('S1', 27), ('S2', 44),
                   ('I0', 8), ('R0', 14), ('C0', 15),
                   ('I1', 19), ('R1', 25), ('C1', 26),
                   ('I2', 36), ('R2', 42), ('C2', 43),
                   ('A1', 45), ('F1', 46))
THEORETICAL_SIGNALS = ('g1', 'J_rev1', 'J_rev2')
EXPERIMENTAL_OBS = (('S0', 16), ('S1', 27), ('S2', 44), ('I0', 8), ('A1', 45))


def param_table(with_translation=False):
    return tuple(PARAMS_CORE) + (tuple(PARAMS_EXTRA) if with_translation else ())


def build_jobs(params):
    jobs = [('baseline', {})]
    for name, _, _ in params:
        for h in STEPS:
            jobs.append((f'{name}|plus|{h:g}', {name: NOMINAL[name] * (1.0 + h)}))
            jobs.append((f'{name}|minus|{h:g}', {name: NOMINAL[name] * (1.0 - h)}))
    return jobs


def job_file(label):
    return JOBS / (label.replace('|', '__') + '.npz')


def build_with(overrides):
    """Build the frozen working point with the given overrides.

    No global patch is left behind: the only table entry that must be patched
    (ZENG['K_D_comp']) is read in Model.__init__, so it is restored immediately
    after construction.
    """
    import model as M
    kinds = {n: k for n, k, _ in param_table(True)}
    attrs = {n: a for n, _, a in param_table(True)}
    ext, carry1, n_gate = frozen_extension(), frozen_carry0(), None
    for name, value in overrides.items():
        k = kinds[name]
        if k == 'extension':
            ext = replace(ext, **{attrs[name]: value})
        elif k == 'carry1':
            carry1 = replace(carry1, **{attrs[name]: value})
        elif k == 'ctor':
            n_gate = float(value)
    kd_scale = float(overrides.get('K_D_comp_scale', 1.0))
    saved = M.ZENG['K_D_comp']
    M.ZENG['K_D_comp'] = saved * kd_scale
    try:
        model = ThreeBit51Model(
            extension=ext,
            carry=ThreeBitCarryParameters(carry0=frozen_carry0(), carry1=carry1),
            n_A1_gate=n_gate)
    finally:
        M.ZENG['K_D_comp'] = saved
    return dict(model=model, k_D_comp_used=float(saved * kd_scale),
                q=float(model.base.q), K_complex=float(model.base.K_complex))


def extract_observables(sol, sig, which, signals, t_min, every_min):
    t = sol.t
    mask = t >= t_min
    cols, names = [], []
    for name, idx in which:
        cols.append(np.asarray(sol.y[idx])[mask])
        names.append(name)
    for s in signals:
        cols.append(np.asarray(sig[s])[mask])
        names.append(s)
    M = np.vstack(cols).T
    tt = t[mask]
    if every_min > SAMPLE_MIN:
        stride = max(1, int(round(every_min / SAMPLE_MIN)))
        M, tt = M[::stride], tt[::stride]
    return dict(Y=M, names=names, t=tt, sample_min=every_min)


def run_one(overrides, hours, t_min):
    """ONE integration; both observation models are extracted from it."""
    r = build_with(overrides)
    sol = r['model'].simulate(hours=hours, sample_min=SAMPLE_MIN, max_step_min=MAX_STEP_MIN)
    sig = r['model'].diagnostic_signals(sol.y)
    return dict(
        theoretical=extract_observables(sol, sig, THEORETICAL_OBS, THEORETICAL_SIGNALS,
                                        t_min, SAMPLE_MIN),
        experimental=extract_observables(sol, sig, EXPERIMENTAL_OBS, (), t_min,
                                         EXPERIMENTAL_SAMPLE_MIN),
        meta=dict(q=r['q'], K_complex=r['K_complex'],
                  k_D_comp_used=r['k_D_comp_used']))


# ------------------------------------------------------------------ run/merge
def run(args):
    params = select_params(args)
    jobs = build_jobs(params)
    mine = jobs[args.shard::args.nshards]
    JOBS.mkdir(parents=True, exist_ok=True)
    tag = f'shard{args.shard:02d}of{args.nshards:02d}'
    started = time.perf_counter()
    for label, ov in mine:
        res = run_one(ov, args.hours, args.analysis_start_h)
        np.savez_compressed(job_file(label), theo_Y=res['theoretical']['Y'],
                            exp_Y=res['experimental']['Y'])
        job_file(label).with_suffix('.json').write_text(json.dumps(dict(
            label=label, overrides={k: float(v) for k, v in ov.items()},
            theo_names=res['theoretical']['names'],
            exp_names=res['experimental']['names'],
            theo_n=int(res['theoretical']['Y'].shape[0]),
            exp_n=int(res['experimental']['Y'].shape[0]),
            meta=res['meta']), ensure_ascii=False, indent=2), encoding='utf-8')
        print(f'[{tag}] {label:<44} done ({time.perf_counter()-started:.0f}s)', flush=True)
    (OUT / f'meta_li_{tag}.json').write_text(json.dumps(dict(
        shard=args.shard, nshards=args.nshards, jobs=len(mine), hours=args.hours,
        analysis_start_h=args.analysis_start_h, steps=list(STEPS),
        parameters=[n for n, _, _ in params],
        model_sha256=sha256(MODEL_FILE), source_sha256=source_hashes()),
        ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'[{tag}] wrote {len(mine)} jobs')


def load_jobs(params):
    theo, exp = {}, {}
    for label, _ in build_jobs(params):
        f = job_file(label)
        if not f.exists():
            raise SystemExit(f'missing job file {f}')
        z = np.load(f)
        meta = json.loads(f.with_suffix('.json').read_text(encoding='utf-8'))
        theo[label] = dict(Y=z['theo_Y'], names=meta['theo_names'])
        exp[label] = dict(Y=z['exp_Y'], names=meta['exp_names'])
    return theo, exp


def select_params(args):
    params = param_table(args.with_translation)
    if args.params:
        want = {x.strip() for x in args.params.split(',') if x.strip()}
        unknown = want - {n for n, _, _ in param_table(True)}
        if unknown:
            raise SystemExit(f'unknown parameter(s): {sorted(unknown)}')
        params = tuple(p for p in params if p[0] in want)
        if not params:
            raise SystemExit('no parameters selected')
    return params


# ------------------------------------------------------------------- analysis
def central_sensitivity(plus, minus, base_scale, h):
    """(n_time, n_obs) normalised central-difference sensitivity for one parameter.

    base_scale has one entry per OBSERVABLE and broadcasts along axis 1; the time
    axis is not scaled.
    """
    dln = np.log(1.0 + h) - np.log(1.0 - h)
    n = min(plus.shape[0], minus.shape[0])
    return (plus[:n] - minus[:n]) / (base_scale[None, :] * dln)


def stacked_sensitivity(results, params, h):
    """(n_time * n_obs, n_params).  Rows are observable-minor, so the rows of
    observable i are exactly `i::n_obs`."""
    base = results['baseline']
    scale = np.max(np.abs(base['Y']), axis=0)
    scale[scale <= 0] = 1.0
    cols = []
    for p in params:
        Sj = central_sensitivity(results[f'{p}|plus|{h:g}']['Y'],
                                 results[f'{p}|minus|{h:g}']['Y'], scale, h)
        cols.append(Sj.reshape(-1))
    S = np.column_stack(cols)
    if not np.all(np.isfinite(S)):
        raise SystemExit('non-finite sensitivity entry')
    return S, scale


def analyse(results, params, names, h, rank_tol=RANK_TOL):
    S, scale = stacked_sensitivity(results, params, h)
    n_obs = len(names)
    F = S.T @ S
    U, sv, Vt = np.linalg.svd(S, full_matrices=False)
    rank = int((sv > rank_tol * sv[0]).sum())
    cond = float(sv[0] / sv[-1]) if sv[-1] > 0 else float('inf')
    C = np.linalg.pinv(F)
    d = np.sqrt(np.maximum(np.diag(C), 1e-300))
    corr = C / np.outer(d, d)
    return dict(S=S, sv=sv, rank=rank, cond=cond, corr=corr, weakest=Vt[-1],
                names=list(names), params=list(params), scale=scale, n_obs=n_obs)


def rms_summary(S, names, params):
    n_obs = len(names)
    n_t = S.shape[0] // n_obs
    out = np.zeros((n_obs, len(params)))
    for i in range(n_obs):
        out[i, :] = np.sqrt(np.mean(S[i::n_obs, :] ** 2, axis=0))
    df = pd.DataFrame(out, columns=list(params))
    df.insert(0, 'observable', list(names))
    df.insert(1, 'n_timepoints', n_t)
    return df


def observable_contribution(S, names):
    n_obs = len(names)
    sv_full = np.linalg.svd(S, full_matrices=False)[1]
    rows = []
    for i, nm in enumerate(names):
        keep = np.arange(S.shape[0]) % n_obs != i
        sv = np.linalg.svd(S[keep, :], full_matrices=False)[1]
        rows.append(dict(observable=nm, row_norm=float(np.linalg.norm(S[i::n_obs, :])),
                         sigma_min_full=float(sv_full[-1]),
                         sigma_min_without=float(sv[-1]),
                         sigma_min_gain=float(sv[-1] - sv_full[-1]),
                         rank_without=int((sv > RANK_TOL * sv[0]).sum()),
                         variance_fraction_lost=float(
                             1.0 - np.sum(sv ** 2) / max(np.sum(sv_full ** 2), 1e-300))))
    return rows


def merge(args):
    params = select_params(args)
    names_all = [n for n, _, _ in params]
    theo, exp = load_jobs(params)
    verdict = dict(
        hours=args.hours, analysis_start_h=args.analysis_start_h,
        sample_min=SAMPLE_MIN, experimental_sample_min=EXPERIMENTAL_SAMPLE_MIN,
        steps=list(STEPS), rank_tol=RANK_TOL, parameters=names_all,
        nominal_values={n: NOMINAL[n] for n in names_all},
        model_sha256=sha256(MODEL_FILE), source_sha256=source_hashes(),
        observation_models=dict(
            theoretical=dict(observables=[n for n, _ in THEORETICAL_OBS] +
                             list(THEORETICAL_SIGNALS), sample_min=SAMPLE_MIN),
            experimental=dict(observables=[n for n, _ in EXPERIMENTAL_OBS],
                              sample_min=EXPERIMENTAL_SAMPLE_MIN,
                              excluded=['RDF pools', 'complexes', 'immature and mRNA states',
                                        'internal gate g1', 'all fluxes'])),
        conditional_on=dict(
            perturbation_scale_pct=[s * 100 for s in STEPS],
            observables='listed under observation_models; the experimental set is a subset',
            sampling='2 min theoretical / 30 min experimental',
            noise='IDEALIZED NOISELESS observation; no measurement-error model',
            quantitation=('states treated as absolutely quantified in model a.u.; scale_i is '
                          'the baseline trajectory maximum of each observable'),
            initial_condition=('the frozen model default initial state, not fitted; a different '
                               'or unknown initial condition changes the sensitivities'),
            note=('these are assumptions of the EXPERIMENT, not of the model; identifiability '
                  'holds only within them')),
        results={})

    for tag, res in (('theoretical', theo), ('experimental', exp)):
        names = res['baseline']['names']
        entry = dict(observables=list(names), n_observables=len(names),
                     n_timepoints=int(res['baseline']['Y'].shape[0]), per_step={})
        sens_frames, corr_frames = {}, {}
        contrib_rows = []
        for h in STEPS:
            a = analyse(res, names_all, names, h)
            entry['per_step'][f'{h:g}'] = dict(
                singular_values=[float(x) for x in a['sv']],
                effective_rank=a['rank'], condition_number=a['cond'],
                weakest_direction=dict(zip(a['params'], [float(x) for x in a['weakest']])),
                normalised_singular_values=[float(x / a['sv'][0]) for x in a['sv']])
            sens_frames[f'{h:g}'] = rms_summary(a['S'], names, a['params'])
            corr_frames[f'{h:g}'] = pd.DataFrame(a['corr'], index=a['params'],
                                                columns=a['params'])
            if h == STEPS[0]:
                contrib_rows = observable_contribution(a['S'], names)
                entry['pairwise_correlations'] = a['corr'].tolist()
        entry['finite_difference_convergence'] = []
        for h1, h2 in zip(STEPS, STEPS[1:]):
            S1 = sens_frames[f'{h1:g}'].drop(columns=['observable', 'n_timepoints']).to_numpy()
            S2 = sens_frames[f'{h2:g}'].drop(columns=['observable', 'n_timepoints']).to_numpy()
            entry['finite_difference_convergence'].append(dict(
                from_pct=h1 * 100, to_pct=h2 * 100,
                relative_matrix_change=float(np.linalg.norm(S1 - S2) /
                                             max(np.linalg.norm(S1), 1e-300))))
        entry['observable_contribution'] = contrib_rows
        verdict['results'][tag] = entry
        for h in STEPS:
            sens_frames[f'{h:g}'].to_csv(OUT / f'sensitivity_{tag}_h{h:g}.csv', index=False,
                                         encoding='utf-8')
            corr_frames[f'{h:g}'].to_csv(OUT / f'correlations_{tag}_h{h:g}.csv',
                                         encoding='utf-8')
        pd.DataFrame(contrib_rows).to_csv(OUT / f'observable_contribution_{tag}.csv',
                                          index=False, encoding='utf-8')

    th = verdict['results']['theoretical']
    corr_th = np.asarray(th['pairwise_correlations'])
    idx = {n: i for i, n in enumerate(names_all)}

    def rho(a, b):
        return None if a not in idx or b not in idx else float(corr_th[idx[a], idx[b]])

    q_terms = ['complex_on', 'complex_off', 'complex_decay']
    verdict['predicted_confusions'] = dict(
        q_vs_kd_comp=dict(
            note=('model.py: q = kon/(koff+complex_decay+growth) and K_complex = q*K_D_comp; '
                  'the reverse Hill reads C/K_complex, so q and K_D_comp enter ONLY as a '
                  'product. Only K_complex is identifiable; the individual factors are not.'),
            rho_K_D_comp_scale_vs_complex_on=rho('K_D_comp_scale', 'complex_on'),
            rho_K_D_comp_scale_vs_complex_off=rho('K_D_comp_scale', 'complex_off'),
            rho_K_D_comp_scale_vs_complex_decay=rho('K_D_comp_scale', 'complex_decay')),
        kon_koff_decay=dict(
            correlations={f'{a}~{b}': rho(a, b) for a in q_terms for b in q_terms if a < b},
            note='these three reach the reverse Hill only through q'),
        uM_per_au_vs_translation=(dict(rho=rho('uM_per_au', 'translation_h'),
                                       tested=True)
                                  if 'translation_h' in idx else
                                  dict(tested=False, reason=(
                                      'translation_h is not in the 9-parameter set; rerun with '
                                      '--with-translation to test this confusion'))),
        mRNA_vs_maturation=dict(
            rho_mRNA_vs_A1_mat=rho('carry1_mrna_half_life_min', 'A1_maturation_half_life_min'),
            rho_mRNA_vs_F1_mat=rho('carry1_mrna_half_life_min', 'F1_maturation_half_life_min'),
            rho_A1mat_vs_F1mat=rho('A1_maturation_half_life_min', 'F1_maturation_half_life_min')),
        n_A1_gate_vs_scale=dict(
            rho_vs_uM_per_au=rho('n_A1_gate', 'uM_per_au'),
            rho_vs_complex_on=rho('n_A1_gate', 'complex_on'),
            note='n_A1_gate is an exponent; a scale confusion shows as large |rho| with a scale'),
    )
    h0 = f'{STEPS[0]:g}'
    verdict['headline'] = dict(
        theoretical=dict(effective_rank=th['per_step'][h0]['effective_rank'],
                         condition_number=th['per_step'][h0]['condition_number'],
                         singular_values=th['per_step'][h0]['singular_values'],
                         weakest_direction=th['per_step'][h0]['weakest_direction']),
        experimental=dict(effective_rank=verdict['results']['experimental']
                          ['per_step'][h0]['effective_rank'],
                          condition_number=verdict['results']['experimental']
                          ['per_step'][h0]['condition_number'],
                          singular_values=verdict['results']['experimental']
                          ['per_step'][h0]['singular_values'],
                          weakest_direction=verdict['results']['experimental']
                          ['per_step'][h0]['weakest_direction']),
        rank_lost_by_restricting_to_experimental=(
            th['per_step'][h0]['effective_rank'] -
            verdict['results']['experimental']['per_step'][h0]['effective_rank']))
    (OUT / 'local_identifiability_verdict.json').write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding='utf-8')
    write_manifest()
    print(json.dumps(dict(headline=verdict['headline'],
                          predicted_confusions=verdict['predicted_confusions']),
                     ensure_ascii=False, indent=2))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('--shard', type=int, required=True)
    r.add_argument('--nshards', type=int, required=True)
    r.add_argument('--hours', type=float, default=HOURS)
    r.add_argument('--analysis-start-h', type=float, default=ANALYSIS_START_H)
    r.add_argument('--with-translation', dest='with_translation', action='store_true',
                   default=False)
    r.add_argument('--params', type=str, default='')
    m = sub.add_parser('merge')
    m.add_argument('--nshards', type=int, required=True)
    m.add_argument('--hours', type=float, default=HOURS)
    m.add_argument('--analysis-start-h', type=float, default=ANALYSIS_START_H)
    m.add_argument('--with-translation', dest='with_translation', action='store_true',
                   default=False)
    m.add_argument('--params', type=str, default='')
    args = ap.parse_args()
    (run if args.cmd == 'run' else merge)(args)


if __name__ == '__main__':
    main()
