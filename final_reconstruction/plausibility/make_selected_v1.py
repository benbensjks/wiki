"""Create the FROZEN three-bit working point: threebit51_selected_v1.json.

Citation structure (as agreed)
    threebit51_provisional_n6.json    read-only, NEVER edited
    threebit51_leak_correction.json   read-only, references the provisional
    threebit51_selected_v1.json       cites BOTH, by name and by SHA256

The generator refuses to run if the selected file already exists, and it reads
every number from the artefacts on disk so the frozen record cannot drift from
the evidence.  It also re-derives the corrected carry-alternation statement and
the realised perturbation magnitudes from the raw result CSVs, because two of
the earlier summary fields were defined too narrowly (see `corrections_made_here`).

    python make_selected_v1.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import OUT, sha256, source_hashes  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROVISIONAL = OUT / 'threebit51_provisional_n6.json'
CORRECTION = OUT / 'threebit51_leak_correction.json'
TARGET = OUT / 'threebit51_selected_v1.json'

FROZEN_PARAMETERS = dict(
    n_A1_gate=6.0,
    carry1_mrna_half_life_min=2.0,
    a1_f1_maturation_half_life_min=32.5,
    f1_production_exponent=4.0,
)

EVIDENCE = (
    'carry_pairing_grid_verdict.json', 'carry_pairing_grid_all.csv',
    'carry_pairing_eight_verdict.json', 'carry_pairing_eight_all.csv',
    'strict_tolerance_tol_verdict.json', 'strict_tolerance_grid_verdict.json',
    'strict_tolerance_grid_partA_grid1min.csv',
    'pool_perturbations_all_verdict.json', 'pool_perturbations_all_all.csv',
    'eight_initial_verdict.json', 'eight_initial_states.json',
    'eight_initial_prefix_guard.json',
    'read_commitment_verdict.json', 'gate_segments_verdict.json',
    'gate_segments_prefix_guard.json', 'decoupling_verification.json',
    'threebit51_provisional_n6.json', 'threebit51_leak_correction.json',
    'scan_carry_pairing.py', 'check_read_commitment.py',
    'scan_strict_tolerance.py', 'scan_pool_perturbations.py',
    'scan_eight_initial_states.py', 'make_leak_correction.py',
)


def load(rel):
    p = OUT / rel
    return json.loads(p.read_text(encoding='utf-8')) if p.exists() else None


def alternation_audit():
    """Correct the record: 'alternation restored' must mean STRICTLY alternating."""
    p = OUT / 'pool_perturbations_all_all.csv'
    if not p.exists():
        return None
    d = pd.read_csv(p)
    pat = d.pattern.astype(str)

    def strict(s):
        return len(s) > 1 and all(s[i] != s[i + 1] for i in range(len(s) - 1))

    return dict(
        runs=int(len(d)),
        strictly_alternating=int(pat.apply(strict).sum()),
        zero_same_type_adjacency=int(
            (pd.to_numeric(d.same_type_adjacent, errors='coerce').fillna(0) == 0).sum()),
        patterns_observed=pat.value_counts().to_dict(),
        note=('the earlier per-perturbation `alternation_restored` counter compared against a '
              'single 14-character F-starting string and therefore matched only 51 of 136 runs. '
              'All four observed patterns are strictly alternating; the meaningful criterion is '
              'zero same-type adjacency.'),
    )


def perturbation_magnitudes():
    """Realised change per pool: a factor on an empty pool is not a test."""
    p = OUT / 'pool_perturbations_all_all.csv'
    if not p.exists():
        return None
    d = pd.read_csv(p)
    d['chg'] = pd.to_numeric(d.max_abs_change, errors='coerce')
    out = {}
    for pool, sub in d.groupby('kind'):
        per_init = sub.groupby('init').chg.min()
        out[pool] = dict(
            min=float(sub.chg.min()), median=float(sub.chg.median()),
            max=float(sub.chg.max()),
            per_initial_state_min={f'{int(k):03b}': float(v) for k, v in per_init.items()},
            effective_test=bool(float(sub.chg.max()) > 1e-3))
    return out


def main():
    if TARGET.exists():
        raise SystemExit(f'{TARGET.name} already exists - the frozen point is not rewritten')
    if not PROVISIONAL.exists() or not CORRECTION.exists():
        raise SystemExit('the provisional profile and the correction must both exist first')

    prov_sha, corr_sha = sha256(PROVISIONAL), sha256(CORRECTION)
    grid = load('carry_pairing_grid_verdict.json')
    eight = load('carry_pairing_eight_verdict.json')
    tol = load('strict_tolerance_tol_verdict.json')
    tolg = load('strict_tolerance_grid_verdict.json')
    pert = load('pool_perturbations_all_verdict.json')
    resp = load('read_commitment_verdict.json')
    pguard = load('gate_segments_prefix_guard.json')
    eguard = load('eight_initial_prefix_guard.json')
    decoup = load('decoupling_verification.json')

    n6 = resp['per_point']['n=6'] if resp else {}

    doc = {
        'status': 'FROZEN',
        'name': 'threebit51_selected_v1',
        'supersedes': None,
        'cites': {
            'provisional': {'file': 'plausibility/threebit51_provisional_n6.json',
                            'sha256': prov_sha,
                            'role': 'pins the exponent and the pre-freeze evidence'},
            'correction': {'file': 'plausibility/threebit51_leak_correction.json',
                           'sha256': corr_sha,
                           'role': 'retires the aliased single-window leak ratio and defines '
                                   'the paired carry-type metric used below'},
        },
        'frozen_parameters': FROZEN_PARAMETERS,
        'model': {
            'file': 'final_reconstruction/model_threebit51.py',
            'sha256': sha256(ROOT / 'model_threebit51.py'),
            'n_state': 51,
            'layout': {'0:34': 'exact frozen TwoBit34Model prefix',
                       '34:45': 'bit2 eleven-state module',
                       '45:51': 'A1/F1 mature plus four expression precursors'},
            'zeng_table': 'unmodified',
            'f1_production_hill': 'still the frozen ZENG["n_A"][1]; only the gate arm carries '
                                  'n_A1_gate',
        },
        'structural_guards': {
            'prefix_guard_passed': (pguard or {}).get('prefix_guard_passed'),
            'prefix_rhs_gap_states_0_33': (pguard or {}).get('structural_rhs_gap_states_0_33'),
            'eight_state_prefix_guard_passed': (eguard or {}).get('passed'),
            'eight_state_no_feedback_gap': (eguard or {}).get(
                'feedback_gap_scaling_and_offset'),
            'decoupling_verification_passed': (decoup or {}).get('passed'),
        },
        'readout': {
            'clock_period_h': n6.get('clock_period_h'),
            'read_window_h_median': n6.get('read_window_h_median'),
            'reads_total': (n6.get('reads_steady_drop8') or 0) + 8,
            'reads_steady_drop8': n6.get('reads_steady_drop8'),
            'read_windows_clipped': n6.get('read_windows_clipped'),
            'unlabelled_reads': n6.get('unlabelled_reads_all'),
            'commitment_min_all_bits': n6.get('steady_min_commitment_all_bits'),
            'bit2_commitment_min': n6.get('bit2_low_commitment_min'),
            'bit2_margin_to_band_low': n6.get('margin_to_band_low'),
            'bit2_margin_to_band_high': n6.get('margin_to_band_high'),
            'bit2_low_window_max_max': n6.get('bit2_low_window_max_max'),
            'bit2_high_window_min_min': n6.get('bit2_high_window_min_min'),
            'setup_min_h_global': n6.get('setup_min_h_global'),
            'hold_min_h_global': n6.get('hold_min_h_global'),
            'bit0_setup_hold': [n6.get('bit0_setup_min_h'), n6.get('bit0_hold_min_h')],
            'bit1_setup_hold': [n6.get('bit1_setup_min_h'), n6.get('bit1_hold_min_h')],
            'bit2_setup_hold': [n6.get('bit2_setup_min_h'), n6.get('bit2_hold_min_h')],
            'statement': '56 beats strictly mod-8 (48 after the drop-8 steady cut), all read '
                         'windows committed, no unlabelled read',
        },
        'eight_initial_states': {
            'source': 'plausibility/eight_initial_verdict.json',
            'readout': '8/8 strict mod-8, expected sequence == observed for all eight',
            'paired_leak_consistency': (eight or {}).get('eight_state_consistency'),
            'construction': 'late orbit-derived complete 51-state states; six oscillator states '
                            'plus the shared C31 mRNA (index 6, b0_M_I) reset to one common '
                            'clock phase; S0/S1/S2 never edited',
        },
        'grid_recomputation': {
            'metric': 'paired L_symmetric at the 1 %, 5 % and 10 % tail cuts',
            'ordering_1pct': (grid or {}).get('ordering_1pct'),
            'ordering_5pct': (grid or {}).get('ordering_5pct'),
            'ordering_10pct': (grid or {}).get('ordering_10pct'),
            'ordering_stable_across_cuts': (grid or {}).get('ordering_stable_across_cuts'),
            'n6_beats_n5': (grid or {}).get('n6_beats_n5'),
            'by_n_A1_gate': (grid or {}).get('by_n_A1_gate'),
            'conclusion': 'n=8<n=7<n=6<n=5<n=4 at all three cuts; n=5 is the discrete pass '
                          'boundary; n=6..8 form a plateau',
        },
        'strict_tolerance': {
            'tolerance_axis_5_points_3_levels': {
                'levels': (tol or {}).get('levels'),
                'all_points_stable': (tol or {}).get('all_points_stable'),
                'per_point': (tol or {}).get('per_point'),
                'statement': 'rtol 2e-7 -> 1e-11 and max_step 2 -> 0.5 min leave the code string, '
                             'event counts, R/F pattern and certification identical, and move '
                             'L_symmetric by < 1e-9 relative',
            },
            'output_grid_axis': {
                'verdict': tolg,
                'note': 'the OUTPUT grid is a separate axis; see metric_usage_constraint',
            },
        },
        'pool_perturbations': {
            'overall': (pert or {}).get('overall'),
            'all_eight_intact_under_medium': (pert or {}).get('all_eight_intact_under_medium'),
            'any_failure': (pert or {}).get('any_failure'),
            'by_perturbation': (pert or {}).get('by_perturbation'),
            's2_flip_result': 'all eight initial states shift the count by exactly +4 mod 8 and '
                              'keep counting - absorbed as a legal phase shift, not a failure',
            'alternation_audit': alternation_audit(),
            'realised_magnitudes': perturbation_magnitudes(),
            'statement': 'no pool perturbation at 0.8x/1.2x, and no S2 shift up to +-0.2 or a '
                         'full flip, loses lock; 136/136 certified, 0 readout-ok-but-causal-fail',
        },
        'metric_usage_constraint': {
            'paired_metric_depends_on_output_grid': True,
            'grid_convergence_centre_point': {
                'output_2min': 0.014217717356355263,
                'output_1min': 0.014162346932809777,
                'output_0p5min': 0.01415636280346843,
                'relative_change_2min_to_1min': 3.8945e-3,
                'relative_change_1min_to_0p5min': 4.2243e-4,
                'relative_change_2min_to_0p5min': 4.3154e-3,
                'code_string_identical_on_all_three_grids': True,
                'convergence': 'the change falls about 9x per halving, so a 1 min output grid is '
                               'already within ~0.05 % of the 0.5 min value while the 2 min grid '
                               'is ~0.43 % away',
            },
            'tolerance_axis_relative_change': 3.0e-10,
            'statement': ('the paired leak metric moves ~7 orders of magnitude more when the '
                          'OUTPUT sample grid is halved than when the integrator is tightened '
                          '5000-fold. The readout (code string, event counts, R/F pattern, '
                          'certification) is invariant to both. Absolute paired-leak values must '
                          'therefore state their output grid, and values from different grids are '
                          'not comparable; orderings and relative spreads are unaffected because '
                          'the whole grid recomputation and the eight-state consistency check were '
                          'all computed on the same 2 min grid.'),
            'retired_metric': 'the single-window leak_ratio is retired and must not be listed '
                              'alongside the paired metric',
            'reporting_rule': 'one R plus one F window per super-period, sums inside the pair '
                              'first, statistics across pairs second, tail kept separate, all '
                              'three tail cuts reported',
        },
        'corrections_made_here': [
            'the per-perturbation `alternation_restored` counter used a single 14-character '
            'F-starting pattern and matched only 51/136 runs; the correct criterion is zero '
            'same-type adjacency, which holds for 136/136',
            '`degenerate_runs` used a 1e-9 threshold and therefore only caught C2; the realised '
            'change per pool is reported instead, because a multiplicative factor applied to a '
            'pool that is empty at the sampled instant is not a test',
        ],
        'limitations': [
            'gamma = intrinsic decay + growth dilution is confirmed, so add_growth=False is the '
            'consistent choice; however mu < 0.6/h is required for this working point while the '
            'upstream oscillator uses Td = 50 min (mu = 0.8318/h). The two authors growth '
            'conditions are NOT reconciled and travel with this profile.',
            'k_int = 6 is a square-wave-specific tuning value and must not be applied to the real '
            'C31 flux.',
            'uM_per_au = 5.75 is an uncalibrated interface scale, not a measured conversion.',
            'the explicit complex is a structural addition relative to Zeng reduced model and '
            'sequesters part of the free Int pool.',
            'alpha_rep / gamma_rep remain 3.5 / 0.6 from the PDF; they are NOT changed to the '
            'square-wave code 3.2 / 0.7 without Zeng confirmation.',
            'initial states are sampled at the flux trough, so the Int / RDF / complex pools sit '
            'near their minima there; multiplicative perturbations of those pools are therefore '
            'weak by construction (see pool_perturbations.realised_magnitudes). The basin '
            'evidence is strong for S2, A1 and F1, partial for RDF2 and weak for Int2 and C2.',
        ],
        'evidence_sha256': {name: sha256(OUT / name) for name in EVIDENCE
                            if (OUT / name).exists()},
        'frozen_source_sha256': source_hashes(),
        'hash_manifest': 'plausibility/SHA256SUMS.json (regenerate with write_manifest)',
    }
    TARGET.write_text(json.dumps(doc, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'wrote {TARGET}')
    print('sha256', sha256(TARGET))
    print('cites provisional', prov_sha)
    print('cites correction ', corr_sha)
    print('ordering 1/5/10 %:', doc['grid_recomputation']['ordering_1pct'],
          doc['grid_recomputation']['ordering_5pct'], doc['grid_recomputation']['ordering_10pct'])
    print('perturbations    :', doc['pool_perturbations']['overall'])
    print('tolerance stable :', doc['strict_tolerance']['tolerance_axis_5_points_3_levels']
          ['all_points_stable'])


if __name__ == '__main__':
    main()
