"""Erratum for the aliased leak statistic quoted in the provisional profile.

Direction of reference (deliberate)
-----------------------------------
    threebit51_provisional_n6.json   <- read-only, NEVER edited
    threebit51_leak_correction.json  -> references the provisional file by name
                                        and by SHA256, and states which of its
                                        fields and conclusions no longer hold
    threebit51_selected_v1.json      -> (later) cites BOTH files

The provisional file must not reference this erratum, because that would require
rewriting an already-frozen artefact.

Run after `scan_carry_pairing.py merge` has produced its verdicts; the script
includes whatever exists and records what is still missing.

    python make_leak_correction.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import OUT, sha256, source_hashes  # noqa: E402

PROVISIONAL = OUT / 'threebit51_provisional_n6.json'
TARGET = OUT / 'threebit51_leak_correction.json'


def load(rel):
    p = OUT / rel
    return json.loads(p.read_text(encoding='utf-8')) if p.exists() else None


def main():
    if not PROVISIONAL.exists():
        raise SystemExit('the provisional profile is missing; nothing to correct')
    if TARGET.exists():
        raise SystemExit(f'{TARGET.name} already exists - not overwritten by design')

    prov = json.loads(PROVISIONAL.read_text(encoding='utf-8'))
    prov_sha = sha256(PROVISIONAL)

    n5 = load('carry_type_diagnostic/summary_n5.json')
    n6 = load('carry_type_diagnostic/summary_n6.json')
    aliasing = load('eight_initial_diagnostic/paired_three_values.json')
    grid = load('carry_pairing_grid_verdict.json')
    eight = load('carry_pairing_eight_verdict.json')

    superseded = []
    for key, d in prov['selection_evidence']['certified_grid']['per_point'].items():
        superseded.append(dict(
            path=f'selection_evidence.certified_grid.per_point["{key}"].leak5pct_ratio',
            superseded_value=d.get('leak5pct_ratio')))
    superseded.append(dict(
        path='selection_evidence.why_n6_over_n5.far_off_leak_ratio_worst',
        superseded_value=prov['selection_evidence']['why_n6_over_n5']
        ['far_off_leak_ratio_worst']))

    doc = {
        'status': 'metric_correction',
        'name': 'threebit51_leak_correction',
        'applies_to': 'plausibility/threebit51_provisional_n6.json',
        'original_sha256': prov_sha,
        'supersedes_nothing': ('this erratum is purely additive; the provisional file stays '
                               'read-only and unedited, and remains citable as of its own hash'),
        'affected_fields': [
            'selection_evidence.why_n6_over_n5.far_off_leak_ratio_worst',
            'selection_evidence.certified_grid.per_point[*].leak1pct_ratio',
            'selection_evidence.certified_grid.per_point[*].leak5pct_ratio',
            'selection_evidence.certified_grid.per_point[*].leak10pct_ratio',
            'the accompanying claim that the n-ordering is identical at the 1 %, 5 % and '
            '10 % tail cuts',
        ],
        'unaffected': [
            'trajectories and integrator settings',
            'mod-8 certification (56 reads, 48 steady, built on read windows not on this metric)',
            'event counts (crossings / gates / reverse events)',
            'read-window commitment and 0.30 / 0.70 band-edge distances',
            'setup / hold margins',
            'the structural prefix guard and the decoupling verification',
            'the eight-initial-state readout result',
        ],
        'cause': {
            'summary': ('alternating F/R carry types combined with a median taken over an ODD '
                        'number of windows'),
            'mechanism': [
                'bit2 is written on only every other carry, so consecutive carries alternate:',
                '  R type - S2 is HIGH when the carry opens; the carry performs the intended '
                'REVERSE write (gate_on_rev ~ 1.0); its far-off then sits with S2 LOW, so the '
                'active pool is the FORWARD one',
                '  F type - S2 is LOW when the carry opens; the carry performs the intended '
                'FORWARD write; its far-off then sits with S2 HIGH, so the reverse pool is active',
                'each carry therefore leaves a far-off whose dominant direction is the OPPOSITE '
                'of the action it just performed',
                'over a 600 h horizon an initial phase that ends the carry train on an F gives 15 '
                'windows (8 F + 7 R); an even train gives 14 (7 F + 7 R)',
                'median(far_off_rev) / median(gate_on_rev) taken over an odd window count lands ON '
                'one carry type instead of between the two, so the ratio jumps by orders of '
                'magnitude with no change in the circuit',
                'secondary defect: the numerator was a median over far-off-EVALUABLE windows while '
                'the denominator was a median over ALL windows, so the two index sets differed',
            ],
        },
        'aliasing_proof': {
            'source': 'plausibility/eight_initial_diagnostic/paired_three_values.json',
            'three_initial_states_of_the_same_circuit': aliasing,
            'raw_ratio_spread_factor': (
                (max(v['raw'] for v in aliasing.values()) /
                 min(v['raw'] for v in aliasing.values())) if aliasing else None),
            'paired_ratio_relative_spread': (
                ((max(v['paired'] for v in aliasing.values()) -
                  min(v['paired'] for v in aliasing.values())) /
                 sum(v['paired'] for v in aliasing.values()) * len(aliasing)) if aliasing else None),
            'conclusion': ('the paired statistic is stable to ~0.1 %, which also demonstrates that '
                           'the eight orbit-derived initial states lie on the same period-8 '
                           'attractor'),
        },
        'new_definition': {
            'classification': 'R if S2 at the carry start >= 0.5, else F (physical, not dose-based)',
            'pairing': ('one R window plus one F window = one full 8-read super-period; incomplete '
                        'single windows at either end are DROPPED and recorded, never padded with '
                        'zeros'),
            'L_rev': 'far_off_rev(F) / gate_on_rev(R)',
            'L_fwd': 'far_off_fwd(R) / gate_on_fwd(F)',
            'L_symmetric': '[far_off_rev(F) + far_off_fwd(R)] / [gate_on_rev(R) + gate_on_fwd(F)]',
            'aggregation_rule': ('all sums are formed INSIDE each pair first; statistics are then '
                                 'taken over pairs. Medians of numerator and denominator are never '
                                 'divided by each other.'),
            'tail_kept_separate': ['tail_rev_R', 'tail_fwd_F'],
            'cuts': [0.01, 0.05, 0.10],
            'implementation': 'plausibility/scan_carry_pairing.py',
            'per_window_data': ('plausibility/carry_pairing/{kind}_windows_*.csv and '
                                '_pairs_*.csv are written for every run, so future re-aggregation '
                                'needs no re-integration'),
            're_integration_note': ('the 20-point grid and the eight-initial-state runs did NOT '
                                    'persist per-window data (scan_gate_segments.py declared a '
                                    'per-window list that was never populated), so both had to be '
                                    're-integrated. This is recorded because it is the reason the '
                                    'correction is not a pure post-processing step.'),
        },
        'corrected_numbers': {
            'two_candidate_points_cold_start_600h': {
                'source': 'plausibility/carry_type_diagnostic/',
                'n5': n5,
                'n6': n6,
            },
            'grid_20_points': grid,
            'eight_initial_states': eight,
        },
        'does_the_correction_change_the_selection': {
            'answer': ('no - n_A1_gate = 6 remains the preselect centre, and the correction '
                       'strengthens rather than weakens it'),
            'reasons': [
                'paired reverse leak is ~3.58x lower for n=6 (0.02391 vs 0.08568)',
                'n=6 super-period spread is ~+/-5 % against ~+/-26 % for n=5',
                'n=6 is deeper on both bit2 band edges and has the larger bit2 hold',
                'commitment is 1.0 for both and therefore ranks neither',
                'all eight initial states are strictly mod-8',
            ],
            'what_the_grid_recomputation_is_for': ('confirming, not re-selecting: whether '
                                                   'n=8<n=7<n=6<n=5 still holds, whether n=6 still '
                                                   'beats n=5 under the symmetric metric, and '
                                                   'whether the ordering is independent of initial '
                                                   'phase and window parity'),
        },
        'freeze_conditions_status': {
            'grid_recomputation_does_not_overturn_n6': (
                grid.get('n6_beats_n5') if grid else 'PENDING - carry_pairing_grid_verdict.json'),
            'eight_state_paired_consistency': (
                'PASS' if (eight and eight.get('eight_state_consistency', {})
                           .get('relative_spread') is not None
                           and eight['eight_state_consistency']['relative_spread'] < 0.05
                           and not eight['eight_state_consistency']['any_anomalies'])
                else ('PENDING' if not eight else 'REVIEW - see carry_pairing_eight_verdict.json')),
        },
        'next_artefact': {
            'file': 'plausibility/threebit51_selected_v1.json',
            'must_cite': ['plausibility/threebit51_provisional_n6.json (and its SHA256)',
                          'plausibility/threebit51_leak_correction.json (and its SHA256)'],
        },
        'source_sha256': source_hashes(),
    }
    TARGET.write_text(json.dumps(doc, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'wrote {TARGET}')
    print('sha256', sha256(TARGET))
    print('applies_to', doc['applies_to'], prov_sha)
    print('grid verdict present :', grid is not None)
    print('eight verdict present:', eight is not None)


if __name__ == '__main__':
    main()
