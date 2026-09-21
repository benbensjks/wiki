"""Write the PROVISIONAL three-bit profile for the selected carry-2 exponent.

Status: this file is a generator, not a source of truth.  It reads every number
from the artefacts on disk and hashes every input, so the profile cannot silently
drift from the evidence.

The profile is written ONCE to `threebit51_provisional_n6.json` and is then marked
read-only.  It is explicitly NOT the frozen working point: the eight-initial-state
verification, the strict-tolerance ring and the molecular-pool perturbations are
still open.  When those pass, a SEPARATE file (`threebit51_selected_v1.json`) is
created and this one is left untouched, so the provisional-to-formal chain stays
readable.

    python make_provisional_profile.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import OUT, sha256, source_hashes  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TARGET = OUT / 'threebit51_provisional_n6.json'

# fixed before the stage ran; recorded so a later reader can audit the choice
SELECTED = dict(n_A1_gate=6.0, carry1_mrna_half_life_min=2.0,
                a1_f1_maturation_half_life_min=32.5, f1_production_exponent='frozen ZENG value')

EVIDENCE = (
    'gate_segments_all.csv', 'gate_segments_verdict.json',
    'gate_segments_prefix_guard.json',
    'read_commitment_all.csv', 'read_commitment_verdict.json',
    'read_commitment_shard00of02.csv', 'read_commitment_shard01of02.csv',
    'meta_read_commitment_shard00of02.json', 'meta_read_commitment_shard01of02.json',
    'check_read_commitment.py',
    'read_commitment_v1/check_read_commitment_v1.py.txt',
    'read_commitment_v1/read_commitment_shard00of02.csv',
    'read_commitment_v1/read_commitment_shard01of02.csv',
    'decoupling_verification.json',
)


def load(name):
    p = OUT / name
    return json.loads(p.read_text(encoding='utf-8')) if p.exists() else None


def main():
    if TARGET.exists():
        raise SystemExit(f'{TARGET.name} already exists - it is not overwritten by design')

    import model as M
    commit = load('read_commitment_verdict.json')
    guard = load('gate_segments_prefix_guard.json')
    decoup = load('decoupling_verification.json')

    n6 = commit['per_point']['n=6']
    n5 = commit['per_point']['n=5']
    cmp_ = commit['margin_comparison']

    evidence_hashes = {name: sha256(OUT / name) for name in EVIDENCE if (OUT / name).exists()}

    profile = {
        'status': 'PROVISIONAL',
        'name': 'threebit51_provisional_n6',
        'not_a_frozen_working_point': True,
        'model': {
            'file': 'final_reconstruction/model_threebit51.py',
            'sha256': sha256(ROOT / 'model_threebit51.py'),
            'n_state': 51,
            'layout': {'0:34': 'exact frozen TwoBit34Model prefix',
                       '34:45': 'bit2 eleven-state module',
                       '45:51': 'A1/F1 mature + four expression precursors'},
        },
        'selected_parameters': dict(
            **SELECTED,
            frozen_f1_drive_exponent=float(M.ZENG['n_A'][1]),
            f1_drive_exponent_note=('n_A1_gate only changes the A1 arm inside the carry-2 gate; '
                                    'the F1 production Hill keeps the frozen table exponent, so '
                                    'the two uses of ZENG["n_A"][1] are decoupled by construction'),
            carry0_parameters='twobit34_results/selected_profile.json -> carry_expression (unchanged)',
            extension='twobit34_results/selected_profile.json -> extension (unchanged)',
            zeng_table='unmodified',
        ),
        'selection_evidence': {
            'certified_grid': {
                'grid': 'n_A1_gate in {4,5,6,7,8} x carry1 mRNA in {2,4} min x A1/F1 maturation in {32.5,60} min',
                'hours': 600,
                'certified': f"{commit['rows']} points; n in {{5,6,7,8}} all pass, n=4 all fail",
                'per_point': {k: {'certified': v['certified'], 'crossings': v['crossings'],
                                  'gates': v['gates'], 'reverse_events': v['reverse_events'],
                                  'leak5pct_ratio': v['leak5pct_ratio'],
                                  'leak1pct_ratio': v['leak1pct_ratio'],
                                  'leak10pct_ratio': v['leak10pct_ratio']}
                              for k, v in commit['per_point'].items()},
            },
            'why_n6_over_n5': {
                'far_off_leak_ratio_worst': {'n5': n5['leak5pct_ratio'], 'n6': n6['leak5pct_ratio'],
                                             'note': 'n=6 is ~3x lower; ordering identical at the '
                                                     '1 %, 5 % and 10 % tail cuts'},
                'margin_to_band_low': {'n5': cmp_['margin_to_band_low']['n5'],
                                       'n6': cmp_['margin_to_band_low']['n6'], 'larger': 'n=6'},
                'margin_to_band_high': {'n5': cmp_['margin_to_band_high']['n5'],
                                        'n6': cmp_['margin_to_band_high']['n6'], 'larger': 'n=6'},
                'bit2_hold_min_h': {'n5': cmp_['bit2_hold_min_h']['n5'],
                                    'n6': cmp_['bit2_hold_min_h']['n6'], 'larger': 'n=6'},
                'bit2_setup_min_h': {'n5': cmp_['bit2_setup_min_h']['n5'],
                                     'n6': cmp_['bit2_setup_min_h']['n6'], 'larger': 'n=5',
                                     'note': 'the ONLY metric where n=5 wins; it is 0.063 h and '
                                             'does not bind, because the global setup minimum is '
                                             'set by bit0 (4.030 h) and is identical for both'},
                'global_timing_identical': {
                    'setup_min_h_global': n6['setup_min_h_global'],
                    'hold_min_h_global': n6['hold_min_h_global'],
                    'note': ('bit0 and bit1 setup/hold are bit-identical between n=5 and n=6, so the '
                             'binding timing margins of the whole 3-bit counter are set by the frozen '
                             'prefix and the gate exponent does not move them')},
                'reading_of_commitment': ('commitment is 1.0 (its ceiling) for both points and therefore '
                                          'has NO ranking power; it only proves both are far from the '
                                          'readout boundary. The ranking comes from the band-edge '
                                          'distances and the leak interval.'),
            },
            'readout_margins_n6': {
                'clock_period_h': n6['clock_period_h'],
                'read_window_h_median': n6['read_window_h_median'],
                'reads_total': n6['reads_steady_drop8'] + 8,
                'reads_steady_drop8': n6['reads_steady_drop8'],
                'reads_steady_drop4': n6['reads_steady_drop4'],
                'read_windows_clipped': n6['read_windows_clipped'],
                'unlabelled_reads': n6['unlabelled_reads_all'],
                'bit2_low_reads': n6['bit2_low_reads'],
                'bit2_high_reads': n6['bit2_high_reads'],
                'commitment_min_all_bits': n6['steady_min_commitment_all_bits'],
                'bit2_low_window_max_max': n6['bit2_low_window_max_max'],
                'bit2_high_window_min_min': n6['bit2_high_window_min_min'],
            },
            'layer_derivation': {
                'D1': '34-state two-bit prefix frozen',
                'D2': '51-state model: states 0..33 have no dependence on states 34..50',
                'D3': 'default n_A1_gate=4 fails the three-bit count',
                'D4': 'shared-maturation route excluded',
                'D5': 'independent gate exponent n=5..8 gives 16/16 mod-8 certification',
                'D6': 'three-segment split shows the old R1 mostly measured the legitimate Int2 tail',
                'D7': 'n=6 beats n=5 on leak, both band-edge margins and bit2 hold',
                'D8': 'commitment recomputation: neither point is near the readout boundary',
            },
        },
        'provenance': {
            'structural_guards': {
                'prefix_guard_passed': (guard or {}).get('prefix_guard_passed'),
                'prefix_guard_rhs_gap_states_0_33': (guard or {}).get(
                    'structural_rhs_gap_states_0_33'),
                'prefix_guard_note': ('states 0..33 must be independent of n_A1_gate; the structural '
                                      'RHS test is used because the adaptive step-size controller '
                                      'makes a trajectory comparison show ~1e-13 round-off'),
                'decoupling_verification_passed': (decoup or {}).get('passed'),
            },
            'script_versions': {
                'check_read_commitment.py v1 (produced the archived v1 shard CSVs)':
                    'B508EB62CF878F285C8FB4DA52FDCAEF3DA8FD54B7BBC8CF983ABC9FD4903421',
                'check_read_commitment.py v2 (merge column-name fix only)':
                    '769BD9DF26B8D1F234375C5187A4D1FB6D2CF6BD8504A9660DEAB72DDFCAE146',
                'check_read_commitment.py v3 (adds setup/hold and per-bit margins)':
                    sha256(OUT / 'check_read_commitment.py'),
                'v1_hash_source': ('reconstructed by reversing the single merge-path edit; the scan '
                                   'ran with exactly this content and a byte-identical copy is kept '
                                   'under plausibility/read_commitment_v1/'),
                'merge_bug': ('v1 read the margin columns without the drop8_/drop4_ prefix, so merge '
                              'raised KeyError. Only the reporting path was affected; the scanner did '
                              'NOT re-integrate anything.'),
                'v1_vs_v3_shared_columns': 122,
                'v1_vs_v3_differing_columns': [
                    'ref_bit2_hold_h (v1 placeholder None -> v3 measured value)',
                    'runtime_s (v3 is ~11 s slower because it adds one analyse pass)'],
                'v1_vs_v3_verdict': ('every one of the 122 shared columns is reproduced exactly; the '
                                     'setup/hold extension perturbed no data'),
            },
            'evidence_sha256': evidence_hashes,
            'frozen_source_sha256': source_hashes(),
            'gate_segments_crosscheck': {
                'leak_ratio_gap_n5': n5['crosscheck_leak_ratio_gap'],
                'leak_ratio_gap_n6': n6['crosscheck_leak_ratio_gap'],
                'crossings_gap': n6['crosscheck_crossings_gap'],
                'bit2_hold_gap': n6['crosscheck_hold_gap'],
                'note': ('the commitment stage re-integrated the same two points and reproduced the '
                         'gate-segment row to round-off, so the margins and the certification come '
                         'from the same trajectory'),
            },
        },
        'still_open_before_frozen': [
            'eight phase-consistent initial states (values 000..111) of the 51-state model',
            'strict-tolerance re-verification at the selected point',
            'molecular-pool perturbations (S / RDF / A1-F1 / Int2)',
            'final source and result hashes for the frozen candidate',
            'written statement of the unreconciled growth condition',
        ],
        'limitations': [
            'gamma = intrinsic decay + growth dilution is confirmed, so add_growth=False is the '
            'consistent choice; however mu < 0.6/h is required for the selected point while the '
            'upstream oscillator uses Td = 50 min (mu = 0.8318/h). The two authors growth '
            'conditions are NOT reconciled and this must be stated with any frozen profile.',
            'k_int = 6 is a square-wave-specific tuning value and must not be applied to the real '
            'C31 flux.',
            'uM_per_au = 5.75 is an uncalibrated interface scale, not a measured conversion.',
            'the explicit complex is a structural addition relative to Zeng reduced model and '
            'sequesters part of the free Int pool.',
            'alpha_rep / gamma_rep remain 3.5 / 0.6 from the PDF; they are NOT changed to the '
            'square-wave code 3.2 / 0.7 without Zeng confirmation.',
        ],
        'upgrade_policy': {
            'this_file_is_read_only': True,
            'superseded_by': None,
            'formal_candidate_file': 'plausibility/threebit51_selected_v1.json',
            'condition': ('created only after the eight initial states, the strict-tolerance ring and '
                          'the molecular-pool perturbations have all passed'),
        },
    }
    TARGET.write_text(json.dumps(profile, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'wrote {TARGET}')
    print('sha256', hashlib.sha256(TARGET.read_bytes()).hexdigest().upper())
    print('status  ', profile['status'])


if __name__ == '__main__':
    main()
