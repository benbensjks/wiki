"""Root-level SHA256 manifest for the project's AUTHORITY documents.

Why this exists
---------------
`plausibility/SHA256SUMS.json` covers `plausibility/**` plus a few named source
files, and nothing outside `final_reconstruction/`.  The project's most
authoritative document - `CURRENT_BASELINE.md`, which carries section 11
("判据与验证纪律", the delivery-facing form of HARD_CONSTRAINTS H1-H7) - sits at
the wiki ROOT and was therefore covered by NO manifest at all, while its peer
`N6_FREEZE_JUSTIFICATION.md` was covered.  Two documents of equal standing must
not differ in traceability, so this script closes that gap.

This is NOT a whole-tree hash.  The wiki also holds unrelated material (other
groups' packages, delivery drafts, working notes); hashing all of it would make
the manifest churn on files that have nothing to do with the model.  The scope is
an explicit NAMED LIST, and a listed path that is missing is a hard error - the
manifest never silently shrinks.

Two subcommands
---------------
    write   regenerate the root manifest (and its .txt mirror)
    check   READ-ONLY verification, exit 1 on any problem; checks
              (a) every listed document still has the recorded hash,
              (b) the plausibility manifest is not stale, i.e. every one of its
                  entries still matches the file on disk,
              (c) the plausibility manifest's .json and .txt agree with each other,
              (d) the authority manifest does not list itself.

The staleness check is deliberately strict: any edit anywhere inside
`plausibility/` invalidates the root manifest and `check` fails loudly until
`write` is run again.  A stale manifest that passes silently is worse than none.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

# wiki root: <wiki>/final_reconstruction/plausibility/make_authority_manifest.py
WIKI = Path(__file__).resolve().parents[2]
FR = WIKI / 'final_reconstruction'
PLAUS = FR / 'plausibility'

OUT_JSON = WIKI / 'SHA256SUMS_authority.json'
OUT_TXT = WIKI / 'SHA256SUMS_authority.txt'

# Named authority documents, relative to the wiki root.  Order is for humans.
AUTHORITY = (
    # --- the authoritative baseline and its discipline section ---------------
    'CURRENT_BASELINE.md',
    # --- frozen model sources (also pinned inside the plausibility manifest) --
    'final_reconstruction/model.py',
    'final_reconstruction/model_twobit34.py',
    'final_reconstruction/model_threebit51.py',
    'final_reconstruction/verify_threebit51.py',
    'final_reconstruction/verify_twobit_causal.py',
    # --- project-level readmes ----------------------------------------------
    'final_reconstruction/README.md',
    'final_reconstruction/plausibility/README.md',
    # --- the two discipline documents ---------------------------------------
    'final_reconstruction/plausibility/HARD_CONSTRAINTS.md',
    'final_reconstruction/plausibility/N6_FREEZE_JUSTIFICATION.md',
    # --- frozen parameter artefacts -----------------------------------------
    'final_reconstruction/plausibility/threebit51_selected_v1.json',
    'final_reconstruction/plausibility/threebit51_provisional_n6.json',
    'final_reconstruction/plausibility/threebit51_leak_correction.json',
    # --- verdicts a reader is most likely to quote ---------------------------
    'final_reconstruction/plausibility/decoupling_verification.json',
    'final_reconstruction/plausibility/carry_pairing_grid_verdict.json',
    'final_reconstruction/plausibility/nA1_vs_maturation_verdict.json',
    'final_reconstruction/plausibility/translation_unobservable.json',
    'final_reconstruction/plausibility/local_identifiability/ANALYSIS_STATUS.json',
    # the pre-audit report and the per-item table it is now read through; the table is
    # the only citable per-component statement, so it must not sit outside the manifest
    'final_reconstruction/plausibility/PREAUDIT_STOCHASTICITY_REPORT.md',
    'final_reconstruction/plausibility/preaudit_acceptance_breakdown.md',
    # --- retraction / invalidation records (evidence of process failures) ----
    'final_reconstruction/plausibility/INVALIDATED_ARTIFACTS.json',
    'final_reconstruction/plausibility/local_identifiability/RETRACTION_n4_run.md',
    'final_reconstruction/plausibility/invalidated_n4/ARCHIVED_HASHES.json',
    'final_reconstruction/plausibility/preaudit_report_archive/ARCHIVED_HASHES.json',
    # --- the manifest this one pins -----------------------------------------
    'final_reconstruction/plausibility/SHA256SUMS.json',
    'final_reconstruction/plausibility/SHA256SUMS.txt',
    # --- wiki-facing notes --------------------------------------------------
    'final_reconstruction/threebit51_results/wiki_n6_20260924/README_图注与参数.md',
    'final_reconstruction/twobit34_results/selected_profile.md',
    # --- scan figures (provenance of every plotted number) -------------------
    'final_reconstruction/scan_figures/README_图注与参数.md',
    'final_reconstruction/scan_figures/figure_data_manifest.json',
)

SELF = ('SHA256SUMS_authority.json', 'SHA256SUMS_authority.txt')


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest().upper()


def build_entries() -> dict:
    """Hash every listed document.  A missing listed document is a hard error."""
    entries, missing = {}, []
    for rel in AUTHORITY:
        p = WIKI / rel
        if not p.is_file():
            missing.append(rel)
            continue
        entries[rel] = sha256_file(p)
    if missing:
        raise SystemExit(
            'FATAL: %d listed authority document(s) are missing. The manifest is\n'
            'never allowed to shrink silently; either restore the file or remove it\n'
            'from AUTHORITY in this script with a recorded reason.\n  %s'
            % (len(missing), '\n  '.join(missing)))
    for name in SELF:
        if name in entries:
            raise SystemExit(f'FATAL: {name} must not list itself')
    return entries


def check_plausibility_manifest() -> tuple:
    """Recompute every entry of plausibility/SHA256SUMS.json from disk.

    Returns (problems, n_entries, n_checked).
    """
    problems = []
    jpath = PLAUS / 'SHA256SUMS.json'
    tpath = PLAUS / 'SHA256SUMS.txt'
    if not jpath.is_file():
        return ['plausibility/SHA256SUMS.json is missing'], 0, 0
    manifest = json.loads(jpath.read_text(encoding='utf-8'))
    checked = 0
    for rel, recorded in manifest.items():
        p = FR / rel
        if not p.is_file():
            problems.append(f'listed in plausibility manifest but absent: {rel}')
            continue
        actual = sha256_file(p)
        checked += 1
        if actual != recorded:
            problems.append(f'STALE plausibility entry: {rel}\n'
                            f'    recorded {recorded}\n    actual   {actual}')
    # the .txt mirror must describe exactly the same content
    if tpath.is_file():
        mirror = {}
        for line in tpath.read_text(encoding='utf-8').splitlines():
            if not line.strip():
                continue
            digest, _, rel = line.partition('  ')
            mirror[rel] = digest
        only_txt = sorted(set(mirror) - set(manifest))
        only_json = sorted(set(manifest) - set(mirror))
        disagree = sorted(k for k in set(mirror) & set(manifest)
                          if mirror[k] != manifest[k])
        if only_txt:
            problems.append(f'SHA256SUMS.txt lists {len(only_txt)} path(s) absent '
                            f'from the .json: {only_txt[:5]}')
        if only_json:
            problems.append(f'SHA256SUMS.json lists {len(only_json)} path(s) absent '
                            f'from the .txt: {only_json[:5]}')
        if disagree:
            problems.append(f'{len(disagree)} path(s) disagree between .json and '
                            f'.txt: {disagree[:5]}')
    else:
        problems.append('plausibility/SHA256SUMS.txt is missing')
    return problems, len(manifest), checked


def write_cmd(_args) -> int:
    entries = build_entries()
    problems, n_entries, n_checked = check_plausibility_manifest()
    doc = {
        'kind': 'authority_manifest',
        'scope': ('named authority documents at the wiki root and inside '
                  'final_reconstruction; NOT a whole-tree hash'),
        'wiki_root': str(WIKI),
        'entries': entries,
        'entry_count': len(entries),
        'plausibility_manifest': {
            'path': 'final_reconstruction/plausibility/SHA256SUMS.json',
            'sha256': entries.get('final_reconstruction/plausibility/SHA256SUMS.json'),
            'entries': n_entries,
            'entries_reverified': n_checked,
            'stale_entries': problems,
            'stale': bool(problems),
        },
        'notes': [
            'excludes itself: ' + ', '.join(SELF),
            'run `check` for a read-only verification of every entry',
            'any edit inside plausibility/ makes this manifest stale by design',
        ],
    }
    OUT_JSON.write_text(json.dumps(doc, ensure_ascii=False, indent=2, sort_keys=True),
                        encoding='utf-8')
    OUT_TXT.write_text(''.join(f'{v}  {k}\n' for k, v in sorted(entries.items())),
                       encoding='utf-8')
    print(f'wrote {OUT_JSON}')
    print(f'wrote {OUT_TXT}')
    print(f'entries: {len(entries)}')
    print(f'plausibility manifest: {n_entries} entries, {n_checked} reverified, '
          f'{len(problems)} problem(s)')
    for p in problems:
        print('  PROBLEM:', p.replace('\n', '\n  '))
    return 1 if problems else 0


def check_cmd(_args) -> int:
    problems = []
    if not OUT_JSON.is_file():
        print(f'FATAL: {OUT_JSON} does not exist; run `write` first')
        return 1
    doc = json.loads(OUT_JSON.read_text(encoding='utf-8'))
    recorded = doc.get('entries', {})
    # (a) every listed document still matches
    for rel, digest in sorted(recorded.items()):
        p = WIKI / rel
        if not p.is_file():
            problems.append(f'authority document is missing: {rel}')
            continue
        actual = sha256_file(p)
        if actual != digest:
            problems.append(f'AUTHORITY DOC CHANGED: {rel}\n'
                            f'    recorded {digest}\n    actual   {actual}')
    # (b) the listed set must equal the script's current named list
    listed = set(recorded)
    wanted = set(AUTHORITY)
    for rel in sorted(wanted - listed):
        problems.append(f'named in this script but absent from the manifest: {rel}')
    for rel in sorted(listed - wanted):
        problems.append(f'in the manifest but no longer named in this script: {rel}')
    # (c) plausibility manifest not stale
    sub, n_entries, n_checked = check_plausibility_manifest()
    problems.extend(sub)
    # (d) the .txt mirror
    if OUT_TXT.is_file():
        mirror = {}
        for line in OUT_TXT.read_text(encoding='utf-8').splitlines():
            if not line.strip():
                continue
            digest, _, rel = line.partition('  ')
            mirror[rel] = digest
        if mirror != recorded:
            problems.append('SHA256SUMS_authority.txt does not mirror the .json')
    else:
        problems.append(f'{OUT_TXT} is missing')
    # (e) self-exclusion
    for name in SELF:
        if name in recorded:
            problems.append(f'the manifest lists itself: {name}')

    print(f'authority entries checked : {len(recorded)}')
    print(f'plausibility entries      : {n_entries} recorded, {n_checked} reverified')
    if problems:
        print(f'FAIL: {len(problems)} problem(s)')
        for p in problems:
            print('  -', p.replace('\n', '\n    '))
        return 1
    print('OK: every authority document and every plausibility entry matches on disk')
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('write', help='regenerate the root-level authority manifest')
    sub.add_parser('check', help='read-only verification; exit 1 on any problem')
    args = ap.parse_args()
    return (write_cmd if args.cmd == 'write' else check_cmd)(args)


if __name__ == '__main__':
    sys.exit(main())
