r"""Archive the pre-revision pre-audit report and both hash manifests.  (H3/H7)

Why
---
The pre-audit report is about to be revised (the "20 molecules execute the write"
framing, the `1/sqrt(N)` CV claim, the `conc_scale` naming and the "a stochastic study
is REQUIRED" conclusion all change).  Under H3/H7 an invalidated or superseded artefact
is ARCHIVED FIRST, never overwritten in place, so the previous wording stays quotable and
the previous manifests stay verifiable.

What is archived, verbatim, with the hash each file had at archive time:
  PREAUDIT_STOCHASTICITY_REPORT.md   -> ARCHIVED_PREAUDIT_REPORT_v1_20260924.md
  plausibility/SHA256SUMS.json       -> ARCHIVED_plausibility_SHA256SUMS.v1.json
  plausibility/SHA256SUMS.txt        -> ARCHIVED_plausibility_SHA256SUMS.v1.txt
  SHA256SUMS_authority.json (root)   -> ARCHIVED_authority_SHA256SUMS.v1.json
  SHA256SUMS_authority.txt  (root)   -> ARCHIVED_authority_SHA256SUMS.v1.txt

The archive file names deliberately do NOT start with `SHA256SUMS`: `write_manifest()`
in plausibility_common.py skips such names, and an archived manifest that the manifest
cannot see would be the very kind of untraceable artefact this round is fixing.

`ARCHIVED_HASHES.json` records, for each pair, the hash BEFORE and AFTER copying, so a
reader can confirm the copy is byte-identical to what was on disk at archive time.

Run:  & 'D:\aconade\python.exe' .\archive_preaudit_report.py
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WIKI = HERE.parents[1]
ARCHIVE = HERE / 'preaudit_report_archive'

PAIRS = (
    ('PREAUDIT_STOCHASTICITY_REPORT.md',
     'ARCHIVED_PREAUDIT_REPORT_v1_20260924.md'),
    ('SHA256SUMS.json', 'ARCHIVED_plausibility_SHA256SUMS.v1.json'),
    ('SHA256SUMS.txt', 'ARCHIVED_plausibility_SHA256SUMS.v1.txt'),
    (WIKI / 'SHA256SUMS_authority.json', 'ARCHIVED_authority_SHA256SUMS.v1.json'),
    (WIKI / 'SHA256SUMS_authority.txt', 'ARCHIVED_authority_SHA256SUMS.v1.txt'),
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for blk in iter(lambda: fh.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest().upper()


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:                                                   # noqa: BLE001
        pass
    ARCHIVE.mkdir(parents=True, exist_ok=True)
    records = []
    for src, dst_name in PAIRS:
        src_path = src if isinstance(src, Path) else (HERE / src)
        dst = ARCHIVE / dst_name
        if not src_path.is_file():
            raise SystemExit(f'FATAL: {src_path} does not exist; nothing to archive')
        before = sha256(src_path)
        if dst.exists() and sha256(dst) != before:
            raise SystemExit(
                f'FATAL: {dst} already exists with DIFFERENT content. Refusing to '
                f'overwrite an archive (that is exactly the H3/H7 failure).')
        shutil.copy2(src_path, dst)
        after = sha256(dst)
        if after != before:
            raise SystemExit(f'FATAL: copy of {src_path} is not byte-identical')
        try:
            rel_src = str(src_path.relative_to(WIKI))
        except ValueError:
            rel_src = str(src_path)
        records.append(dict(source=rel_src, archived=str(dst.relative_to(WIKI)),
                            sha256_before=before, sha256_after_copy=after,
                            size_bytes=dst.stat().st_size))
        print(f'archived {rel_src}\n      -> {dst.name}  {before}')

    doc = dict(
        kind='preaudit_report_archive',
        created_for=('the 2026-09-26 pre-audit wrap-up round: the report and the two '
                     'hash manifests as they stood before the revision'),
        reason=('H3/H7 - artefacts are archived before being superseded, never '
                'overwritten in place. The revised report changes the wording of the '
                'small-number claim, the CV statement, the conc_scale naming and the '
                '"stochastic study REQUIRED" conclusion, so the previous wording must '
                'stay quotable and its manifests verifiable.'),
        naming_note=('archive file names deliberately do not start with "SHA256SUMS", '
                     'because write_manifest() skips such names and an archived '
                     'manifest invisible to the manifest would defeat the purpose'),
        files=records,
    )
    (ARCHIVE / 'ARCHIVED_HASHES.json').write_text(
        json.dumps(doc, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'wrote {ARCHIVE / "ARCHIVED_HASHES.json"}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
