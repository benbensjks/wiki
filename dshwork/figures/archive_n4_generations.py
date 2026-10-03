r"""Durable archive of the two superseded n=4 matched-trajectory data files. (H3/H7)

Both files are superseded; neither is deleted.  The point is that a reader must be able
to verify the Figure 4-3 incident from the bytes themselves rather than from a claim
about them:

  generation 1  data/n4_matched_600h.csv as first written.
                Columns were `S0 = y[44]`, `S1 = y[45]`, `S2 = y[46]`, so the file held
                bit2's **S2, A1 and F1** under the names S0, S1, S2.  Figure 4-3 Panel A
                therefore plotted F1 in its n=4 panel and the real S2 in its n=6 panel.

  generation 2  after the index bug was fixed (columns renamed to the model's own
                `b0_S/b1_S/b2_S/A1/F1`, indices resolved by `state_names`), but before
                the bit-2 pools and fluxes were added.  Its b2_S column is bit-identical
                to generation 1's S0 column, which is how "the fix was a pure relabel"
                is established.

Run:  & 'D:\aconade\python.exe' .\archive_n4_generations.py
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / 'archive'

# (source on disk now, name in the archive, what generation it is)
SOURCES = (
    (Path(r'C:\Users\18633\AppData\Local\Temp\dsh-bd092P\n4_matched_600h_PREBUG.csv'),
     'n4_matched_600h_gen1_mislabelled_S0S1S2.csv',
     'generation 1: columns S0/S1/S2 actually held b2_S/A1/F1'),
    (Path(r'C:\Users\18633\AppData\Local\Temp\dsh-bd092P\n4_pre_extend.csv'),
     'n4_matched_600h_gen2_relabelled_9col.csv',
     'generation 2: correct column names, before b2_I/b2_R/b2_C/J_fwd2/J_rev2 added'),
)
EXPECTED = {
    'n4_matched_600h_gen1_mislabelled_S0S1S2.csv':
        'FBE88CF7FFE9568E47BDF7C4F2D8D31894C64664405B00AE79D45665043F1692',
    'n4_matched_600h_gen2_relabelled_9col.csv':
        'F2DB0AA432325491F385C4EEFC2A411F74398418EA4D89ED8A8623DBE8CEBC6F',
}
COLUMNS = {
    'n4_matched_600h_gen1_mislabelled_S0S1S2.csv':
        ['time_h', 'S0', 'S1', 'S2', 'g1', 'clock_gate', 'C31_flux'],
    'n4_matched_600h_gen2_relabelled_9col.csv':
        ['time_h', 'g1', 'clock_gate', 'C31_flux', 'b0_S', 'b1_S', 'b2_S', 'A1', 'F1'],
}


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
    for src, name, generation in SOURCES:
        if not src.is_file():
            raise SystemExit(f'FATAL: {src} is gone; the archive cannot be rebuilt from '
                             f'what is on disk. Do not invent it.')
        before = sha256(src)
        want = EXPECTED[name]
        if before != want:
            raise SystemExit(f'FATAL: {src}\n  hash {before}\n  expected {want}\n'
                             f'The file on disk is not the generation this archive '
                             f'claims to hold.')
        dst = ARCHIVE / name
        if dst.exists() and sha256(dst) != before:
            raise SystemExit(f'FATAL: {dst} exists with different content; refusing to '
                             f'overwrite an archive')
        shutil.copy2(src, dst)
        after = sha256(dst)
        if after != before:
            raise SystemExit(f'FATAL: copy of {src} is not byte-identical')
        records.append(dict(file=name, generation=generation,
                            columns=COLUMNS[name],
                            sha256=before, size_bytes=dst.stat().st_size,
                            superseded_by='data/n4_matched_600h.csv'))
        print(f'archived {name}\n      {generation}\n      {before}')
    doc = dict(
        kind='n4_matched_trajectory_archive',
        reason=('H3/H7 - superseded artefacts are archived, never overwritten in place. '
                'The Figure 4-3 column-index incident must be verifiable from the bytes '
                'and not only from a written claim.'),
        how_to_read=('gen1 is the file that produced the wrong figure; gen2 is the '
                     'relabelled file. gen2[b2_S] == gen1[S0] and gen2[A1] == gen1[S1] '
                     'and gen2[F1] == gen1[S2], entry-wise identical, which is what '
                     'establishes that the fix changed labels only and no numbers.'),
        files=records,
    )
    (ARCHIVE / 'ARCHIVED_HASHES.json').write_text(
        json.dumps(doc, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'wrote {ARCHIVE / "ARCHIVED_HASHES.json"}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
