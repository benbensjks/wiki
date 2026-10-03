"""重建某个 run 目录的 `SHA256SUMS.json`（含事后分析产物），并留下重建记录。

起因（本轮自己踩的坑，记录在案）
--------------------------------
`scan_clock_axis.py` 的 `Run.save()` 在整轮跑完时对目录内**当时存在**的文件做哈希清单。
我在整轮还在跑的时候，为了早一点检查汇总脚本，就拿 `partial.json` 跑了一次
`analyze_clock_axis.py` —— 于是它在 run 目录里留下了 `summary.json`。
等整轮结束时，`SHA256SUMS.json` 就把这份**中间态**的 `summary.json` 的哈希写了进去。
随后正式分析又覆盖了 `summary.json`，清单与其内容不再对应：

    SHA256SUMS.json 里   summary.json = 5384AD4E…
    实际                  summary.json = B12626C9…

这正是哈希清单本来要防的那类漂移，只不过这次是**做清单的人自己造成的**。
教训：事后分析产物要么写进 `analysis/` 子目录，要么在分析之后重建清单。
本脚本做后者。

用法
----
    python refresh_manifest.py                      # 取 results/ 下最新的 clock_axis_*
    python refresh_manifest.py --dir <路径>
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common as C                                                  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default=None)
    args = ap.parse_args()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')

    d = (Path(args.dir) if args.dir
         else sorted((HERE / 'results').glob('clock_axis_*'))[-1])
    old_path = d / 'SHA256SUMS.json'
    old = json.loads(old_path.read_text(encoding='utf-8')) if old_path.exists() else {}

    new = {}
    for p in sorted(d.iterdir()):
        if p.is_file() and p.name != 'SHA256SUMS.json':
            new[p.name] = C.sha256(p)

    print(f'dir = {d.name}')
    for k in sorted(set(old) | set(new)):
        a, b = old.get(k), new.get(k)
        flag = '  ' if a == b else ('+ ' if a is None else '! ')
        print(f'  {flag}{k:22s} {str(a)[:12]:14s} -> {str(b)[:12]}')
    new['_note'] = ('rebuilt after the post-hoc analysis files were written; see '
                    'refresh_manifest.py for why the run-time manifest held a stale '
                    'summary.json')
    new['_rebuilt_at'] = datetime.now().isoformat(timespec='seconds')
    C.dump(old_path, new)
    print(f'\nwrote {old_path}')


if __name__ == '__main__':
    main()
