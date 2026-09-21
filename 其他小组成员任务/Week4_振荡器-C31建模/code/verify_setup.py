"""Static delivery-package check. This script does not run any simulation."""

from __future__ import annotations

import ast
import importlib.util
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

PROGRAMS = (
    "Mechanistic_Repressilator_C31_ODE.py",
    "Oscillator_Engineering_Sensitivity_Analysis.py",
    "Flux_Driven_Translation_Burden_Model.py",
    "Fixed_PLtetO1_Interface_Design_Scan.py",
    "Self_Consistent_AB_System_Robustness_Validation.py",
    "Stochastic_Sponge_Period_and_Noise_Analysis.py",
)

PACKAGES = ("numpy", "pandas", "scipy", "matplotlib")
REQUIRED_FILES = (
    SCRIPT_DIR / "zhao_core.py",
    PROJECT_ROOT / "data" / "reference_inputs" /
    "unloaded_C31_translation_trajectory.csv",
)


def main() -> int:
    failures: list[str] = []

    print(f"Project root: {PROJECT_ROOT}")
    print("\nPython dependencies:")
    for package in PACKAGES:
        available = importlib.util.find_spec(package) is not None
        print(f"  {package:<12} {'OK' if available else 'MISSING'}")
        if not available:
            failures.append(f"missing package: {package}")

    print("\nRequired files:")
    for path in REQUIRED_FILES:
        available = path.is_file()
        print(f"  {path.relative_to(PROJECT_ROOT)}: "
              f"{'OK' if available else 'MISSING'}")
        if not available:
            failures.append(f"missing file: {path}")

    print("\nProgram syntax and portability:")
    for name in PROGRAMS:
        path = SCRIPT_DIR / name
        if not path.is_file():
            print(f"  {name}: MISSING")
            failures.append(f"missing program: {name}")
            continue
        text = path.read_text(encoding="utf-8-sig")
        try:
            ast.parse(text, filename=str(path))
            syntax = "OK"
        except SyntaxError as exc:
            syntax = f"ERROR ({exc})"
            failures.append(f"syntax error: {name}: {exc}")
        portable = (
            'if __name__ == "__main__":' in text
            and "Path.cwd()" not in text
            and "D:\\" not in text
            and "E:\\" not in text
        )
        if not portable:
            failures.append(f"portability check failed: {name}")
        print(f"  {name}: syntax={syntax}; portable={'OK' if portable else 'ERROR'}")

    if failures:
        print("\nSetup check failed:")
        for failure in failures:
            print(f"  - {failure}")
        return 1

    print("\nSetup check passed. No simulation was run.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

