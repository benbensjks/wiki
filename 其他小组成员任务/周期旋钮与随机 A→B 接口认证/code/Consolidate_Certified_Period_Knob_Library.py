"""Consolidate the deterministic A+B-certified period-knob library.

This post-processing script does not run module A or B.  It combines the
previous 11-13 generation certification with the alternative 10-generation
candidate that passed both PB and LR initial states.  Failed and pending
candidates are excluded from the certified library.

Run from the project root:

    python ./code/Consolidate_Certified_Period_Knob_Library.py
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


VERSION = "Certified_Period_Knob_Library_20260819"
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_ROOT = PROJECT_ROOT / "outputs"
FIGURE_ROOT = PROJECT_ROOT / "figures"
TD_MIN = 50.0
EXPECTED_TARGETS = (10.0, 11.0, 12.0, 13.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-certification-dir", type=Path, default=None)
    parser.add_argument("--fast-certification-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def newest_directory(pattern: str) -> Path:
    matches = sorted(
        (path for path in OUTPUT_ROOT.glob(pattern) if path.is_dir()),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not matches:
        raise FileNotFoundError(f"No directory matched: {OUTPUT_ROOT / pattern}")
    return matches[0]


def resolve_directory(explicit: Path | None, pattern: str) -> Path:
    if explicit is None:
        return newest_directory(pattern)
    path = explicit.expanduser().resolve()
    if not path.is_dir():
        raise NotADirectoryError(path)
    return path


def normalise_source(root: Path, source_name: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    merge = pd.read_csv(root / "period_knob_B_validation_for_merge.csv")
    conditions = pd.read_csv(root / "AB_initial_state_condition_summary.csv")
    tags = pd.read_csv(root / "AB_tag_specific_certification.csv")
    merge["source_certification"] = source_name
    conditions["source_certification"] = source_name
    tags["source_certification"] = source_name
    if "candidate_id" not in merge:
        merge["candidate_id"] = merge.apply(
            lambda row: (
                f"T{row['target_period_generations']:g}__RBS{row['design_RBS_scale']:g}__"
                f"mRNA{row['design_oscillator_mRNA_half_life_min']:g}min"
            ).replace(".", "p"),
            axis=1,
        )
    return merge, conditions, tags


def matching_rows(frame: pd.DataFrame, design: pd.Series) -> pd.DataFrame:
    mask = (
        np.isclose(frame["target_period_generations"], float(design["target_period_generations"]))
        & np.isclose(frame["design_RBS_scale"], float(design["design_RBS_scale"]))
        & np.isclose(
            frame["design_oscillator_mRNA_half_life_min"],
            float(design["design_oscillator_mRNA_half_life_min"]),
        )
    )
    return frame.loc[mask]


def select_certified_designs(merge: pd.DataFrame) -> pd.DataFrame:
    certified = merge.loc[
        merge["B_pass"].astype(bool)
        & merge["source_period_tolerance_pass"].astype(bool)
        & merge["complete_deterministic_knob_certification"].astype(bool)
        & merge["target_period_generations"].isin(EXPECTED_TARGETS)
    ].copy()
    if "realised_A_period_generations" in certified:
        certified["selection_error"] = (
            certified["realised_A_period_generations"]
            - certified["target_period_generations"]
        ).abs()
    else:
        certified["selection_error"] = np.nan
    certified = (
        certified.sort_values(
            ["target_period_generations", "selection_error", "best_worst_one_crossing_complete_fraction"],
            ascending=[True, True, False], na_position="last",
        )
        .groupby("target_period_generations", as_index=False, group_keys=False)
        .head(1)
        .sort_values("target_period_generations")
        .reset_index(drop=True)
    )
    found = tuple(certified["target_period_generations"].astype(float))
    if found != EXPECTED_TARGETS:
        raise RuntimeError(f"Expected certified targets {EXPECTED_TARGETS}, found {found}.")
    return certified


def build_library(
    certified: pd.DataFrame,
    condition_tables: dict[str, pd.DataFrame],
    tag_tables: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    rows = []
    for _, design in certified.iterrows():
        source = str(design["source_certification"])
        conditions = matching_rows(condition_tables[source], design)
        tags = matching_rows(tag_tables[source], design)
        if conditions.empty or tags.empty:
            raise RuntimeError(f"Missing detailed certification rows for {design['candidate_id']}.")
        passing_tags = tags.loc[tags["both_initial_states_B_pass"].astype(bool)]
        if passing_tags.empty:
            raise RuntimeError(f"No passing tag rows for {design['candidate_id']}.")

        realised = float(conditions["realised_A_period_generations"].median())
        rows.append(
            {
                "knob_setting_id": f"Certified_T{float(design['target_period_generations']):g}",
                "target_period_generations": float(design["target_period_generations"]),
                "realised_period_generations": realised,
                "realised_period_h_at_Td50": realised * TD_MIN / 60.0,
                "period_error_generations": abs(
                    realised - float(design["target_period_generations"])
                ),
                "synchronous_oscillator_RBS_scale": float(design["design_RBS_scale"]),
                "synchronous_oscillator_mRNA_half_life_min": float(
                    design["design_oscillator_mRNA_half_life_min"]
                ),
                "shared_PLtetO1": True,
                "C31_RBS_scale": 0.45,
                "C31_mRNA_total_half_life_min": 2.0,
                "sponge_plasmid_copies": 40.0,
                "sponge_accessible_fraction": 1.0,
                "peak_resource_load_fraction": 0.0,
                "C31_flux_peak_uM_h": float(conditions["A_flux_peak_uM_h_median"].median()),
                "C31_flux_FWHM_h": float(conditions["A_flux_FWHM_h_median"].median()),
                "C31_flux_dose_uM_per_cycle": float(
                    conditions["A_flux_dose_uM_cycle_median"].median()
                ),
                "B_passing_Int_tag_rates_h_inv": ";".join(
                    f"{value:g}" for value in sorted(passing_tags["Int_tag_rate_h_inv"].unique())
                ),
                "B_worst_one_crossing_complete_fraction": float(
                    passing_tags["worst_one_crossing_complete_fraction"].min()
                ),
                "B_worst_alternation_fidelity": float(
                    passing_tags["worst_alternation_fidelity"].min()
                ),
                "max_DNA_conservation_residual": float(
                    passing_tags["max_DNA_conservation_residual"].max()
                ),
                "deterministic_A_pass": True,
                "deterministic_B_pass": True,
                "stochastic_validation_status": "pending",
                "source_certification": source,
                "source_candidate_id": str(design["candidate_id"]),
            }
        )
    return pd.DataFrame(rows)


def make_figure(library: pd.DataFrame, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    x = library["target_period_generations"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    axes[0].plot(x, library["realised_period_generations"], "o-", label="realised")
    axes[0].plot(x, x, "k--", label="target")
    axes[0].set(xlabel="target period (gen)", ylabel="realised period (gen)")
    axes[0].legend()

    axes[1].plot(x, library["C31_flux_peak_uM_h"], "o-", label="peak (uM/h)")
    axes[1].plot(x, library["C31_flux_FWHM_h"], "s-", label="FWHM (h)")
    axes[1].set_xlabel("target period (gen)")
    axes[1].legend()

    axes[2].plot(x, library["C31_flux_dose_uM_per_cycle"], "o-", color="tab:purple")
    axes[2].set(xlabel="target period (gen)", ylabel="C31 dose (uM/cycle)")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.suptitle("Deterministic A+B-certified period-knob library")
    fig.tight_layout()
    fig.savefig(figure_dir / "Certified_period_knob_library.png", dpi=190, bbox_inches="tight")
    plt.close(fig)


def write_markdown(library: pd.DataFrame, output_path: Path) -> None:
    lines = [
        "# Deterministic A+B-certified period-knob library",
        "",
        "Only settings that passed module A, both PB/LR initial states in module B,",
        "and the period-tolerance criterion are included.",
        "",
        "| Target (gen) | Realised (gen) | Hours at Td=50 min | Oscillator RBS | mRNA half-life (min) | Passing Int tags (h^-1) |",
        "|---:|---:|---:|---:|---:|---|",
    ]
    for row in library.itertuples(index=False):
        lines.append(
            f"| {row.target_period_generations:.0f} | {row.realised_period_generations:.3f} | "
            f"{row.realised_period_h_at_Td50:.2f} | {row.synchronous_oscillator_RBS_scale:.2f} | "
            f"{row.synchronous_oscillator_mRNA_half_life_min:.1f} | "
            f"{row.B_passing_Int_tag_rates_h_inv} |"
        )
    lines.extend(
        [
            "",
            "The current deterministic system-certified range is approximately 10-13 generations.",
            "The 8-9 generation A-module settings remain outside the present B recovery window.",
            "Stochastic CV and phase-drift certification remain pending.",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    legacy_root = resolve_directory(
        args.legacy_certification_dir, "Period_Knob_AB_Certification_full_*"
    )
    fast_root = resolve_directory(
        args.fast_certification_dir, "Fast_Period_AB_Candidates_full_*"
    )
    legacy_merge, legacy_conditions, legacy_tags = normalise_source(legacy_root, "legacy_11_to_13")
    fast_merge, fast_conditions, fast_tags = normalise_source(fast_root, "fast_alternative_candidates")
    combined_merge = pd.concat([legacy_merge, fast_merge], ignore_index=True, sort=False)
    certified = select_certified_designs(combined_merge)
    library = build_library(
        certified,
        {
            "legacy_11_to_13": legacy_conditions,
            "fast_alternative_candidates": fast_conditions,
        },
        {
            "legacy_11_to_13": legacy_tags,
            "fast_alternative_candidates": fast_tags,
        },
    )

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = args.output_dir or OUTPUT_ROOT / f"Certified_Period_Knob_Library_{stamp}"
    figure_dir = FIGURE_ROOT / output_dir.name
    output_dir.mkdir(parents=True, exist_ok=False)
    figure_dir.mkdir(parents=True, exist_ok=False)
    combined_merge.to_csv(output_dir / "all_available_AB_certifications.csv", index=False)
    library.to_csv(output_dir / "certified_period_knob_library.csv", index=False)
    library.to_csv(output_dir / "period_knob_library_for_stochastic_validation.csv", index=False)
    make_figure(library, figure_dir)
    write_markdown(library, output_dir / "CERTIFIED_PERIOD_KNOB_LIBRARY.md")
    metadata = {
        "version": VERSION,
        "Td_min": TD_MIN,
        "legacy_certification_dir": str(legacy_root),
        "fast_certification_dir": str(fast_root),
        "included_targets_generations": EXPECTED_TARGETS,
        "exclusion_rule": "Exclude failed, pending, or period-out-of-tolerance settings.",
        "stochastic_validation": "pending",
    }
    (output_dir / "scope_and_provenance.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(f"Completed: {VERSION}")
    print(f"Raw tables: {output_dir}")
    print(f"Figures: {figure_dir}")
    print("\nCertified deterministic library:")
    for row in library.itertuples(index=False):
        print(
            f"target={row.target_period_generations:.0f} gen | "
            f"realised={row.realised_period_generations:.3f} gen/"
            f"{row.realised_period_h_at_Td50:.2f} h | "
            f"RBS={row.synchronous_oscillator_RBS_scale:.2f} | "
            f"mRNA t1/2={row.synchronous_oscillator_mRNA_half_life_min:.1f} min | "
            f"B tags={row.B_passing_Int_tag_rates_h_inv} | SSA=pending"
        )


if __name__ == "__main__":
    main()
