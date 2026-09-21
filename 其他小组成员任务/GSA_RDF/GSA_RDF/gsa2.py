"""
============================================================
gsa2.py

RDF system Global Sensitivity Analysis framework.

Same main steps as the oscillator GSA:
1. Saltelli sampling
2. Run model wrapper
3. Sobol analysis
4. Plot and save results
============================================================
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from SALib.analyze import sobol
from SALib.sample import saltelli

from model_wrapper2 import model_wrapper
from parameters2 import MODE, N, N_JOBS, OUTPUT, SOURCE_MODE, problem


def result_output_dir(output=OUTPUT):
    return os.path.join(
        os.path.dirname(__file__),
        "results",
        MODE,
        SOURCE_MODE,
        output,
    )


def result_stem(output=OUTPUT):
    return f"sobol_{MODE}_{SOURCE_MODE}_{output}"


def main():
    print("Generating Sobol samples...")

    param_values = saltelli.sample(
        problem,
        N,
        calc_second_order=False,
    )

    print(param_values[:5])
    print(f"Total simulations : {len(param_values)}")
    print("Running simulations...")

    Y = Parallel(
        n_jobs=N_JOBS,
        verbose=10,
    )(
        delayed(model_wrapper)(params)
        for params in param_values
    )

    Y = np.nan_to_num(np.array(Y))

    print("Y =", Y[:10])

    if np.var(Y) == 0:
        raise RuntimeError(
            "Output variance is zero.\n"
            "Please check parameter bounds, SOURCE_MODE, OUTPUT, or the model."
        )

    print("Calculating Sobol indices...")

    Si = sobol.analyze(
        problem,
        Y,
        calc_second_order=False,
        print_to_console=False,
    )

    S1 = Si["S1"]
    ST = Si["ST"]

    print("\n========== Sobol Results ==========\n")

    for name, s1, st in zip(problem["names"], S1, ST):
        print(
            f"{name:14s}"
            f"S1 = {s1:.4f}"
            f"    ST = {st:.4f}"
        )

    out_dir = result_output_dir()
    os.makedirs(out_dir, exist_ok=True)
    stem = result_stem()

    np.savez(
        os.path.join(out_dir, f"{stem}.npz"),
        Y=Y,
        S1=S1,
        ST=ST,
        names=np.array(problem["names"]),
        mode=MODE,
        source_mode=SOURCE_MODE,
        output=OUTPUT,
    )

    with open(os.path.join(out_dir, f"{stem}.txt"), "w") as f:
        f.write("========== Sobol Results ==========\n\n")
        f.write(f"MODE = {MODE}\n")
        f.write(f"SOURCE_MODE = {SOURCE_MODE}\n")
        f.write(f"OUTPUT = {OUTPUT}\n\n")
        for name, s1, st in zip(problem["names"], S1, ST):
            f.write(
                f"{name:14s}"
                f"S1 = {s1:.4f}"
                f"    ST = {st:.4f}\n"
            )

    x = np.arange(problem["num_vars"])
    width = 0.35

    plt.figure(figsize=(9, 5))

    plt.bar(
        x - width / 2,
        S1,
        width,
        label="First-order (S1)",
    )

    plt.bar(
        x + width / 2,
        ST,
        width,
        label="Total-effect (ST)",
    )

    plt.xticks(
        x,
        problem["names"],
        rotation=30,
        ha="right",
    )

    plt.ylabel("Sobol Index")
    plt.title(f"RDF GSA: {MODE} / {SOURCE_MODE} / {OUTPUT}")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    fig_path = os.path.join(out_dir, f"{stem}.png")
    plt.savefig(fig_path, dpi=150)
    print(f"Saved figure : {fig_path}")
    plt.show()


if __name__ == "__main__":
    main()
