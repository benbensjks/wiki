"""
============================================================
Interaction2.py

RDF parameter interaction contribution.

This mirrors GSA_Oscillator-C31/Interaction.py:
interaction = ST - S1 for one selected output.
============================================================
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from parameters2 import MODE, OUTPUT, SOURCE_MODE


result_dir = os.path.join(os.path.dirname(__file__), "results", MODE, SOURCE_MODE, OUTPUT)
path = os.path.join(result_dir, f"sobol_{MODE}_{SOURCE_MODE}_{OUTPUT}.npz")

if not os.path.exists(path):
    raise RuntimeError(
        f"Missing Sobol result file: {path}\n"
        f"Run gsa2.py first with OUTPUT = {OUTPUT!r}."
    )

data = np.load(path, allow_pickle=True)

params = [str(x) for x in data["names"]]
S1 = data["S1"]
ST = data["ST"]
interaction = ST - S1


print("Interaction")

for p, i in zip(params, interaction):
    print(f"{p}: {i:.4f}")


plt.figure(figsize=(8, 4))

sns.barplot(
    x=params,
    y=interaction,
)

plt.ylabel("ST - S1")
plt.title(f"RDF Parameter Interaction Contribution\n{MODE} / {SOURCE_MODE} / {OUTPUT}")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

fig_path = os.path.join(result_dir, f"interaction_{MODE}_{SOURCE_MODE}_{OUTPUT}.png")
plt.savefig(fig_path, dpi=150)
print(f"Saved figure : {fig_path}")
plt.show()
