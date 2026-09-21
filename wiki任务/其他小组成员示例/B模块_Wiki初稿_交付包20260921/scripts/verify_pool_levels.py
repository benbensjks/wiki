# -*- coding: utf-8 -*-
"""Provenance check for the pool levels shown in the mechanism figure:
  RDF pool ~4.5 uM (LR) vs ~0.002 uM (PB) at pulse arrival.

Method: for the primary non-zero-leak scheme (beta_B=2, beta_R=8, scale=0.30,
k_I=12, leak=0.008), find each clock-pulse onset (flux crosses 0.3 uM/h upward,
skipping the first two pulses) and sample the FREE RDF 0.2 h before the onset
(i.e. just before the pulse consumes it). Classify the state by the LR fraction
at that time. Report medians over all pulses of both initial-state runs.

Run:  python verify_pool_levels.py
"""
import csv
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(
    HERE, "..", "..",
    "B模块_非零泄漏方案_20260919", "deliverables", "reproduction", "results",
)


def analyse(fname):
    rows = list(csv.DictReader(open(os.path.join(DATA, fname))))
    t = [float(r["time_h"]) for r in rows]
    flux = [float(r["Int_flux_uM_h"]) for r in rows]
    rdf = [float(r["RDF_free_uM"]) for r in rows]
    lr = [float(r["LR_fraction"]) for r in rows]

    onsets, below = [], True
    for i in range(len(t)):
        if t[i] < 40:
            continue
        if below and flux[i] > 0.3:
            onsets.append(t[i])
            below = False
        elif not below and flux[i] < 0.1:
            below = True

    pb, lr_state = [], []
    for on in onsets:
        j = min(range(len(t)), key=lambda k: abs(t[k] - (on - 0.2)))
        (lr_state if lr[j] >= 0.5 else pb).append(rdf[j])
    return pb, lr_state, len(onsets)


for fname in ["point_B2.0_R8.0_s0.3_I12.0_f1.0_PB.csv",
              "point_B2.0_R8.0_s0.3_I12.0_f1.0_LR.csv"]:
    pb, lr, n = analyse(fname)
    print(f"{fname}")
    print(f"  pulses analysed: {n}")
    print(f"  PB state: median free RDF = {st.median(pb):.5f} uM (n={len(pb)})")
    print(f"  LR state: median free RDF = {st.median(lr):.4f} uM (n={len(lr)})")
    print(f"  ratio of medians = {st.median(lr) / st.median(pb):.0f}x")
    print()

print("Figure rounding used: LR -> 4.5 uM, PB -> 0.002 uM, contrast ~2500x.")
