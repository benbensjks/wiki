"""Machine-readable working-point record for any artefact.

Why a separate module
---------------------
Six analyses were silently produced at n_A1_gate = 4 (the ZENG table value)
instead of the frozen 6.0, and three of them recorded nothing at all about the
model that produced them, so the mistake was invisible in the artefacts
themselves.  Every artefact must therefore carry the point it was taken at.

This lives in its own file on purpose: plausibility_common.py is covered by
decoupling_verification.json's implementation_sha256, so adding the helper there
would invalidate that verification and force it to be re-run.

    from working_point import working_point_block
    report['working_point'] = working_point_block(model)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plausibility_common import (ROOT, SELECTED_N_A1_GATE, ZENG_N_A1_TABLE,  # noqa: E402
                                 sha256, source_hashes)

MODEL_FILE = ROOT / 'model_threebit51.py'


def working_point_block(model=None, extra: dict | None = None) -> dict:
    """Record the working point an artefact was produced at.

    `model` is optional; when supplied the realised gate exponent and the
    interface scale are read back from the constructed model, which is what makes
    the record a measurement rather than a restatement of intent.
    """
    import model as M

    block = dict(
        model_file='final_reconstruction/model_threebit51.py',
        model_sha256=sha256(MODEL_FILE),
        zeng_n_A1_table=ZENG_N_A1_TABLE,
        selected_n_A1_gate=SELECTED_N_A1_GATE,
        live_zeng_n_A1=float(M.ZENG['n_A'][1]),
        source_sha256=source_hashes(),
        note=('n_A1_gate is a constructor parameter that overrides ONLY the carry-2 gate '
              'arm; leaving it unset falls back to ZENG["n_A"][1] = 4.0, which is NOT the '
              'frozen working point. Every artefact must show which value was in force.'),
    )
    if model is not None:
        block['n_A1_gate_effective'] = float(model.n_A1_gate_effective)
        block['n_A1_gate_is_frozen_point'] = bool(
            float(model.n_A1_gate_effective) == float(SELECTED_N_A1_GATE))
        block['uM_per_au'] = float(model.e.uM_per_au)
        try:
            block['carry1_mrna_half_life_min'] = float(model.carry.carry1.mrna_half_life_min)
            block['A1_maturation_half_life_min'] = float(
                model.carry.carry1.activator_maturation_half_life_min)
            block['F1_maturation_half_life_min'] = float(
                model.carry.carry1.repressor_maturation_half_life_min)
        except AttributeError:
            pass
    if extra:
        block.update(extra)
    return block
