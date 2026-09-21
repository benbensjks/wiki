# Chapter 3 Figure Brief — Concept Figures and Art Hand-off

> Companion to `01_Chapter3_Wiki_Draft.md`. Chinese reading aid: `zh_aux/02_机制图需求与草图_中文辅助.md`.
> Scope: the concept figures handed to the art team, plus one optional concept asset.
>
> **Delivered this round**
> - Figure 3-1 mechanism sketch: `figures/fig03_1_single_bit_mechanism_sketch.png`
> - Figure 3-2 design-window sketch: `figures/fig03_2_design_window_sketch.png`
> - Figure 3-7 parameter table image: `figures/fig03_7_parameter_table.png`
> - optional generative concept draft: `figures/optional_concept_cell_context.png`
>
> **Status: draft**; may be adjusted together with the text until the 9.30 freeze.

---

## 1. Figure 3-1 — State-dependent switch (mechanism)

### 1.1 What the reader must understand

One sentence:

> **Why the single bit remembers 0/1: the DNA holds the state, and the RDF pool holds the
> direction of the next flip; at each clock pulse, the pool level decides which reaction wins.**

Two design conclusions that must come across:

1. **One pulse → one flip** holds only if the **RDF accumulation delay exceeds the Int pulse
   width**; otherwise the DNA flips back within the same pulse;
2. **BM3R1's role is delay, not "stronger is better"**: it keeps RDF production off in the PB
   state, and decays after the flip so that RDF can accumulate.

### 1.2 Layout

Two state cards side by side, with two pulse arcs between them:

- **PB card (left)**: DNA icon `attP × attB`; BM3R1 expressed (arrow from the DNA) and
  repressing the RDF gene (T-bar); RDF pool gauge nearly empty; note "pool emptied between
  pulses".
- **LR card (right)**: DNA icon `attL × attR`; BM3R1 faded ("decays after the flip");
  RDF gene de-repressed (arrow from the DNA); pool gauge nearly full; note "pool fills during
  the LR dwell".
- **Top arc (PB → LR)**: labelled `clock pulse 1: free Int drives forward recombination`;
  **bottom arc (LR → PB)**: labelled `clock pulse 2: Int–RDF complex drives reverse
  recombination`. A small square-pulse glyph precedes each label.
- **Footer rules**: "Memory = the RDF pool: ≈4.5 µM (LR) vs ≈0.002 µM (PB) at pulse arrival,
  ≈2500× contrast" and "Constraint: RDF delay > pulse width — otherwise RDF unlocks during
  the pulse and the DNA flips back within the same cycle".

### 1.3 What must not be drawn wrong

1. **Memory is stored in the RDF pool, not in Int.** Int is the "read/write head" arriving
   once per clock pulse;
2. **Reverse recombination is driven by the Int–RDF complex**, not by free Int or free RDF alone;
3. **State-dependent production must be correct**: PB → BM3R1 ON, RDF OFF; LR → BM3R1 OFF
   (decaying), RDF ON;
4. **The pool levels are the mechanism**: nearly empty in PB, nearly full in LR at pulse arrival;
5. **Do not draw RDF as an inhibitor** (e.g., inhibiting or degrading Int) — it is a
   directionality factor forming a complex with Int.

### 1.4 Copy-ready English labels

```text
Single-bit switch: state-dependent memory
PB state - bit = 0
LR state - bit = 1
DNA: attP x attB
DNA: attL x attR
BM3R1
expressed
RDF gene
pBM3R1 repressed -> RDF OFF
de-repressed -> RDF ON
decays after the flip
RDF pool (memory)
~0.002 uM at pulse arrival
~4.5 uM at pulse arrival
Pool emptied between pulses (BM3R1 keeps RDF production off)
Pool fills during the LR dwell (no BM3R1 left to repress it)
clock pulse 1: free Int drives forward recombination (PB -> LR)
clock pulse 2: Int-RDF complex drives reverse recombination (LR -> PB)
Memory = the RDF pool: ~4.5 uM (LR) vs ~0.002 uM (PB) at pulse arrival, ~2500x contrast
Constraint: RDF delay > pulse width
```

### 1.5 Caption draft

> **Figure 3-1. State-dependent switch (concept sketch).**
> The DNA element adopts two configurations, PB and LR. In the PB state BM3R1 is expressed and
> represses RDF production, so the RDF pool is nearly empty at pulse arrival (~0.002 µM); the
> arriving clock pulse therefore drives the forward reaction (PB→LR). In the LR state BM3R1
> decays and RDF accumulates (~4.5 µM at pulse arrival, ≈2500× contrast), so the next pulse
> forms the Int–RDF complex and drives the reverse reaction (LR→PB). Memory resides in the
> RDF pool; the delay must exceed the pulse width.

---

## 2. Figure 3-2 — Design window and failure modes (concept)

### 2.1 What the reader must understand

> **Why the design window has edges**: each side fails for a distinct mechanistic reason,
> so the window is a joint property of the tag rate and the expression ratio, not a single
> safe number.

### 2.2 Layout

Two panels:

- **Panel A — tag window**: a horizontal band over $k_{tag,int}$ with three zones
  (too weak / deterministic pass [3, 24] h⁻¹ with the stochastically stable sub-band
  [6, 18] / too strong); below it, three schematic $S(t)$ mini-traces with captions
  ("inter-pulse Int → extra flips", "one pulse = one flip", "peak clipped → stuck in LR").
- **Panel B — joint expression window**: schematic $(\beta_R,\beta_B)$ plane with a green
  pass plateau ($\beta_R$ 4–8 × $\beta_B$ 0.75–6, measured), failure labels on each side
  (RDF too weak → reverse flip incomplete; RDF too strong → leak amplified → spontaneous
  flips; BM3R1 too weak → pool not emptied → multi-flip), a dashed "beyond tested range"
  box above, and a dashed arrow "higher promoter leak shrinks the window toward lower RDF
  strength".

### 2.3 What must not be drawn wrong

1. The figure is **schematic**; the measured pass/fail map is the feasible-region data figure
   (Figure 3-4B) — do not present the schematic as measured data;
2. The window numbers: deterministic **[3, 24] h⁻¹**, stochastically stable **[6, 18] h⁻¹**;
   the plateau **β_R 4–8 × β_B 0.75–6**;
3. The failure mechanisms must stay on the correct sides (weak tag → extra flips; strong tag →
   stuck in LR; weak RDF → incomplete reverse flip; strong RDF → leak-amplified flips);
4. The leak arrow direction: higher leak → window moves toward **lower** RDF strength.

### 2.4 Caption draft

> **Figure 3-2. Design window and failure modes (concept sketch).**
> (A) The Int degradation-tag window has edges: too weak leaves inter-pulse Int and causes
> extra flips; too strong clips the pulse peak and the reverse flip fails (stuck in LR);
> the deterministic pass band is [3, 24] h⁻¹, the stochastically stable band [6, 18] h⁻¹.
> (B) The joint expression window is a plateau (β_R 4–8 × β_B 0.75–6, measured under the
> non-zero-leak scheme); each side fails for a distinct mechanistic reason, and a higher
> promoter leak shifts the window toward lower RDF strength. Schematic; the measured
> pass/fail map is the feasible-region data figure.

---

## 3. Optional concept asset — cell context (generative draft)

`figures/optional_concept_cell_context.png` — a flat-style illustration of the switch in its
cellular context (E. coli outline, plasmid, highlighted DNA segment, Integrase, BM3R1, RDF
pool, clock pulse), with English labels added.

- **Use**: chapter opener / hero image, or a context illustration near the mechanism section.
- **Status**: generative concept draft; the art team should unify style, palette and line
  weights before publication. Keep or drop as the page design requires.

---

## 4. Reserved companion figure — one-pulse-one-flip timing schematic

A timing schematic (Int pulse / $S$ / BM3R1 / RDF pool over one cycle plus the next pulse,
annotated `delay > pulse width`) is **reserved for later refinement** and is not part of this
hand-off. When refined, it can use real simulation traces instead of schematic curves.

---

## 5. Visual conventions (project-wide)

| Use | Color | Hex |
|---|---|---|
| Structure boxes / DNA states / stable states | TEMPO deep blue | `#304B53` |
| Time evolution / thresholds / risks and emphasis | TEMPO orange | `#D29144` |
| Recombinase and model results | Dry Lab teal-green | `#4F9194` |
| Module light fill | light blue | `#DCE5E7` |
| Threshold / memory region background | light orange | `#F2E2CF` |
| Passing states / complexes | light teal | `#DCEBEC` |
| Pending parameters / secondary lines | blue-gray | `#8CA0A5` |
| Body and axis text | charcoal | `#33383A` |

Line conventions: repression = solid line with T-bar; activation = solid arrow;
reversible binding = double-headed arrow; module boundary = dashed box.
Suggested component colors: Int = teal; RDF = orange; BM3R1 = deep blue;
Int–RDF complex = charcoal outline + light-teal fill; DNA PB = light blue, DNA LR = light teal.

---

## 6. Interaction / animation suggestions (for the web team)

- **Figure 3-1 flip animation**: play one clock cycle; pause at the flip and annotate
  "direction set by the RDF pool"; hover definitions for Int / RDF / BM3R1 / complex / PB / LR;
- **Figure 3-2 hover**: hover a failure mode to expand the mechanism and a mini trace;
- **Optional**: replace the reserved timing figure's schematic curves with real simulation
  traces after the 9.30 freeze.

---

## 7. Deliverables

| File | Content | Status |
|---|---|---|
| `figures/fig03_1_single_bit_mechanism_sketch.png` | Figure 3-1 state-dependent switch | delivered |
| `figures/fig03_2_design_window_sketch.png` | Figure 3-2 design window and failure modes | delivered |
| `figures/fig03_7_parameter_table.png` | Figure 3-7 parameter table (image) | delivered |
| `figures/optional_concept_cell_context.png` | optional generative concept draft (cell context) | delivered |
| `scripts/make_mechanism_v2.py` | Figure 3-1 generator (editable, reproducible) | delivered |
| `scripts/make_design_window_v1.py` | Figure 3-2 generator | delivered |
| `scripts/make_parameter_table.py` | Figure 3-7 generator | delivered |
| `scripts/verify_pool_levels.py` | provenance check for the ≈0.002 / ≈4.5 µM pool levels | delivered |
| reserved timing schematic | see §4 | deferred |

> Candidate data figures (`candidate_*.png`) and their source mapping are listed in
> `figures/README.md`.
