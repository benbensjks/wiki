# Engineering — Model-guided Design of the Single-bit Counter

> Companion to `01_Chapter3_Wiki_Draft.md` and `02_Figure_Brief.md`.
> Chinese reading aid: `zh_aux/03_Engineering构想_中文辅助.md`.
> Purpose: source material for the Engineering page and the Modeling-Guided Design panel.
> Engineering logic to be fixed by 9.30; final text by 10.7.
> Evidence strength is labeled throughout: **[certified]** strictly verified inside the model /
> **[mechanism-inferred]** supported by model mechanism and scans, not yet experimentally
> verified / **[pending calibration]** requires measurement before a value can be fixed.

---

## 1. Engineering Question and Narrative

The engineering question of the single-bit module:

> How can a continuously oscillating upstream signal be converted into a reliable digital
> memory element that advances exactly one step per cycle? Which biological parts must be
> chosen, over what ranges must their expression and clearance be tuned, and how can those
> ranges be turned into measurable targets for the wet-lab team?

The narrative follows the project-wide **Question → Model → Result → Design decision** and
**Design → Build → Test → Learn → Redesign** double loop.

---

## 2. Engineering Cycles (three rounds + one cross-module iteration)

### Round 1 — From the original TetR to BM3R1, exposing the input-interface problem

| Step | Content |
|---|---|
| **Design** | Replace the TetR delay circuit of the Zhao 2019 single-input counting module with **BM3R1**, which is orthogonal to the oscillator; keep the ϕC31 Int/RDF recombination core unchanged |
| **Build** | 38-state ODE: the original 35 recombination equations + 3 equations for the BM3R1 delay circuit; the A-module v36 waveform (period 6.77 h, baseline 0.28 µM) as the Int source |
| **Test** | Porting validation (reproducing the original LR 0.31/0.88 alternation), default-parameter failure diagnosis, multi-dimensional scans (input baseline / delay strength / tag / K·n) |
| **Learn** | 1) default parameters fail under the new input; 2) the failure mechanism is "baseline leak causes spontaneous flipping + insufficient delay causes within-pulse back-flipping"; 3) **the only hard interface requirement is int0 promoter leak ≤1.5%**; 4) counter-side tuning alone cannot rescue the old waveform |
| **Decision** | Request a tighter PLtetO1 leak from the A module / wet-lab team; publish a first BM3R1 expression range ($k_{BM3,tsl}$ 15–30 h⁻¹) |

**Round takeaway** — the model did not prove "it works"; it **defined the interface
boundary**, turning "why the counter fails" into a concrete quantity the A module could change.

> After this round the B report was revised to v3 (2026-08-01) following the group's internal
> audit report (corrected statements; rebuilt the second-round scan). This is a
> "results audited → delivery revised" record and can be mentioned briefly on the wiki.

### Cross-module iteration — A↔B interface feedback and rework (between Rounds 1 and 2)

This is the most complete **design-feedback loop** in the project: the interface between two
models was not defined once and for all — it was corrected jointly by failure diagnosis on
both sides. Evidence from both sides is documented in their Week 4 reports:

```text
B (Week 3): counting fails under the real waveform
    ↓ feedback: leak / baseline / pulse width are interface requirements
A (Week 4): explicit-mRNA model + flux interface v31(t)=ρ(t)β31m31(t)
    ↓ A-side check: old waveform + new period still flips; new waveform + old period fails
    ↓ diagnosis: the interface problem comes from the flux peak / dose / width combination, not the period
B (Week 4): re-run (rework) with A's frozen input
    ↓ cross-validation: independent agreement with A's v53g/v53h (RBS 0.45 / tag 8), 80/80
B additions: tag 12 gives a wider window; krdf,tsl = 100 usable (tag 12); K 12–50 nM all pass
```

| Side | Action | Key quantitative result | Source |
|---|---|---|---|
| B (Week 3) | Drove the counter with A's v36 waveform; found no solution | only hard requirement: int0 leak ≤1.5%; baseline ≤48 copies needed to count | B Week 3 §7.2–7.4 |
| A (Week 4) | Switched the interface to production flux; checked old/new waveform combinations | flux pulse width 2.91 h → 3.45–3.62 h; "old waveform + new period flips, new waveform + old period fails" | A Week 4 §6.1 |
| A (Week 4) | Independently scanned the B-side RBS×tag and ran an A/B self-consistent load validation | RBS 0.45 stable for tag 4–20 h⁻¹; RBS 0.45/tag 8 passes 80/80 over 4 load levels × dual initial states | A Week 4 §7.1–7.2 |
| B (Week 4) | Rework: re-ran all conclusions with A's frozen input | old working point fails → new working point; K-sensitive band disappears (12–50 nM all pass) | B Week 4 §5, §9.1 |
| B (9/19) | Promoter audit overturns the working point again | under 0.8% non-zero leak the Week 4 point fails → new scheme (β_B=2, β_R=8, scale=0.30) | B 9/19 package |

**Writing suggestion (Engineering page)**: present this as "the interface was designed, not
assumed" — the two models were aligned at the interface before each side ran its own scans;
every rework was triggered by a **reproducible failure** and produced a new design constraint.
Do not write it as "we talked and fixed a bug" routine troubleshooting (see §6).

### Round 2 — Re-validation under the frozen clock and Int degradation-tag design

| Step | Content |
|---|---|
| **Design** | Accept A's Week 4 frozen input (production-flux interface, period 10.59 h, measured leak 0.5%); make the "Int degradation tag" a formal design variable |
| **Build** | Interface becomes $J_I(t)=s_I\cdot v_A(t)$ (µM/h); Int binding / dilution / tag clearance handled inside the B module; the tag acts on all Int-containing species |
| **Test** | Working-point re-scan (expression × input × tag), RBS×tag joint workspace, K/n robustness, heterogeneity Monte Carlo (30 instances per level) |
| **Learn** | 1) the working point holds under the new input (dual-initial-state score 0.995–0.998); 2) **the K-sensitive band disappears** (12–50 nM all pass) — the longer period gives RDF enough recovery time; 3) the tag is **required** (tag = 0 always fails) and the window is bounded on both sides; 4) the upper edge is not a monotonic "peak clipping" process but a fragmented band of failed LR→PB flips; 5) the stochastically stable window (6–18 h⁻¹) is narrower than the deterministic window (3–24 h⁻¹) |
| **Decision** | Give the wet-lab team **RBS 0.45 / expression ratio / tag gradient ruler (native ssrA / LAA+4 / SsrA2X / LAA-LAA)**; state that "tag rates must be measured under matched conditions" as a prerequisite |

**Round takeaway** — **the tag is not an optional optimization but a functional necessity**,
and its rate must be paired with the input strength; the sequence alone is not enough.

### Round 3 — Promoter audit and scheme reconstruction under non-zero leak

| Step | Content |
|---|---|
| **Design** | Audit the RDF promoter about to be synthesized by the wet-lab team; incorporate an explicit non-zero leak into the model (Cello B1 gate $y_{min}/y_{max}=0.8\%$) |
| **Build** | Add the leak floor $\ell_R$ that cannot be pushed lower by infinite repression to the RDF mRNA equation; introduce the fraction $\eta$ of RDF returned upon Int clearance |
| **Test** | Joint search over expression × input × clearance; 430 h dual-initial-state long runs; feasible-region grids ($\beta_B\times\beta_R$, $\ell_R\times\beta_R$); start-phase and parameter perturbations; provenance audit of the 36 bp sequence |
| **Learn** | 1) **the Week 4 working point fails under 0.8% leak**; 2) a re-matched scheme exists: $\beta_B=2$, $\beta_R=8$, scale 0.30, $k_I=12$; 3) the feasible region is a **wide plateau** ($\beta_R$ 4–8 × $\beta_B$ 0.75–6); 4) **leak tolerance ≈ 1%**, and the leak floor scales with the maximum promoter strength — "a stronger RDF promoter can be worse"; 5) the wet-lab 36 bp variant is AI-generated with missing bases and should not be synthesized as is |
| **Decision** | 1) switch the recommendation to the characterized **Cello pBM3R1 66 bp**; 2) publish the expression-ratio target (BM3R1:RDF ≈ 1:4 maximum capacity); 3) write "trough start" into the experimental conditions; 4) make the η test the next critical experiment |

**Round takeaway** — **promoter "strength" is not better-is-better, because leak scales with
strength**; the design goal moves from "pick a strong promoter" to "calibrate expression ratio
and leak together".

---

## 3. Design Rules (Modeling-Guided Design panel material)

| Design quantity | Recommended range / target | Basis | Evidence strength |
|---|---|---|---|
| RDF-side promoter | Cello pBM3R1 66 bp (reference) | characterized; 36 bp variant audited and rejected | [certified] (sequence) / [pending calibration] (MC4100 strength) |
| RDF relative leak $\ell_R$ | ≤0.8% (design cap 1%) | feasible-region scan: at $\beta_R=8$, $\ell_R\le0.01$ | [mechanism-inferred] |
| BM3R1 : RDF maximum capacity | ≈ 1 : 4 (0.25) | center of the non-zero-leak feasible region | [mechanism-inferred] |
| BM3R1 expression | max flux 1.02 µM/h; PB steady state ≈740 copies/cell | model-target conversion | [pending calibration] |
| RDF expression | max flux 4.08 µM/h; de-repressed steady state ≈2950 copies/cell | model-target conversion | [pending calibration] |
| Int input scale | 0.30 (paired with $k_I=12$); adjust via the pairing rule once $k_I$ is measured | $s_I(k_I)=0.30(\mu+k_I)/(\mu+12)$ | [mechanism-inferred] |
| Int degradation tag | deterministic [3, 24] h⁻¹; stochastically stable 6–18 h⁻¹; operating point 8–16 | dual-initial-state scans + heterogeneity MC | [mechanism-inferred] (rate) / [pending calibration] (fusion protein) |
| Start-up condition | known initial state + trough entry into the first complete pulse | 20 phases × dual initial states 37/40 | [mechanism-inferred] |
| Shared interface | BM3R1 ($T$) pool visible to the shutdown module; no extra BM3R1/RDF tags in this module | Chapter 5 BM3R1–sRNA mechanism | [mechanism-inferred] |

---

## 4. Learn List: Where the Model Changed Design Intuition

1. **"A stronger repressor gives a longer delay" is wrong.** Too strong → RDF unlocks too late
   → reverse flip fails; too weak → within-pulse double flipping. The delay must fall inside a
   window that "covers the pulse width but does not exceed the dwell".
2. **"A faster tag clears more cleanly" is wrong.** The tag upper edge is a fragmented band
   (failed LR→PB flips, stalled cycles), and it moves up with RBS.
3. **"A stronger RDF promoter is better" is wrong.** Under non-zero leak, the leak floor
   scales with maximum strength, and an over-strong promoter destroys the feasible region
   (high leak requires low $\beta_R$).
4. **The working point is not a constant.** Changing the input waveform (period 6.77 → 10.59 h)
   moves the whole operating band (K-sensitive band disappears, tag window widens); changing
   the promoter leak (0 → 0.8%) overturns the entire working point. A design delivery must be
   "window + conditions", not "a set of numbers".
5. **Initialization is a design condition.** "Arbitrary phase entry" does not pass (37/40);
   trough start plus a known initial state is part of the current scheme and must be arranged
   explicitly in the experimental plan.
6. **Interface semantics decide conclusions.** A→B must transmit production flux (µM/h), not
   concentration; otherwise tag/dilution are double-counted, and such interface errors do not
   surface in single-module debugging.
7. **Rework is part of the design process, not a failure.** Three working-point replacements
   (v36 → Week 4 → non-zero leak) each came from a new, reproducible failure and produced a
   sharper design constraint.

---

## 5. Interfaces to Experiment (inputs to the next Redesign)

By priority:

1. **Expression calibration**: measure the steady-state expression and leak of the current
   constructs for BM3R1 / RDF, convert to µM / copy number, compare with §3 targets, and
   adjust RBS (or promoter choice);
2. **Int clearance rate**: measure the de-induction decay curve of the Int–tag fusion and
   select the input scale via the pairing rule; if measured $k<3\ \mathrm{h^{-1}}$, consider
   ClpX co-expression;
3. **η determination**: use same-condition RDF stability / content changes to test whether Int
   clearance additionally consumes RDF; if significant, switch to the $\eta=0$ alternative
   parameter set and re-scan;
4. **Stepwise pBAD-INT characterization**: first characterize the switch with a controlled
   input (square-wave comparison: peak 2 µM/h, width 2–3 h passes), then couple to the A
   oscillator and check initialization.

---

## 6. Judging Criteria Alignment (2026 Judge Handbook / Best Model)

### 6.1 The four Best Model rubric aspects (handbook wording)

1. How impressive is the modeling? (span / integration / engineering use)
2. Did the model help the team understand a part, device, or system?
3. Did the team use measurements of a part, device, or system to develop the model?
4. Does the modeling approach provide a good example for others?

**B-module material**: the single bit sits at the **device level** (switch mechanism + flip
window). "Using measurements to develop the model" currently relies on published measurements
(Zhao/Pokhilko parameters, Cello gate characterization); the team's own measurements are
still missing — the text must state this honestly and must not fabricate a validation loop.
Reusability material: the window-scan scripts, criteria scripts and parameter tables can be
reused by other teams with replaced parameters.

### 6.2 iGEM engineering methodology checklist (handbook "On Engineering")

- Identify the problem and demonstrate understanding → gather data (and cite sources) and
  recognize unknowns and constraints → select applicable guiding principles and theories →
  **list assumptions, approximations and simplifications** → **establish quantifiable measures
  of success** → show how the problem was solved → validate the results → communicate the solution.

**B-module mapping**: assumptions and evidence labels are explicit (§3; wiki draft §12); the
success criteria are endpoint $Q$ + one midpoint crossing + counting parity (publicly defined);
"validation" is dual-initial-state long runs + perturbations + heterogeneity MC (all in-model).

### 6.3 Two warnings from the handbook (must be avoided in writing)

1. **Do not write engineering as routine troubleshooting**: the handbook quotes a judge
   feedback — "much of your work under the engineering criterion reads more like routine
   troubleshooting than actual design-driven engineering."
   → Every rework here must be written as "design question → model criterion → design decision",
   not "found an error → fixed it".
2. **"Failure is central to engineering"**: the handbook states that success and failure are
   equally valuable.
   → The level-three cascade failure, the fragmented tag upper edge, and the non-zero-leak
   overturn should all remain on the wiki and show how they changed the design.

### 6.4 Group-level writing requirements already issued

| Source | Requirement | Implementation here |
|---|---|---|
| 《注意事项.pdf》 (9/13 group file) | **Prepare your own engineering as Design → Build → Test → Learn**; coordinate the mechanism figure with the art team; unified colors (blue/orange, see Heidelberg 2025 Model); large fonts, thick lines, English; parameter tables as images; collect references; single-page anchor TOC (Waseda-Tokyo 2024 Model pattern) | this document is the D/B/T/L material; figures in `02_Figure_Brief.md`; parameter table to be rendered as an image; anchor TOC planned |
| 《Best Model思路.pptx》 (9/13, prepared by the group lead) | Judges will ask: what modeling was done, assumptions and rationale, what data were used, how results affected design; includes a "Design question → Modeling insight → Consequence" table | every section of the wiki draft follows "question → model → result → design decision" |
| Cui Xiyan (senior, 9/13 chat) | Read other teams' wikis (especially foreign teams); be rigorous, concise, understandable; coordinate page splits with the web team | benchmark analysis in §7 |
| Cao Ruoxi (senior, 9/20 chat) | No external website links on the wiki; last year's team hosted files in the web git | internal-hosting plan for the Software page (already designed by the group lead) |

---

## 7. External Benchmarks: Engineering and Wiki Writing (critical review)

> Observations from public wikis accessed on 2026-09-21. Cited only for writing-method comparison.
> **Method note**: some iGEM wikis render their content client-side; raw HTTP fetches can show an
> empty page. All pages below were checked in a browser render (headless Chrome) before analysis.

### 7.1 IZJU-China 2025 (Engineering page)

**Form**: the whole page is organized as DBTL1 → DBTL2 → DBTL3 → DBTL4 …, each containing
fixed **Design / Build / Test / Learn** blocks; long processes are split into sub-numbers
(e.g., DBTL3.1–3.4), each with figures and captions, data tables, and references at the end.

**Worth adopting**:

1. **One DBTL = one design decision** (target selection, cleavage strategy, delivery choice,
   membrane disruption); the headings themselves are decision names;
2. **Failures and abandoned routes are reported honestly**: Cas9-FITC labeling had poor
   recovery → abandoned → switched to sgRNA-FITC; gel readout was confounded by the marker →
   abandoned → switched to UV-Vis. Each switch states why and what came next;
3. **High numerical density**: AUC, recovery rates, labeling efficiency, significance tests,
   all tabulated;
4. **Captions carry the argument**: not "as shown", but conditions, controls, and conclusions.

**Critique / not transferable**:

- Extremely wet-lab-heavy (gels, AFM, UV-Vis); **a dry-lab project has no such material
  density**, and imitating it would mean padding with experiments;
- Some Learn blocks merely restate results ("signal increases with kinase") without
  extracting a transferable design rule;
- The page is very long, managed by sub-numbering and anchors; without equivalent navigation
  it would overwhelm readers.

### 7.2 Heidelberg 2025 (Engineering page)

**Form**: overview cards for Wet Lab / Dry Lab / HP, then engineering phases
("Processing → Responding → …"), each containing Iteration 1…n, each iteration again
**DESIGN / BUILD / TEST / LEARN**; dense figures with statistical tests marked on the plots.

**Worth adopting**:

1. **Iteration numbering + a one-line takeaway**: readers see at a glance which round solved what;
2. **Strong controls and counter-evidence** (catalytically dead mutants, coiled-coil-free
   variants) proving the function truly comes from the designed mechanism;
3. **Captions contain method details** (transfection amounts, timing, normalization),
   supporting reproducibility.

**Critique / not transferable**:

- Also wet-lab driven; its "Learn" sometimes only states the next step, with less on how the
  round changed the overall design — exactly where a modeling-led project can do better.

### 7.3 USTC 2025 (our previous team)

**Form**: the Engineering page is organized as **Iteration 1–4**, each containing several
cycles, each cycle with **Design / Build / Test / Learn** (≈35k characters). The Model page
presents five models, each with Objective / Challenges / Equations / Parameters / Results /
**Guidance to Wet Lab** / Discussion / References, plus a self-hosted interactive platform
(static HTML inside the team's own web path).

**Worth adopting**:

1. **Iteration + cycle numbering with a one-line title** (e.g., "Iteration 3: Optimize the
   strain ratio (Model)") — the reader can trace every design move;
2. **Honest, specific failure records**: the two-dimensional Gaussian fitting model could not
   be validated because experimental data were too scarce → the team learned to prioritize
   mechanistic models; a "self-updating" algorithm was abandoned after generating large
   computational errors → it led to the new "flux" algorithm. The failure and the pivot are
   both documented;
3. **"Guidance to Wet Lab" as an explicit closing block of each model** — exactly the
   "model → design decision" link the judges look for; we should mirror this in Chapter 3;
4. **Self-hosted interactive content** in the team's own web directory, consistent with the
   group rule that the wiki must not link to external sites.

**Take-away for our page**: the iteration numbering, the honest failure records and the
"Guidance to Wet Lab" closing block are directly reusable; keep the model-driven D/B/T/L
blocks on one page with explicit links.

### 7.4 Writing decisions derived from the benchmarks (for the B-module Engineering page)

| Decision | Approach |
|---|---|
| Structure | D/B/T/L blocks, but each round named as a model-driven design cycle (no padding with experiments) |
| Granularity | one round = one reproducible failure + one new design constraint; three rounds + one cross-module iteration |
| Captions | every figure states condition / control / conclusion and carries the argument |
| Failures | level-three failure, fragmented tag edge, non-zero-leak overturn all retained with design consequences |
| Numbers | all ranges carry units and criteria definitions; evidence labels retained in text |
| Navigation | single-page anchor TOC + one-line TL;DR per round; confirm page split with the web team |
| Closing block | end the chapter with an explicit **"Guidance to Wet Lab"** block (pattern adopted from USTC 2025), listing model → experiment decisions |
| Interactive content | self-host any interactive component inside the wiki path (no external links, per the group rule) |
| Reuse | emphasize parameterized scripts and criteria separated from the model, reusable by other teams |

**Status of adoption in this delivery (traceability)**:

| Borrowed element | Source | Implemented in |
|---|---|---|
| Single-page anchor TOC | group rule / Waseda-Tokyo pattern | `01_Chapter3_Wiki_Draft.md` §Contents |
| "Guidance to Wet Lab" closing block | USTC 2025 | `01_Chapter3_Wiki_Draft.md` §12 |
| Data-figure caption drafts (condition / control / conclusion) | IZJU-China 2025, Heidelberg 2025 | `01_Chapter3_Wiki_Draft.md` §15 |
| Iteration titles + one-line takeaways per round | Heidelberg 2025 | `03_Engineering.md` §2 ("Round takeaway") |
| Parameter table delivered as an image | group rule (《注意事项.pdf》) | `figures/fig03_7_parameter_table.png` (+ `scripts/make_parameter_table.py`) |
| Failure and pivot records retained | IZJU-China 2025, USTC 2025 | §2 rounds, §4 Learn list, §7.3 |
| Self-hosted interactive content | USTC 2025 / group rule (no external links) | planned with the web team (§9) |

---

## 8. Boundaries and Writing Conventions

- **Do not write "validated"**: all conclusions are model predictions; use "candidate working
  point" / "model target" until experimental calibration is complete;
- **Do not mix parameter sets**: Week 4 (zero-leak assumption) and the non-zero-leak scheme
  are two different model conditions; parameters and conclusions must not be cross-cited;
- **Do not present engineering criteria as experimental standards**: $Q$, midpoint crossing
  and counting parity are engineering acceptance criteria; definitions and raw traces must be
  given together;
- **Do not write rework as troubleshooting**: each rework is "new failure → new constraint →
  new working point";
- **Do not claim a complete suicide / clearance system**: this chapter covers single-bit
  counting and the output-duration interface only; effectors and cell clearance are out of scope.

---

## 9. Web Component Suggestions (for the web team)

| Page location | Suggested component | Interaction |
|---|---|---|
| Mechanism (Figure 3-1) | mechanism figure + flip animation | hover definitions; play one cycle |
| Design window (Figure 3-2) | concept figure + hover | hover a failure mode to expand its mechanism |
| Operating window (Figure 3-4) | heatmaps + sliders | drag (β_B, β_R) or ℓ_R and see pass/fail live |
| Tag window (Figure 3-5) | window strip chart | hover a tag level to see its failure mode (multi-flip / failed reverse flip) |
| Start-up conditions (wiki draft §5.4) | pass-rate bars | toggle "fixed trough / random phase" |
| Cross-module iteration (§2) | timeline component | click a node to expand the A / B failures and corrections |
| Engineering panel | Design → Build → Test → Learn → Decision cards | consistent with the Modeling-Guided Design panel |
