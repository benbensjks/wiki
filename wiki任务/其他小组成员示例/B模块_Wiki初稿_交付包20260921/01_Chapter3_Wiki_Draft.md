# Coupling the Clock to a Reliable Single-bit Counter

> **Chapter 3 wiki draft (EN, primary).** Companion Chinese reading aid: `zh_aux/01_第3章wiki初稿_中文辅助.md`.
>
> This chapter corresponds to Chapter 3 of the modeling structure agreed at the 9.15 meeting
> and answers the B-module questions listed on slide 3 of that deck: the recommended BM3R1
> promoter range, the BM3R1–RDF relationship, where a degradation tag is needed,
> whether earlier conclusions still hold under the final clock waveform, and how the
> one-pulse-one-flip operating window is defined. Per the meeting, single-bit / coupling
> robustness is placed in this chapter.
>
> **Status: draft (not frozen).** Intended to let the web team build the page structure and
> the art team start the mechanism figure; content will be revised until the 9.30 freeze and
> updated as wet-lab calibration arrives.
> Numerical sources: B-module Week 3 / Week 4 delivery packages and the two 2026-09-19
> packages (non-zero-leak scheme; promoter audit). Evidence labels follow §12; no model
> result is presented as an experimental result.

---

## Contents

1. [Page Flow](#1-page-flow)
2. [Modeling Status](#2-modeling-status)
3. [The Problem](#3-the-problem)
4. [How the Switch Works](#4-how-the-switch-works)
5. [One Pulse, One Flip — Criteria and Operating Window](#5-one-pulse-one-flip--criteria-and-operating-window)
6. [Re-validation under the Final Clock](#6-re-validation-under-the-final-clock)
7. [The BM3R1–RDF Delay Circuit](#7-the-bm3r1rdf-delay-circuit)
8. [Int Degradation Tag](#8-int-degradation-tag)
9. [Promoter Recommendation](#9-promoter-recommendation)
10. [Robustness and Boundaries](#10-robustness-and-boundaries)
11. [Engineering — Model-guided Decisions (summary)](#11-engineering--model-guided-decisions-summary)
12. [Guidance to Wet Lab](#12-guidance-to-wet-lab)
13. [Scope, Limitations and Reproducibility](#13-scope-limitations-and-reproducibility)
14. [Chapter Conclusions](#14-chapter-conclusions)
15. [Figure List and Caption Drafts](#15-figure-list-and-caption-drafts)
16. [References](#16-references)

---

## 1. Page Flow

```text
Real upstream ϕC31 waveform (period 10.59 h, production-flux interface, µM/h)
          ↓ drives
Single-bit switch: Int / RDF / BM3R1 + one PB ⇄ LR DNA state
          ↓ mechanism
One pulse = one flip: memory lives in the RDF pool; the pool at pulse arrival sets the direction
          ↓ criteria
State-oriented acceptance criteria (endpoint Q + one midpoint crossing per cycle + counting parity)
          ↓ design
Expression-ratio window (BM3R1 / RDF) × Int degradation-tag window × RDF promoter leak tolerance
          ↓ stress tests
Start-up phase, parameter perturbations, heterogeneity Monte Carlo, fate of RDF upon Int clearance
          ↓ design decisions
66 bp promoter recommendation, expression-ratio targets, tag gradient ruler, initialization requirement
```

---

## 2. Modeling Status

| Item | Status | Current evidence boundary |
|---|---|---|
| ϕC31–RDF–BM3R1 single-bit model (38-state ODE) | Complete | Recombination core ported line-by-line from Zhao 2019 (35 equations); BM3R1 delay circuit replaces 3 equations; DNA conservation to 5 significant figures |
| One-pulse-one-flip mechanism | Complete | Failure diagnosis + working-point scans; flip direction set by the RDF pool at pulse arrival |
| Re-validation under the final clock (Week 4) | Complete | Default parameters fail; new working point holds; K-sensitive band disappears |
| Int degradation tag: window and mechanism | Complete | Deterministic window, stochastic-stable window, fragmented upper edge, literature rate check |
| Scheme under non-zero RDF promoter leak (0.8%) | Complete (model) | Dual initial states, 430 h, 35 steady cycles; fixed-phase perturbations pass; arbitrary phase still fails |
| Feasible region (expression ratio × leak) | Complete (discrete grid) | Wide plateau β_R 4–8 × β_B 0.75–6; leak tolerance ≈ ≤1% |
| Promoter sequence and wet-lab hand-off | Complete | Cello pBM3R1 66 bp recommended; the 36 bp variant was audited and rejected |
| RDF fate upon Int clearance (η) | Not done (experiment) | η = 1 and η = 0 parameter sets must not be mixed; η ≥ 0.9 passes, η ≤ 0.75 fails |
| Expression calibration (α, β, copy number) | Not done | β is a model translation coefficient, not a promoter RPU, and does not map to a unique RBS sequence |
| Molecular noise / single-cell success rate | Not done | Current Monte Carlo perturbs parameter scenarios; it is not a Gillespie single-cell simulation |

---

## 3. The Problem

> Can a bistable switch built from site-specific recombination still act as a **digital memory
> element** when driven by the **real upstream waveform** rather than an ideal square wave?
> What is the criterion for saying "yes"? How wide is its operating window? Which design
> quantities actually determine success?

Four sub-questions:

1. **Mechanism** — why does it flip once per pulse instead of following the input continuously?
2. **Interface** — what quantity is actually transmitted from the clock module, and do earlier
   conclusions survive the waveform update?
3. **Criteria and window** — $S$ is a continuous variable; on what basis is it read as 0 or 1,
   and how is the one-pulse-one-flip window defined?
4. **Design quantities** — what ranges are required for the BM3R1/RDF expression ratio,
   the Int degradation tag, and the RDF promoter leak?

---

## 4. How the Switch Works

### 4.1 State variables and the two DNA states

The counting unit is a flippable DNA segment flanked by ϕC31 `attP`/`attB` sites:

- **PB state (bit 0)**: `attP × attB`; with the LR fraction denoted $S$, the PB fraction is $1-S$;
- **LR state (bit 1)**: `attL × attR`.

$S$ is a continuous fraction, not a logic bit: the model has no threshold events, no event
queue, and no forced freeze during idle periods. Every digital behavior must emerge from
continuous dynamics. The single-bit model has 38 ODE states: mRNA–immature–mature chains for
Int, BM3R1 and RDF; free Int and its dimer; RDF; the Int–RDF complex; and all DNA–Int(±RDF)
intermediates and synaptic complexes.

### 4.2 Recombination in two directions

$$v_f = k_{fwd}\,H(I;K_{D,int},2)\,\frac{K_{inh}}{K_{inh}+R},
\qquad
v_r = k_{rev}\,H(C;K_C,2)$$

$$\frac{dS}{dt} = v_f\,(1-S) - v_r\,S$$

where $H(x;K,n)=x^n/(K^n+x^n)$; $I$ and $R$ are free Int and free RDF, and $C$ is the
Int–RDF complex. **Forward (PB→LR) is driven by free Int; reverse (LR→PB) is driven by the
Int–RDF complex.**

> Key mechanism: the reverse reaction depends on the **instantaneous product** $I\cdot R$
> (under complex quasi-steady state, $C \approx q\,I\,R$, equivalent to
> $H(I\cdot R;K_{D,comp},2)$), not on $R$ alone. Reverse recombination does **not consume**
> RDF — the drop in free RDF comes from binding sequestration and complex clearance.

### 4.3 The BM3R1 delay circuit

$$ \dot m_B = \alpha_B D_{PB} - \gamma_B m_B, \qquad \dot B = \beta_B m_B - \mu B $$

$$ \dot m_R = \alpha_R D_{LR}\left[\ell_R + \frac{1-\ell_R}{1+(B/K_B)^n}\right] - \gamma_R m_R $$

- BM3R1 (denoted $T$) is expressed **only in the PB state**; after switching, its mRNA and
  protein decay;
- RDF is transcribed **only in the LR state** and is repressed by BM3R1 through a Hill
  function; $\ell_R$ is the leak floor that cannot be pushed lower by arbitrarily strong repression.

The delay is approximately $\tau \approx \ln(B_{ss}/K_B)/\mu$ when only growth dilution clears
the protein. **The delay must exceed the entire pulse width**; otherwise RDF is unlocked
within the same pulse, the DNA flips back and forth, and counting fails — the central design
constraint of this module.

### 4.4 Why one pulse produces exactly one flip

1. **The RDF pool is the memory.** RDF accumulates during the LR dwell; in the PB state
   BM3R1 pushes RDF production close to zero and the pool is emptied. In the primary
   non-zero-leak scheme (430 h, dual initial states), the free RDF pool **just before pulse
   arrival** is ≈ **4.46 µM** in the LR state versus ≈ **0.0018 µM** in the PB state —
   a ≈ **2500-fold** (three orders of magnitude) difference.
2. **The pool at pulse arrival sets the direction.** If the previous state was LR, the
   complex is large and the reverse reaction wins; if the previous state was PB, the pool is
   empty and the forward reaction has no competitor.
3. **The flip is a monotonic single-direction transition**: LR→PB and PB→LR each complete
   monotonically within ≈ 2 h (LR falls from ~1.0 to ~0.003, or the reverse), with no
   back-and-forth inside one pulse — the direction is already determined before the pulse arrives.

> In the model, $dS/dt$ depends only on the instantaneous $I, R, C$; there is no delay term
> and no bistability term. All memory resides in the RDF pool, which is preserved between
> clock pulses.

### Figure 3-1 (mechanism figure, for the art team)

State-dependent switch: the PB and LR states, their state-dependent production, the RDF pool
level at pulse arrival, and the two pulse-driven flip directions. Layout, arrow semantics,
the points that must not be drawn incorrectly, and copy-ready English labels are given in
`02_Figure_Brief.md`; this round's sketch is
`figures/fig03_1_single_bit_mechanism_sketch.png`. A companion one-pulse-one-flip timing
schematic is reserved for later refinement.

---

## 5. One Pulse, One Flip — Criteria and Operating Window

### 5.1 Strict counting criteria

$S$ is continuous, so the readout must be stated explicitly. This chapter uses:

1. **Endpoint standard** — the high and low endpoints must be $\ge 0.95$ and $\le 0.05$;
   the worst per-point fraction of correct states, $Q$, is taken over all points, so
   intermediate states cannot be skipped;
2. **Exactly one midpoint crossing per cycle** — the $LR=0.5$ crossing is detected by
   **root-finding events on the integrator's internal steps**, not by sparse plotting points;
3. **Counting parity** — with a defined initial state and start phase, the counting parity and
   the total number of crossings are checked from the first complete input onward.

$Q$ is a DNA-endpoint metric, **not a single-cell counting success rate**.

> Statement: these are engineering acceptance criteria, not an experimental standard.
> Changing the criteria changes conclusions, so the definitions are published.

### 5.2 Operating window under the final clock input (zero RDF leak assumption)

Driven by the Week 4 frozen input (period 10.59 h, production-flux interface, measured
PLtetO1 leak 0.5%), the untuned default parameters fail. The re-scanned working point:

| Component / parameter | Recommended value | Window |
|---|---|---|
| ϕC31 RBS scale (A input interface) | 0.45 | — |
| $k_{BM3,tsl}$ (BM3R1 translation) | 15 h⁻¹ | 6–30 h⁻¹ |
| $k_{rdf,tsl}$ (RDF translation) | 200 h⁻¹ | ≥200 h⁻¹ robust (100 usable) |
| $K_{BM3}$ / $n$ | 18.6 nM / 3.4 | 12–50 nM × 2.9–3.4 all pass |
| $k_{tag,int}$ (Int degradation tag) | 8–16 h⁻¹ | deterministic [3, 24]; stochastically stable 6–18 |

Dual-initial-state scores are 0.995–0.998, and the steady-state valley sequence alternates
cleanly between 0.00 and 1.00 with fidelity 1.00.

> **The disappearance of the K-sensitive band** is one of the most important robustness
> results of this re-run: with the old v36 input the working band was only 18–35 nM
> ($K\le12$ collapsed, $K\ge50$ failed), whereas under the new waveform (period 10.59 h,
> inter-pulse gap ≈ 8 h) any $K$ between 12 and 50 nM and both Hill configurations count
> reliably — the wet-lab team no longer needs to hit $K_{BM3}$ precisely.

### 5.3 Current scheme under non-zero RDF promoter leak (0.8%)

The promoter audit revealed a non-negligible leak floor (Cello B1 gate,
$y_{min}/y_{max}=0.8\%$). Under this explicit leak the Week 4 working point fails; matching
expression, input and clearance again yields a scheme:

| Quantity | Week 4 point | Non-zero-leak primary scheme | Nature |
|---|---:|---:|---|
| BM3R1 translation coefficient $\beta_B$ | 15 h⁻¹ | **2 h⁻¹** | model target |
| RDF translation coefficient $\beta_R$ | 200 h⁻¹ | **8 h⁻¹** | model target |
| Int synthesis flux scale | 0.45 | **0.30** | relative to the frozen A input |
| Int extra clearance $k_I$ | 12 h⁻¹ | **12 h⁻¹** | effective model parameter |
| RDF promoter leak floor $\ell_R$ | 0 (assumed) | **0.008** | Cello B1 literature prior |
| Fraction of RDF returned after Int clearance $\eta$ | 0 (implicit) | **1** | working hypothesis, to be tested |

**Results**: dual initial states, 430 h, 35 steady full cycles pass strictly, worst
$Q=0.9966$; with a known initial state and trough start, all 40 complete pulses count
correctly; fixed-phase perturbations pass 200/200 (small) and 199/200 (wide); arbitrary
phase gives 37/40 — **the initialization requirement cannot be omitted**.

**Feasible region** (primary-scheme conditions, fixed trough start, 430 h):

- $(\beta_B,\beta_R)$ wide plateau: β_R 4–8 × β_B 0.75–6 all pass (60 cells); including
  extensions, 76/130 cells pass; β_R ≥ 16 has no solution under these conditions;
- leak tolerance: at β_R = 8, ℓ_R ≤ 0.01 passes (0.0125 starts to fail); for ℓ_R ≤ 0.006,
  β_R 4–16 all pass; higher leak requires lower β_R;
- boundary: the grid is discrete; random phase narrows the region; scale, tag and mechanism
  η are fixed and require re-scanning if changed.

**Figure 3-2** (concept sketch) explains *why* the window has edges: the failure mode on each
side (multi-flip / failed reverse flip / leak-amplified spontaneous flipping) and the shift of
the window toward lower RDF strength as the promoter leak rises.

### 5.4 Start-up conditions (initialization is a design condition)

| Scenario | Passed / tested | Note |
|---|---:|---|
| 20 start phases × PB/LR | **37/40** | arbitrary phase can still fail |
| Small perturbations, random phase | 193/200 | ±10% scenario set |
| Small perturbations, fixed trough | **200/200** | worst $Q=0.9933$ |
| Wide perturbations, random phase | 194/200 | ±20% scenario set |
| Wide perturbations, fixed trough | 199/200 | one case at $Q=0.9249$ |

Conclusion: **a known initial state plus entry from the flux trough into the first complete
pulse** is the start-up requirement of the current scheme; "insensitive to any access phase"
has not been achieved.

---

## 6. Re-validation under the Final Clock

### 6.1 Interface semantics: production flux, not concentration

The A→B interface is the **ϕC31 production flux** $J_I(t)=s_I \cdot v_A(t)$ (µM/h, where
$s_I$ is the relative RBS translation strength):

$$ \dot I_{tot} = J_I(t) - (\mu + k_I)\,I_{tot} $$

Int binding, dilution and tag clearance are handled inside the B module; C31 concentration is
**no longer imposed as an external state** (which would double-count tag and dilution). This
implements the project-wide interface requirement (concentration vs flux vs normalized
waveform) on the B side.

### 6.2 Waveform specification (Week 4 frozen input, RBS = 1 reference)

| Metric | v36 (old, used in Week 3) | Week 4 unloaded (new) |
|---|---|---|
| Period | 6.77 h | **10.59 h** |
| C31 protein peak | 3.99 µM | 7.01 µM |
| C31 protein trough | 0.28 µM | 0.104 µM |
| Flux peak / trough (RBS = 1) | — (inferred from protein) | 6.64 / 0.0360 µM/h (peak/trough 185) |
| Per-cycle flux dose (RBS = 1) | — | 22.8 µM (steady valley-to-valley integral) |
| PLtetO1 leak (model) | 5% (assumed) | **0.5%** (measured at the trough) |

The non-zero-leak scheme uses scale = 0.30: flux peak ≈ 1.992 µM/h, trough ≈ 0.0108 µM/h,
per-cycle flux integral ≈ 6.853 µM.

### 6.3 Findings

1. **Old conclusions fail under the old input and are re-established under the new one**:
   untuned combinations still fail with the new input; the re-scanned working point (§5.2)
   gives dual-initial-state scores of 0.995–0.998.
2. **The failure mechanism decomposes cleanly** (identified in Week 3, alleviated by the new
   input): baseline leak causes spontaneous flipping between pulses; a pulse that is too wide
   with too short a delay causes back-and-forth flipping within one pulse.
3. **Interface requirement**: the old waveform required PLtetO1 leak ≤1.5%; the new waveform's
   measured 0.5% satisfies it.
4. **But the non-zero RDF leak then overturns the Week 4 point** — this is why the current
   scheme (§5.3) exists. Two rounds of "failure → re-matching" show that the single-bit
   working point is not a fixed set of constants but a joint window that moves with the input
   waveform and the promoter leak.

---

## 7. The BM3R1–RDF Delay Circuit

### 7.1 What the BM3R1/RDF relationship actually is

- BM3R1's role is to **create the correct delay in RDF expression**; stronger is **not**
  better:
  - too strong → RDF unlocks too late → the reverse flip fails (LR stuck high);
  - too weak → delay too short → back-and-forth flipping within one pulse.
- RDF can interfere with the next forward flip while also being required for the reverse
  flip; it must **not** be treated as a simple "inhibitor" to be cleared.
- The delay is set by BM3R1 clearance: $\tau \approx \ln(B_{ss}/K_B)/\mu$.

### 7.2 Expression-ratio window (translated into experimentally meaningful targets)

The non-zero-leak primary scheme, converted under the model conditions
($\alpha_B=\alpha_R=120\ \mathrm{h^{-1}}$, $\gamma=4\ \mathrm{h^{-1}}$,
$\mu=0.832\ \mathrm{h^{-1}}$, $D_{tot}=0.017\ \mu M$, 1 fL assumption):

| Metric | Model target | Suggested measurable form |
|---|---:|---|
| $\alpha_B\beta_B$ | 240 h⁻² | — |
| $\alpha_R\beta_R$ | 960 h⁻² | — |
| BM3R1 maximum synthesis flux | 1.02 µM/h | PB steady state 1.226 µM ≈ **740 copies/cell** |
| RDF maximum synthesis flux | 4.08 µM/h | de-repressed steady state 4.905 µM ≈ **2950 copies/cell** |
| BM3R1 : RDF maximum capacity | 1 : 4 | compare in the same host / vector / growth condition |
| RDF relative leak | ≤0.8% | ratio of saturated-repressed to fully de-repressed fluorescence |

**Boundary statement**: $\beta$ is a translation coefficient, not a promoter strength (RPU),
and does not map to a unique RBS sequence; "15→2, 200→8" compares model working points, not
measured expression of the team's constructs. If the measured α, γ or copy number differ,
re-fit to the measured production rates rather than copying β.

### 7.3 Promoter and K/n evidence level

- The reference promoter is the characterized Cello **pBM3R1** (66 bp), see §9;
- $K_{BM3}=18.6$ nM is an **anchored estimate**: back-calculated from the measured Cello leak
  $y_{min}/y_{max}=0.008$ and the original model scale ($S/K=4.12$, with $S=0.0765\ \mu M$),
  falling in the typical in vivo half-repression range of TetR-family regulators; the Cello
  $K=0.04$ RPU cannot be converted directly to nM;
- the width of the working band in $K$ is given in §5.2/§5.3; if the construct behaves as the
  B2 configuration ($n=2.9$), the optimal $K$ must be re-selected.

---

## 8. Int Degradation Tag

### 8.1 Target and roles

**The tag is placed on Int** (all Int-containing species: free Int, dimer, Int–RDF complex,
and all DNA–Int complexes); the current scheme needs **no** degradation tag on BM3R1 or RDF.

Two roles:

1. **Suppressing inter-pulse free Int**: tag 0→4 lowers the inter-pulse free Int by ≈ 56-fold
   (0.0560 → 0.0010 µM), preventing spontaneous back-flipping after RDF unlocks and multiple
   flips — the failure mode at low tag;
2. **Clipping the pulse Int peak**: the peak falls monotonically with tag; clipping too hard
   causes the LR→PB flip to fail — the failure mode at high tag.

### 8.2 Functional window and failure modes

At RBS 0.45:

- **deterministic clean window $k_{tag,int}\in[3,24]$ h⁻¹** (dual initial states, no
  intermediate-state samples; dual-initial-state passing begins at 2.8, fragmentation from 25);
- **stochastically stable window 6–18 h⁻¹** (30 instances per level; tag 6–18 pass 30/30;
  tag 5 and 22 marginal; tag 3 and 25–26 clearly degraded);
- **RBS × tag joint workspace** is a "plateau with narrowing edges": the fully stable plateau
  tag 7–18 passes over the whole RBS 0.30–0.60 range; low-speed tags (≤4 h⁻¹) cannot be
  paired with high RBS, and high-speed tags (≥19 h⁻¹) require correspondingly higher RBS;
- **the upper edge is not a monotonic failure**: tag 25–28 is a fragmented degradation band
  (adjacent levels can jump), dominated by failed LR→PB flips and stalled cycles; at tag ≥30
  most instances are stuck in the LR high state.

### 8.3 Tag rates must be measured

Literature rates for the same ssrA tag vary widely with the measurement method:

- native ssrA: E. coli direct measurements give apparent $t_{1/2}\approx6$–14 min
  ($k\approx3$–8 h⁻¹), covering the lower edge of the window;
- high-speed variants such as LAA-LAA: ≈14–38 h⁻¹ after SI residual conversion, partially
  beyond the RBS 0.45 upper edge (can be compensated by higher RBS, but instability under
  heterogeneity appears above ≈22 h⁻¹);
- if Int fusion slows the rate to 1–2.7 h⁻¹ (pessimistic case), the system fails — changing
  the tag sequence will not help; ClpX co-expression must be considered.

**Experimental suggestion (gradient ruler)**: build sfGFP fusions with
{no tag, native ssrA, LAA+4, SsrA2X, LAA-LAA} under the same promoter/RBS; measure $k$ in
MC4100 at 37 °C in exponential phase at low induction, using both de-induction and
translation-inhibition methods; select the tag by where the measured $k$ lands.

### 8.4 Pairing with input strength

The tag rate need not be exactly 12 h⁻¹. The non-zero-leak scheme tested an interface family
that approximately preserves the effective Int amplitude:

$$ s_I(k_I) = 0.30\,\frac{\mu+k_I}{\mu+12} $$

| Measured $k_I$ / h⁻¹ | Corresponding input scale |
|---:|---:|
| 2 | 0.0662 |
| 4 | 0.1130 |
| 8 | 0.2065 |
| 12 | 0.3000 |
| 20 | 0.4870 |

This is a **pairing rule for the tested discrete points, not a continuous-interval proof**;
the tag must be matched to the measured clearance rate of the Int fusion protein, and a tag
name must not be equated with 12 h⁻¹.

---

## 9. Promoter Recommendation

### 9.1 Recommended part

The traceable reference part is the Cello **pBM3R1 (66 bp)**, 5′→3′ along transcription:

```text
AATCCGCGTGATAGGTCTGATTCGTTACCAATTGACGGAATGAACGTTCATTCCGATAATGCTAGC
```

| Reference gate | De-repressed output | Repressed output | Relative floor |
|---|---:|---:|---:|
| Cello B1_BM3R1 | 0.5 RPU | 0.004 RPU | **0.8%** |

These are measured scales in a specific Cello host/expression cassette, relative to the
J23101 reference cassette; they are **not** absolute strengths in MC4100. B1/B2/B3 share the
same promoter with different RBSs and must not be treated as three independent low-leak promoters.

### 9.2 Audit of the 36 bp variant

The 36 bp sequence provided by the wet-lab team was audited:

```text
DeepSeek original output (37 bp): TTGACACGGAATGAACGTTCATTCCGATAATGCTAGC
As transcribed (36 bp):           TTGACACGGAATGAACGTTCATCCGATAATGCTAGC
Corrected (38 bp):                TTGACACGGAATGAACGTTCATTCCGTATAATGCTAGC
```

- The source is AI-generated ("J23119 + BM3R1 operator") and is **not** a literature or
  Registry part;
- the 36 bp version breaks the 20 bp perfect-palindrome operator into a 19 bp non-palindrome;
  both 36 bp and 37 bp versions are missing one T relative to the J23119 -10 box (`TATAAT`);
- Cello/Stanton strength, leak, or $K/n$ parameters must not be applied to it.

**Recommendation: do not synthesize the 36/37 bp version as is.** Prefer the Cello pBM3R1
66 bp; if the J23119 backbone is required, use the corrected 38 bp and characterize it as a
new part.

> Directional reminder: the non-zero-leak scheme requires **moderate, calibrated** RDF
> expression; "pick the strongest promoter" comes from the old zero-leak working point and
> points the opposite way.

---

## 10. Robustness and Boundaries

### 10.1 Numerical and model boundaries

- DNA conservation: 0.01700 µM throughout (5 significant figures);
- conclusions unchanged under tightened integration tolerances and maximum step
  (see `data/nominal_comparison.csv` in the delivery package);
- all results are **deterministic mean-field**: no molecular noise, plasmid segregation, or
  single-cell success rate is simulated; the heterogeneity Monte Carlo perturbs parameter
  scenarios, it is not a Gillespie single-cell simulation.

### 10.2 Mechanistic uncertainty: η (RDF fate upon Int clearance)

The primary scheme assumes that when Int is cleared by the tag, non-covalently bound RDF is
returned completely ($\eta=1$); the Week 4 code implicitly co-degrades RDF with Int-containing
complexes ($\eta=0$). The two parameter sets must not be mixed:

- $\eta=1$: primary scheme passes;
- $\eta=0.9$: still passes ($Q\approx0.983$);
- $\eta\le0.75$: tested levels fail;
- the original $\eta=0$ treatment also has a solution: $\beta_B=1.5$, $\beta_R=20$,
  scale = 0.085, $k_I=2.2$, worst $Q=0.9879$.

**Which treatment applies must be decided by experiment** — whether Int clearance
additionally consumes RDF (same-condition RDF stability / content measurements).

### 10.3 Controlled-input comparison and pBAD-INT

For the primary B parameters, a smooth square wave (flux peak 2 µM/h, relative floor 0.5%)
was also tested: at period 10.6 h and width 2 or 3 h both initial states pass; widths of
0.25/0.5/1 h do not meet the strict endpoint standard. This is a **flux-condition comparison
for planning stepwise pBAD-INT characterization, not an operating specification for
arabinose concentration or induction time**. The induction-to-Int-flux mapping of the real
pBAD-INT construct has not been measured.

### 10.4 Not yet done / not modeled

1. Expression calibration (α, β, copy number) and the Int fusion clearance rate;
2. η determination (RDF co-loss);
3. Molecular noise and single-cell success rates (Gillespie / plasmid segregation);
4. Stochastic robustness re-test of high-tag levels (≥22) in the RBS-compensated configuration;
5. Self-consistent A+B resource-load audit (this module's results use the unloaded A input).

---

## 11. Engineering — Model-guided Decisions (summary)

| Step | Content |
|---|---|
| **Design** | A one-bit store from Int / RDF / BM3R1 plus one DNA two-state element; BM3R1 provides the delay, the RDF pool carries the memory; no logic bits or event-driven flips |
| **Build** | An explicit 38-state ODE (original recombination core + BM3R1 delay circuit + explicit leak + Int tag clearance); the real upstream flux is the only driver |
| **Test** | Strict counting criteria, dual-initial-state long runs (430 h), start-phase and parameter perturbations, heterogeneity Monte Carlo, promoter sequence audit |
| **Learn** | Memory lives in the RDF pool; BM3R1 is not "stronger is better"; the tag window is bounded on both sides and its upper edge is non-monotonic; initialization and expression calibration are the two critical experimental prerequisites |
| **Decision** | Use the 66 bp reference promoter; publish the BM3R1/RDF expression-ratio and Int-clearance pairing rules; write "trough start" into the experimental conditions; make the η test the next experiment |

> The full Design → Build → Test → Learn → Redesign narrative, including the A↔B
> feedback/rework loop, is in `03_Engineering.md`.

---

## 12. Guidance to Wet Lab

The model-to-experiment decisions of this chapter, in the order the wet-lab team should
execute them:

| # | Decision | What to do in the lab | Status |
|---|---|---|---|
| 1 | **RDF promoter** | Use the characterized Cello pBM3R1 66 bp; do **not** synthesize the 36/37 bp variant as is; if the J23119 backbone is required, use the corrected 38 bp as a new, uncharacterized part | decided from the audit |
| 2 | **Expression calibration** | Measure steady-state expression and leak of the current BM3R1 / RDF constructs in the same host/vector; convert to µM and copies/cell; match the targets (BM3R1 ≈ 740 copies/cell, RDF ≈ 2950 copies/cell, ratio ≈ 1:4) by adjusting RBS or promoter strength | first experimental priority |
| 3 | **Int degradation tag** | Build the sfGFP gradient ruler (no tag / native ssrA / LAA+4 / SsrA2X / LAA-LAA); measure $k$ in MC4100 at 37 °C in exponential phase at low induction, using both de-induction and translation-inhibition methods; select the tag by where $k$ falls in the window; if $k<3\ \mathrm{h^{-1}}$, consider ClpX co-expression | tag is required; rate must be measured |
| 4 | **RDF fate upon Int clearance (η)** | Compare RDF stability/content under matched conditions with and without Int tag clearance; if significant co-loss exists, switch to the η = 0 alternative parameter set and re-scan | resolves the largest mechanistic uncertainty |
| 5 | **Start-up arrangement** | Provide a known initial state and enter the first complete pulse from the flux trough; do not assume arbitrary-phase entry (37/40) | part of the current scheme |
| 6 | **Stepwise pBAD-INT characterization** | Use the square-wave flux comparison (peak 2 µM/h, width 2–3 h passes) to characterize the switch before coupling to the clock; model flux is **not** an arabinose concentration or induction-time specification | planned in the experimental schedule |
| 7 | **Re-calibration loop** | Feed measured values back into the parameterized scan scripts and re-check the operating window before freezing construct choices | continuous |

## 13. Scope, Limitations and Reproducibility

### In scope

- Continuous dynamics of Int / RDF / BM3R1 and the PB ⇄ LR two-state element;
- the real upstream ϕC31 production flux as the driver;
- strict counting criteria, operating window, start-up conditions and perturbations;
- Int degradation tag, BM3R1/RDF expression ratio, RDF promoter leak and sequence.

### Model boundary

- This module covers **one bit only**; multi-bit carry is covered in Chapter 4;
- output-duration programming belongs to the shutdown module; this module only registers the
  shared interface (the BM3R1 pool);
- cell-to-cell variability, molecular noise and population phase are out of scope here.

### Parameter evidence labels

| Label | Objects in this chapter |
|---|---|
| literature value | recombination core rate constants (Zhao 2019 / Pokhilko 2016); Cello B1 gate ymax/ymin |
| literature-informed range | rate ranges for native ssrA and other tags |
| model-effective parameter | $\beta_B$, $\beta_R$, $K_B$, $n$, $\gamma$, $k_I$, $\ell_R$, $\eta$ |
| design variable | Int input scale, RBS, promoter choice |
| measurement pending | expression calibration, Int fusion clearance rate, η, single-cell noise |

### Reproducibility

- Model and scan scripts ship with the delivery packages (Week 4 package + non-zero-leak package);
- computing environment: srv2026 conda `igem-tempo-2026` (Python 3.12; numpy/scipy/matplotlib);
- outputs are saved as CSV / JSON with plotting data; scans keep parameter grids and metadata;
- run levels: `smoke` for pipeline checks, `scan` for statistics, `confirmation` for reported results.

---

## 14. Chapter Conclusions

1. The single bit can act as a **digital memory element** under the real upstream waveform:
   memory resides in the RDF pool, and the Int–RDF complex sets the flip direction at pulse
   arrival. One-pulse-one-flip is achievable by design, provided the delay exceeds the pulse width.
2. The criteria for "yes" must be public and state-oriented (endpoint $Q$ + midpoint crossing +
   counting parity); the operating window is not a fixed set of constants but a joint window
   that moves with the input waveform and the promoter leak.
3. The current candidate scheme (non-zero leak ℓ_R = 0.008): $\beta_B=2$, $\beta_R=8$,
   scale = 0.30, $k_I=12\ \mathrm{h^{-1}}$; dual initial states, 430 h / 35 steady cycles pass
   with worst $Q=0.9966$; the start-up requirement is "known initial state + trough entry
   into the first complete pulse".
4. Two direct experimental interfaces are provided: the **BM3R1/RDF expression ratio
   (1:4 maximum capacity)** and the **Int tag gradient ruler (native ssrA / LAA+4 / SsrA2X /
   LAA-LAA)**; the RDF promoter recommendation is Cello pBM3R1 66 bp, and the 36 bp variant
   should not be synthesized as is.
5. Main open items: expression calibration, η (RDF co-loss), and single-cell noise. Until
   these are resolved, this is a **candidate working point under explicit model conditions,
   not an experimental success rate**.

---

## 15. Figure List and Caption Drafts

| Figure | Type | Content | Status |
|---|---|---|---|
| Figure 3-1 | Mechanism (concept) | State-dependent switch: PB/LR states, state-dependent production, RDF pool levels at pulse arrival, the two pulse-driven flip directions | sketch delivered (`figures/fig03_1_single_bit_mechanism_sketch.png`) |
| Figure 3-2 | Design (concept) | Why the window has edges: tag window with three regimes (multi-flip / pass / stuck); joint expression window with failure modes | sketch delivered (`figures/fig03_2_design_window_sketch.png`) |
| Figure 3-3 | Data | Clock input flux + steady-state trajectory (Int pulses / $S$ / RDF / BM3R1) | candidate figures available |
| Figure 3-4 | Data | Operating window triptych: (A) RBS×tag; (B) $(\beta_B,\beta_R)$ feasible region; (C) leak×$\beta_R$ | candidate figures available |
| Figure 3-5 | Data | Tag mechanism (inter-pulse / peak Int vs $k_{tag}$) + heterogeneity MC | candidate figures available |
| Figure 3-6 | Data | Start-up and perturbation pass counts + long-run comparison + input–clearance pairing | candidate figures available |
| Figure 3-7 | Table image | Parameter table: conditions, primary scheme, design targets, operating windows | delivered (`figures/fig03_7_parameter_table.png`) |

> A one-pulse-one-flip timing schematic is reserved as a companion concept figure (refinement
> pending). An optional generative concept draft (cell context) is provided as
> `figures/optional_concept_cell_context.png` for the art team.

### Caption drafts

**Figure 3-1. State-dependent switch (concept sketch).**
The DNA element adopts two configurations, PB and LR. In the PB state BM3R1 is expressed and
represses RDF production, so the RDF pool is nearly empty at pulse arrival (~0.002 µM); the
arriving clock pulse therefore drives the forward reaction (PB→LR). In the LR state BM3R1
decays and RDF accumulates (~4.5 µM at pulse arrival, ≈2500× contrast), so the next pulse
forms the Int–RDF complex and drives the reverse reaction (LR→PB). Memory resides in the RDF
pool; the delay must exceed the pulse width.

**Figure 3-2. Design window and failure modes (concept sketch).**
(A) The Int degradation-tag window has edges: too weak leaves inter-pulse Int and causes extra
flips; too strong clips the pulse peak and the reverse flip fails (stuck in LR); the
deterministic pass band is [3, 24] h⁻¹, the stochastically stable band [6, 18] h⁻¹.
(B) The joint expression window is a plateau (β_R 4–8 × β_B 0.75–6, measured under the
non-zero-leak scheme); each side fails for a distinct mechanistic reason, and a higher
promoter leak shifts the window toward lower RDF strength. Schematic; the measured pass/fail
map is the feasible-region data figure.

**Figure 3-3. Clock input and steady-state counting.**
(A) ϕC31 production flux delivered by the clock module (Week 4 frozen input, RBS scale 1);
the shaded band marks the trough used as the start-up reference.
(B) Steady-state trajectory of the primary non-zero-leak scheme ($\beta_B=2$, $\beta_R=8$,
scale 0.30): each Int pulse (top) triggers exactly one LR transition (middle), while BM3R1
and the RDF pool (bottom) set the flip direction. Model output; no experimental replicates.

**Figure 3-4. Operating windows.**
(A) RBS × tag joint workspace at the Week 4 zero-leak point; circles mark strict
dual-initial-state passes; the fully stable plateau is tag 7–18 h⁻¹.
(B) $(\beta_B,\beta_R)$ feasible region under 0.8% RDF leak (430 h, fixed trough start);
the red outline marks strict passes; the wide plateau is $\beta_R$ 4–8 × $\beta_B$ 0.75–6.
(C) Leak × $\beta_R$ tolerance: higher leak requires lower $\beta_R$.
Colour: worst dual-initial-state $Q$.

**Figure 3-5. Degradation-tag window and heterogeneity.**
(A) Inter-pulse and peak free Int versus tag rate (log scale); the grey dotted line is the
~0.1 µM flipping threshold; the green band is the RBS 0.45 working window.
(B) Heterogeneity Monte Carlo (30 instances per tag level): tag 6–18 is fully stable,
5 and 22 marginal, 3 and 25–26 clearly degraded. These are model parameter perturbations,
not an experimental noise distribution.

**Figure 3-6. Start-up, long-run behaviour and input–clearance pairing.**
(A) Strict-criteria pass counts for start-phase and parameter perturbations (fixed trough vs
random phase).
(B) Long-run comparison: the Week 4 point fails under non-zero leak, while both the
RDF-return scheme and the retuned original treatment recover counting.
(C) Input–clearance pairing rule used to match a measured Int tag rate to the input scale.

**Figure 3-7. Parameter table.** Model conditions, the primary non-zero-leak scheme, design
targets in measurable form, and operating windows; every row carries its evidence label.

> Candidate figures are in `figures/` (`candidate_*.png`); the mapping to source packages is in
> `figures/README.md`. Figure numbering overlaps with the Chapter 3 draft delivered by another
> member and must be unified when the drafts are merged (see delivery `README.md`).

---

## 16. References

1. Zhao, J., Pokhilko, A., Ebenhöh, O., Rosser, S. J., & Colloms, S. D. (2019). A single-input
   binary counting module based on serine integrase site-specific recombination.
   *Nucleic Acids Research, 47*(9), 4896–4909. DOI: 10.1093/nar/gkz245.
2. Pokhilko, A., Zhao, J., Ebenhöh, O., Colloms, S. D., et al. (2016). The mechanism of ϕC31
   integrase directionality: experimental analysis and computational modelling.
   *Nucleic Acids Research, 44*(15), 7360–7372. DOI: 10.1093/nar/gkw616.
3. Potvin-Trottier, L., Lord, N. D., Vinnicombe, G., & Paulsson, J. (2016). Synchronous
   long-term oscillations in a synthetic gene circuit. *Nature, 538*, 514–517.
   DOI: 10.1038/nature19841.
4. Stanton, B. C., Nielsen, A. A. K., Tamsir, A., Clancy, K., Peterson, T., & Voigt, C. A.
   (2014). Genomic mining of prokaryotic repressors for orthogonal logic gates.
   *Nature Chemical Biology, 10*, 99–105. DOI: 10.1038/nchembio.1411.
5. Nielsen, A. A. K., Der, B. S., Shin, J., Vaidyanathan, P., Paralanov, V., Strychalski, E. A.,
   Ross, D., Densmore, D., & Voigt, C. A. (2016). Genetic circuit design automation.
   *Science, 352*, aac7341. DOI: 10.1126/science.aac7341.
6. Liang, Q., & Fulco, A. J. (1995). Transcriptional regulation of the genes encoding
   cytochromes P450BM-1 and P450BM-3 in *Bacillus megaterium* by the binding of Bm3R1 repressor
   to Barbie box elements and operator sites. *Journal of Biological Chemistry, 270*(31),
   18606–18614. DOI: 10.1074/jbc.270.31.18606.
7. Andersen, J. B., Sternberg, C., Poulsen, L. K., Bjorn, S. P., Givskov, M., & Molin, S.
   (1998). New unstable variants of green fluorescent protein for studies of transient gene
   expression in bacteria. *Applied and Environmental Microbiology, 64*(6), 2240–2246.
   DOI: 10.1128/AEM.64.6.2240-2246.1998.
8. Lies, M., & Maurizi, M. R. (2008). Turnover of endogenous SsrA-tagged proteins mediated by
   ATP-dependent proteases in *Escherichia coli*. *Journal of Biological Chemistry, 283*(32),
   22918–22929. DOI: 10.1074/jbc.M801692200.
9. Szydlo, K., Ignatova, Z., & Gorochowski, T. E. (2022). Improving the robustness of
   engineered bacteria to nutrient stress using programmed proteolysis.
   *ACS Synthetic Biology, 11*(3), 1049–1059. DOI: 10.1021/acssynbio.1c00490.
10. Jadhav, P., Roy, S., Butzin, X. Y., & Butzin, N. C. (2025). Engineering a new SsrA-based
    degradation tag (LAA-LAA) and a bacterial synthetic oscillator.
    *ACS Synthetic Biology, 14*(4), 1062–1071. DOI: 10.1021/acssynbio.4c00612.
11. Klimecka, M. M., Antosiewicz, A., Izert, M. A., et al. (2021). A uniform benchmark for
    testing SsrA-derived degrons in *Escherichia coli*. *Molecules, 26*(19), 5936.
    DOI: 10.3390/molecules26195936.

> Reference style: full author list, year, title, journal, volume(issue), pages, DOI.
> Sources cited here are literature values or literature-informed ranges; model-effective
> parameters and pending measurements are labeled in the main text.
