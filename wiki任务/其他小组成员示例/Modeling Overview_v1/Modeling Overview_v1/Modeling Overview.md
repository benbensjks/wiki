# Modeling Overview

## 1. Why Modeling Matters for TEMPO

TEMPO is designed to transform biological oscillations into programmable cellular decisions.

At the system level, its logic appears simple:

**generate a periodic signal → count each pulse → propagate the count → trigger an output → shut the output down**

大概这个位置，放 overview 的流程图


However, none of these steps is determined by circuit topology alone.

Whether TEMPO functions correctly depends on dynamic quantities such as pulse amplitude, pulse width, protein accumulation and clearance, DNA-state transitions, regulatory thresholds, and the timing relationships between modules.

For example, an oscillator may generate a stable waveform but still fail to drive the counter correctly. A recombinase switch may flip successfully once, yet fail when the next pulse arrives before RDF or Integrase has sufficiently cleared. Similarly, a stable DNA state may determine **when** downstream expression begins, but does not automatically determine **how long** that expression should persist.

Therefore, mathematical modeling was not used as an isolated descriptive step in TEMPO. It became an engineering tool for asking whether individually plausible biological modules could operate together as a coherent timing system.

Our modeling work was built around five central questions:

* Can the oscillator generate a stable and tunable timing signal?
* Can each oscillator cycle be translated into exactly one counter transition?
* What regulatory conditions allow the recombinase switch to operate reliably?
* How can a single-bit switch be extended into a multi-bit binary counter?
* Once the desired count is reached, how can downstream expression be limited to a finite time window?

These questions define the structure of our modeling framework.

---

## 2. From Biological Functions to Modeling Questions

TEMPO processes temporal information through a sequence of biological operations.

Although these operations form one continuous system, each introduces a different dynamic constraint. We therefore divided the overall problem into several modeling questions rather than constructing a single monolithic model.

| Biological function          | Modeling question                                                                      | Why it matters                                                                                                      |
| ---------------------------- | -------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| Oscillator                   | What determines the period, amplitude, and pulse width of the clock signal?            | The downstream counter can only operate reliably if the timing signal itself is stable and tunable.                 |
| Oscillator–counter interface | What signal is actually transmitted from the oscillator to the recombinase system?     | A biologically generated Integrase pulse may behave very differently from an idealized square-wave input.           |
| Single-bit counter           | Under what conditions does one pulse produce exactly one PB↔LR transition?             | Reliable counting requires each oscillatory event to correspond to one and only one state transition.               |
| Regulatory design            | How should BM3R1, RDF, and protein degradation be configured?                          | Their relative dynamics determine whether switching, recovery, and repression occur within the correct time window. |
| Multi-bit counter            | How can a completed transition in one bit generate a carry signal for the next bit?    | Scaling from one bit to binary counting requires reliable information transfer between stages.                      |
| Timed output                 | How can a stable counter state generate only a finite period of downstream expression? | Reaching the target state determines when expression begins, but not automatically when it should end.              |

These questions define the structure of the modeling work described in the following sections.

Each model focuses on a specific biological function, but none is treated as an isolated system. Instead, the output or constraint identified in one module is passed forward to the next stage of the TEMPO architecture.

---

## 3. Connecting Models to Engineering Decisions

The central purpose of our modeling framework was not only to describe each TEMPO module independently, but to determine whether the modules could operate together under compatible dynamic conditions.

The models were therefore linked through explicit biological interfaces.

The oscillator model provides the time-dependent signal delivered to the recombinase system. The single-bit counter model tests whether this input can produce reliable PB↔LR transitions and whether Integrase, RDF, and regulatory proteins recover before the following pulse. The resulting timing constraints then become requirements for multi-bit carry propagation. Finally, the downstream model asks whether the counter state can be converted into a controlled expression window rather than persistent output.

In this way, information flows through the modeling framework in the same direction as it flows through the biological system:

**oscillator dynamics → Integrase input → DNA-state switching → carry propagation → timed downstream expression**

This connection between modules was important because a design that performs well in isolation does not necessarily remain functional when coupled to the rest of the system. Modeling therefore allowed us to identify not only feasible operating regions, but also incompatibilities that were not obvious from circuit diagrams alone.

Throughout the project, simulation results were repeatedly translated into engineering consequences.

For example, our models were used to:

* identify parameters capable of tuning oscillator dynamics;
* replace idealized counter inputs with biologically generated Integrase waveforms;
* determine conditions compatible with one-pulse-one-flip behavior;
* constrain BM3R1, RDF, and protein-clearance requirements;
* reveal timing and leakage limitations during multi-bit scaling;
* and identify persistent downstream expression as a separate design problem, motivating the introduction of a delayed sRNA-based shutdown mechanism.

Importantly, not every simulation confirmed the original design.

Some results defined a successful operating region, whereas others exposed failure modes or showed that an initial assumption was too optimistic. We treated these outcomes as part of the engineering process: when a model revealed a limitation, the result was used to refine the mechanism, redefine the feasible design space, or motivate an alternative architecture.

**Biological question → model → simulation → insight or limitation → design decision**

This workflow connects the individual modeling sections into a single engineering framework and provides the basis for the design choices discussed throughout the rest of the page.

---

## 4. Modeling as an Engineering Cycle

Our models were not developed only to reproduce expected biological behavior.

Throughout the project, they were repeatedly used to test assumptions, expose failure modes, and guide redesign.

The recurring workflow was:

**Design question → mechanistic model → simulation → failure or constraint → interpretation → design revision**

This process contributed to several engineering decisions, including:

* identifying parameters capable of tuning oscillator dynamics;
* replacing an overly idealized oscillator–counter connection with a biologically meaningful time-dependent input;
* defining conditions for reliable one-pulse-one-flip behavior;
* constraining suitable BM3R1 and RDF expression regimes;
* evaluating whether protein degradation is required for reliable recovery;
* comparing alternative carry mechanisms for multi-bit counting;
* identifying leakage and timing limitations that emerge during cascade scaling;
* recognizing that the original counter-output architecture could lead to persistent downstream expression;
* redesigning the output module around delayed sRNA-mediated shutdown;
* and testing how robust these conclusions remain under parameter variation.

Not every simulation produced an ideal result.

Some models confirmed that a proposed mechanism was feasible. Others revealed that a design was fragile, poorly scalable, or incompatible with the timing constraints imposed by upstream modules.

We consider both outcomes important.

A model that identifies **why a design fails** can provide as much engineering value as one that successfully reproduces the desired behavior, because it defines the boundary between plausible and implausible designs.

For this reason, the following sections present not only our successful simulations, but also the limitations, failed assumptions, and redesigns that shaped the final TEMPO architecture.

---

## 5. From Models to an Interactive Design Tool

As our modeling framework expanded, another problem became apparent.

The information needed to design TEMPO was distributed across multiple simulations:

* oscillator parameters determined pulse timing;
* Integrase and RDF dynamics determined whether counting succeeded;
* carry dynamics determined whether multi-bit propagation remained reliable;
* and sRNA parameters determined whether downstream output closed within the desired time window.

Inspecting these modules through separate scripts made it difficult to evaluate a complete design under one consistent set of parameters.

We therefore developed **TEMPO Design Explorer**, an interactive interface that connects our existing models into a common design workflow.

Instead of asking users to directly modify abstract ODE parameters, the Explorer exposes biologically meaningful design controls and allows users to inspect their effects across the coupled system.

A typical workflow is:

**select a timing design
→ adjust biological parameters
→ simulate the coupled dynamics
→ inspect oscillator, counter, and shutdown behavior
→ identify the failed requirement
→ modify the design**

The platform does not introduce a separate biological model.

Rather, it transforms the models described on this page into a reusable design and decision-support tool.

Through the Explorer, users can examine how changes in parameters such as Integrase input, degradation, RDF/BM3R1 regulation, and sRNA production affect:

* pulse timing,
* counter transitions,
* recovery before subsequent pulses,
* flip fidelity,
* residual downstream expression,
* and overall design feasibility.

This extends the role of modeling beyond our own project.

Instead of providing only static figures describing one set of simulations, TEMPO Design Explorer allows future teams to interact with the modeling framework and explore alternative timing designs for themselves.

**→ Explore TEMPO Design Explorer on our Software page.**

---

## 6. How to Read the Modeling Section

The following sections follow the same logic as the biological information flow through TEMPO:

**Oscillator**
↓
**Oscillator–Counter Coupling**
↓
**Single-Bit Recombinase Switch**
↓
**Multi-Bit Carry Architectures**
↓
**Timed Shutdown**
↓
**Sensitivity, Validation, and Design Guidance**

For each model, we focus on four questions:

1. **What engineering problem were we trying to solve?**
2. **How did we represent the biology mathematically?**
3. **What did the simulations reveal?**
4. **How did those results affect the design of TEMPO?**

Our goal is therefore not only to show mathematical descriptions of the system, but to demonstrate how modeling helped transform TEMPO from a conceptual circuit architecture into a quantitatively testable engineering framework.
