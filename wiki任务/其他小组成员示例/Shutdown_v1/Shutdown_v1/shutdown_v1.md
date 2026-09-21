# BM3R1-Driven sRNA Shutdown Module

> Draft for the Model Wiki. The wording, figure layout, and numerical presentation can be refined during final webpage design.

## 1. Programming the Output Duration

The TEMPO counter uses stable DNA states to record the number of input cycles. After the counter switches from the PB state to the LR state, the LR state remains stable and continues to drive the downstream output. This is useful for storing the counting result, but it also creates a problem: the output may remain active long after the required response period.

We therefore separated two timing functions in the system:

- the counter determines **when output begins**;
- the shutdown module determines **how long output lasts**.

The design goal was to create a finite expression window without disrupting the PB–LR state transition or erasing the stored counting result.

## 2. Design Concept: Delayed sRNA Shutdown

Our design uses residual BM3R1 as a molecular timer. Before the PB-to-LR transition, BM3R1 represses the shutdown promoter. After the transition, new BM3R1 production stops, but the existing BM3R1 pool does not disappear immediately. Its gradual clearance creates a delay between counter switching and shutdown activation.

When BM3R1 falls below the effective repression threshold, the shutdown sRNA begins to accumulate. The sRNA binds the output mRNA and promotes its removal, reducing new protein synthesis even though the LR DNA state remains unchanged.

The expected sequence is therefore:

**PB-to-LR switching → BM3R1 decline → delayed sRNA activation → output shutdown**

> **Figure 1 placement — immediately after this paragraph**  
> **Figure:** *Delayed sRNA Shutdown Mechanism*  
> **Type:** mechanism diagram to be redrawn by the art team.  
> **Content:** three stages arranged from left to right: (1) BM3R1 represses sRNA transcription before switching; (2) residual BM3R1 gradually declines after PB-to-LR switching; (3) sRNA accumulates and suppresses output mRNA.  
> **Existing reference for the art team:** `wiki_report/assets/00_shutdown_mechanism.png`.

## 3. Mathematical Model

We translated the proposed circuit into a mechanistic ODE model. The model describes free BM3R1, shutdown sRNA, output mRNA, and output protein. The counter trajectory provides the upstream PB/LR state, while BM3R1 links the counter to the shutdown layer.

A simplified form of the shutdown equations is:

$$
\frac{dS}{dt}=\alpha_S f_{\mathrm{rep}}(B)-\delta_S S-kSM
$$

$$
\frac{dM}{dt}=\alpha_M L(t)-\delta_M M-kSM
$$

$$
\frac{dP}{dt}=\beta M-\delta_P P
$$

where $B$ is free BM3R1, $S$ is the shutdown sRNA, $M$ is output mRNA, $P$ is output protein, and $L(t)$ represents activation by the LR state. BM3R1-dependent repression is described using a Hill function. Growth dilution is included in the effective loss rates.

The model uses RyhB-like sRNA regulation as a coarse-grained reference. It does not attempt to reproduce every molecular step of sRNA–mRNA interaction.

## 4. Baseline Dynamics

We first tested whether the model could generate the intended temporal sequence after a PB-to-LR transition. Under the nominal 50-minute doubling condition, BM3R1 decreased after switching, sRNA accumulated after a delay, and output translation first increased and then declined.

The shutdown model produced a new-translation expression window with a full width at half maximum of **6.98 h**, compared with **10.61 h** in the matched no-shutdown control. At the late LR stage, translation remained at **17.0%** of its peak value with shutdown, compared with **97.9%** without shutdown. The peak translation rate was still **98.6%** of the no-shutdown control, indicating that the module shortened the output period without preventing initial activation.

These results support the main design idea: a stable LR state can preserve the counting result while the downstream regulatory layer limits expression duration.

> **Figure 2 placement — directly below the baseline-results paragraph**  
> **Figure:** *Baseline Dynamics of the Shutdown Module*  
> **Use this existing result:** `wiki_report/assets/01_baseline_timecourse.png`.  
> **What it should show:** BM3R1 decline, delayed sRNA accumulation, and output activation followed by shutdown. The final Wiki caption should point out the 6.98 h translation FWHM and compare the shutdown trajectory with the no-shutdown control. The art team may adjust fonts, labels, and colors, but the plotted data should remain unchanged.

## 5. Tuning the Expression Window

We next varied shutdown-related parameters to examine whether the expression window occurred only at one parameter set. In the two-dimensional scan of sRNA production strength and the relative BM3R1 repression threshold, **91 of 169** parameter combinations met the predefined shutdown criteria. Their translation-window widths ranged from **5.42 to 7.36 h**.

The one-at-a-time analysis suggested that the sRNA maximum transcription rate had the largest effect on window width in the tested range. The relative promoter threshold and the copy-number or transcription-supply ratio also affected the timing and completeness of shutdown.

These results suggest that the duration can be adjusted through the regulatory layer. However, they do not show that every parameter combination produces successful shutdown. In the current model, sufficient sRNA supply relative to output-mRNA production is required.

> **Figure 3 placement — after the tunability discussion**  
> **Figure 3A:** *Shutdown Tunability Map* — use `wiki_report/assets/02_tunability_heatmap.png`. This is the main evidence for the 91/169 passing combinations and the 5.42–7.36 h window range.  
> **Figure 3B:** *One-at-a-Time Parameter Sensitivity* — use `wiki_report/assets/06_shutdown_oat_sensitivity.png`. This panel identifies which tested parameters have the largest influence on expression-window width.  
> These two existing plots can be placed side by side as one combined Wiki figure.

## 6. Robustness Tests

We performed a joint parameter perturbation analysis to test whether finite output was retained outside the nominal parameter set. Of 512 Latin-hypercube samples, **277 samples (54.1%)** passed the criteria for output activation, finite pulse duration, late shutdown, and peak preservation. The median window width among the passing samples was **6.86 h**.

> **Figure 4 placement — after the joint-perturbation paragraph**  
> **Figure:** *Joint Parameter Robustness Analysis*  
> **Use this existing result:** `wiki_report/assets/07_shutdown_joint_uncertainty.png`.  
> The caption should state that 277/512 samples passed the predefined criteria and that this represents a feasible parameter region rather than a predicted experimental success rate.

We also tested growth-dependent timing. At doubling times of 40, 50, and 60 minutes, the predicted translation-window widths were **5.48, 6.98, and 8.54 h**, respectively. All three conditions retained the required order of counter switching, delayed sRNA activation, and shutdown, although the exact timing changed.

> **Figure 5 placement — after the growth-effects paragraph**  
> **Figure:** *Growth-Dependent Shutdown Dynamics*  
> **Use this existing result:** `wiki_report/assets/05_self_consistent_growth.png`.  
> The figure should be used to compare the 40-, 50-, and 60-minute doubling conditions and show that growth changes the predicted window length without removing the qualitative shutdown sequence.

The robustness results therefore indicate the presence of a feasible operating region rather than universal robustness under all possible parameter changes.

## 7. Design Development

The shutdown architecture was selected after comparing several possible ways to limit persistent output. Fixed counter-state windows would couple output duration directly to the timing of later input cycles. A second counter or an additional DNA-state mechanism would increase circuit complexity. Protein degradation tags could accelerate output removal, but they would not stop continued production from a stable active state. RNA-level regulation offered a more direct way to stop new expression while leaving the counter state intact.

This comparison led us to reuse a molecule already present in the counter. Residual BM3R1 provides the delay, while the sRNA supplies the shutdown action. The resulting design assigns DNA state and expression duration to different regulatory layers.

These alternatives were evaluated conceptually during early design. We did not construct or simulate a recombination-intermediate shutdown mechanism, so it should not be described as an experimentally failed design.

> **Figure 6 placement — at the end of the design-development section**  
> **Figure:** *Development of the Shutdown Design*  
> **Type:** flowchart to be drawn by the art team.  
> **Suggested sequence:** persistent output from the stable LR state → comparison of fixed state windows, a second counter, STOP elements, degradation tags, and RNA-level regulation → limitations identified → residual BM3R1 selected as the delay signal → BM3R1-gated sRNA selected as the shutdown layer.  
> This is a design-history figure. It should not present the recombination-intermediate route as a constructed or experimentally tested failure.

## 8. Current Limitations

The present model is deterministic and does not include cell-to-cell variability. Effective degradation terms combine molecular degradation with growth dilution, and the sRNA interaction is represented by a simplified mass-action model. Several parameters are based on literature-scale estimates rather than time-resolved measurements from the final construct.

Therefore, the model currently supports the feasibility and design logic of delayed shutdown. Experimental measurements will still be required to determine the real BM3R1 clearance rate, sRNA expression strength, repression threshold, and output decay dynamics.

The main conclusion at this stage is that the counter and shutdown module can perform different temporal functions: **the counter determines when expression begins, whereas the BM3R1–sRNA module determines how long it persists.**
