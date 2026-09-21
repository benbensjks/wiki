# Gamma coupling test: 1.20x

This case changes the repressilator gamma term in the oscillator equations:

`dTetR/dt`, `dCI/dt`, and `dLacI/dt` use `gamma = gamma_base * factor`.

The generated C31 waveform is then used as the time-dependent input for the RDF counter model.

- gamma_factor: 1.2
- gamma_value_per_min: 0.02772588722239781
- period_min: 338.7487658199443
- period_generations: 11.291625527331478
- c31_peak_copies: 2280.424908835341
- c31_peak_uM: 3.788081243912526
- c31_trough_copies: 182.44333125413857
- duty_above_half_peak: 0.314
- toggle_score: 0.7532513152514871
- H_mean: 0.960271813710204
- L_mean: 0.24674868474851283
- alternation_fidelity: 1.0
- Int_peak_uM: 0.5793883317229942
- Int_final_uM: 0.09747404924546692
- RDF_peak_uM: 34.789881317719264
- BM3R1_peak_uM: 3.825
- sample_count: 10
