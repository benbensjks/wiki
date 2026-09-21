# Gamma coupling test: 0.80x

This case changes the repressilator gamma term in the oscillator equations:

`dTetR/dt`, `dCI/dt`, and `dLacI/dt` use `gamma = gamma_base * factor`.

The generated C31 waveform is then used as the time-dependent input for the RDF counter model.

- gamma_factor: 0.8
- gamma_value_per_min: 0.018483924814931874
- period_min: 508.2097016169362
- period_generations: 16.940323387231206
- c31_peak_copies: 2532.045703548422
- c31_peak_uM: 4.206055985960834
- c31_trough_copies: 162.26809080304872
- duty_above_half_peak: 0.2797777777777778
- toggle_score: 0.15991288685033186
- H_mean: 0.8037302934591135
- L_mean: 0.36034845259867254
- alternation_fidelity: 0.25
- Int_peak_uM: 0.6042363960549271
- Int_final_uM: 0.02275443746034377
- RDF_peak_uM: 13.766931032548719
- BM3R1_peak_uM: 3.825
- sample_count: 5
