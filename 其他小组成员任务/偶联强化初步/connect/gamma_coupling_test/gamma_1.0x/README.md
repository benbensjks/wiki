# Gamma coupling test: 1.00x

This case changes the repressilator gamma term in the oscillator equations:

`dTetR/dt`, `dCI/dt`, and `dLacI/dt` use `gamma = gamma_base * factor`.

The generated C31 waveform is then used as the time-dependent input for the RDF counter model.

- gamma_factor: 1.0
- gamma_value_per_min: 0.023104906018664842
- period_min: 406.5677612935489
- period_generations: 13.552258709784963
- c31_peak_copies: 2404.228483379289
- c31_peak_uM: 3.9937350222247328
- c31_trough_copies: 170.2836817011985
- duty_above_half_peak: 0.2853333333333333
- toggle_score: 0.8542011391556787
- H_mean: 0.8817441462882764
- L_mean: 0.14579886084432123
- alternation_fidelity: 1.0
- Int_peak_uM: 0.5913740231440374
- Int_final_uM: 0.05399442262560622
- RDF_peak_uM: 5.0837713951475765
- BM3R1_peak_uM: 3.825
- sample_count: 8
