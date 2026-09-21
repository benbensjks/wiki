# eta coupling test: 1.20x

This case changes the TetR sponge accessibility parameter in the oscillator model:

`S_eff = S_phys * eta`, which changes TetR buffering and therefore PLtetO1/C31 timing.

The generated C31 waveform is used as the time-dependent input for the RDF counter model.

- eta_factor: 1.2
- eta_value: 0.6
- period_min: 403.4672445407567
- period_generations: 13.448908151358557
- c31_peak_uM: 4.050772164659662
- duty_above_half_peak: 0.2911111111111111
- toggle_score: 0.0
- H_mean: 0.6148944274372862
- L_mean: nan
- alternation_fidelity: 0.0
- Int_peak_uM: 0.5973696726191751
- Int_final_uM: 0.029270949266707823
- RDF_peak_uM: 4.021291690289487
- BM3R1_peak_uM: 3.825
- sample_count: 8
