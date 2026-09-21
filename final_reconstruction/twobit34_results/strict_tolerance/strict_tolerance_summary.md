# 中心与边界严格容差复核

- 轨迹：27/27
- 积分失败：0
- 全部点结论收敛：**True**

| 点 | 类型 | uM/mat/mRNA | ultra结论 | 收敛 | baseline相对ultra最大数值差 |
|---|---|---|---|:---:|---:|
| selected_centre | centre | 5.75/32.5/2 | PASS | True | 9.73e-10 |
| pass_low_uM_edge | pass_boundary | 5.5/30/2 | PASS | True | 9.74e-10 |
| pass_low_mat_edge | pass_boundary | 6/25/2 | PASS | True | 1.01e-09 |
| pass_notch_lip | pass_boundary | 6.25/35/2 | PASS | True | 3.58e-08 |
| pass_high_corner | pass_boundary | 6.5/25/4 | PASS | True | 9.74e-10 |
| fail_low_uM | fail_boundary | 5.5/27.5/2 | readout_unlabelled_window | True | 9.87e-10 |
| fail_low_mat | fail_boundary | 5.75/25/2 | readout_unlabelled_window | True | 9.74e-10 |
| fail_notch_mid | fail_boundary | 6.25/32.5/2 | readout_unlabelled_window | True | 2.34e-08 |
| fail_notch_low | fail_boundary | 6.25/30/2 | readout_unlabelled_window | True | 1.05e-08 |
