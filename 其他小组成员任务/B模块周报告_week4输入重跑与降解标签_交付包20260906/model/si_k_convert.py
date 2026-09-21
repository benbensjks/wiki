"""si_k_convert.py — Jadhav 2025 SI (Table S2/S3) 残留% → 有效降解速率 k 的换算表.
假定: 稳态 CFP = 产率/(k_dil + k_tag) (一阶清除); 对照 No-Tag 仅稀释.
R_tag = CFP_tag/CFP_notag → k_tag = k_dil·(1/R_tag − 1).
输出 S3 全诱导点×三标签的 k(k_dil=0.83 h⁻¹ 主表 + 0.83/1.39 敏感性),
以及 S2 高诱导档的对照值, 供报告与实验组校准.
"""
import numpy as np

KD = 0.8317766  # 50-min 倍增时的稀释率, h^-1 (A 模块 week4 口径)

# ---- Table S3 (SI, IPTG 0.1 mM, 诱导 Ara 0.1%→1%; CFP 残留%, vs No-Tag=100%) ----
S3 = {
    0.1: {'SsrA': 16.65362437, 'LAA-LAA': 4.018802, 'SsrA2X': 3.363274},
    0.2: {'SsrA': 33.51584962, 'LAA-LAA': 7.379777, 'SsrA2X': 13.30575},
    0.4: {'SsrA': 49.74359939, 'LAA-LAA': 21.89051, 'SsrA2X': 39.99999},
    0.6: {'SsrA': 58.81850522, 'LAA-LAA': 32.71032, 'SsrA2X': 56.82141},
    0.8: {'SsrA': 74.8854229, 'LAA-LAA': 41.02415, 'SsrA2X': 73.24916},
    1.0: {'SsrA': 76.47357525, 'LAA-LAA': 47.53178, 'SsrA2X': 80.74372},
}
# ---- Table S2 (SI, 1mM IPTG + 1% Ara, 高诱导矩阵) ----
S2 = {'SsrA': 37.18906, 'LAA-LAA': 13.67582, 'SsrA2X': 17.81488,
      'LAA+4': 37.38304, 'LAA': 51.30243, 'LAA-LAA-LAA': 53.41405,
      'AANDENY': 43.23287, 'AANDENY-AANDENY': 80.84964}

def k_of(Rpct, kd=KD):
    R = Rpct / 100.0
    return kd * (1.0 / R - 1.0) if R > 0 else float('inf')

print('=== Table S3: k_tag (h^-1), k_dil=0.83 h^-1 ===')
print('Ara%   SsrA    LAA-LAA SsrA2X   | 相对速率比( vs SsrA )')
for ara, vals in S3.items():
    ks = {t: k_of(v) for t, v in vals.items()}
    print(f'{ara:5.1f} {ks["SsrA"]:7.2f} {ks["LAA-LAA"]:7.2f} {ks["SsrA2X"]:7.2f} | '
          f'{ks["LAA-LAA"]/ks["SsrA"]:6.2f} {ks["SsrA2X"]/ks["SsrA"]:6.2f}')

print('\n=== Table S3 @k_dil=1.39 h^-1 (30-min Td 敏感性) ===')
for ara in (0.1, 0.2, 0.4):
    ks = {t: k_of(v, kd=1.3852) for t, v in S3[ara].items()}
    print(f'{ara:5.1f} {ks["SsrA"]:7.2f} {ks["LAA-LAA"]:7.2f} {ks["SsrA2X"]:7.2f}')

print('\n=== Table S2 (高诱导 1mM IPTG/1% Ara): k_tag = k_dil(1/R-1) ===')
for t, v in S2.items():
    print(f'{t:14s} {v:8.2f}%  ->  {k_of(v):7.2f} h^-1')

print('\n=== 跨表相对速率(速率比, 非残留比) ===')
print('S2: LAA-LAA/SsrA =', round(k_of(S2['LAA-LAA']) / k_of(S2['SsrA']), 2),
      ' SsrA2X/SsrA =', round(k_of(S2['SsrA2X']) / k_of(S2['SsrA']), 2))
print('S3@0.1%: LAA-LAA/SsrA =', round(k_of(S3[0.1]['LAA-LAA']) / k_of(S3[0.1]['SsrA']), 2),
      ' SsrA2X/SsrA =', round(k_of(S3[0.1]['SsrA2X']) / k_of(S3[0.1]['SsrA']), 2))
print('\n=== 若以文献直测 SsrA k=3-8 h^-1 锚定 S3@0.1% 速率比, LAA-LAA 推算区间 ===')
ratio = k_of(S3[0.1]['LAA-LAA']) / k_of(S3[0.1]['SsrA'])
print('ratio@0.1%% = %.2f -> LAA-LAA 区间 = %.1f-%.1f h^-1' % (ratio, 3*ratio, 8*ratio))
