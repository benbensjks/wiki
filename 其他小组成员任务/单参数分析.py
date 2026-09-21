import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd

# =====================================================================
# 1. 定义基础参数字典 p (全局基准参数)
# =====================================================================
p = {
    # Pokhilko 单 Bit 拓扑与重组动力学参数
    'k_int': 6.0, 
    'gamma_int': 2.0,       
    'alpha_rep': 3.2, 
    'gamma_rep': 0.7,       
    'alpha_rdf': 6.0, 
    'gamma_rdf': 0.8,       
    'K_rep': 0.85,          # 抬高 K_rep 防止 R0 过度抑制 RDF0
    'n': 3.9,               
    'k_fwd': 7.0, 
    'k_rev': 5.0,           # 降低逆向重组速率常数，防止暂态回退
    'K_D_int0': 1.0, 
    'K_D_comp': 1.2, 
    'Kinh': 0.1,            
    
    # I1-FFL 0 (Bit 0 -> Bit 1 进位调控参数)
    'alpha_A0': 8.0,  
    'gamma_A0': 1.9,
    'alpha_R0': 3.8,        
    'gamma_R0': 0.6, 
    'K_A0': 0.5,      
    'n_A0': 2.0,
    'K_R0': 0.6,     
    'n_R0': 4.0, 
    'alpha_Int1': 18.0,
    'gamma_int1': 1.4,
    'K_D_int1': 1.0,

    # I1-FFL 1 (Bit 1 -> Bit 2 进位调控参数：A1 负自馈 + 双因子门控)
    'alpha_A1': 16.0, 
    'gamma_A1': 2.5,
    'K_auto1': 0.6,         
    'n_auto1': 2.0,         
    'alpha_R1': 5.0,  
    'gamma_R1': 1.1, 
    'K_A1': 1.2,            
    'n_A1': 4.0,            
    'K_R1': 0.4,      
    'n_R1': 4.0, 
    'alpha_Int2': 38.0,     
    'gamma_int2': 2.2,      
    'K_D_int2': 0.6         
}

# 脉冲时间定义
pulse_starts = [1.0, 9.0, 17.0, 25.0, 33.0, 41.0, 49.0, 57.0]
pulse_duration = 0.8

def get_u_in(t):
    for t_s in pulse_starts:
        if t_s <= t <= t_s + pulse_duration:
            return 1.0
    return 0.0

# =====================================================================
# 2. 定义 ODE 方程系统
# =====================================================================
def counter_3bit_ode(t, y, p_in):
    pb0, int0, rep0, rdf0, a0, r0, \
    pb1, int1, rep1, rdf1, a1, r1, \
    pb2, int2, rep2, rdf2 = y

    u_in = get_u_in(t)
    lr0, lr1, lr2 = 1.0 - pb0, 1.0 - pb1, 1.0 - pb2

    def calc_rates(pb, lr, int_c, rdf_c, k_d_int):
        vf = p_in['k_fwd'] * pb * ((int_c**2) / (k_d_int**2 + int_c**2)) * (p_in['Kinh'] / (p_in['Kinh'] + rdf_c))
        vr = p_in['k_rev'] * lr * (((int_c * rdf_c)**2) / (p_in['K_D_comp']**2 + (int_c * rdf_c)**2))
        return vf, vr

    vf0, vr0 = calc_rates(pb0, lr0, int0, rdf0, p_in['K_D_int0'])
    vf1, vr1 = calc_rates(pb1, lr1, int1, rdf1, p_in['K_D_int1'])
    vf2, vr2 = calc_rates(pb2, lr2, int2, rdf2, p_in['K_D_int2'])

    dpb0 = -vf0 + vr0
    dint0 = p_in['k_int'] * u_in - p_in['gamma_int'] * int0
    drep0 = p_in['alpha_rep'] * pb0 - p_in['gamma_rep'] * rep0
    drdf0 = p_in['alpha_rdf'] * lr0 * (1.0 / (1.0 + (rep0 / p_in['K_rep'])**p_in['n'])) - p_in['gamma_rdf'] * rdf0

    da0 = p_in['alpha_A0'] * pb0 - p_in['gamma_A0'] * a0
    act0 = (a0**p_in['n_A0']) / (p_in['K_A0']**p_in['n_A0'] + a0**p_in['n_A0'])
    dr0 = p_in['alpha_R0'] * act0 - p_in['gamma_R0'] * r0
    rep_gate0 = (p_in['K_R0']**p_in['n_R0']) / (p_in['K_R0']**p_in['n_R0'] + r0**p_in['n_R0'])
    dint1 = p_in['alpha_Int1'] * act0 * rep_gate0 - p_in['gamma_int1'] * int1

    dpb1 = -vf1 + vr1
    drep1 = p_in['alpha_rep'] * pb1 - p_in['gamma_rep'] * rep1
    drdf1 = p_in['alpha_rdf'] * lr1 * (1.0 / (1.0 + (rep1 / p_in['K_rep'])**p_in['n'])) - p_in['gamma_rdf'] * rdf1

    auto_gate1 = (p_in['K_auto1']**p_in['n_auto1']) / (p_in['K_auto1']**p_in['n_auto1'] + a1**p_in['n_auto1'])
    da1 = p_in['alpha_A1'] * pb1 * auto_gate1 - p_in['gamma_A1'] * a1 
    
    act1 = (a1**p_in['n_A1']) / (p_in['K_A1']**p_in['n_A1'] + a1**p_in['n_A1'])
    dr1 = p_in['alpha_R1'] * act1 - p_in['gamma_R1'] * r1
    rep_gate1 = (p_in['K_R1']**p_in['n_R1']) / (p_in['K_R1']**p_in['n_R1'] + r1**p_in['n_R1'])
    pulse_gate = (int0**3.0) / (0.4**3.0 + int0**3.0)
    dint2 = p_in['alpha_Int2'] * act1 * rep_gate1 * pulse_gate - p_in['gamma_int2'] * int2

    dpb2 = -vf2 + vr2
    drep2 = p_in['alpha_rep'] * pb2 - p_in['gamma_rep'] * rep2
    drdf2 = p_in['alpha_rdf'] * lr2 * (1.0 / (1.0 + (rep2 / p_in['K_rep'])**p_in['n'])) - p_in['gamma_rdf'] * rdf2

    return [dpb0, dint0, drep0, drdf0, da0, dr0, 
            dpb1, dint1, drep1, drdf1, da1, dr1, 
            dpb2, dint2, drep2, drdf2]

# =====================================================================
# 2. 核心评估函数 (计算系统的综合表现得分 Score)
# =====================================================================
def evaluate_circuit_performance(p_test):
    t_span = (0, 35)
    t_eval = np.linspace(0, 35, 1000)

    rep_ss = p_test['alpha_rep'] / p_test['gamma_rep']
    a0_ss = p_test['alpha_A0'] / p_test['gamma_A0']
    
    a1_val = 0.5
    for _ in range(50):
        auto = (p_test['K_auto1']**p_test['n_auto1']) / (p_test['K_auto1']**p_test['n_auto1'] + a1_val**p_test['n_auto1'])
        a1_val = p_test['alpha_A1'] * auto / p_test['gamma_A1']

    y0 = [1.0, 0.0, rep_ss, 0.0, a0_ss, rep_ss,  
          1.0, 0.0, rep_ss, 0.0, a1_val, rep_ss,  
          1.0, 0.0, rep_ss, 0.0]

    try:
        sol = solve_ivp(counter_3bit_ode, t_span, y0, args=(p_test,), t_eval=t_eval, method='RK45')
        
        LR1 = 1.0 - sol.y[6]
        LR2 = 1.0 - sol.y[12]
        
        # 1. P3 阶段 Bit 1 的维持表现 (期望 >= 0.9)
        mask_p3 = (sol.t >= 17.5) & (sol.t <= 24.5)
        min_LR1_P3 = np.min(LR1[mask_p3]) if np.any(mask_p3) else 0.0
        
        # 2. P2/P3 静态期 Bit 2 的泄漏情况 (期望 <= 0.05)
        mask_leak = (sol.t >= 9.0) & (sol.t <= 24.5)
        max_leak_P23 = np.max(LR2[mask_leak]) if np.any(mask_leak) else 1.0
        
        # 3. P4 阶段 Bit 2 的翻转深度 (期望 >= 0.9)
        mask_p4 = (sol.t >= 25.5) & (sol.t <= 32.5)
        max_LR2_P4 = np.max(LR2[mask_p4]) if np.any(mask_p4) else 0.0

        # 综合得分惩罚机制：翻转度 - 泄漏量 - Bit1塌陷惩罚
        score = max_LR2_P4 - max_leak_P23 - (1.0 - min_LR1_P3)
        return max(0.0, score)
    except:
        return 0.0

# =====================================================================
# 3. 全参数扫描：计算每个参数的“敏感度指数”
# =====================================================================
print("正在扫描所有参数的影响程度，请稍候...")
base_score = evaluate_circuit_performance(p)
sensitivity_results = {}

PERTURBATION = 0.50  # 对每个参数施加 +/- 20% 的微扰

for param_name, base_val in p.items():
    # 上浮 20%
    p_up = p.copy()
    p_up[param_name] = base_val * (1 + PERTURBATION)
    score_up = evaluate_circuit_performance(p_up)
    
    # 下浮 20%
    p_down = p.copy()
    p_down[param_name] = base_val * (1 - PERTURBATION)
    score_down = evaluate_circuit_performance(p_down)
    
    # 相对得分波动绝对值作为“敏感度”
    sensitivity = (abs(score_up - base_score) + abs(score_down - base_score)) / 2.0
    sensitivity_results[param_name] = sensitivity

# 转为 DataFrame 排序
df_sens = pd.DataFrame(list(sensitivity_results.items()), columns=['Parameter', 'Sensitivity'])
df_sens = df_sens.sort_values(by='Sensitivity', ascending=False).reset_index(drop=True)
df_sens_active = df_sens[df_sens['Sensitivity'] > 0.001]

# 打印排名 Top 10 关键参数
print("\n=== 对电路影响最大的前 10 个关键参数 ===")
print(df_sens.head(10).to_string(index=False))

# =====================================================================
# 4. 绘图：参数敏感度排行榜 (Tornado Plot)
# =====================================================================
plt.figure(figsize=(9, 6), dpi=150)
top_n = 9
plt.barh(df_sens['Parameter'][:top_n][::-1], df_sens['Sensitivity'][:top_n][::-1], color='crimson', alpha=0.8)
plt.xlabel('Sensitivity Index (Impact on Circuit Performance)', fontsize=10, fontweight='bold')
plt.title(f'Global Parameter Sensitivity Ranking (Top {top_n})', fontsize=11, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()