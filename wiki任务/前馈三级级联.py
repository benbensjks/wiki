import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd

# 1. 模型参数字典：引入 A 蛋白负自馈 (Autoregulation) 参数
p = {
    # --- Pokhilko 单 Bit 拓扑与重组动力学参数 ---
    'k_int': 6.0, 
    'gamma_int': 2.0,       # 输入整合酶 Int0 参数
    'alpha_rep': 3.5, 
    'gamma_rep': 0.6,       # 通用阻遏蛋白 Rep 参数
    'alpha_rdf': 6.0, 
    'gamma_rdf': 0.8,       # 表达方向性因子 RDF 参数
    'K_rep': 0.45, 
    'n': 3.9,               # Rep 蛋白的 Hill 协同抑制
    'k_fwd': 7.0, 
    'k_rev': 7.0,           # Pokhilko 重组速率常数
    'K_D_int0': 1.0, 
    'K_D_comp': 0.8, 
    'Kinh': 0.1,            # 解离与竞争抑制常数
    
    # --- I1-FFL 0 (Bit 0 -> Bit 1 进位调控参数) ---
    'alpha_A0': 8.0,  
    'gamma_A0': 1.9,
    'alpha_R0': 4.3, 
    'gamma_R0': 0.6, 
    'K_A0': 0.5,      
    'n_A0': 2.0,
    'K_R0': 0.6,     
    'n_R0': 4.0, 
    'alpha_Int1': 18.0,
    'gamma_int1': 1.4,
    'K_D_int1': 1.0,

    # --- I1-FFL 1 (Bit 1 -> Bit 2 进位调控参数：A1 负自馈 + 超敏激活) ---
    'alpha_A1': 16.0,        # 提高初始爆发合成速率，保证脉冲峰值
    'gamma_A1': 2.5,
    'K_auto1': 0.6,         # A1 自抑制门槛 (控制 A1 稳态在低位)
    'n_auto1': 2.0,         # A1 自抑制 Hill 系数
    
    'alpha_R1': 5.0,  
    'gamma_R1': 1.1, 
    'K_A1': 1.2,            # 抬高 downstream 激活门槛 (确保稳态低于门槛)
    'n_A1': 4.0,            # 提高激活 Hill 协同阶数 (增强 Sigmoidal 效果)
    'K_R1': 0.4,      
    'n_R1': 4.0, 
    'alpha_Int2': 28.0,
    'gamma_int2': 2.0,      # 加快 Int2 消除率，防止残留
    'K_D_int2': 1.0         # 抬高 Bit 2 对 Int2 的重组门槛
}

# 2. 外源驱动脉冲 u_in(t)
pulse_starts = [1.0, 9.0, 17.0, 25.0, 33.0, 41.0, 49.0, 57.0]
pulse_duration = 0.8

def get_u_in(t):
    for t_s in pulse_starts:
        if t_s <= t <= t_s + pulse_duration:
            return 1.0
    return 0.0

# 3. 3-Bit 级联 ODE 系统
def counter_3bit_ode(t, y, p):
    pb0, int0, rep0, rdf0, a0, r0, \
    pb1, int1, rep1, rdf1, a1, r1, \
    pb2, int2, rep2, rdf2 = y

    u_in = get_u_in(t)
    lr0, lr1, lr2 = 1.0 - pb0, 1.0 - pb1, 1.0 - pb2

    # Pokhilko 重组速率计算
    def calc_rates(pb, lr, int_c, rdf_c, k_d_int):
        vf = p['k_fwd'] * pb * ((int_c**2) / (k_d_int**2 + int_c**2)) * (p['Kinh'] / (p['Kinh'] + rdf_c))
        vr = p['k_rev'] * lr * (((int_c * rdf_c)**2) / (p['K_D_comp']**2 + (int_c * rdf_c)**2))
        return vf, vr

    vf0, vr0 = calc_rates(pb0, lr0, int0, rdf0, p['K_D_int0'])
    vf1, vr1 = calc_rates(pb1, lr1, int1, rdf1, p['K_D_int1'])
    vf2, vr2 = calc_rates(pb2, lr2, int2, rdf2, p['K_D_int2'])

    # --- Bit 0 动态 ---
    dpb0 = -vf0 + vr0
    dint0 = p['k_int'] * u_in - p['gamma_int'] * int0
    drep0 = p['alpha_rep'] * pb0 - p['gamma_rep'] * rep0
    drdf0 = p['alpha_rdf'] * lr0 * (1.0 / (1.0 + (rep0 / p['K_rep'])**p['n'])) - p['gamma_rdf'] * rdf0

    # Bit 0 -> Bit 1 进位 (I1-FFL 0)
    da0 = p['alpha_A0'] * pb0 - p['gamma_A0'] * a0
    act0 = (a0**p['n_A0']) / (p['K_A0']**p['n_A0'] + a0**p['n_A0'])
    dr0 = p['alpha_R0'] * act0 - p['gamma_R0'] * r0
    rep_gate0 = (p['K_R0']**p['n_R0']) / (p['K_R0']**p['n_R0'] + r0**p['n_R0'])
    dint1 = p['alpha_Int1'] * act0 * rep_gate0 - p['gamma_int1'] * int1

    # --- Bit 1 动态 ---
    dpb1 = -vf1 + vr1
    drep1 = p['alpha_rep'] * pb1 - p['gamma_rep'] * rep1
    drdf1 = p['alpha_rdf'] * lr1 * (1.0 / (1.0 + (rep1 / p['K_rep'])**p['n'])) - p['gamma_rdf'] * rdf1

    # Bit 1 -> Bit 2 进位 (引入 A1 负自馈负反馈)
    auto_gate1 = (p['K_auto1']**p['n_auto1']) / (p['K_auto1']**p['n_auto1'] + a1**p['n_auto1'])
    #da1 = p['alpha_A1'] * pb1 * auto_gate1 - p['gamma_A1'] * a1  # A1 自抑制
    da1 = p['alpha_A1'] * pb1 - p['gamma_A1'] * a1
    act1 = (a1**p['n_A1']) / (p['K_A1']**p['n_A1'] + a1**p['n_A1'])
    dr1 = p['alpha_R1'] * act1 - p['gamma_R1'] * r1
    rep_gate1 = (p['K_R1']**p['n_R1']) / (p['K_R1']**p['n_R1'] + r1**p['n_R1'])
    pulse_gate = (int0**3.0) / (0.4**3.0 + int0**3.0)
    dint2 = p['alpha_Int2'] * act1 * rep_gate1 *pulse_gate- p['gamma_int2'] * int2

    # --- Bit 2 动态 ---
    dpb2 = -vf2 + vr2
    drep2 = p['alpha_rep'] * pb2 - p['gamma_rep'] * rep2
    drdf2 = p['alpha_rdf'] * lr2 * (1.0 / (1.0 + (rep2 / p['K_rep'])**p['n'])) - p['gamma_rdf'] * rdf2

    return [dpb0, dint0, drep0, drdf0, da0, dr0, 
            dpb1, dint1, drep1, drdf1, da1, dr1, 
            dpb2, dint2, drep2, drdf2]

# 4. 数值求解与初值计算
t_span = (0, 65)
t_eval = np.linspace(0, 65, 3000)

a0_ss = p['alpha_A0'] / p['gamma_A0']
act0_ss = (a0_ss**p['n_A0']) / (p['K_A0']**p['n_A0'] + a0_ss**p['n_A0'])
r0_ss = p['alpha_R0'] * act0_ss / p['gamma_R0']

# A1 稳态通过包含自抑制的非线性方程计算
def get_a1_ss(p):
    a1_val = 0.5
    for _ in range(100):
        auto = (p['K_auto1']**p['n_auto1']) / (p['K_auto1']**p['n_auto1'] + a1_val**p['n_auto1'])
        a1_val = p['alpha_A1'] * auto / p['gamma_A1']
    return a1_val

a1_ss = get_a1_ss(p)
act1_ss = (a1_ss**p['n_A1']) / (p['K_A1']**p['n_A1'] + a1_ss**p['n_A1'])
r1_ss = p['alpha_R1'] * act1_ss / p['gamma_R1']

rep_ss = p['alpha_rep'] / p['gamma_rep']

y0 = [1.0, 0.0, rep_ss, 0.0, a0_ss, r0_ss,  # R0 初值设为 rep_ss
      1.0, 0.0, rep_ss, 0.0, a1_ss, r1_ss,  # R1 初值设为 rep_ss (关键！压住 t=0 的 Int2)
      1.0, 0.0, rep_ss, 0.0]

sol = solve_ivp(counter_3bit_ode, t_span, y0, args=(p,), t_eval=t_eval, method='RK45')

# 5. 提取数据
LR0, LR1, LR2 = 1.0 - sol.y[0], 1.0 - sol.y[6], 1.0 - sol.y[12]
A0, R0, Int1 = sol.y[4], sol.y[5], sol.y[7]
A1, R1, Int2 = sol.y[10], sol.y[11], sol.y[13]
Int0 = sol.y[1]

# 6. 三面板绘图
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 9), sharex=True, dpi=150)

# --- Subplot 1: DNA 拓扑状态 ---
for t_s in pulse_starts:
    ax1.axvspan(t_s, t_s + pulse_duration, color='lightgray', alpha=0.5)

ax1.plot(sol.t, LR0, 'darkslategray', linewidth=2.0, label=r'Bit 0 ($LR_0$)')
ax1.plot(sol.t, LR1, 'blue', linewidth=2.0, label=r'Bit 1 ($LR_1$)')
ax1.plot(sol.t, LR2, 'darkorange', linewidth=2.0, label=r'Bit 2 ($LR_2$)')

labels = [
    "P1: (0,0,1)=1", "P2: (0,1,0)=2", "P3: (0,1,1)=3", "P4: (1,0,0)=4",
    "P5: (1,0,1)=5", "P6: (1,1,0)=6", "P7: (1,1,1)=7", "P8: (0,0,0)=0"
]
text_pos = [4.5, 12.5, 20.5, 28.5, 36.5, 44.5, 52.5, 60.5]
bbox_props = dict(boxstyle="round,pad=0.25", fc="yellow", ec="k", lw=0.8, alpha=0.8)
for pos, label in zip(text_pos, labels):
    ax1.text(pos, 0.5, label, fontsize=8, fontweight='bold', ha='center', va='center', bbox=bbox_props)

ax1.set_title("1. DNA Topology Recombination States (Pokhilko Model + Autoregulated I1-FFL)", fontsize=10, fontweight='bold')
ax1.set_ylabel("DNA LR Fraction", fontsize=8.5)
ax1.set_ylim(-0.05, 1.35)
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.legend(loc='upper right', ncol=3, fontsize=8)

# --- Subplot 2: I1-FFL 调控因子 (A & R) ---
for t_s in pulse_starts:
    ax2.axvspan(t_s, t_s + pulse_duration, color='lightgray', alpha=0.5)

ax2.plot(sol.t, A0, color='steelblue', linestyle='-', linewidth=1.5, label=r'$A_0$ (Activator 0)')
ax2.plot(sol.t, R0, color='darkorange', linestyle='--', linewidth=1.5, label=r'$R_0$ (Repressor 0)')
ax2.plot(sol.t, A1, color='saddlebrown', linestyle='-', linewidth=1.5, label=r'$A_1$ (Autoregulated Activator 1)')
ax2.plot(sol.t, R1, color='brown', linestyle='--', linewidth=1.5, label=r'$R_1$ (Repressor 1)')
ax2.set_ylim(-0.2, 8)
ax2.set_title("2. I1-FFL Regulation Factors with Negative Autoregulation on $A_1$", fontsize=10, fontweight='bold')
ax2.set_ylabel("Concentration (a.u.)", fontsize=8.5)
ax2.grid(True, linestyle='--', alpha=0.5)
ax2.legend(loc='upper right', ncol=4, fontsize=8)

# --- Subplot 3: 级联整合酶脉冲 (Int) ---
for t_s in pulse_starts:
    ax3.axvspan(t_s, t_s + pulse_duration, color='lightgray', alpha=0.5)

ax3.plot(sol.t, Int0, 'k:', linewidth=1, label=r'$Int_0$ (External Input)')
ax3.plot(sol.t, Int1, 'g-', linewidth=1, label=r'$Int_1$ (Carry Pulse to Bit 1)')
ax3.plot(sol.t, Int2, 'darkorange', linewidth=1, label=r'$Int_2$ (Carry Pulse to Bit 2)')

ax3.set_title("3. Cascaded Integrase Pulses ($Int_i$)", fontsize=10, fontweight='bold')
ax3.set_xlabel("Time (Hours)", fontsize=9)
ax3.set_ylabel("Concentration (a.u.)", fontsize=8.5)
ax3.set_xlim(0, 65)
ax3.grid(True, linestyle='--', alpha=0.5)
ax3.legend(loc='upper right', ncol=3, fontsize=8)

plt.tight_layout()
plt.show()