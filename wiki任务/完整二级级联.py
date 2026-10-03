import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# =====================================================
# 1. 外部输入信号 get_u_in(t) (Repressilator 驱动 Bit 0)
# =====================================================
class P:
    lambd, K, n, Td = 1000, 13, 3, 30
    gamma = math.log(2) / Td
    S_eff, Kd, n_s = 40, 6, 1.5
    K_R, n_R, H_leak = 6, 2, 0.05
    N_C, k_tx, gamma_m = 10, 5, math.log(2) / 2
    beta, mu = 0.5, math.log(2) / 30

def hill_rep(x): return (P.K**P.n) / (P.K**P.n + x**P.n + 1e-12)
def PLtetO(Tf): return P.H_leak + (1 - P.H_leak) / (1 + (Tf / P.K_R)**P.n_R)

_Ttot_grid = np.linspace(0, 3000, 3000)
def _exact_free(Ttot, S):
    if Ttot <= 0: return 0.0
    low, high = 0.0, float(Ttot)
    while high - low > 1e-5:
        mid = (low + high) / 2.0
        bound = S * (mid**P.n_s) / (mid**P.n_s + P.Kd**P.n_s)
        if (mid + bound) > Ttot: high = mid
        else: low = mid
    return (low + high) / 2.0
_Tfree_grid = np.array([_exact_free(x, P.S_eff) for x in _Ttot_grid])
def free_TetR(T): return float(np.interp(T, _Ttot_grid, _Tfree_grid))

def repressilator_ode(t, y):
    TetR, CI, LacI, mRNA, C31 = y
    Tf = free_TetR(TetR)
    return [
        P.gamma * (P.lambd * hill_rep(LacI) - TetR),
        P.gamma * (P.lambd * hill_rep(Tf) - CI),
        P.gamma * (P.lambd * hill_rep(CI) - LacI),
        P.N_C * P.k_tx * PLtetO(Tf) - P.gamma_m * mRNA,
        P.beta * mRNA - P.mu * C31
    ]

sol_rep = solve_ivp(repressilator_ode, (0, 3000), [0.5, 0.2, 1, 0, 0], 
                    t_eval=np.linspace(0, 3000, 3000), method="BDF")
C31_max = np.max(sol_rep.y[4][1000:])
def get_u_in(t): 
    val = np.interp(t, sol_rep.t, sol_rep.y[4])
    return float(np.clip(val / C31_max, 0, 2.0))  

'''


pulse_starts = [300.0, 900.0, 1500.0, 2100.0]
pulse_duration = 100

def get_u_in(t):
    for t_s in pulse_starts:
        if t_s <= t <= t_s + pulse_duration:
            return 1.0
    return 0.0

'''


# =====================================================
# 2. 你的结构精简版：二级计数器 (counter_2bit_ode)
# =====================================================
p = {
    # --- 1. 重组核心动力学 (速率平滑化，适应 ~200min 脉冲宽度) ---
    'k_fwd': 0.118,        # 正向重组速率 (min^-1)
    'k_rev': 0.08,        # 逆向重组速率 (min^-1)
    'Kinh': 0.1026,          # RDF 对正向反应的阻断门槛
    'K_D_comp': 3.2,      # (Int * RDF) 逆向重组复合物激活门槛
    'k_int': 6.0, 
    'gamma_int': 2.0,     # Int0 响应速率

    # --- 2. Bit 0 专属参数 ---
    'K_D_int0': 1.8,      # 避开 ~0.24 的输入底噪，在波峰(~2.5)处翻转
    'n_int0': 4.0,        # 协同结合阶数
    
    'alpha_rep0': 0.1821,
    'gamma_rep0': 0.0119,
    
    'alpha_rdf0': 0.0325, 
    'gamma_rdf0': 0.0055, 

    # --- 3. Bit 1 专属参数 ---
    'K_D_int1': 1.8, 
    'n_int1': 4.0,
    'alpha_rep1': 0.2, 
    'gamma_rep1': 0.005,
    'alpha_rdf1': 0.0265, 
    'gamma_rdf1': 0.005,

    # --- 4. 共享抑制门槛 ---
    'K_rep': 1.8, 
    'n': 3.0,

    # --- 5. Bit 0 -> Bit 1 进位逻辑 (I1-FFL 时序强匹配) ---
    'alpha_A0': 0.1413,    'gamma_A0': 0.0731,    # A0 快响应 (半衰期 ~14 min)
    'n_A0': 2.0,        'K_A0': 1.0,
    'alpha_R0': 0.0979,   'gamma_R0': 0.0116,   # R0 极慢响应 (半衰期 ~350 min)，成功挤出进位窗口
    'n_R0': 5.4332,        'K_R0': 2.6454,
    
    'alpha_Int1': 0.3788,  'gamma_int1': 0.08   # 产生幅值适中的 Int1 脉冲
}

def counter_2bit_ode(t, y, p):
    pb0, int0, rep0, rdf0, a0, r0, \
    pb1, int1, rep1, rdf1 = y

    u_in = get_u_in(t)
    lr0, lr1 = 1.0 - pb0, 1.0 - pb1

    # 引入高阶协同结合 (n_int0)，使电路在波峰到达时才翻转
    int0_act = (int0**p['n_int0']) / (p['K_D_int0']**p['n_int0'] + int0**p['n_int0'] )
    int1_act = (int1**p['n_int1']) / (p['K_D_int1']**p['n_int1'] + int1**p['n_int1'] )

    vf0 = p['k_fwd'] * pb0 * int0_act * (p['Kinh'] / (p['Kinh'] + rdf0 ))
    vr0 = p['k_rev'] * lr0 * (((int0 * rdf0)**2) / (p['K_D_comp']**2 + (int0 * rdf0)**2 ))

    vf1 = p['k_fwd'] * pb1 * int1_act * (p['Kinh'] / (p['Kinh'] + rdf1 ))
    vr1 = p['k_rev'] * lr1 * (((int1 * rdf1)**2) / (p['K_D_comp']**2 + (int1 * rdf1)**2 ))

    # --- Bit 0 动态 ---
    dpb0 = -vf0 + vr0
    dint0 = p['k_int'] * u_in - p['gamma_int'] * int0
    drep0 = p['alpha_rep0'] * pb0 - p['gamma_rep0'] * rep0
    drdf0 = p['alpha_rdf0'] * lr0 * (1.0 / (1.0 + (rep0 / p['K_rep'])**p['n'])) - p['gamma_rdf0'] * rdf0

    # --- Bit 0 -> Bit 1 进位 ---
    da0 = p['alpha_A0'] * pb0 - p['gamma_A0'] * a0
    act0 = (a0**p['n_A0']) / (p['K_A0']**p['n_A0'] + a0**p['n_A0'])
    dr0 = p['alpha_R0'] * act0 - p['gamma_R0'] * r0
    rep_gate0 = (p['K_R0']**p['n_R0']) / (p['K_R0']**p['n_R0'] + r0**p['n_R0'])

    dint1 = p['alpha_Int1'] * act0 * rep_gate0 - p['gamma_int1'] * int1

    # --- Bit 1 动态 (使用 Bit 1 独立参数) ---
    dpb1 = -vf1 + vr1
    drep1 = p['alpha_rep1'] * pb1 - p['gamma_rep1'] * rep1
    drdf1 = p['alpha_rdf1'] * lr1 * (1.0 / (1.0 + (rep1 / p['K_rep'])**p['n'])) - p['gamma_rdf1'] * rdf1

    return [dpb0, dint0, drep0, drdf0, da0, dr0, 
            dpb1, dint1, drep1, drdf1]

# =====================================================
# 4. 仿真与绘图
# =====================================================
a0_ss = p['alpha_A0'] / p['gamma_A0']
act0_ss = (a0_ss**p['n_A0']) / (p['K_A0']**p['n_A0'] + a0_ss**p['n_A0'])
r0_ss = p['alpha_R0'] * act0_ss / p['gamma_R0']

rep0_ss = p['alpha_rep0'] / p['gamma_rep0']
rep1_ss = p['alpha_rep1'] / p['gamma_rep1']

# 初始条件 (按你原始代码的默认初始设置，初始 PB0=1, PB1=1)
y0 = [1.0, 0.0, rep0_ss, 0.0, a0_ss, r0_ss,
      1.0, 0.0, rep1_ss, 0.0]

t_eval = np.linspace(0, 2000, 2000)
sol = solve_ivp(counter_2bit_ode, (0, 2000), y0, args=(p,), t_eval=t_eval, method="RK45")
t_hrs = sol.t / 60.0

LR0, LR1= 1.0 - sol.y[0], 1.0 - sol.y[6], 
A0, R0, Int1 = sol.y[4], sol.y[5], sol.y[7]
Int0 = sol.y[1]

# =====================================================
# 3. 绘图 (严格分为 3 个子图)
# =====================================================
fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True, dpi=120)

# 子图 1: LR 状态 (显示重组状态曲线)
axs[0].plot(t_hrs, LR0, 'steelblue', linewidth=2, label=r'Bit 0 ($LR_0$)')
axs[0].plot(t_hrs, LR1, 'goldenrod', linewidth=2, label=r'Bit 1 ($LR_1$)')
axs[0].set_ylim(-0.05, 1.15)
axs[0].set_ylabel('Bit State (LR)', fontsize=9)
axs[0].set_title('1. Recombination States ', fontsize=10, fontweight='bold')
axs[0].legend(loc='upper right', ncol=2)
axs[0].grid(True, linestyle=':', alpha=0.6)

# 子图 2: I1-FFL 调控因子 (A0 与 R0)
axs[1].plot(t_hrs, A0, 'steelblue', linewidth=1.8, label=r'Activator $A_0$')
axs[1].plot(t_hrs, R0, 'goldenrod', linewidth=1.8, label=r'Repressor $R_0$')
axs[1].set_ylabel('Concentration (a.u.)', fontsize=9)
axs[1].set_title(r'2. I1-FFL Carry Regulation Factors ($A_0$ & $R_0$)', fontsize=10, fontweight='bold')
axs[1].legend(loc='upper right', ncol=2)
axs[1].grid(True, linestyle=':', alpha=0.6)

# 子图 3: 整合酶 (Int0 输入与生成的 Int1 脉冲)
axs[2].plot(t_hrs, Int0, 'k:', linewidth=1.3, label=r'$Int_0$ (External Input)')
axs[2].plot(t_hrs, Int1, 'darkslategray', linewidth=2.0, label=r'$Int_1$ (Carry Pulse to Bit 1)')
axs[2].set_xlabel('Time (Hours)', fontsize=9)
axs[2].set_ylabel('Concentration (a.u.)', fontsize=9)
axs[2].set_title(r'3. Integrase Dynamics ($Int_0$ driving Bit 0, $Int_1$ driving Bit 1)', fontsize=10, fontweight='bold')
axs[2].legend(loc='upper right', ncol=2)
axs[2].grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.show()