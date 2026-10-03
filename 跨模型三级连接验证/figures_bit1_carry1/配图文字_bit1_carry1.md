# bit1 / carry1 配图文字与出处

对应两张图：

| 文件 | 内容 | 状态 |
|---|---|---|
| `out/fig_bit1_internal.png` | bit1（ZMH 约化位，6 状态）内部流程 | 已过 4 项数值检查 |
| `out/fig_carry1_wired.png` | bit1 → carry1 → bit2 的完整接线 | 已过 4 项数值检查 |

生成脚本 `make_bit1_carry1_sketches.py`；布局自检 `check_bit1_carry1_layout.py`。
**这是接线草图（wiring draft），不是数据图**：不读任何仿真输出。

---

## 1 可以直接用的图注

> **图 A. bit1（ZMH 约化位）的内部流程机制。**
> 进位输入来自 bit0 的 `b0_S`：`A0` 的源项是 `0.1413·(1 − b0_S)`，这是 A0 唯一的
> 天然驱动（面板 ①）。A0 一边激活整合酶 `I1`，一边通过延迟臂 `F0` 抑制 I1，
> 构成非相干前馈环：`I1 源 = 0.3788·H(A0;1.0,2)·(1 − H(F0;2.6454,5.4332))`。
> 面板 ② 是开关 `pb1`（= `PB1`）：`dpb1/dt = −vf + vr`，正向通量由 I1 驱动并被
> `RDF1` 抑制，反向通量由代数复体 `(I1·RDF1)²` 驱动。`pb1` 驱动 `T1`，
> `1 − pb1`（即 LR1）驱动 `RDF1`，而 `T1` 反过来抑制 `RDF1`。
> **与 bit0 的结构差别**：bit1 没有 C 状态，复体是代数量 `(I1·RDF1)²`；
> bit0/bit2 有显式 C 状态与 `k_complex = q·K_D_comp`。
> **参数版本**：A0/F0/pb1/I1/T1/RDF1 取自 donor 自己的 `p`（`wiki任务/完整二级级联.py`）；
> 对应 `model_hzh.py:65-77` 的 `middle_rhs`。示意图，非数据图。

> **图 B. bit1 接上 carry1 与 bit2 之后的完整流程机制。**
> `carry1` 是 HBY 的 6 状态块（A1、F1 各一条 mRNA→未成熟→成熟 的三级链）。
> 两条链的源分别由 bit1 的 `pb1`（`A1 源 = 16.0·pb1·(1 − H(A1;0.6,2))`，带 A1 负自调）
> 与 A1（`F1 源 = 5.0·H(A1;1.2,4)`）给出。三因子 AND 门
> `g1 = H(A1;1.2,6)·(1 − H(F1;0.4,4))·clock` 产生 `Int2 源 = α_Int[1]·g1 = 38.0·g1`，
> 进入 bit2 的 I 链。`clock = H(b0_I;0.3,2.0)` 读的是 **bit0** 的成熟游离整合酶——
> 它跨过 carry0 直接门控 carry1，是本模型的关键设计。
> **参数版本**：A1/F1/bit2 取自 ZENG 表（`final_reconstruction/model.py`）；
> carry1 见 `hybrid_model.py:200-218`，门与时钟见 `:189-198`。示意图，非数据图。

---

## 2 两套参数表不能混（本图最容易出错的地方）

同一个符号在两张表里存在且**数值不同**，图里必须分开标注：

| 符号 | donor（`self.z`，用于 bit1） | ZENG（`self.tail.p`，用于 carry1/bit2） |
|---|---|---|
| `k_fwd` | **0.118** | **7.0** |
| `k_rev` | **0.08** | **5.0** |
| `K_D_comp` | **3.2** | **1.2** |
| 抑制常数 | `Kinh = 0.1026` | `K_inh = 0.1` |

`model_hzh.py` 里 `self.z = donor_namespace()['p'].copy()`（`:22`）与
`self.tail = HbyReceiver()`（`:20`，其 `self.p = load_zeng_table()`）是两个独立来源。
bit0 与 bit2 用 ZENG，bit1 的 6 个状态用 donor——**看箭头上的数字时先看它属于哪一段**。

---

## 3 每条箭头对应的 RHS 项

| 图 | 箭头 | 源码 |
|---|---|---|
| A | `b0_S → A0` | `model_hzh.py:73` `alpha_A0*(1-b0s)`（`b0s` = `y[10]` = `b0_S`，见 `:89`） |
| A | `A0 → ×g0 → I1` | `:68` `act` + `:69` `gate` + `:75` `alpha_Int1*act*gate` |
| A | `A0 → F0` | `:74` `alpha_R0*act`（`act` 见 `:68`；代码里抑制臂叫 `R0`/`r0`） |
| A | `F0 ⊣ ×g0` | `:69` `K_R0^n_R0/(K_R0^n_R0 + f^n_R0)` |
| A | `I1 →- vf → pb1` | `:71` `k_fwd*pb*ia*(Kinh/(Kinh+r))` + `:74` `-vf+vr` |
| A | `RDF1 ⊣ vf` | `:71` 分母的 `Kinh/(Kinh+r)`（与 bit0 的 ⑨ 同义） |
| A | `pb1 ⇄（vr）` | `:72` `k_rev*(1-pb)*((i*r)^2/(K_D_comp^2+(i*r)^2))` |
| A | `pb1 → T1` | `:76` `alpha_rep1*pb - gamma_rep1*tr` |
| A | `(1−pb1) → RDF1`、`T1 ⊣ RDF1` | `:77` `alpha_rdf1*(1-pb)*(1/(1+(tr/K_rep)^n)) - gamma_rdf1*r` |
| B | `pb1 → A1 链`、`A1 负自调` | `hybrid_model.py:204` `alpha_A[1]*pb1*(1-hill(a,K_auto1,n_auto1))` |
| B | `A1 → F1 链` | `:205` `alpha_F[1]*hill(a,K_A[1],n_A[1])` |
| B | 两条链的三级 | `:184-187` `expression()`（`translation_h=30`、carry mRNA 2.0 min / 成熟 32.5 min） |
| B | `A1 → g1`、`F1 ⊣ g1`、`clock → g1` | `:194` `hill(a,K_A[1],n_A1_gate)*(1-hill(f,K_F[1],n_F[1]))*clock` |
| B | `clock` | `:193` `hill(clock_scale*int0, clock_K_au, clock_n)`，`int0 = y[2] = b0_I`（`:90`） |
| B | `g1 → bit2` | `:197` `u2_target_au_per_h = alpha_Int[1]*gate`；bit2 链见 `:210` |

---

## 4 一处需要注意的"两个 Hill 系数"

`A1` 进门用 `n = 6.0`（`HbyConfig.n_A1_gate`，frozen extension 的值），
但 A1 驱动 F1 时用 `n_A[1] = 4.0`；两处 `K_A[1] = 1.2` 相同。
图里两侧都标注了，引用时别把 6.0 和 4.0 混成一个。

---

## 5 这两张图做过的数值检查

画图的人看不到图，所以全部用数值验证：

| 检查 | 结果 |
|---|---|
| 框-框重叠 | 0 |
| 箭头穿框 | 0 |
| 悬空箭头端（tip 不落在任何框/线上） | 0 |
| **箭头互相交叉**（bit0 的检查里没有这一项，本图新增） | 0 |
| 文字-文字重叠 | 0 |
| 文字-框重叠 | 0 |
| T 形横杠落点：`RDF1 ⊣ vf` | 到 `vf` 折线距离 **0.000e+00** |
| T 形横杠落点：`A1 负自调 ⊣ A1 源` | 到该折线距离 **0.000e+00** |
| 每个框都有入边（除 `b1.blk` 这个最左输入块） | 通过 |
| 图上 27 个常数逐个与模型比对 | **27/27 一致** |

**为什么 bit1 那张分成两个面板**：RDF1 既要抑制 `vf`（一条反向边）、又要喂代数复体，
放在同一面板里任何摆法都会让某两条线交叉。拆成"进位门"与"开关+输出"后各自可平面化，
跨面板的线用 `I1 >` / `< I1`、`PB1 >` / `< PB1` 网络标签衔接（与 bit0 图同一约定）。

---

## 6 复现

```
cd 跨模型三级连接验证/figures_bit1_carry1
& 'D:\aconade\python.exe' -B .\make_bit1_carry1_sketches.py      # 出图
& 'D:\aconade\python.exe' -B .\check_bit1_carry1_layout.py       # 布局自检，应 exit 0
```

脚本复用 `figures_bit0/make_bit0_sketches.py` 里的 `Draft` 与几何检查，
不另抄一份，避免两处风格漂移。注意 SVG 哈希不可复现（matplotlib 写 `<dc:date>`
和随机元素 id），要指纹请用 PNG。
