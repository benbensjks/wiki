# bit0 草图（手绘底稿）

给混合三级计数器 **bit0** 画的两张流程示意图，用作手绘底稿。

| 文件 | 内容 | 尺寸 |
|---|---|---|
| `out/fig_bit0_internal.png` / `.svg` | **bit0 内部**：11 个状态、全部反应、开关 `S/PB` | 20 × 13.6 in |
| `out/fig_bit0_wired.png` / `.svg` | **bit0 接上上游与下游 carry**：上游振荡器 + bit0 + carry0 + bit1 + clock 接口 | 32 × 18.6 in |

SVG 用 `svg.fonttype='none'` 导出，**文字是真 `<text>`**（64 / 90 个），可直接在 Illustrator / Inkscape 里改字，也可以逐字校验。

---

## 怎么生成

```powershell
cd C:\Users\18633\Desktop\wiki\跨模型三级连接验证\figures_bit0
& 'D:\aconade\python.exe' -B .\make_bit0_sketches.py          # 出图
& 'D:\aconade\python.exe' -B .\check_bit0_sketch_layout.py    # 排版自检
```

---

## 抽象序号 → 真实生化元件

### bit0 的 11 个状态（两张图共用同一套小写标签）

| 标号 | 代码名 | 真实生化元件 |
|---|---|---|
| **a** | `b0_M_I` | 整合酶 mRNA（ϕC31 整合酶转录本） |
| **b** | `b0_I_u` | 未成熟整合酶（翻译产物） |
| **c** | `b0_I` | 成熟游离整合酶（重组酶本体） |
| **d** | `b0_M_T` | BM3R1 / Rep 的 mRNA |
| **e** | `b0_T_u` | 未成熟 BM3R1 |
| **f** | `b0_T` | 成熟 BM3R1 / Rep（阻遏蛋白） |
| **g** | `b0_M_R` | RDF 的 mRNA |
| **h** | `b0_R_u` | 未成熟 RDF |
| **i** | `b0_R` | 成熟游离 RDF（决定重组方向） |
| **j** | `b0_C` | Int·RDF 复合物（显式结合态） |
| **k** | `b0_S` | LR 构象比例；`PB = 1 − S` |

（bit2 用同一套 11 状态模板，代码名加 `b2_` 前缀。）

### 反应编号

| 编号 | 反应 | 依据 |
|---|---|---|
| ① | 转录：上游 C31 启动子活性 → `M_I` | `model_hzh.py:49` |
| ② | 翻译：`M_I` → `I_u` | `model_hzh.py:52` |
| ③ | 成熟：`I_u` → `I` | `model_hzh.py:54` |
| ④ | Rep 表达链，源 `α_rep·(1−S)`（PB 驱动） | `model_hzh.py:55` |
| ⑤ | RDF 表达链，源 `α_rdf·S·(1−H(T))`（LR 驱动） | `model_hzh.py:56` |
| ⑥ | `T` 抑制 RDF 表达，`K_rep = 0.85`、`n_rep = 3.9` | 同上 |
| ⑦ | 可逆结合 `I + R` 生成 `C`；`k_on = 0.1`、`k_off = 1`、`δ_C = 1` | `model_hzh.py:57-59` |
| ⑧ | 正向重组 `v_f = k_fwd·H(I;K_D_int,2)·K_inh/(K_inh+R)` | `model_hzh.py:60` |
| ⑨ | 游离 RDF 抑制正向 —— 就是 ⑧ 分母里的 `K_inh/(K_inh+R)` | 同上 |
| ⑩ | 反向重组 `v_r = k_rev·H(C;K_complex,2)` | `model_hzh.py:61` |
| (11) (12) | 开关的两条输出：PB 与 LR 分别驱动两条表达链 | `model_hzh.py:55-56` |
| 开关 | `dS/dt = v_f·(1−S) − v_r·S` | `model_hzh.py:62` |

### 第二张图里的大写标号

| 标号 | 真实元件 | 来源 |
|---|---|---|
| **A** | TetR：振荡器阻遏蛋白 | 韩亚轩 v53d |
| **B** | CI | 同上 |
| **C** | LacI | 同上 |
| **D** | C31 启动子（PLtetO1 型，输出活性 `h31(t)`） | 同上 |
| **E** | `A0_zmh`：进位激活臂 | 曾同学现行 Python 代码 |
| **F** | `F0_zmh`：进位抑制臂（曾代码里变量名叫 `R0`） | 同上 |
| **G** | `I1_zmh`：bit1 的整合酶 | 同上 |
| **H** | bit1 约化模块 `pb1 / I1 / T1 / RDF1`（4 状态） | 同上 |
| **I** | clock 门 `H(b0_I; K=0.3, n=2)` | `hybrid_model.py:193` |

---

## 画图时最容易错的四点（已固化在脚本里）

1. **转录与翻译分开**：上游启动子活性喂的是 `M_I`（①），翻译通量喂的是 `I_u`（②）。两者不是同一根线。
2. **clock 不是总线**：`g0` 只有两项、**没有时钟**；只有 `g1` 才有三因子 AND（`H(A1;1.2,6)·G(F1;0.4,4)·clock`），且 clock 读的是 **bit0 自己的成熟游离整合酶 `c`（`b0_I`）**。
3. **两级 carry 不对称**：`g0` 扇入 2、`g1` 扇入 3；`g1` 的门臂指数 6 与 F1 产生臂的 4 是解耦的。
4. **bit1 是约化位**：4 个状态，**没有** mRNA、**没有**未成熟态、**没有** `C` 复合物，且 `pb1` 是直接积分的（不是 `1 − S`）。

---

## 验证记录

`check_bit0_sketch_layout.py` 做两类检查，**全部通过**：

| 检查 | fig_bit0_internal | fig_bit0_wired |
|---|---|---|
| 方框—方框重叠 | 0 | 0 |
| 连线穿过方框 | 0 | 0 |
| 文字—文字重叠（取真实渲染包围盒） | 0 | 0 |
| 文字—方框重叠 | 0 | 0 |

说明：这三张图**没有经过目视检查**（当前模型无图像输入），上表是数值替代方案。
检查器取的是 matplotlib 实际渲染的字形包围盒，不是估算值。

另外确认：**没有缺字形警告**。微软雅黑不含 `⇄ ⊣ ⊗ ▸ ◂ ⑪ ⑫`，脚本里已替换为
`/`、`抑制`、`×`、`>`/`<`、`(11)`/`(12)`。

---

## 产物哈希（SHA256）

```
fig_bit0_internal.png  26D23632CD5580CFAFD599B9D30B4A98DE19009ADE01786DD8525407B04BD314
fig_bit0_internal.svg  95155608D5C4007083C6F0EDD2CA09BA4D66EF9530B716AF45BF150ED142857A
fig_bit0_wired.png     39985C733C1CE3C01CE00AAAF6B7193E1E0EDB0490090DC32D364D55B8828A37
fig_bit0_wired.svg     7DC79F2E932D8D8E1384F0C27E465FDE2A5E7EEA5A30B31D844E73D7EEDA9A21
```

---

## 文件

```
figures_bit0/
├─ make_bit0_sketches.py          生成两张图（内含几何自检）
├─ check_bit0_sketch_layout.py    文字/方框排版自检（本次为数值替代目视）
├─ tempo_style.py                 项目统一配色与图元（从 dshwork/figures 复制）
├─ 配图文字_bit0.md                图注、逐元素描述、英文对照、正文段落
├─ 生图提示词_bit0.md               交给专用生图 AI 的提示词（含负向提示词与校验清单）
├─ README.md                      本文件
└─ out/                           PNG + SVG
```

`tempo_style.py` 是从 `dshwork/figures/tempo_style.py` 复制的副本
（源文件 SHA256 `70E7F7A5D1C93907DD352E6000563B331052C9865240739E5C5E346F191049E6`），
复制过来是为了本目录自包含；改配色要改源文件再重新复制。

---

## 尚未做的

- **carry1 与 bit2 没有展开**：第二张图只画到 clock 接口，carry1（HBY A1/F1）与 bit2 未画。
- **清除/衰减项未逐条画箭头**：各分子按各自的 λ 或 γ 衰减，只在图例里说明。
- **数值未标注**：图上是结构，`K`、`n`、半衰期等具体数值留给标注轮。
