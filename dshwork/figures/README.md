# 机制图 / 流程图的代码生成

对应 `04_第3章_作图说明.md` 与 `05_第4章_作图说明.md` 中**可以用代码绘制**的示意图。
数据图（读数条带、热图、剂量—响应）不在本目录，它们需要仿真输出。

## 运行

```powershell
$env:MPLBACKEND='Agg'
& 'D:\aconade\python.exe' .\make_schematics.py --out .\out
# 只画其中一张：
& 'D:\aconade\python.exe' .\make_schematics.py --out .\out --only fig04_ffl_carry
```

依赖只有 `matplotlib`（`Agg` 后端，无需图形界面）与 `numpy`。

## 生成的文件

| 图 | 文件名 | 类型 | 来源 |
|---|---|---|---|
| Figure 1 | `fig01_modeling_workflow.png` | 流程图 | `01` 图 1：Design → Build → Test → Learn → Decision 闭环 |
| Figure 3-1 | `fig03_single_bit_mechanism.png` | 机制图 | `04` §3：Int / RDF / BM3R1 开关 |
| Figure 4-1 | `fig04_ffl_carry.png` | 机制图 | `05` §3：I1-FFL 进位门（含两个时间过程小图） |
| Figure 4-6 | `fig04_6_bm3r1_shutdown_interface.png` | 机制图 | `05` §8：BM3R1 共享接口 |
| Figure 7 | `fig07_full_system_coupling.png` | 框架图 | `01` 图 7：全系统耦合 |
| Figure 8 | `fig08_design_panel.png` | 面板 | `01` 图 8：Modeling-Guided Design |

## 两个文件的分工

- **`tempo_style.py`**：调色板、线型约定与绘图原语。所有颜色与线型与 Oscillator 章节一致，
  这样同一个 Wiki 上的机制图属于同一视觉家族。
  - 原语：`box` / `arrow` / `repression`（T 形端点）/ `binding`（双向箭头）/
    `and_node`（× 号 AND 节点）/ `hline`（阈值虚线）/ `scope`（模块边界虚线框）/
    `hollow_arrow`（"尚未建模"的空心虚线箭头）/ `legend_lines`（线型图例）
- **`make_schematics.py`**：每张图一个函数，坐标全部用数据坐标显式写出，便于微调。

## 约定

1. **全部标签只用英文**。作图说明要求"英文栏可直接复制到图片中"，
   而且英文彻底避开中文字体依赖，headless 渲染不会出现缺字方框。
2. **不要改颜色**。若要调整，改 `tempo_style.py` 里的常量，所有图一起变。
3. **AND 必须画成 `and_node`（× 号圆节点）**，不要用加号或并联箭头——
   "相乘"是进位机制成立的前提。
4. **抑制必须用 `repression`（T 形端点）**，与激活箭头区分开。
5. 需要新增图时，在 `make_schematics.py` 里写一个返回 `save(fig, ...)` 的函数，
   并注册到末尾的 `FIGURES` 字典即可。

## 尚未由代码生成的图

| 图 | 原因 |
|---|---|
| 3-2 单周期翻转分解 | 需要冻结模型的真实轨迹 |
| 3-3 读窗几何与驻留带 | Panel 1–2 是定义、可以画；Panel 3 需要真实轨迹做新旧判据对照 |
| 3-4 工作窗口 | 需要 `uM × 进位成熟时间` 二维细化扫描结果 |
| 3-5 扰动与容差 | 需要四组分子池扰动扫描结果 |
| 4-2 门对比度 | 需要单周期的 $g_0$ 曲线 |
| 4-3 进位因果链 | 需要 $J_{rev,0}$ / $g_0$ / 重组酶 / $S_1$ 四行同步轨迹 |
| 4-4 300 h 两位读数 | 需要 300 h 轨迹 |
| 4-5 工作窗口热图 | 同上，需要二维扫描结果 |
| 4-7 三级机制与悬崖 | 需要 63 点扫描的逐点结果 |

这些可以由同一套 `tempo_style.py` 原语加上仿真输出来生成，配色自动保持一致。
