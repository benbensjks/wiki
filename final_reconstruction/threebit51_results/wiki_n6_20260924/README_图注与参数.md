# 51状态成功模8计数图（2026-09-24）

这组图按已归档成功扫描的实际构造复算600 h，未改变模型、ZENG或冻结参数文件。

## 选择哪张图

- `01_threebit_overview_600h.png`：全程56个有效读窗，适合完整结果展示。
- `02_threebit_zoom_90_190h.png`：真实稳态片段，完整显示1→2→…→7→0，建议Wiki正文优先用这张。
- `03_mod8_readout.png`：三位数字条带，适合紧凑展示。

每张同时提供PDF与SVG；所有曲线来自真实ODE轨迹，整数阶梯仅连接读窗解码结果，不表示转变期间的连续数字输出。

## 本次确认的配置差异

成功扫描 `scan_gate_segments.build` / `check_read_commitment.build` 使用 `frozen_extension()`，从34状态selected_profile继承clock_K_au=0.3、clock_n=2。
51状态默认函数和曾同学参考代码中的0.4/3是另一套值。因此本次图使用实际已成功的0.3/2，不声称已复现0.4/3版本。

共同参数：uM_per_au=5.75；bit mRNA半衰期2 min、成熟半衰期20 min；carry0/1 mRNA半衰期2 min、A/F成熟半衰期32.5 min；独立n_A1_gate=6；F1产生指数4；add_growth=False。

## 可直接使用的图注

**中文：** 51状态三级计数器的确定性ODE轨迹。上游C31翻译通量驱动bit0，两级进位依次驱动bit1与bit2。600 h仿真产生56个有效读窗，按模8重复输出1–2–3–4–5–6–7–0。放大图取90–190 h真实片段；0.3/0.7虚线为读带边界，浅色竖带为周期20%的读窗。独立A1门指数为6；本图实际时钟门为K=0.3、n=2。结果为模型预测，尚未经湿实验验证。

**English:** Deterministic simulation of the 51-state three-bit counter. The C31 translation flux drives bit0, while two carry modules drive bit1 and bit2. The 600 h trajectory yields 56 valid finite-window reads following 1–2–3–4–5–6–7–0. The zoom shows the actual 90–190 h segment. The independent A1 gate exponent is 6; the clock gate used here has K=0.3 and n=2. Steps connect decoded reads only.

## 数据与复现

- `trajectory.npz`：全51状态、时间和状态名。
- `trajectories.csv`：全状态及物理信号。
- `read_windows.csv`：56个读窗、三个位和值。
- `parameters.json`：完整实际配置、求解器设置及源文件SHA256。
- `verification.json` / `summary.json`：实际验收结果。
- `make_wiki_threebit_n6_figures.py`（模型根目录）：生成脚本；仅复用参数与源文件哈希匹配的缓存。

C31 flux轴为µM/h；底部carry轴为下一位整合酶的目标产生速率a.u./h，经过现有转录与成熟链后才贡献成熟整合酶，不是游离Int浓度。
