# 生图提示词：bit0 内部级联机制图

给专用生图 AI 用的提示词。严格按 `hby_zmh_hby/model_hzh.py:39-63` 的 bit0 真实级联写。

---

## 0 使用前必读：文字策略

**文生图模型画不准文字，长句一定糊。** 所以这份提示词把图上文字压缩到极限：

| 图上会出现 | 图上**不会**出现 |
|---|---|
| 11 个单字母 `a b c d e f g h i j k` | 任何英文/中文化学名称 |
| 10 个单数字 `1 2 3 4 5 6 7 8 9 10` | 标题、图例、公式、参数值、说明句 |
| 约 4 个纯符号 `PB`、`LR` | 任何长度超过 3 字符的字符串 |

**图例（字母 → 真实元件）由你自己后期贴上去**，中英文对照表在本文 §6。

这样做的第二个好处：生图模型最擅长画"几何 + 箭头"，最不擅长排版文字。把文字拿掉，
它就能把力气全用在拓扑上。

---

## 1 主提示词（英文版 · 推荐）

> 英文提示词对多数生图模型效果更好。直接整段粘贴。

```
A flat vector textbook-style schematic diagram of one genetic toggle switch (a
site-specific recombinase bit), on a pure white background. Hub-and-spokes
composition. Absolutely no photorealism, no 3D, no glow, no cells, no DNA
helices, no molecular surface renderings, no decorative flourishes. Only rounded
rectangles, arrows and one circle.

=== LAYOUT ===

CENTER of the canvas: one circle with a thin dark outline, containing the single
lowercase letter "k". This is the ONLY circle in the whole image.

UPPER LEFT: a vertical stack of three identical rounded rectangles, evenly
spaced, labelled top to bottom with one lowercase letter each: "a", "b", "c".
A short straight arrow points down from "a" into "b", and another from "b" into
"c". A separate short arrow enters the top of "a" coming from above; its tail is
free and does not touch anything.

UPPER RIGHT: a vertical stack of three identical rounded rectangles, top to
bottom: "d", "e", "f". A short straight arrow points down from "d" into "e" and
from "e" into "f".

LOWER RIGHT: a vertical stack of three identical rounded rectangles, top to
bottom: "g", "h", "i". A short straight arrow points down from "g" into "h" and
from "h" into "i".

LOWER LEFT: a single rounded rectangle containing the letter "j".

The three stacks and the single box are arranged symmetrically around the centre
circle, leaving generous white space between them.

=== CONNECTORS (the arrowheads matter most) ===

1) From the circle "k", a smooth curved arrow sweeps up and to the right, ending
   with a solid triangular arrowhead on the LEFT side of box "d".

2) From the circle "k", a smooth curved arrow sweeps down and to the right,
   ending with a solid triangular arrowhead on the LEFT side of box "g".

3) From the RIGHT side of box "c", a smooth curved arrow sweeps right and
   slightly down, ending with a solid triangular arrowhead on the upper-left of
   circle "k".

4) From the upper-right corner of box "j", a smooth curved arrow sweeps up and
   to the right, ending with a solid triangular arrowhead on the lower-left of
   circle "k".

5) A DOUBLE-HEADED arrow — arrowheads at BOTH ends — links the bottom of box "c"
   to the upper-left corner of box "j".

6) A DOUBLE-HEADED arrow — arrowheads at BOTH ends — links the LEFT side of box
   "i" to the RIGHT side of box "j".

7) From the bottom of box "f" a short straight line points straight down and ends
   in a FLAT PERPENDICULAR BAR (a "T" shape, like a blunt wall), positioned just
   above box "g". This is a repression symbol, not an arrowhead.

8) From the upper-left corner of box "i", a short line curves up and to the LEFT
   and ends in a FLAT PERPENDICULAR BAR (a "T" shape, like a blunt wall). The bar
   must sit ON the curved arrow from "c" into "k" (arrow 8), crossing it, so that
   the T visibly blocks that arrow. This is a repression symbol, not an
   arrowhead. Do NOT leave it floating in empty space: the bar has to touch the
   arrow it blocks.

=== SMALL NUMBER LABELS ===

Small plain digits, each inside a tiny open circle, placed beside the
corresponding arrow:

"1"  beside the free incoming arrow above box "a"
"2"  beside the short arrow from "a" to "b"
"3"  beside the short arrow from "b" to "c"
"4"  beside the curved arrow from "k" up to "d"
"5"  beside the curved arrow from "k" down to "g"
"6"  beside the short T-barred line between "f" and "g"
"7"  beside the double-headed arrow between "c" and "j"
"8"  beside the curved arrow from "c" into "k"
"9"  beside the T-bar that blocks the curved arrow from "c" into "k"
"10" beside the curved arrow from "j" into "k"

=== EXTRA SHORT LABELS (only these three-character strings allowed) ===

Just to the LEFT of the arrow numbered 4, write the two letters "PB".
Just to the LEFT of the arrow numbered 5, write the two letters "LR".
Nothing else. No other words anywhere in the image.

=== TEXT POLICY ===

The ONLY text allowed in the entire image is:
  - the eleven single lowercase letters a b c d e f g h i j k inside the boxes
    and the circle;
  - the small numbers 1 to 10 beside the arrows;
  - the two two-letter strings "PB" and "LR".
Do NOT render any title, caption, legend, sentence, formula, axis, scale bar,
watermark or signature. Leave the entire bottom third of the canvas as empty
white space.

=== STYLE ===

Flat vector illustration. Uniform thin line weight, about 2 px. Rounded corners
on all rectangles. Pure white background. Exactly three muted colours:
  - dark slate blue-grey #304B53 for every outline, every arrow stroke and all
    lettering;
  - muted teal #4F9194 for the fill of boxes "c", "f" and "i" only;
  - muted orange #D29144 for the two T-barred repression lines (items 7 and 8)
    and for the free incoming arrow above box "a".
All remaining boxes are very light grey #F3F5F6 with a #304B53 outline and no
gradient. Everything snapped to a grid, perfectly aligned, wide margins, no
overlapping lines, no crossing arrows, no shadows.

Landscape aspect ratio 3:2.
```

---

## 2 负向提示词（Negative prompt）

```
photorealistic, 3D render, CGI, glossy plastic, glass, metal, glowing, neon,
gradient shading, drop shadow, depth of field, bokeh, microscope image, petri
dish, E. coli cells, DNA double helix, protein ribbon, molecular surface,
cartoon mascot, cute character, hands, people, laboratory equipment, test tubes,
clip art, clipart, stock illustration, busy background, textured paper,
hand-drawn sketch, pencil, watercolour, low resolution, blurry, jpeg artifacts,

any title, any caption, any paragraph of text, any legend block, any table,
any formula, any axis, any scale bar, any watermark, any signature, any logo,

misspelled letters, extra letters, random words, garbled text, Chinese
characters, Cyrillic characters,

more than one circle, extra circles, extra boxes, extra arrows, arrows that
cross each other, arrowheads at the wrong end, missing arrowheads, curved lines
that miss their target box, overlapping shapes, cramped layout, elements pushed
to the edge of the frame,

repression bars floating in empty space, T-bars that do not touch the arrow they
block, a T-bar parallel to the line instead of crossing it, any box or circle
with no incoming arrow, any box or circle with no outgoing arrow, an arrow that
stops short of the box it should point at.
```

---

## 3 技术参数

| 项 | 建议值 | 说明 |
|---|---|---|
| 比例 | 3:2 横版（或 16:9） | 留出下部三分之一空白给图例 |
| 分辨率 | ≥ 2048 px 长边 | 便于后期在 Inkscape/PPT 里加字 |
| 风格关键词 | flat vector, schematic, textbook diagram, minimal | 别用 "infographic"（会加装饰） |
| 采样步数 | 模型默认偏高档 | 拓扑类图需要收敛 |
| 生成张数 | **至少 8 张** | 这类图单张成功率低，拓扑经常崩 |
| 后期 | 导出 SVG 或矢量图 | 位图后期加字会糊 |

**如果模型支持"参考图"（image-to-image / reference）**：把我生成的
`out/fig_bit0_internal.png` 当参考图一起喂进去，保真度会明显提高——即使它排版不好看，
拓扑是对的。

---

## 4 精简兜底版（主版画不出来时用）

如果模型处理不了 11 个框 + 10 个编号，退到这一版：**只画 7 个元素、5 条线**，
把两条表达链各压成一个框。

```
A flat vector schematic of a genetic toggle switch, hub and spokes, white
background, no text except single letters.

Centre: one circle labelled "k".
Upper left: two stacked rounded rectangles, "a" above "b", a short arrow from
  "a" down into "b", and a free short arrow entering "a" from above.
Upper right: one rounded rectangle labelled "d".
Lower right: one rounded rectangle labelled "g".
Lower left: one rounded rectangle labelled "j".

Arrows: a curved arrow from "k" up-right into "d"; a curved arrow from "k"
down-right into "g"; a curved arrow from "b" right into "k"; a curved arrow from
"j" up-right into "k".
Double-headed arrow from "b" to "j".
Short line from "d" straight down ending in a flat perpendicular bar near "g".

Style: flat vector, thin dark slate #304B53 outlines, teal #4F9194 fill on "d",
white background, no gradients, no shadows, no words, no title. Landscape 3:2.
```

---

## 5 逐条校验清单（拿到图后对照）

生成出来先别急着用，按这张表逐条看。**任何一条不满足就重生成**——
这类图"看着像"和"拓扑对"是两回事。

| # | 检查项 | 依据（源码） |
|---|---|---|
| 1 | 全图**只有一个圆圈** `k` | `b0_S` 是唯一的构象比例状态 |
| 2 | `a→b→c` 三个框竖排、箭头向下 | `model_hzh.py:51-54`（转录→翻译→成熟） |
| 3 | `a` 上方有个**游离箭头**指入（外源输入） | `:49-51`，转录来自上游启动子 |
| 4 | 上右链 `d→e→f` 竖排、箭头向下 | `:55`，Rep 表达链 |
| 5 | 下右链 `g→h→i` 竖排、箭头向下 | `:56`，RDF 表达链 |
| 6 | `f` 到下右链之间是 **T 形横杠**（不是箭头） | `:56` 里的 `(1 − H(T;K_rep,n_rep))` |
| 7 | `c` 与 `i` **都**有**双向箭头**连到 `j` | `:57-59`，`I + R ⇄ C`，两条都消耗游离池 |
| 8 | `k` 有 **2 条入**（来自 `c`、来自 `j`）**2 条出**（去 `d`、去 `g`） | `:60-62` 加上 `:55-56` 的源项 |
| 9 | `i` 左侧有 T 形短截线，且**横杠压在 `c→k` 那条箭头上**（不能悬空） | `:60` 分母 `K_inh/(K_inh+R)` —— 它抑制的是正向通量本身 |
| 10 | **没有任何**从 `k` 到 `j`、或从 `j` 到 `c`/`i` 的箭头 | `d[9]` 只由 `bind − un − δ_C·C` 决定 |
| 11 | **没有第 12 个框**（不要出现 "E0/E1"、"U"、"Int1" 之类） | bit0 恰好 11 个状态 |
| 12 | 图的下三分之一是**空白** | 留给图例 |

**最容易出现的三个错**（生图模型高频踩坑）：

- 把 `i` 那条 T 形截线画成**悬空**、没有压在 `c→k` 的箭头上（它必须压住那条箭头，
  `RDF` 抑制的就是这条正向通量；悬空的 T 形等于说"抑制某个没名字的东西"）
- 把两条双向箭头画成**单向**（`c↔j`、`i↔j` 都必须双头）
- 多画一条 `k→j` 或 `j→k` 之外的线，或把 `f⊣g` 画成箭头

---

## 6 字母 → 真实生化元件（后期贴的图例）

| 字母 | 代码名 | 中文 | English |
|---|---|---|---|
| a | `b0_M_I` | 整合酶 mRNA（ϕC31） | integrase mRNA |
| b | `b0_I_u` | 未成熟整合酶 | immature integrase |
| c | `b0_I` | 成熟游离整合酶 | mature free integrase |
| d | `b0_M_T` | BM3R1 / Rep 的 mRNA | repressor mRNA |
| e | `b0_T_u` | 未成熟 BM3R1 | immature BM3R1 |
| f | `b0_T` | 成熟 BM3R1 / Rep | mature BM3R1 / Rep |
| g | `b0_M_R` | RDF 的 mRNA | RDF mRNA |
| h | `b0_R_u` | 未成熟 RDF | immature RDF |
| i | `b0_R` | 成熟游离 RDF | mature free RDF |
| j | `b0_C` | Int·RDF 复合物 | Int–RDF complex |
| k | `b0_S` | LR 构象比例（`PB = 1 − S`） | LR fraction |

| 数字 | 反应 | English |
|---|---|---|
| 1 | 转录：上游 C31 启动子活性 → `M_I` | transcription (external) |
| 2 | 翻译：`M_I` → `I_u` | translation |
| 3 | 成熟：`I_u` → `I` | maturation |
| 4 | PB 驱动 Rep 表达链 | PB drives the Rep chain |
| 5 | LR 驱动 RDF 表达链 | LR drives the RDF chain |
| 6 | T 抑制 RDF 表达 | T represses RDF expression |
| 7 | 可逆结合 `I + R ⇄ C` | reversible binding |
| 8 | 正向重组（I 激活、R 抑制） | forward recombination |
| 9 | 游离 RDF 抑制正向 | free RDF represses forward |
| 10 | 反向重组（C 驱动） | reverse recombination |

> **注意编号与 `fig_bit0_internal.png` 不同。** 那张图用的是 12 个编号
> （④⑤ 标表达链、(11)(12) 标 `k` 的两条输出）。这份提示词把它们**合并成 4 和 5**，
> 因为要少给生图模型 2 个标签。引用时别混用两套编号。

---

## 7 报告可直接用的图注

> **图 X. 单个 HBY 位（bit0）的分子级联。**
> 整合酶经 `a → b → c`（转录 → 翻译 → 成熟）生成；位内的 Rep 链 `d → e → f`
> 由 PB 驱动（4），RDF 链 `g → h → i` 由 LR 驱动（5），且 T 抑制 RDF 的表达（6）。
> 成熟整合酶 `c` 与游离 RDF `i` 以 1:1 可逆结合成复合物 `j`（7），结合同时消耗两者的游离池。
> 正向重组由 `c` 激活并被 `i` 抑制（8、9），反向重组由 `j` 驱动（10）；
> 两者共同决定构象比例 `S`（LR）与 `PB = 1 − S`（k）。
> **模型 / 输入 / 参数版本**：HBY 早期基础表参数；对应 `model_hzh.py:39-63`。示意图，非数据图。

---

## 8 如果想要"接上下游"的版本

在 §1 主提示词后面追加这一段（其余不变）：

```
ADDITIONALLY, to the left of the whole composition, outside a large dashed
rounded rectangle that encloses everything else, draw a small vertical
repression ring of three rounded rectangles labelled "A", "B", "C" (top to
bottom), each connected to the next by a short line ending in a FLAT
PERPENDICULAR BAR (repression), with the third one also pointing back to the
first. Below them place one more rounded rectangle labelled "D". A single
straight arrow runs from "D" to the free incoming arrow above box "a".

To the right of the dashed rectangle, place two rounded rectangles labelled "E"
and "F" (E above F), a small circle containing a multiplication cross "x" below
them, and one rounded rectangle labelled "G" at the bottom right.

Both "E" and "F" must have an incoming arrow. Draw all four of these:
 - ONE arrow enters the LEFT side of "E", coming from the dashed rectangle. Its
   tail starts outside "E" and does not attach to any box; place the small label
   "PB" beside that tail. This is the ONLY incoming arrow of "E".
 - ONE arrow runs straight DOWN from the bottom of "E" into the top of "F". This
   is the ONLY incoming arrow of "F". Do not omit it: without it "F" has no
   input at all and the loop reads wrong.
 - ONE arrow runs from the right of "E" down-right into the "x" circle.
 - ONE short line from the bottom of "F" points down into the "x" circle and ends
   in a FLAT PERPENDICULAR BAR.

So "E" has one arrow in and two arrows out (it forks); "F" has one in and one out.
An arrow runs from the "x" circle into "G".

The dashed rectangle is labelled only with the digit "34".
```

追加后的字母对照：`A`=TetR、`B`=CI、`C`=LacI、`D`=C31 启动子、
`E`=`A0_zmh`、`F`=`F0_zmh`、`G`=`I1_zmh`；虚线框内的 34 表示模型状态数。

> **`E` / `F` 的入边不能省。** `E` 的入边来自 bit0 的 `PB0 = 1 − b0_S`
> （`model_hzh.py:73`：`alpha_A0*(1-b0s)`，`b0s` 就是 `y[10]` = `b0_S`，
> 见 `:89` 的调用），`F` 的入边来自 `E`（`:74`：`alpha_R0*act`，其中
> `act` 是 `A0` 的 Hill 函数，见 `:68`；**代码里抑制臂叫 `R0`/`r0`，图上叫 `F0`**，
> 两者是同一个量）。这两条画漏了，`carry0` 就会读成两个互不相关的输入做 AND，
> 而不是"激活臂自建延迟抑制臂"的非相干前馈环。
> 配套的 `PB0` 网络标签见 `out/fig_bit0_wired.png` 里的 `PB0 >` / `< PB0`。

### §8 版逐条校验（拿到图后对照）

| # | 检查项 | 依据（源码） |
|---|---|---|
| 1 | `E` **有入边**：一条箭头从虚线框方向进入 `E` 左侧 | `model_hzh.py:73`（源项 `1 − b0_S`） |
| 2 | `F` **有入边**：一条箭头从 `E` 竖直向下进入 `F` | `:74`（源项 `act(A0)`） |
| 3 | `E` 是**分叉**：入 1 条、出 2 条（一条去 `F`、一条去 `×`） | `:73-76` |
| 4 | `F` 出 1 条，末端是 **T 形横杠**压在 `×` 上（不是箭头） | `:69`（`K_R0^n/(K_R0^n+f^n)`） |
| 5 | `×` → `G`，且 `G` 只有这一个入边 | `:75`（`alpha_Int1*act*gate`） |
| 6 | 上游 `A⊣B⊣C⊣A` 三条都是 **T 形横杠**，不是箭头 | 振荡器阻遏环 |

**最容易踩的坑**：只画 `E→×` 与 `F⊣×`，漏掉 `E` 的入边和 `E→F`。这样
`carry0` 会读成两个独立输入做 AND，而 `F` 完全没有任何来源——
非相干前馈环（激活臂自建延迟抑制臂）就看不出来了。

---

## 附：这份提示词的取材

上面每一条几何与箭头都对应源码里的一项 RHS，不是凭印象写的：

| 提示词里的内容 | 源码位置 |
|---|---|
| 11 个状态名与 `a`–`k` 的对应 | `model_hzh.py:10-11`（`NAMES` 前 11 项） |
| `M_I` 的转录源 | `model_hzh.py:49-51` |
| `M_I → I_u` 翻译、`I_u → I` 成熟 | `model_hzh.py:52-54` |
| Rep 链源 `α_rep·(1−S)` | `model_hzh.py:55` |
| RDF 链源 `α_rdf·S·(1−H(T))` | `model_hzh.py:56` |
| `I + R ⇄ C` 双向、两边都扣游离池 | `model_hzh.py:57-59` |
| 正向 `v_f` 含 `H(I)` 与 `K_inh/(K_inh+R)` | `model_hzh.py:60` |
| 反向 `v_r` 含 `H(C)` | `model_hzh.py:61` |
| `dS/dt = v_f(1−S) − v_r·S` | `model_hzh.py:62` |
