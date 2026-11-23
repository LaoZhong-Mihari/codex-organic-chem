"""Draw the benchmark charts used in README.md and README.zh-CN.md.

    uv run --with matplotlib python scripts/make_readme_figures.py

Writes docs/images/<chart>.<en|zh>.<light|dark>.png in a journal-figure style.
Bars are means, error bars are 95% intervals (a percentile bootstrap over items,
10,000 resamples with seed 0, for scores; Wilson score intervals for
proportions), dots are single items, and the rows under a score panel give the
share of runs that used the skill and the mean cost per run.

The numbers are copied from the results named next to each data block. The
skill-trigger, OrgReason-Real and oMeBench runs come from an internal benchmark
harness that is not part of this repository, so they cannot be rerun here; the
intervals are recomputed below from the per-item scores, with the same bootstrap
as the analysis scripts.
"""

from __future__ import annotations

import math
import random
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "docs" / "images"
DPI = 200  # 1 CSS px = 0.01 in, so a 9 in figure is 900 px wide shown at 2x
PX = 0.01  # inches per CSS px
PT = 0.72  # typographic points per CSS px at this scale

THEMES = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781",
              "grid": "#e1e0d9", "axis": "#c3c2b7", "base": "#9a988f", "s1": "#2a78d6", "s2": "#eb6834",
              "wash": 0.30},
    "dark": {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7", "muted": "#898781",
             "grid": "#2c2c2a", "axis": "#383835", "base": "#77756f", "s1": "#3987e5", "s2": "#d95926",
             "wash": 0.42},
}
FONTS = {"en": ["Helvetica Neue", "Arial"], "zh": ["Helvetica Neue", "PingFang SC", "Hiragino Sans GB"]}
# mathtext also sets the plain text of a string that contains math, and it has no
# per-glyph fallback, so the Chinese strings avoid math altogether.
MATH_FONT = {"en": "Helvetica Neue", "zh": "Helvetica Neue"}

# --- data -------------------------------------------------------------------

# Internal skill-trigger test (2026-11-18, Claude Opus 5.5 xhigh, fresh Claude
# Code sessions, two per item): sessions in which loading the skill was the first
# action. Old description vs new description + MCP server instructions.
TRIGGER = [
    ("real", (3, 18), (16, 18)),
    ("pages", (0, 8), (6, 8)),
    ("other", (21, 32), (25, 32)),
    ("near", None, (0, 24)),
]

# Benchmark runs with tools (conditions D and E together, one run per item each):
# runs in which the skill was used at any point. Olympiad: OrgReason-Real,
# 2026-11-17 (old description) and 2026-11-18. oMeBench: the 17 hard reactions
# whose reference mechanism passes the audit; Claude Opus 5.5 on 2026-11-22,
# GPT-6.2 Astra in Codex on 2026-11-21.
LOADING = [
    ("l_olymp_old", "old", (3, 18)),
    ("l_olymp", "new", (14, 18)),
    ("l_ome_opus", "new", (4, 34)),
    ("l_ome_astra", "new", (34, 34)),
]

# scripts/benchmark_route_check.py (results.md in its output directory): one
# compound corrupted at a time in three curated eight-step routes; detected =
# output differs from the clean run.
# (legacy balance + normalization tools, route-check at blocking or review level)
ROUTE = [
    ("all", (67, 282), (280, 282)),
    ("ch2", (6, 77), (77, 77)),
    ("carbon", (6, 52), (51, 52)),
    ("on", (0, 49), (49, 49)),
    ("invert", (0, 36), (35, 36)),
    ("drop", (36, 36), (36, 36)),
    ("move", (16, 27), (27, 27)),
    ("omit", (3, 3), (3, 3)),
    ("ez", (0, 2), (2, 2)),
]

# OrgReason-Real (internal): auto-graded score per problem for conditions
# (A, D, E), one run per problem and condition, Claude Opus 5.5 xhigh.
# A from 2026-11-17, D and E from 2026-11-18.
OLYMP = {
    "CCHO40-P4": (0.2, 0.4, 1.0),
    "ICHO26-T6": (0.9036, 0.9036, 0.9518),
    "CCHO40-P9": (0.8889, 0.9444, 0.9444),
    "CCHO40-P8": (0.9091, 0.9091, 0.9091),
    "CCHO40-P10": (1.0, 1.0, 1.0),
    "ICHO26-T1": (1.0, 1.0, 1.0),
    "ICHO26-T5": (1.0, 1.0, 1.0),
    "ICHO26-T8a": (1.0, 1.0, 1.0),
    "ICHO26-T9": (1.0, 1.0, 1.0),
}
OLYMP_SKILL = (None, (7, 9), (7, 9))
OLYMP_COST = (0.96, 1.41, 1.08)  # USD per run (Claude Code total_cost_usd)

# oMeBench (internal harness, public data): S_total per reaction on the 17 hard
# reactions whose reference mechanism passes the audit, conditions (A, D, E), one
# run each. First tuple Claude Opus 5.5 xhigh (2026-11-22), second GPT-6.2 Astra
# xhigh in Codex (2026-11-21).
OME = {
    "ART-005": ((0.06, 0.0, 0.0), (0.12, 0.05, 0.12)),
    "NG-012": ((0.84, 0.57, 0.18), (0.84, 0.84, 0.84)),
    "NR-043": ((0.08, 0.7, 0.7), (0.08, 0.08, 0.32)),
    "NR-056": ((1.0, 1.0, 1.0), (1.0, 0.76, 0.76)),
    "NR-071": ((0.62, 0.62, 0.62), (0.62, 1.0, 0.62)),
    "NR-204": ((0.78, 0.56, 1.0), (0.4, 0.4, 0.56)),
    "NR-211": ((0.54, 0.32, 0.32), (0.0, 0.0, 0.0)),
    "NR-220": ((0.38, 0.64, 0.64), (0.38, 0.38, 0.64)),
    "NR-221": ((0.23, 0.23, 0.23), (0.23, 0.23, 0.0)),
    "NR-224": ((0.38, 0.38, 0.38), (0.38, 0.38, 0.38)),
    "NR-227": ((0.32, 0.32, 0.32), (0.32, 0.32, 0.32)),
    "NR-242": ((1.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
    "NR-300": ((0.58, 0.58, 0.58), (0.58, 0.58, 0.58)),
    "NR-307": ((0.65, 0.81, 0.71), (0.77, 0.77, 0.77)),
    "NR-311": ((0.35, 0.35, 0.35), (0.0, 0.0, 0.05)),
    "NR-318": ((0.4, 0.4, 0.22), (0.4, 0.22, 0.4)),
    "NR-320": ((0.25, 0.19, 0.19), (0.0, 0.0, 0.0)),
}
OME_SKILL = {"opus": (None, (3, 17), (1, 17)), "astra": (None, (17, 17), (17, 17))}
OME_COST = {"opus": (0.113, 0.171, 0.158), "astra": None}  # Codex CLI reports no cost

TEXT = {
    "en": {
        "old": "Old skill description", "new": "Current skill description",
        "load_title": "How often the model uses the skill",
        "load_a": "Trigger test, Claude Opus 5.5: skill loaded as the first action",
        "load_b": "Benchmark runs with tools (D and E): skill used at any point",
        "real": "OrgReason-Real olympiad items (9)", "pages": "Olympiad pages written for the test (4)",
        "other": "Other organic questions (16)", "near": "Near misses that should not load it (12)",
        "l_olymp_old": "Olympiad · Claude Opus 5.5, Nov 17", "l_olymp": "Olympiad · Claude Opus 5.5, Nov 18",
        "l_ome_opus": "oMeBench · Claude Opus 5.5", "l_ome_astra": "oMeBench · GPT-6.2 Astra in Codex",
        "not_run": "not run", "lower": "lower is better",
        "load_note": ("Error bars: 95% Wilson intervals. a: two fresh Claude Code sessions per item, 2026-11-18. "
                      "b: one run per item in each of D and E; the old description was replaced on Nov 18."),
        "route_title": "Does route-check notice a single planted error?",
        "route_sub": "282 single-compound corruptions of three curated eight-step routes: share that changes the output",
        "legacy": "Earlier checks (balance, normalization)", "audit": "route-check (blocking or review)",
        "all": "All corruptions (282)", "ch2": "Insert one CH₂ (77)", "carbon": "Delete one carbon (52)",
        "on": "Swap O and N (49)", "invert": "Invert one stereocentre (36)", "drop": "Drop one stereo label (36)",
        "move": "Move a substituent to the next carbon (27)", "omit": "Omit the last step (3)", "ez": "Flip E/Z (2)",
        "route_note": ("Error bars: 95% Wilson intervals. On the three clean routes: no blocking issue and 0-3 review "
                       "prompts per route; the earlier checks called all 24 correct steps unbalanced. "
                       "Measures the tool, not a model."),
        "cA": "A  no tools", "cD": "D  MCP tools + skill", "cE": "E  D + Python tool",
        "mean_row": "Mean", "skill_row": "Skill used", "cost_row": "Cost per run", "na": "—",
        "not_reported": "not reported",
        "olymp_title": "OrgReason-Real olympiad problems (IChO 2026, CChO 2026), Claude Opus 5.5 (xhigh)",
        "olymp_a": "Mean auto-graded score (n = 9 problems)", "olymp_b": "Score per problem",
        "score": "Auto-graded score",
        "olymp_note": ("Bars: mean. Error bars: 95% bootstrap interval over problems. Dots: problems. One run per "
                       "problem and condition (A 2026-11-17, D and E 2026-11-18); gold answers are hand-transcribed "
                       "and not yet signed off by a chemist."),
        "ome_title": "oMeBench mechanisms: 17 hard reactions whose reference mechanism passes the audit",
        "ome_opus": "Claude Opus 5.5 (xhigh)", "ome_astra": "GPT-6.2 Astra (xhigh) in Codex",
        "stotal": "Exact-match score, $S_\\mathrm{total}$",
        "ome_note": ("Bars: mean. Error bars: 95% bootstrap interval over reactions. Dots: reactions. One run per "
                     "reaction and condition; Opus on 2026-11-22, Astra on 2026-11-21. Codex reports no cost."),
        "eff_title": "Change in score with the tools, compared with A on the same items",
        "eff_sub": "Dots: mean paired difference. Bars: 95% bootstrap interval over items.",
        "g_olymp": "Olympiad · Claude Opus 5.5 (n = 9)", "g_ome_opus": "oMeBench · Claude Opus 5.5 (n = 17)",
        "g_ome_astra": "oMeBench · GPT-6.2 Astra (n = 17)",
        "dA": "D − A", "eA": "E − A", "diff_axis": "Difference in mean score (scores run from 0 to 1)",
        "delta_col": "Difference [95% interval]",
        "eff_note": ("Olympiad: auto-graded score; oMeBench: $S_\\mathrm{total}$. One run per item and condition. "
                     "Most of the olympiad E − A gain comes from one problem, CCHO40-P4."),
    },
    "zh": {
        "old": "旧版 skill 描述", "new": "现行 skill 描述",
        "load_title": "模型多常用到 skill",
        "load_a": "触发测试，Claude Opus 5.5：第一步就加载 skill",
        "load_b": "带工具的正式测评（D 组和 E 组）：答题过程中用到 skill",
        "real": "OrgReason-Real 竞赛真题（9 道）", "pages": "为测试自编的竞赛题图（4 道）",
        "other": "其他有机化学问题（16 道）", "near": "不该加载的近似问题（12 道）",
        "l_olymp_old": "竞赛题 · Claude Opus 5.5，11 月 17 日", "l_olymp": "竞赛题 · Claude Opus 5.5，11 月 18 日",
        "l_ome_opus": "oMeBench · Claude Opus 5.5", "l_ome_astra": "oMeBench · Codex 中的 GPT-6.2 Astra",
        "not_run": "未测", "lower": "越低越好",
        "load_note": ("误差线：95% Wilson 区间。a：每题两个全新的 Claude Code 会话，2026-11-18。"
                      "b：D、E 两组每题各运行一次；11 月 18 日起换用现行描述。"),
        "route_title": "路线里埋一个错误，route-check 能不能发现？",
        "route_sub": "三条精选的 8 步路线，每次只改坏一个化合物，共 282 例：输出发生变化的比例",
        "legacy": "原有检查（配平、规范化）", "audit": "route-check（blocking 或 review）",
        "all": "全部（282）", "ch2": "插入一个 CH₂（77）", "carbon": "删掉一个碳（52）",
        "on": "O、N 互换（49）", "invert": "翻转一个立体中心（36）", "drop": "去掉一个立体标记（36）",
        "move": "取代基挪到相邻碳上（27）", "omit": "漏掉最后一步（3）", "ez": "E/Z 翻转（2）",
        "route_note": ("误差线：95% Wilson 区间。三条正确路线上：没有 blocking，每条 0-3 个 review 提示；"
                       "原有检查把 24 个正确步骤全部判为不配平。测的是工具的区分能力，不是模型。"),
        "cA": "A  无工具", "cD": "D  MCP 工具 + skill", "cE": "E  D + Python 工具",
        "mean_row": "均值", "skill_row": "用到 skill", "cost_row": "每次费用", "na": "—", "not_reported": "未报告",
        "olymp_title": "OrgReason-Real 竞赛题（IChO 2026、CChO 2026），Claude Opus 5.5（xhigh）",
        "olymp_a": "自动评分均值（9 道题）", "olymp_b": "每道题得分",
        "score": "自动评分",
        "olymp_note": ("柱高：均值。误差线：按题自助法的 95% 区间。圆点：单道题。每题每组运行一次"
                       "（A 组 2026-11-17，D、E 组 2026-11-18）；标准答案为人工转写，尚未经化学专家确认。"),
        "ome_title": "oMeBench 机理题：参考机理通过核查的 17 道 hard 题",
        "ome_opus": "Claude Opus 5.5（xhigh）", "ome_astra": "Codex 中的 GPT-6.2 Astra（xhigh）",
        "stotal": "完全匹配得分",
        "ome_note": ("柱高：均值。误差线：按反应自助法的 95% 区间。圆点：单个反应。每个反应每组运行一次；"
                     "Opus 运行于 2026-11-22，Astra 运行于 2026-11-21。Codex 不报告费用。"),
        "eff_title": "用上工具后得分的变化（与同一批题的 A 组相比）",
        "eff_sub": "圆点：按题配对的平均差。横线：按题自助法的 95% 区间。",
        "g_olymp": "竞赛题 · Claude Opus 5.5（9 道）", "g_ome_opus": "oMeBench · Claude Opus 5.5（17 道）",
        "g_ome_astra": "oMeBench · GPT-6.2 Astra（17 道）",
        "dA": "D − A", "eA": "E − A", "diff_axis": "平均分之差（得分范围 0 到 1）",
        "delta_col": "差值 [95% 区间]",
        "eff_note": ("竞赛题为自动评分，oMeBench 为完全匹配得分。每题每组运行一次。"
                     "竞赛题 E − A 的增益大部分来自 CCHO40-P4 一道题。"),
    },
}

# --- statistics -----------------------------------------------------------------


def boot_ci(values, reps=10000, seed=0):
    """95% percentile bootstrap interval of the mean (the analysis scripts use the same procedure)."""
    rng = random.Random(seed)
    n = len(values)
    means = sorted(sum(values[rng.randrange(n)] for _ in range(n)) / n for _ in range(reps))
    return means[int(0.025 * reps)], means[int(0.975 * reps) - 1]


def wilson(k, n, z=1.959964):
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def rate_label(pair):
    k, n = pair
    return f"{k}/{n} · {round(100 * k / n)}%"


# --- drawing helpers ----------------------------------------------------------


def setup(lang: str, theme: dict) -> None:
    plt.rcParams.update({
        "font.family": FONTS[lang], "font.size": 9.0,
        "axes.edgecolor": theme["axis"], "axes.linewidth": 0.8, "axes.facecolor": theme["surface"],
        "text.color": theme["ink"], "axes.labelcolor": theme["ink2"],
        "xtick.color": theme["axis"], "ytick.color": theme["axis"],
        "xtick.labelcolor": theme["ink2"], "ytick.labelcolor": theme["ink2"],
        "xtick.major.size": 3, "ytick.major.size": 3, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.direction": "out", "ytick.direction": "out",
        "mathtext.fontset": "custom", "mathtext.rm": MATH_FONT[lang], "mathtext.it": f"{MATH_FONT[lang]}:italic",
        "mathtext.bf": f"{MATH_FONT[lang]}:bold", "svg.fonttype": "none",
    })


def frame(width_px: float, height_px: float, theme: dict):
    return plt.figure(figsize=(width_px * PX, height_px * PX), dpi=DPI, facecolor=theme["surface"])


def axes_px(fig, W, H, left, top, width, height):
    """Axes from a CSS-px box measured from the top-left corner of the figure."""
    return fig.add_axes([left / W, (H - top - height) / H, width / W, height / H])


def title(fig, W, H, text, sub=None):
    fig.text(24 / W, 1 - 16 / H, text, fontsize=11.2, fontweight="semibold", va="top")
    if sub:
        fig.text(24 / W, 1 - 38 / H, sub, fontsize=8.8, color=plt.rcParams["axes.labelcolor"], va="top")


def legend(fig, W, H, items, top):
    handles = [Patch(facecolor=color, edgecolor="none", label=label) for label, color in items]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(18 / W, 1 - top / H), ncol=len(handles),
               frameon=False, fontsize=8.8, handlelength=0.9, handleheight=0.9, handletextpad=0.5,
               columnspacing=1.8, labelcolor=plt.rcParams["text.color"])


def heading(ax, letter, text, theme, dx_px=0.0, dy_px=12.0):
    """Panel label: bold letter, then the panel's title, above the axes' top-left corner."""
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(dx_px * PT, dy_px * PT),
                textcoords="offset points", fontsize=10.5, fontweight="bold", color=theme["ink"], va="bottom")
    ax.annotate(text, xy=(0, 1), xycoords="axes fraction", xytext=((dx_px + 16) * PT, dy_px * PT),
                textcoords="offset points", fontsize=9.2, color=theme["ink"], va="bottom")


def footnote(fig, W, H, theme, text):
    width = int(W / 5.9) if " " in text else int(W / 11.5)
    wrapped = "\n".join(textwrap.wrap(text, width=width))
    fig.text(24 / W, 12 / H, wrapped, fontsize=7.8, color=theme["muted"], va="bottom", linespacing=1.5)


def style_h(ax, theme, xmax=1.0):
    """Horizontal bars: categories on the left, a 0-100% axis at the bottom."""
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_xlim(0, xmax)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25%", "50%", "75%", "100%"], fontsize=8.2)
    ax.tick_params(axis="y", length=0, pad=8, labelsize=8.8, labelcolor=theme["ink"])
    ax.grid(axis="x", color=theme["grid"], linewidth=0.6)
    ax.set_axisbelow(True)


def style_v(ax, theme, ylabel, ymax=1.06):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_ylim(0, ymax)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1.00"], fontsize=8.2)
    ax.set_ylabel(ylabel, fontsize=8.8, labelpad=6)
    ax.grid(axis="y", color=theme["grid"], linewidth=0.6)
    ax.set_axisbelow(True)


def hbar_rate(ax, y, pair, color, thick, theme, suffix=""):
    """A proportion as a horizontal bar with its Wilson interval and a k/n label past the interval."""
    k, n = pair
    rate = k / n
    lo, hi = wilson(k, n)
    if rate > 0:
        ax.barh(y, rate, height=thick, color=color, edgecolor="none", zorder=2)
    ax.errorbar(rate, y, xerr=[[max(rate - lo, 0)], [max(hi - rate, 0)]], fmt="none", ecolor=theme["ink"], elinewidth=0.9,
                capsize=2.4, capthick=0.9, zorder=3)
    label = rate_label(pair) + (f"  ({suffix})" if suffix else "")
    ax.text(hi + 0.014, y, label, va="center", ha="left", fontsize=8.0, color=theme["ink2"])


def paired_rates(ax, rows, tx, theme, colors, suffix=None):
    n = len(rows)
    ax.set_ylim(n - 0.5, -0.5)
    thick, gap = 0.3, 0.07
    for i, (key, first, second) in enumerate(rows):
        for j, (pair, color) in enumerate(((first, colors[0]), (second, colors[1]))):
            y = i - (thick + gap) / 2 if j == 0 else i + (thick + gap) / 2
            if pair is None:
                ax.text(0.01, y, tx["not_run"], va="center", ha="left", fontsize=8.0, color=theme["muted"])
                continue
            hbar_rate(ax, y, pair, color, thick, theme, (suffix or {}).get(key, "") if j == 1 else "")
    ax.set_yticks(range(n))
    ax.set_yticklabels([tx[key] for key, *_ in rows])


def score_bars(ax, groups, theme, seed=11):
    """Mean as a washed bar with a solid cap, 95% bootstrap interval, one dot per item."""
    rng = random.Random(seed)
    means = []
    for i, (values, color) in enumerate(groups):
        mean = sum(values) / len(values)
        lo, hi = boot_ci(list(values))
        ax.bar(i, mean, width=0.62, color=color, alpha=theme["wash"], edgecolor="none", zorder=1)
        ax.plot([i - 0.31, i + 0.31], [mean, mean], color=color, linewidth=2.2, solid_capstyle="butt", zorder=2)
        xs = [i + rng.uniform(-0.2, 0.2) for _ in values]
        ax.scatter(xs, values, s=15, color=color, edgecolors=theme["surface"], linewidths=0.7, zorder=3)
        ax.errorbar(i, mean, yerr=[[max(mean - lo, 0)], [max(hi - mean, 0)]], fmt="none", ecolor=theme["ink"], elinewidth=1.1,
                    capsize=4, capthick=1.1, zorder=4)
        means.append(mean)
    ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(["A", "D", "E"][: len(groups)], fontsize=9.0)
    ax.tick_params(axis="x", length=0, labelcolor=theme["ink"])
    return means


def annotation_rows(ax, theme, rows, row_labels=True, first_px=30, step_px=15):
    """Small table under a bar panel: one line per row, one cell per bar."""
    for r, (label, cells) in enumerate(rows):
        dy = -(first_px + r * step_px) * PT
        if row_labels:
            ax.annotate(label, xy=(0, 0), xycoords="axes fraction", xytext=(-8 * PT, dy), textcoords="offset points",
                        ha="right", va="center", fontsize=7.9, color=theme["muted"])
        for i, cell in enumerate(cells):
            ax.annotate(cell, xy=(i, 0), xycoords=("data", "axes fraction"), xytext=(0, dy),
                        textcoords="offset points", ha="center", va="center", fontsize=7.9, color=theme["ink2"])


def skill_cells(tx, skills):
    return [tx["na"] if pair is None else f"{pair[0]}/{pair[1]} · {round(100 * pair[0] / pair[1])}%"
            for pair in skills]


def cost_cells(tx, costs, n):
    return [tx["not_reported"]] * n if costs is None else [f"${c:.2f}" for c in costs]


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, dpi=DPI, facecolor=fig.get_facecolor())
    plt.close(fig)


# --- charts -------------------------------------------------------------------


def chart_loading(lang, tname):
    theme, tx = THEMES[tname], TEXT[lang]
    setup(lang, theme)
    W, H = 900, 640
    fig = frame(W, H, theme)
    title(fig, W, H, tx["load_title"])
    legend(fig, W, H, [(tx["old"], theme["base"]), (tx["new"], theme["s1"])], 40)
    left, width = 320, 470
    ax_a = axes_px(fig, W, H, left, 118, width, 228)
    paired_rates(ax_a, TRIGGER, tx, theme, (theme["base"], theme["s1"]), suffix={"near": tx["lower"]})
    style_h(ax_a, theme)
    heading(ax_a, "a", tx["load_a"], theme, dx_px=-(left - 24))

    ax_b = axes_px(fig, W, H, left, 420, width, 132)
    n = len(LOADING)
    ax_b.set_ylim(n - 0.5, -0.5)
    for i, (key, config, pair) in enumerate(LOADING):
        hbar_rate(ax_b, i, pair, theme["base"] if config == "old" else theme["s1"], 0.42, theme)
    ax_b.set_yticks(range(n))
    ax_b.set_yticklabels([tx[key] for key, *_ in LOADING])
    style_h(ax_b, theme)
    heading(ax_b, "b", tx["load_b"], theme, dx_px=-(left - 24))
    footnote(fig, W, H, theme, tx["load_note"])
    save(fig, f"skill-loading.{lang}.{tname}.png")


def chart_route(lang, tname):
    theme, tx = THEMES[tname], TEXT[lang]
    setup(lang, theme)
    W, H = 900, 660
    fig = frame(W, H, theme)
    title(fig, W, H, tx["route_title"], tx["route_sub"])
    legend(fig, W, H, [(tx["legacy"], theme["base"]), (tx["audit"], theme["s1"])], 60)
    ax = axes_px(fig, W, H, 320, 104, 470, 470)
    paired_rates(ax, ROUTE, tx, theme, (theme["base"], theme["s1"]))
    style_h(ax, theme)
    ax.axhline(0.5, color=theme["axis"], linewidth=0.7)
    ax.get_yticklabels()[0].set_fontweight("semibold")
    footnote(fig, W, H, theme, tx["route_note"])
    save(fig, f"route-audit.{lang}.{tname}.png")


def chart_olymp(lang, tname):
    theme, tx = THEMES[tname], TEXT[lang]
    setup(lang, theme)
    W, H = 900, 580
    fig = frame(W, H, theme)
    colors = (theme["base"], theme["s1"], theme["s2"])
    title(fig, W, H, tx["olymp_title"])
    legend(fig, W, H, list(zip((tx["cA"], tx["cD"], tx["cE"]), colors)), 40)

    ax = axes_px(fig, W, H, 112, 116, 290, 300)
    groups = [([v[c] for v in OLYMP.values()], colors[c]) for c in range(3)]
    means = score_bars(ax, groups, theme)
    style_v(ax, theme, tx["score"])
    annotation_rows(ax, theme, [(tx["mean_row"], [f"{m:.3f}" for m in means]),
                                (tx["skill_row"], skill_cells(tx, OLYMP_SKILL)),
                                (tx["cost_row"], cost_cells(tx, OLYMP_COST, 3))])
    heading(ax, "a", tx["olymp_a"], theme, dx_px=-88)

    axd = axes_px(fig, W, H, 560, 116, 300, 330)
    items = list(OLYMP)
    axd.set_ylim(len(items) - 0.5, -0.5)
    for side in ("top", "right", "left"):
        axd.spines[side].set_visible(False)
    axd.set_xlim(0, 1.04)
    axd.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    axd.set_xticklabels(["0", "0.25", "0.50", "0.75", "1.00"], fontsize=8.2)
    axd.set_yticks(range(len(items)))
    axd.set_yticklabels(items, fontsize=8.4)
    axd.tick_params(axis="y", length=0, pad=8, labelcolor=theme["ink"])
    axd.grid(axis="both", color=theme["grid"], linewidth=0.6)
    axd.set_axisbelow(True)
    for i, key in enumerate(items):
        for j, (v, c) in enumerate(zip(OLYMP[key], colors)):
            axd.plot([v], [i + (j - 1) * 0.24], marker="o", markersize=5.6, markerfacecolor=c,
                     markeredgecolor=theme["surface"], markeredgewidth=1.1, linestyle="", zorder=3)
    axd.set_xlabel(tx["score"], fontsize=8.8, labelpad=6)
    heading(axd, "b", tx["olymp_b"], theme, dx_px=-96)
    footnote(fig, W, H, theme, tx["olymp_note"])
    save(fig, f"orgreason-real.{lang}.{tname}.png")


def chart_omebench(lang, tname):
    theme, tx = THEMES[tname], TEXT[lang]
    setup(lang, theme)
    W, H = 900, 560
    fig = frame(W, H, theme)
    colors = (theme["base"], theme["s1"], theme["s2"])
    title(fig, W, H, tx["ome_title"])
    legend(fig, W, H, list(zip((tx["cA"], tx["cD"], tx["cE"]), colors)), 40)
    for p, (model, left) in enumerate((("opus", 132), ("astra", 540))):
        ax = axes_px(fig, W, H, left, 116, 300, 290)
        groups = [([v[p][c] for v in OME.values()], colors[c]) for c in range(3)]
        means = score_bars(ax, groups, theme, seed=13 + p)
        style_v(ax, theme, tx["stotal"] if p == 0 else "")
        annotation_rows(ax, theme, [(tx["mean_row"], [f"{m:.3f}" for m in means]),
                                    (tx["skill_row"], skill_cells(tx, OME_SKILL[model])),
                                    (tx["cost_row"], cost_cells(tx, OME_COST[model], 3))])
        heading(ax, "ab"[p], tx[f"ome_{model}"], theme, dx_px=-108 if p == 0 else -96)
    footnote(fig, W, H, theme, tx["ome_note"])
    save(fig, f"omebench.{lang}.{tname}.png")


def signed(x: float, digits: int) -> str:
    return f"{x:+.{digits}f}".replace("-", "−")


def chart_effect(lang, tname):
    theme, tx = THEMES[tname], TEXT[lang]
    setup(lang, theme)
    W, H = 900, 520
    fig = frame(W, H, theme)
    title(fig, W, H, tx["eff_title"], tx["eff_sub"])
    groups = [
        ("g_olymp", list(OLYMP.values()), OLYMP_SKILL, OLYMP_COST),
        ("g_ome_opus", [v[0] for v in OME.values()], OME_SKILL["opus"], OME_COST["opus"]),
        ("g_ome_astra", [v[1] for v in OME.values()], OME_SKILL["astra"], OME_COST["astra"]),
    ]
    left, width, top, height = 250, 250, 118, 290
    columns = (528, 688, 800)  # CSS px from the figure's left edge: interval, skill, cost
    ax = axes_px(fig, W, H, left, top, width, height)

    def at(x_px, y, text, **kw):
        ax.annotate(text, xy=(0, y), xycoords=("axes fraction", "data"), xytext=((x_px - left) * PT, 0),
                    textcoords="offset points", va="center", **kw)

    y = 0.0
    for key, rows, skills, costs in groups:
        at(24, y, tx[key], ha="left", fontsize=8.8, fontweight="semibold", color=theme["ink"])
        for cond, label_key, color in ((1, "dA", theme["s1"]), (2, "eA", theme["s2"])):
            y += 1
            d = [r[cond] - r[0] for r in rows]
            mean = sum(d) / len(d)
            lo, hi = boot_ci(d)
            ax.plot([lo, hi], [y, y], color=color, linewidth=1.6, solid_capstyle="butt", zorder=2)
            for x in (lo, hi):
                ax.plot([x, x], [y - 0.16, y + 0.16], color=color, linewidth=1.6, zorder=2)
            ax.plot([mean], [y], marker="o", markersize=6.4, markerfacecolor=color, markeredgecolor=theme["surface"],
                    markeredgewidth=1.2, linestyle="", zorder=3)
            at(44, y, tx[label_key], ha="left", fontsize=8.6, color=theme["ink"])
            at(columns[0], y, f"{signed(mean, 3)}  [{signed(lo, 2)}, {signed(hi, 2)}]", ha="left", fontsize=7.9,
               color=theme["ink2"])
            pair = skills[cond]
            at(columns[1], y, f"{pair[0]}/{pair[1]} · {round(100 * pair[0] / pair[1])}%", ha="left", fontsize=7.9,
               color=theme["ink2"])
            at(columns[2], y, tx["not_reported"] if costs is None else f"${costs[cond]:.2f}", ha="left", fontsize=7.9,
               color=theme["ink2"])
        y += 1.0
    ax.set_ylim(y - 0.4, -0.6)
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_xlim(-0.3, 0.4)
    ticks = [-0.2, -0.1, 0, 0.1, 0.2, 0.3]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["0" if t == 0 else signed(t, 1) for t in ticks], fontsize=8.2)
    ax.axvline(0, color=theme["ink2"], linewidth=0.9, zorder=1)
    ax.grid(axis="x", color=theme["grid"], linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel(tx["diff_axis"], fontsize=8.8, labelpad=6)
    for x_px, text in zip(columns, (tx["delta_col"], tx["skill_row"], tx["cost_row"])):
        ax.annotate(text, xy=(0, 1), xycoords="axes fraction", xytext=((x_px - left) * PT, 10 * PT),
                    textcoords="offset points", ha="left", va="bottom", fontsize=7.9, color=theme["muted"])
    footnote(fig, W, H, theme, tx["eff_note"])
    save(fig, f"tool-effect.{lang}.{tname}.png")


def main() -> None:
    for lang in ("en", "zh"):
        for tname in ("light", "dark"):
            chart_loading(lang, tname)
            chart_route(lang, tname)
            chart_olymp(lang, tname)
            chart_omebench(lang, tname)
            chart_effect(lang, tname)
    print(f"wrote charts to {OUT}")


if __name__ == "__main__":
    main()
