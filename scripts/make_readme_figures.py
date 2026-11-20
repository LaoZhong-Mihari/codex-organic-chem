"""Draw the benchmark charts used in README.md and README.zh-CN.md.

    uv run --with matplotlib python scripts/make_readme_figures.py

Writes docs/images/<chart>.<en|zh>.<light|dark>.png. The numbers are copied
from the results named next to each data block. The skill-trigger and
OrgReason-Real numbers come from an internal benchmark harness that is not part
of this repository, so they cannot be recomputed here.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "docs" / "images"
DPI = 200  # 1 CSS px = 0.01 in, so a 9 in figure is 900 px wide shown at 2x
PX = 0.01  # inches per CSS px

THEMES = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781",
              "grid": "#e1e0d9", "axis": "#c3c2b7", "base": "#9a988f", "s1": "#2a78d6", "s2": "#eb6834"},
    "dark": {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7", "muted": "#898781",
             "grid": "#2c2c2a", "axis": "#383835", "base": "#77756f", "s1": "#3987e5", "s2": "#d95926"},
}
FONTS = {"en": ["Helvetica Neue", "Arial"], "zh": ["Helvetica Neue", "PingFang SC", "Hiragino Sans GB"]}

# --- data -------------------------------------------------------------------

# Internal skill-trigger test (2026-11-18, Opus 5.5, xhigh): sessions in which
# loading the skill was the first action; "new" = new description + MCP server
# instructions. Last row: OrgReason-Real runs of 2026-11-17 (old) and
# 2026-11-18 (new), conditions D and E, skill used at any point.
TRIGGER = [
    ("real", (3, 18), (16, 18)),
    ("pages", (0, 8), (6, 8)),
    ("other", (21, 32), (25, 32)),
    ("near", None, (0, 24)),
    ("bench", (3, 18), (14, 18)),
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

# OrgReason-Real (internal): auto-graded score per item, one run per item and
# condition, Claude Opus 5.5 xhigh. A from 2026-11-17, D and E from 2026-11-18.
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
OLYMP_MEAN = (0.878, 0.906, 0.978)
OLYMP_SKILL = ("0/9", "7/9", "7/9")

TEXT = {
    "en": {
        "trigger_title": "How often Claude Opus 5.5 loads the skill on its own",
        "trigger_sub": "Fresh Claude Code sessions with the skill and the MCP tools available",
        "old": "Old description", "new": "New description + MCP server instructions",
        "real": "OrgReason-Real olympiad items (9)", "pages": "Olympiad pages written for the test (4)",
        "other": "Other organic questions (16)", "near": "Near misses that should not load it (12)",
        "bench": "Benchmark runs, conditions D + E",
        "not_run": "not run", "lower": "lower is better",
        "trigger_note": ("Top four rows: skill loaded as the first action (trigger test, 2026-11-18). "
                         "Last row: skill used at any point in full benchmark runs, 2026-11-17 vs 2026-11-18."),
        "route_title": "Does the check notice a planted error in a route?",
        "route_sub": "282 single-compound corruptions of three curated eight-step routes",
        "legacy": "Earlier checks (balance, normalization)", "audit": "route-check (blocking or review)",
        "all": "All corruptions (282)", "ch2": "Insert one CH₂ (77)", "carbon": "Delete one carbon (52)",
        "on": "Swap O and N (49)", "invert": "Invert one stereocentre (36)", "drop": "Drop one stereo label (36)",
        "move": "Move a substituent to the next carbon (27)", "omit": "Omit the last step (3)", "ez": "Flip E/Z (2)",
        "route_note": ("On the three clean routes: no blocking issue and 0-3 review prompts per route; the earlier "
                       "checks called all 24 correct steps unbalanced. Measures the tool, not a model."),
        "olymp_title": "OrgReason-Real olympiad track: IChO 2026 + CChO 2026, 9 problems",
        "olymp_sub": "Auto-graded score, one run per problem and condition, Claude Opus 5.5 (xhigh)",
        "cA": "A  no tools", "cD": "D  MCP tools + skill", "cE": "E  D + Python tool",
        "mean": "Mean score", "per_item": "Score per problem", "skill_used": "skill used {}",
        "olymp_note": ("Preliminary: one run each (A 2026-11-17, D/E 2026-11-18); gold answers are hand-transcribed "
                       "and not yet signed off by a chemist. A scored 0.923 in an earlier run."),
    },
    "zh": {
        "trigger_title": "Claude Opus 5.5 自己加载 skill 的频率",
        "trigger_sub": "全新的 Claude Code 会话，skill 和 MCP 工具都可用",
        "old": "旧描述", "new": "新描述 + MCP 服务器说明",
        "real": "OrgReason-Real 竞赛真题（9 道）", "pages": "为测试自编的竞赛题图（4 道）",
        "other": "其他有机化学问题（16 道）", "near": "不该加载的近似问题（12 道）",
        "bench": "正式测评 D + E 组",
        "not_run": "未测", "lower": "越低越好",
        "trigger_note": ("前四行：第一步就加载 skill 的会话（触发测试，2026-11-18）。"
                         "最后一行：完整测评中是否用到 skill，2026-11-17 对比 2026-11-18。"),
        "route_title": "路线里埋一个错误，核查能不能发现？",
        "route_sub": "三条精选的 8 步路线，每次只改坏一个化合物，共 282 例",
        "legacy": "原有检查（配平、规范化）", "audit": "route-check（blocking 或 review）",
        "all": "全部（282）", "ch2": "插入一个 CH₂（77）", "carbon": "删掉一个碳（52）",
        "on": "O、N 互换（49）", "invert": "翻转一个立体中心（36）", "drop": "去掉一个立体标记（36）",
        "move": "取代基挪到相邻碳上（27）", "omit": "漏掉最后一步（3）", "ez": "E/Z 翻转（2）",
        "route_note": ("三条正确路线上：没有 blocking，每条 0-3 个 review 提示；原有检查把 24 个正确步骤"
                       "全部判为不配平。测的是工具的区分能力，不是模型。"),
        "olymp_title": "OrgReason-Real 竞赛题：IChO 2026 + CChO 2026，共 9 道",
        "olymp_sub": "自动评分，每道题每个条件跑 1 次，Claude Opus 5.5（xhigh）",
        "cA": "A  无工具", "cD": "D  MCP 工具 + skill", "cE": "E  D + Python 工具",
        "mean": "平均分", "per_item": "每道题得分", "skill_used": "用到 skill {}",
        "olymp_note": ("初步结果：每组只跑 1 次（A 为 2026-11-17，D/E 为 2026-11-18）；标准答案为人工转写，"
                       "尚未经化学专家签字确认。A 组在更早的一次运行中得 0.923。"),
    },
}

# --- drawing helpers ----------------------------------------------------------


def setup(lang: str, theme: dict) -> None:
    plt.rcParams.update({
        "font.family": FONTS[lang], "font.size": 9.4, "axes.edgecolor": theme["axis"],
        "text.color": theme["ink"], "axes.labelcolor": theme["ink2"],
        "xtick.color": theme["muted"], "ytick.color": theme["ink"], "svg.fonttype": "none",
    })


def frame(width_px: float, height_px: float, theme: dict):
    fig = plt.figure(figsize=(width_px * PX, height_px * PX), dpi=DPI, facecolor=theme["surface"])
    return fig


def header(fig, theme, title, subtitle, legend, top_px, width_px, height_px):
    y = 1 - top_px / height_px
    fig.text(24 / width_px, y, title, fontsize=11.6, fontweight="semibold", color=theme["ink"], va="top")
    fig.text(24 / width_px, y - 24 / height_px, subtitle, fontsize=9.4, color=theme["ink2"], va="top")
    if legend:
        handles = [Line2D([], [], marker="s", linestyle="", markersize=7, markerfacecolor=c,
                          markeredgecolor=c, label=lab) for lab, c in legend]
        fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(16 / width_px, y - 44 / height_px),
                   ncol=len(handles), frameon=False, fontsize=9.4, handletextpad=0.4, columnspacing=1.6,
                   labelcolor=theme["ink"])


def footnote(fig, theme, text, width_px, height_px):
    import textwrap
    wrapped = "\n".join(textwrap.wrap(text, width=int(width_px / 6.1))) if " " in text else "\n".join(
        textwrap.wrap(text, width=int(width_px / 11.5)))
    fig.text(24 / width_px, 14 / height_px, wrapped, fontsize=8.2, color=theme["muted"], va="bottom",
             linespacing=1.5)


def hbar(ax, y, length, thick, color, x_per_in, y_per_in):
    """Horizontal bar from x=0 with a 4px rounded data end and a square baseline end."""
    if length <= 0:
        return
    r_x = 4 * PX * x_per_in
    r_x = min(r_x, length / 2)
    aspect = y_per_in / x_per_in
    ax.add_patch(FancyBboxPatch((0, y - thick / 2), length, thick, boxstyle=f"round,pad=0,rounding_size={r_x}",
                                mutation_aspect=aspect, facecolor=color, edgecolor="none", linewidth=0))
    ax.add_patch(Rectangle((0, y - thick / 2), max(length - r_x, 0), thick, facecolor=color, edgecolor="none"))


def style_axis(ax, theme, xmax, ticks, fmt):
    ax.set_facecolor(theme["surface"])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(theme["axis"])
    ax.spines["bottom"].set_linewidth(0.6)
    ax.set_xlim(0, xmax)
    ax.set_xticks(ticks)
    ax.set_xticklabels([fmt(t) for t in ticks], fontsize=8.6, color=theme["muted"])
    ax.tick_params(axis="x", length=0, pad=6)
    ax.tick_params(axis="y", length=0, pad=10)
    ax.grid(axis="x", color=theme["grid"], linewidth=0.6, linestyle="-")
    ax.set_axisbelow(True)


def axes_scale(fig, ax):
    pos = ax.get_position()
    w_in = pos.width * fig.get_figwidth()
    h_in = pos.height * fig.get_figheight()
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    return abs(x1 - x0) / w_in, abs(y1 - y0) / h_in


def paired_bars(fig, ax, rows, labels, theme, colors, fmt_label, suffix=None):
    """Two bars per row (first above second), 2px surface gap, rate on a 0-1 axis."""
    n = len(rows)
    ax.set_ylim(n - 0.5, -0.5)
    style_axis(ax, theme, 1.0, [0, 0.25, 0.5, 0.75, 1.0], lambda t: f"{int(t * 100)}%")
    x_per_in, y_per_in = axes_scale(fig, ax)
    thick = 13 * PX * y_per_in
    gap = 2 * PX * y_per_in
    ax.set_yticks(range(n))
    ax.set_yticklabels([labels[key] for key, *_ in rows], fontsize=9.4, color=theme["ink"])
    for i, (key, a, b) in enumerate(rows):
        for j, (pair, color) in enumerate(((a, colors[0]), (b, colors[1]))):
            yc = i - (thick + gap) / 2 if j == 0 else i + (thick + gap) / 2
            if pair is None:
                ax.text(0.008, yc, labels["not_run"], va="center", ha="left", fontsize=8.2, color=theme["muted"])
                continue
            k, m = pair
            rate = k / m
            hbar(ax, yc, rate, thick, color, x_per_in, y_per_in)
            label = fmt_label(k, m) + (f"  ({suffix[key]})" if suffix and key in suffix and j == 1 else "")
            ax.text(rate + 0.012, yc, label, va="center", ha="left", fontsize=8.4, color=theme["ink2"])


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, dpi=DPI, facecolor=fig.get_facecolor())
    plt.close(fig)


# --- charts -------------------------------------------------------------------


def chart_trigger(lang, tname):
    theme, tx = THEMES[tname], TEXT[lang]
    setup(lang, theme)
    W, H = 900, 490
    fig = frame(W, H, theme)
    header(fig, theme, tx["trigger_title"], tx["trigger_sub"],
           [(tx["old"], theme["base"]), (tx["new"], theme["s1"])], 20, W, H)
    ax = fig.add_axes([300 / W, 96 / H, (W - 300 - 60) / W, (H - 96 - 118) / H])
    paired_bars(fig, ax, TRIGGER, tx, theme, (theme["base"], theme["s1"]),
                lambda k, m: f"{k}/{m} · {round(100 * k / m)}%", suffix={"near": tx["lower"]})
    ax.axhline(3.5, color=theme["grid"], linewidth=0.8)
    footnote(fig, theme, tx["trigger_note"], W, H)
    save(fig, f"skill-loading.{lang}.{tname}.png")


def chart_route(lang, tname):
    theme, tx = THEMES[tname], TEXT[lang]
    setup(lang, theme)
    W, H = 900, 660
    fig = frame(W, H, theme)
    header(fig, theme, tx["route_title"], tx["route_sub"],
           [(tx["legacy"], theme["base"]), (tx["audit"], theme["s1"])], 20, W, H)
    ax = fig.add_axes([300 / W, 96 / H, (W - 300 - 60) / W, (H - 96 - 118) / H])
    paired_bars(fig, ax, ROUTE, tx, theme, (theme["base"], theme["s1"]), lambda k, m: f"{k}/{m}")
    ax.axhline(0.5, color=theme["grid"], linewidth=0.8)
    ax.get_yticklabels()[0].set_fontweight("semibold")
    footnote(fig, theme, tx["route_note"], W, H)
    save(fig, f"route-audit.{lang}.{tname}.png")


def chart_olymp(lang, tname):
    theme, tx = THEMES[tname], TEXT[lang]
    setup(lang, theme)
    W, H = 900, 620
    fig = frame(W, H, theme)
    colors = (theme["base"], theme["s1"], theme["s2"])
    names = (tx["cA"], tx["cD"], tx["cE"])
    header(fig, theme, tx["olymp_title"], tx["olymp_sub"], list(zip(names, colors)), 20, W, H)

    # left: mean score per condition
    axm = fig.add_axes([60 / W, (H - 140 - 120) / H, 300 / W, 120 / H])
    axm.set_ylim(2.6, -0.6)
    style_axis(axm, theme, 1.0, [0, 0.5, 1.0], lambda t: f"{t:.1f}")
    x_per_in, y_per_in = axes_scale(fig, axm)
    thick = 18 * PX * y_per_in
    axm.set_yticks(range(3))
    axm.set_yticklabels(["A", "D", "E"], fontsize=9.4, color=theme["ink"])
    for i, (v, c, used) in enumerate(zip(OLYMP_MEAN, colors, OLYMP_SKILL)):
        hbar(axm, i, v, thick, c, x_per_in, y_per_in)
        axm.text(v + 0.02, i, f"{v:.3f}", va="center", ha="left", fontsize=8.6, color=theme["ink2"])
    axm.set_title(tx["mean"], loc="left", fontsize=9.4, color=theme["ink2"], pad=8)
    for i, used in enumerate(OLYMP_SKILL):
        fig.text(60 / W, (H - 140 - 120) / H - (36 + 17 * i) / H, f"{'ADE'[i]}: {tx['skill_used'].format(used)}",
                 fontsize=8.2, color=theme["muted"], va="top")

    # right: per-problem dot plot
    axd = fig.add_axes([520 / W, 100 / H, (W - 520 - 40) / W, (H - 140 - 100) / H])
    items = list(OLYMP)
    axd.set_ylim(len(items) - 0.5, -0.5)
    style_axis(axd, theme, 1.04, [0, 0.25, 0.5, 0.75, 1.0], lambda t: f"{t:.2f}".rstrip("0").rstrip("."))
    axd.set_xlim(0, 1.04)
    axd.set_yticks(range(len(items)))
    axd.set_yticklabels(items, fontsize=8.8, color=theme["ink"])
    axd.grid(axis="y", color=theme["grid"], linewidth=0.6)
    for i, key in enumerate(items):
        for j, (v, c) in enumerate(zip(OLYMP[key], colors)):
            axd.plot([v], [i + (j - 1) * 0.24], marker="o", markersize=7.2, markerfacecolor=c,
                     markeredgecolor=theme["surface"], markeredgewidth=1.44, linestyle="", zorder=3)
    axd.set_title(tx["per_item"], loc="left", fontsize=9.4, color=theme["ink2"], pad=8)
    footnote(fig, theme, tx["olymp_note"], W, H)
    save(fig, f"orgreason-real.{lang}.{tname}.png")


def main() -> None:
    for lang in ("en", "zh"):
        for tname in ("light", "dark"):
            chart_trigger(lang, tname)
            chart_route(lang, tname)
            chart_olymp(lang, tname)
    print(f"wrote charts to {OUT}")


if __name__ == "__main__":
    main()
