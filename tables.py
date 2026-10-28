"""tables — 页附的变异明细表（A4 横向，可分页，列可自定义）。

`columns` 是 (key, header, fmt) 的有序列表；某个 key 在所有行都取不到值时该列自动省略。
可用 key 见 `FIELDS`，用 `extra_getters` 还能挂任意派生列。
"""
from __future__ import annotations

import math

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

PAGE_SIZES = {"A4L": (297 / 25.4, 210 / 25.4), "A4": (210 / 25.4, 297 / 25.4),
              "A3L": (420 / 25.4, 297 / 25.4)}

DEFAULT_COLUMNS = [
    ("id", "rsid", None),
    ("pos", "position", None),
    ("alleles", "ea/oa", None),
    ("maf", "maf", 3),
    ("beta", "BETA", 4),
    ("se", "SE", 4),
    ("p", "P", None),
    ("r2", "r²", 3),
]


def _fmt(v, nd):
    if v is None:
        return None
    if isinstance(v, float):
        if math.isnan(v):
            return None
        if nd is None:
            return f"{v:.3g}"
        return f"{v:.{nd}f}"
    return str(v)


def _cell(row, key, fmt, extra_getters):
    if key == "alleles":
        a1, a2 = row.get("a1") or "", row.get("a2") or ""
        return f"{a1}/{a2}" if (a1 or a2) else None
    if key == "pos":
        return f"{row.get('chr','')}:{int(row['pos']):,}" if row.get("pos") else None
    if extra_getters and key in extra_getters:
        return extra_getters[key](row)
    return _fmt(row.get(key), fmt)


def detail_pages(rows, *, columns=None, title="", lines=(), page="A4L",
                 fontsize=5.8, rows_per_page=None, source_note="",
                 extra_getters=None, row_filter=None):
    """返回 figure 列表（每个元素是一页）。

    rows        : dict 列表（含 id/pos/a1/a2/p/... 及任意附加字段）
    columns     : [(key, header, fmt), ...]；缺省用 DEFAULT_COLUMNS
    lines       : 页眉下方的说明行（str 列表）
    extra_getters: {key: callable(row) -> str}
    row_filter  : callable(row) -> bool，在成表前过滤
    """
    if row_filter:
        rows = [r for r in rows if row_filter(r)]
    if not rows:
        return []
    columns = columns or DEFAULT_COLUMNS
    usable = []
    for c in columns:
        key, header, fmt = c[0], c[1], (c[2] if len(c) > 2 else None)
        if any(_cell(r, key, fmt, extra_getters) not in (None, "") for r in rows):
            usable.append((key, header, fmt))
    if not usable:
        return []

    size = PAGE_SIZES.get(str(page).upper(), PAGE_SIZES["A4L"])
    xs = [0.030 + i * (0.940 / len(usable)) for i in range(len(usable))]
    top, bot, dy = 0.866, 0.070, 0.0166
    per = int(rows_per_page or max(8, (top - bot) / dy - 2))
    pages = [rows[i:i + per] for i in range(0, len(rows), per)]
    figs = []
    for pi, pg in enumerate(pages, 1):
        f = plt.figure(figsize=size)
        f.text(0.030, 0.968, title, fontsize=8.5, fontweight="bold", va="top")
        y_line = 0.941
        for ln in list(lines)[:3]:
            f.text(0.030, y_line, ln, fontsize=6.0, va="top")
            y_line -= 0.017
        if source_note:
            f.text(0.030, y_line - 0.002, source_note, fontsize=5.4, va="top", color="#555555")
        y = top
        for xi, (_, header, _f) in zip(xs, usable):
            f.text(xi, y, header, fontsize=6.0, fontweight="bold", va="center")
        f.lines.append(plt.Line2D([0.030, 0.970], [y - 0.012, y - 0.012],
                                  color="#333333", lw=0.6, transform=f.transFigure))
        for ri, row in enumerate(pg):
            yy = y - 0.030 - ri * dy
            if ri % 2 == 0:
                f.patches.append(Rectangle((0.028, yy - dy * 0.46), 0.944, dy * 0.92,
                                           transform=f.transFigure, facecolor="#f2f2f2",
                                           edgecolor="none", zorder=0))
            for xi, (key, _h, fmt) in zip(xs, usable):
                v = _cell(row, key, fmt, extra_getters)
                if v not in (None, ""):
                    f.text(xi, yy, v, fontsize=fontsize, va="center", zorder=1)
        f.text(0.030, 0.060, f"page {pi} of {len(pages)}   ·   "
                             f"{len(pg)} of {len(rows)} rows",
               fontsize=5.4, color="#555555", va="top")
        figs.append(f)
    return figs
