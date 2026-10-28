"""render — 位点图渲染：关联 / 基因模型 / Ensembl 调控元件 / ENCODE 支持度。

一条 90 mm 宽的图最多 4 条轨道，全部共享 x 轴，可选内容自动决定版面。
"""
from __future__ import annotations

import math

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt            # noqa: E402
from matplotlib.lines import Line2D        # noqa: E402
from matplotlib.patches import Rectangle   # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

# ---- 样式（默认 = 自动模式；plot_locus 可覆盖）
from .style import Style, DEFAULTS as _STYLE_DEFAULTS

_ST = Style()

# ---- 配色（保留旧名，供外部引用；实际绘制走 _ST）
GREY = _STYLE_DEFAULTS["background_color"]
GENE = "#4d4d4d"
HIGENE = "#b2182b"
IDX = "#6a3d9a"
RING = "#000000"
LD_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "locuszoom", _STYLE_DEFAULTS["ld_cmap"])
REG_COLORS = list(_STYLE_DEFAULTS["reg_colors"])
ENC_COLORS = list(_STYLE_DEFAULTS["encode_colors"])
ENC_COLOR = dict(ENC_COLORS)
REG_COLORS = [("promoter", "#e6550d"), ("enhancer", "#31a354"),
              ("CTCF_binding_site", "#3182bd"), ("open_chromatin_region", "#756bb1")]
REG_LABEL = {"promoter": "promoter", "enhancer": "enhancer",
             "CTCF_binding_site": "CTCF", "open_chromatin_region": "open chromatin"}
ENC_COLORS = [("H3K27ac", "#D95F02"), ("H3K4me3", "#1B9E77"), ("H3K27me3", "#7570B3"),
              ("CTCF", "#4D4D4D"), ("DNase", "#E7298A")]
ENC_COLOR = dict(ENC_COLORS)


def style(font_size=None, label_size=None, **over):
    """设置 matplotlib rcParams；参数缺省时用 Style 默认值。"""
    st = Style(**{k: v for k, v in over.items() if v is not None})
    if font_size is not None:
        st["font_size"] = font_size
    if label_size is not None:
        st["axis_label_size"] = label_size
    st.apply_rcparams()
    return


def _legacy_rcparams(font_size=7, label_size=6.8):
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": font_size, "axes.labelsize": label_size, "axes.titlesize": 7.5,
        "xtick.labelsize": 6, "ytick.labelsize": 6,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "axes.linewidth": 0.6, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.0, "ytick.major.size": 2.0,
    })


# =============================================================================
# 基因轨道
# =============================================================================
def _base_points(ax, x, y):
    """所有变异先铺一层（灰）底；关闭时只画有着色的点。"""
    if _ST["draw_background_points"]:
        ax.scatter(x, y, s=_ST["background_size"], c=_ST["background_color"],
                   marker=_ST["point_marker"], linewidths=0,
                   alpha=_ST["background_alpha"], zorder=2)


def _display_name(g, style_name="short"):
    if g.get("symbol"):
        return g["symbol"]
    if style_name == "full":
        return g["id"]
    if style_name == "biotype":
        return g.get("biotype") or g["id"]
    return "\u2026" + g["id"][-6:]


def _layout_genes(axg, genes, lo, hi, highlight, rows=5, name_style="short", fontsize=None):
    """把基因排进行，并预留"基因 + 标签"的空间，返回带 row/label_ok 的列表。"""
    fontsize = _ST["gene_label_size"] if fontsize is None else fontsize
    gl = []
    for g in genes:
        d = dict(g)
        d["name"] = _display_name(g, name_style)
        merged = []
        for s, e in sorted(g["exons"]):
            if merged and s <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        d["merged"] = merged
        gl.append(d)
    gl.sort(key=lambda z: z["start"])
    rend = axg.figure.canvas.get_renderer()
    probe = axg.text(0, 0, "MMMMMMMMMM", fontsize=fontsize + 0.4)
    charw = (probe.get_window_extent(renderer=rend).width / 10.0
             / axg.get_window_extent(renderer=rend).width * (hi - lo))
    probe.remove()
    pad = (hi - lo) * 0.013
    taken = [[] for _ in range(rows)]

    def fits(row, a0, b0):
        return all(b0 < x0 - pad or a0 > x1 + pad for x0, x1 in row)

    for v in gl:
        w = len(v["name"]) * charw
        if (v["end"] - v["start"]) >= w + 2 * pad:
            v["lmode"], a0, a1 = "center", v["start"], v["end"]
        elif v["end"] + pad + w <= hi:
            v["lmode"], a0, a1 = "right", v["start"], v["end"] + pad + w
        elif v["start"] - pad - w >= lo:
            v["lmode"], a0, a1 = "left", v["start"] - pad - w, v["end"]
        else:
            v["lmode"], a0, a1 = "right", v["start"], v["end"] + pad + w
        v["label_ok"] = False
        for i in range(rows):
            if fits(taken[i], a0, a1):
                taken[i].append((a0, a1))
                taken[i].sort()
                v["row"], v["label_ok"] = i, True
                break
        if not v["label_ok"]:
            for i in range(rows):
                if fits(taken[i], v["start"], v["end"]):
                    taken[i].append((v["start"], v["end"]))
                    taken[i].sort()
                    v["row"] = i
                    break
            else:
                v["row"] = rows - 1
    for v in gl:
        if v["name"] in highlight or v.get("symbol") in highlight:
            v["label_ok"] = True
    return gl


def _draw_genes(axg, gl, lo, hi, highlight, rows=5):
    axg.set_ylim(-0.55, rows - 0.40)
    axg.set_yticks([])
    axg.set_xlim(lo, hi)
    for v in gl:
        r = v["row"]
        hi_ = v["name"] in highlight or v.get("symbol") in highlight
        col = _ST["gene_hl_color"] if hi_ else _ST["gene_color"]
        lw = _ST["gene_lw_hl"] if hi_ else _ST["gene_lw"]
        axg.plot([v["start"], v["end"]], [r, r], color=col, lw=lw,
                 solid_capstyle="butt", zorder=3)
        for s0, e0 in v["merged"]:
            _eh = _ST["gene_exon_h"]
            axg.add_patch(Rectangle((s0, r - _eh / 2), max(e0 - s0, (hi - lo) * 0.0011),
                                    _eh, facecolor=col, edgecolor="none", zorder=4))
        vs, ve = max(v["start"], lo), min(v["end"], hi)
        if _ST["gene_arrows"] and ve - vs > (hi - lo) * 0.015:
            xe = ve - (hi - lo) * 0.004 if v["strand"] == "+" else vs + (hi - lo) * 0.004
            axg.plot([xe], [r], marker=">" if v["strand"] == "+" else "<",
                     ms=_ST["gene_arrow_size"],
                     color=col, markeredgewidth=0, zorder=5)

    rend = axg.figure.canvas.get_renderer()
    axb = axg.get_window_extent(renderer=rend)

    def put(v, dy):
        hi_ = v["name"] in highlight or v.get("symbol") in highlight
        fs = _ST["gene_label_size_hl"] if hi_ else _ST["gene_label_size"]
        col = _ST["gene_hl_color"] if hi_ else _ST["gene_color"]
        fw = "bold" if hi_ else "normal"
        y = v["row"] + 0.30 + dy
        if v["lmode"] == "center":
            t = axg.text((v["start"] + v["end"]) / 2, y, v["name"], fontsize=fs, va="center",
                         ha="center", color=col, fontweight=fw, clip_on=True, zorder=6)
        elif v["lmode"] == "left":
            t = axg.text(v["start"] - (hi - lo) * 0.006, y, v["name"], fontsize=fs, va="center",
                         ha="right", color=col, fontweight=fw, clip_on=True, zorder=6)
        else:
            t = axg.text(v["end"] + (hi - lo) * 0.006, y, v["name"], fontsize=fs, va="center",
                         ha="left", color=col, fontweight=fw, clip_on=True, zorder=6)
        bb = t.get_window_extent(renderer=rend)
        if bb.x1 > axb.x1 or bb.x0 < axb.x0:       # 越界钳回轴内
            x0d, y0d = axg.transData.transform(t.get_position())
            shift = (axb.x1 - bb.x1) if bb.x1 > axb.x1 else (axb.x0 - bb.x0)
            t.set_position((axg.transData.inverted().transform((x0d + shift, y0d))[0],
                            t.get_position()[1]))
        return t

    placed, n_lab = [], 0
    for v in sorted(gl, key=lambda z: (z["name"] not in highlight
                                       and z.get("symbol") not in highlight,
                                       z["row"], z["start"])):
        hi_ = v["name"] in highlight or v.get("symbol") in highlight
        if not v["label_ok"] and not hi_:
            continue
        done = False
        for dy in (0, 0.20, 0.40):
            t = put(v, dy)
            bb = t.get_window_extent(renderer=rend).expanded(1.12, 1.30)
            if not any(bb.overlaps(q) for q in placed):
                placed.append(bb)
                n_lab += 1
                done = True
                break
            t.remove()
        if not done and hi_:
            t = put(v, 0)
            bb = t.get_window_extent(renderer=rend).expanded(1.12, 1.30)
            if not any(bb.overlaps(q) for q in placed):
                placed.append(bb)
                n_lab += 1
            else:
                t.remove()
    return n_lab, len(gl)


# =============================================================================
# 调控元件轨道
# =============================================================================
def _draw_reg(axr, elems, lo, hi, types=None):
    types = types or _ST["reg_colors"]
    order = [t for t, _ in types]
    color = dict(types)
    axr.set_ylim(0, max(1, len(order)))
    axr.set_yticks([])
    axr.set_xlim(lo, hi)
    minw = (hi - lo) * _ST["reg_min_width_frac"]
    counts = {k: 0 for k in order}
    for s, e, tp in elems:
        if tp not in color:
            continue
        counts[tp] += 1
        y = len(order) - order.index(tp) - 0.5
        axr.add_patch(Rectangle((max(s, lo), y - 0.34), max(min(e, hi) - max(s, lo), minw),
                                0.68, facecolor=color[tp], edgecolor="none", zorder=3))
    return counts


# =============================================================================
# ENCODE 支持度轨道
# =============================================================================
def _draw_encode(axe, segments, lo, hi, marks, n_datasets=None):
    """segments: {mark: [(start, end, support), ...]}"""
    n_datasets = _ST["encode_datasets"] if n_datasets is None else n_datasets
    _ecol = _ST.encode_color_map()
    axe.set_ylim(0, max(1, len(marks)))
    axe.set_yticks([])
    axe.set_xlim(lo, hi)
    n = len(marks)
    for i, mark in enumerate(marks):
        y = n - 1 - i
        axe.hlines(y + 0.16, lo, hi, color="#DDDDDD", lw=0.5, zorder=0)
        for s, e, sup in segments.get(mark, []):
            if e <= lo or s >= hi:
                continue
            a = 0.35 + 0.65 * min(1.0, sup / float(n_datasets))
            axe.add_patch(Rectangle((max(s, lo), y + 0.16),
                                    max(min(e, hi) - max(s, lo), (hi - lo) * 0.0006),
                                    0.68, facecolor=_ecol.get(mark, "#888888"),
                                    alpha=a, lw=0, zorder=2))
        axe.text(-0.012, y + 0.5, mark, transform=axe.get_yaxis_transform(),
                 fontsize=_ST["tick_size"], ha="right", va="center",
                 color=_ST["gene_color"])


# =============================================================================
# 主渲染
# =============================================================================
def plot_locus(records, *, chrom=None, index_snp=None, index_pos=None,
               credset=None, genes=None, reg=None, reg_types=None,
               encode=None, encode_marks=None,
               r2=None, color_by="auto", highlight_genes=(),
               gene_rows=5, gene_biotypes=None, gene_name_style="short",
               build_label="GRCh37", figsize_mm=None,
               title=None, subtitle=None, note=None, red_line=None,
               show_reg=True, show_encode=False, ld_label=None,
               region=None, style=None, dpi=None):
    """画一张位点图，返回 (fig, info)。

    只有 assoc 轨道是必需的；genes / reg / encode 给数据才画。
    color_by: auto | ld | pip | none；r2 传入时优先按 LD 上色。
    """
    if not records:
        raise ValueError("records 为空")
    global _ST
    _ST = Style.of(style)
    w_mm = float(figsize_mm[0]) if figsize_mm else float(_ST["width_mm"])
    h_mm = float(figsize_mm[1]) if figsize_mm else float(_ST["height_mm"])
    dpi = int(dpi or _ST["dpi"])
    red_line = _ST["vline"] if red_line is None else red_line
    ld_label = ld_label or _ST["ld_label"]
    chrom = str(chrom or records[0]["chr"]).replace("chr", "")
    pos = np.array([r["pos"] for r in records], dtype=np.int64)
    lo, hi = int(pos.min()), int(pos.max())
    if hi - lo <= 0:
        hi = lo + 600_000
        lo = max(1, lo - 300_000)
    pv = np.array([r["p"] for r in records], dtype=float)
    good = (pv > 0) & (pv <= 1)
    pos, pv = pos[good], pv[good]
    recs = [r for r, g in zip(records, good) if g]
    y = -np.log10(pv)
    if index_pos is None:
        index_pos = next((r["pos"] for r in recs if r["id"] == index_snp), None)

    # ---------- 是否会有色标（决定右侧是否留位）
    _has_pip = bool(credset) and any(
        not math.isnan(c.get("pip", float("nan"))) for c in (credset or []))
    will_cbar = ((r2 is not None and color_by in ("auto", "ld"))
                 or (r2 is None and _has_pip and color_by in ("auto", "pip")))

    # ---------- 轨道与版面（按毫米分配，随画布大小自适应）
    # 页眉高度 = max(样式给定, 内容实际所需)，避免字号变大后图例压到轴上沿
    _t_mm = _ST["title_size"] * 0.36 + 0.2
    _s_mm = _ST["subtitle_size"] * 0.36 + 0.2
    _n_mm = _ST["note_size"] * 0.36 + 0.2
    _l_mm = _ST["legend_size"] * 0.36 + 1.6
    _need_mm = 1.6 + _t_mm + _s_mm + _n_mm + 2 * _l_mm + 0.4
    header_mm = max(float(_ST["header_mm"] if (genes and show_reg and reg is not None)
                          else _ST["header_mm_small"]), _need_mm)
    footer_mm = float(_ST["footer_mm"])
    gap_mm = float(_ST["gap_mm"])
    rows_mm = []
    names = []
    names.append("assoc");  rows_mm.append(0.0)      # 占位，稍后按剩余高度算
    if genes:
        names.append("genes"); rows_mm.append(gene_rows * float(_ST["gene_row_mm"]))
    if genes and show_reg and reg is not None:
        names.append("reg"); rows_mm.append(float(_ST["reg_track_mm"]))
    if show_encode and encode is not None:
        n_mk = len(encode_marks or ENC_COLORS)
        names.append("encode"); rows_mm.append(
            max(11.0, n_mk * float(_ST["encode_row_mm"])))

    fixed = sum(rows_mm[1:]) + gap_mm * max(0, len(rows_mm) - 1)
    avail = h_mm - header_mm - footer_mm
    assoc_h = max(float(_ST["assoc_min_mm"]), avail - fixed)
    rows_mm[0] = assoc_h
    total = sum(rows_mm)
    height_ratios = [r / total for r in rows_mm]
    bottom = footer_mm / h_mm
    top = 1.0 - header_mm / h_mm
    fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4),
                     facecolor=_ST["facecolor"])
    gs = fig.add_gridspec(len(names), 2, height_ratios=height_ratios,
                          width_ratios=[1, float(_ST["cbar_width_ratio"])]
                          if will_cbar else [1, 0.001],
                          hspace=float(_ST["hspace"]), wspace=float(_ST["wspace"]),
                          left=float(_ST["left"]),
                          right=float(_ST["right"]) if will_cbar else float(_ST["right_no_cbar"]),
                          top=top, bottom=bottom)
    axes = {}
    order = names
    for i, n in enumerate(order):
        axes[n] = fig.add_subplot(gs[i, 0], sharex=axes[order[0]] if i else None)
    cax = fig.add_subplot(gs[0, 1].subgridspec(10, 1)[0:7, 0])

    ax = axes["assoc"]
    drew_cbar = False

    # ---------- 关联
    ax.axhline(0, color="none")
    mode = color_by
    if mode == "auto":
        mode = "ld" if r2 is not None else ("pip" if credset else "none")
    if mode == "ld" and r2 is not None:
        r2 = np.asarray(r2, dtype=float)
        ok = ~np.isnan(r2)
        _base_points(ax, pos, y)
        if ok.any():
            m = ax.scatter(pos[ok], y[ok], s=_ST["point_size"], c=r2[ok],
                           cmap=_ST.ld_colormap(), vmin=_ST["ld_vmin"], vmax=_ST["ld_vmax"],
                           marker=_ST["point_marker"], linewidths=_ST["point_edge_width"],
                           edgecolors=_ST["point_edge"], alpha=_ST["point_alpha"], zorder=3)
            cb = fig.colorbar(m, cax=cax)
            drew_cbar = True
            cb.set_label(ld_label, fontsize=_ST["cbar_label_size"])
            cb.set_ticks(list(_ST["ld_ticks"]))
            cb.ax.tick_params(labelsize=_ST["cbar_tick_size"], length=1.5)
            cb.outline.set_linewidth(_ST["spine_width"])
    elif mode == "pip" and credset:
        pip = {}
        for r in credset:
            pip[(r.get("id") or "").strip()] = r.get("pip", float("nan"))
            if r.get("pos"):
                pip[r["pos"]] = r.get("pip", float("nan"))
        vals = np.array([pip.get(r["id"], pip.get(r["pos"], np.nan)) for r in recs], float)
        ok = ~np.isnan(vals)
        _base_points(ax, pos, y)
        if ok.any():
            m = ax.scatter(pos[ok], y[ok], s=_ST["point_size"], c=vals[ok],
                           cmap=_ST.pip_colormap(), vmin=_ST["ld_vmin"], vmax=_ST["ld_vmax"],
                           marker=_ST["point_marker"], linewidths=_ST["point_edge_width"],
                           edgecolors=_ST["point_edge"], alpha=_ST["point_alpha"], zorder=3)
            cb = fig.colorbar(m, cax=cax)
            drew_cbar = True
            cb.set_label(_ST["pip_label"], fontsize=_ST["cbar_label_size"])
            cb.ax.tick_params(labelsize=_ST["cbar_tick_size"], length=1.5)
            cb.outline.set_linewidth(_ST["spine_width"])
    else:
        _base_points(ax, pos, y)

    # credible 圈
    if credset:
        cids = {(r.get("id") or "").strip() for r in credset}
        cpos = {r["pos"] for r in credset if r.get("pos")}
        ring = np.array([(r["id"] in cids) or (r["pos"] in cpos) for r in recs])
        if ring.any():
            ax.scatter(pos[ring], y[ring], s=_ST["cred_size"],
                       facecolors=_ST["cred_face"], edgecolors=_ST["cred_edge"],
                       marker=_ST["cred_marker"], linewidths=_ST["cred_edge_width"], zorder=5)
    # index SNP
    if index_pos is not None:
        k = int(np.argmin(np.abs(pos - index_pos)))
        ax.scatter([index_pos], [y[k]], marker=_ST["index_marker"], s=_ST["index_size"],
                   c=_ST["index_color"], edgecolors=_ST["index_edge"],
                   linewidths=_ST["index_edge_width"], zorder=6)
        if _ST["index_label"] and index_snp:
            # 靠近轴顶时改画在点下方，避免与页眉/图例相撞
            _ymax = float(y.max()) if len(y) else 1.0
            _above = y[k] < _ymax * 0.90
            ax.annotate(str(index_snp), (index_pos, y[k]),
                        xytext=(0, _ST["index_label_dy"] if _above else -_ST["index_label_dy"]),
                        textcoords="offset points", ha="center",
                        va="bottom" if _above else "top",
                        fontsize=_ST["index_label_size"], color=_ST["index_label_color"],
                        zorder=7)
    ax.set_xlim(lo, hi)
    ax.set_ylim(bottom=0)
    if region and _ST["region_span"]:
        ax.axvspan(max(region[0], lo), min(region[1], hi), color=_ST["region_color"],
                   alpha=_ST["region_alpha"], zorder=0, lw=0)
    ax.set_ylabel("\u2212log\u2081\u2080(P)", fontsize=_ST["axis_label_size"])
    ax.tick_params(labelbottom=False)

    # ---------- 基因
    n_lab = n_gene = 0
    hl = set(highlight_genes or ())
    if genes:
        gl = _layout_genes(axes["genes"], genes, lo, hi, hl, rows=gene_rows,
                           name_style=gene_name_style)
        n_lab, n_gene = _draw_genes(axes["genes"], gl, lo, hi, hl, rows=gene_rows)
        axes["genes"].tick_params(labelbottom=False)
        axes["genes"].spines["left"].set_visible(True)
    # ---------- 调控
    reg_counts = {}
    if genes and show_reg and reg is not None:
        reg_counts = _draw_reg(axes["reg"], reg, lo, hi, reg_types)
        axes["reg"].tick_params(labelbottom=False)
    # ---------- ENCODE
    enc_segments = {}
    if show_encode and encode is not None:
        marks = encode_marks or [m for m, _ in ENC_COLORS]
        enc_segments = {m: encode.get(m, []) for m in marks}
        _draw_encode(axes["encode"], enc_segments, lo, hi, marks)
        axes["encode"].tick_params(labelbottom=False)

    # ---------- 红线 + x 轴
    if red_line and index_pos is not None:
        for n in order:
            axes[n].axvline(index_pos, color=_ST["vline_color"], lw=_ST["vline_width"],
                            ls=_ST["vline_style"], alpha=_ST["vline_alpha"], zorder=7,
                            solid_capstyle="butt")
    bottom = axes[order[-1]]
    bottom.set_xlabel(f"Chromosome {chrom} position (Mb, {build_label})",
                      fontsize=_ST["axis_label_size"])
    bottom.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v/1e6:.2f}"))

    # ---------- 页眉与图例（位置按毫米反推，字号取样式）
    if title:
        ttl = title
    else:
        names_hl = sorted({n for n in hl if n})
        ttl = str(index_snp) + (("  |  " + "-".join(names_hl)) if names_hl else "")
    # 页眉行距按字号推算（0.36 mm/pt ≈ 1.2 倍行高），保证改字号时不重叠
    def _mm(pt, pad=0.0):
        return pt * 0.36 + pad
    _y = 1.6
    y_title = 1.0 - _y / h_mm
    _y += _mm(_ST["title_size"], 0.2)
    y_sub = 1.0 - _y / h_mm
    _y += _mm(_ST["subtitle_size"], 0.2)
    y_note = 1.0 - _y / h_mm
    _y += _mm(_ST["note_size"], 0.2)
    y_leg1 = 1.0 - _y / h_mm
    _y += _mm(_ST["legend_size"], 1.2)
    y_leg2 = 1.0 - _y / h_mm
    fig.text(float(_ST["left"]), y_title, ttl, fontsize=_ST["title_size"],
             fontweight=_ST["title_weight"], va="top", color=_ST["title_color"])
    if subtitle:
        sub = subtitle
    elif len(recs):
        _k = int(np.nanargmin(pv))
        sub = f"{len(recs):,} variants   ·   min P {pv[_k]:.2g} ({recs[_k]['id']})"
    else:
        sub = "0 variants"
    fig.text(float(_ST["left"]), y_sub, sub, fontsize=_ST["subtitle_size"], va="top",
             color=_ST["subtitle_color"])
    if note:
        fig.text(float(_ST["left"]), y_note, note, fontsize=_ST["note_size"], va="top",
                 color=_ST["note_color"])
    handles = []
    if credset:
        handles.append(Line2D([], [], marker=_ST["cred_marker"], ls="",
                              ms=math.sqrt(_ST["cred_size"]) * 0.9,
                              mfc=_ST["cred_face"], mec=_ST["cred_edge"],
                              mew=_ST["cred_edge_width"], label="credible-set member"))
    handles.append(Line2D([], [], marker=_ST["index_marker"], ls="",
                          ms=math.sqrt(_ST["index_size"]) * 0.82,
                          mfc=_ST["index_color"], mec=_ST["index_edge"],
                          mew=_ST["index_edge_width"], label="index SNP"))
    fig.legend(handles=handles, loc="upper left",
               bbox_to_anchor=(float(_ST["left"]), y_leg1),
               ncol=3, fontsize=_ST["legend_size"], handletextpad=.35,
               columnspacing=1.5, borderpad=0, borderaxespad=0)
    if reg_counts:
        _rc = _ST["reg_colors"]
        rh = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="none") for _, c in _rc]
        fig.legend(handles=rh, labels=[REG_LABEL.get(t, t) for t, _ in _rc],
                   loc="upper left", bbox_to_anchor=(float(_ST["left"]), y_leg2),
                   ncol=min(4, len(_rc)), fontsize=_ST["legend_size"],
                   handlelength=1.1, handleheight=0.55, handletextpad=.35,
                   columnspacing=1.2, borderpad=0, borderaxespad=0)

    # ---------- 脊线 / 网格 / 背景
    for name, a in axes.items():
        a.set_facecolor(_ST["axes_facecolor"])
        if _ST["grid"]:
            a.grid(True, color=_ST["grid_color"], linewidth=_ST["grid_width"])
        else:
            a.grid(False)
        for sp in ("left", "right", "top", "bottom"):
            a.spines[sp].set_visible(sp in _ST["spines"])
            a.spines[sp].set_linewidth(_ST["spine_width"])
        if "left" not in _ST["spines"]:
            a.tick_params(axis="y", length=0)

    if not drew_cbar:
        cax.set_visible(False)
    info = {"lo": lo, "hi": hi, "chrom": chrom, "index_pos": index_pos,
            "n_variants": len(recs), "n_genes": n_gene, "n_labelled": n_lab,
            "reg_counts": reg_counts, "color_mode": mode,
            "n_encode": sum(len(v) for v in enc_segments.values())}
    return fig, info
