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

# ---- 配色
GREY = "#d9d9d9"
GENE = "#4d4d4d"
HIGENE = "#b2182b"
IDX = "#6a3d9a"
RING = "#000000"
LD_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("locuszoom", [
    (0.00, "#d9d9d9"), (0.20, "#d9d9d9"), (0.20, "#1f78b4"),
    (0.40, "#33a02c"), (0.60, "#fdb863"), (0.80, "#e31a1c"), (1.00, "#8c0d0d")])
REG_COLORS = [("promoter", "#e6550d"), ("enhancer", "#31a354"),
              ("CTCF_binding_site", "#3182bd"), ("open_chromatin_region", "#756bb1")]
REG_LABEL = {"promoter": "promoter", "enhancer": "enhancer",
             "CTCF_binding_site": "CTCF", "open_chromatin_region": "open chromatin"}
ENC_COLORS = [("H3K27ac", "#D95F02"), ("H3K4me3", "#1B9E77"), ("H3K27me3", "#7570B3"),
              ("CTCF", "#4D4D4D"), ("DNase", "#E7298A")]
ENC_COLOR = dict(ENC_COLORS)


def style(font_size=7, label_size=6.8):
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
def _display_name(g, style_name="short"):
    if g.get("symbol"):
        return g["symbol"]
    if style_name == "full":
        return g["id"]
    if style_name == "biotype":
        return g.get("biotype") or g["id"]
    return "\u2026" + g["id"][-6:]


def _layout_genes(axg, genes, lo, hi, highlight, rows=5, name_style="short", fontsize=5.1):
    """把基因排进行，并预留"基因 + 标签"的空间，返回带 row/label_ok 的列表。"""
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
        col, lw = (HIGENE, 0.75) if hi_ else (GENE, 0.5)
        axg.plot([v["start"], v["end"]], [r, r], color=col, lw=lw,
                 solid_capstyle="butt", zorder=3)
        for s0, e0 in v["merged"]:
            axg.add_patch(Rectangle((s0, r - 0.072), max(e0 - s0, (hi - lo) * 0.0011),
                                    0.144, facecolor=col, edgecolor="none", zorder=4))
        vs, ve = max(v["start"], lo), min(v["end"], hi)
        if ve - vs > (hi - lo) * 0.015:
            xe = ve - (hi - lo) * 0.004 if v["strand"] == "+" else vs + (hi - lo) * 0.004
            axg.plot([xe], [r], marker=">" if v["strand"] == "+" else "<", ms=2.2,
                     color=col, markeredgewidth=0, zorder=5)

    rend = axg.figure.canvas.get_renderer()
    axb = axg.get_window_extent(renderer=rend)

    def put(v, dy):
        hi_ = v["name"] in highlight or v.get("symbol") in highlight
        fs, col = (5.5, HIGENE) if hi_ else (5.1, GENE)
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
    types = types or REG_COLORS
    order = [t for t, _ in types]
    color = dict(types)
    axr.set_ylim(0, max(1, len(order)))
    axr.set_yticks([])
    axr.set_xlim(lo, hi)
    minw = (hi - lo) * 0.0016
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
def _draw_encode(axe, segments, lo, hi, marks, n_datasets=20):
    """segments: {mark: [(start, end, support), ...]}"""
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
            axe.add_patch(Rectangle((max(s, lo), y + 0.16), max(min(e, hi) - max(s, lo), (hi - lo) * 0.0006),
                                    0.68, facecolor=ENC_COLOR.get(mark, "#888888"),
                                    alpha=a, lw=0, zorder=2))
        axe.text(-0.012, y + 0.5, mark, transform=axe.get_yaxis_transform(),
                 fontsize=6.5, ha="right", va="center", color="#222222")


# =============================================================================
# 主渲染
# =============================================================================
def plot_locus(records, *, chrom=None, index_snp=None, index_pos=None,
               credset=None, genes=None, reg=None, reg_types=None,
               encode=None, encode_marks=None,
               r2=None, color_by="auto", highlight_genes=(),
               gene_rows=5, gene_biotypes=None, gene_name_style="short",
               build_label="GRCh37", figsize_mm=(90, 92),
               title=None, subtitle=None, note=None, red_line=True,
               show_reg=True, show_encode=False, ld_label="LD r² with index SNP",
               dpi=300):
    """画一张位点图，返回 (fig, info)。

    只有 assoc 轨道是必需的；genes / reg / encode 给数据才画。
    color_by: auto | ld | pip | none；r2 传入时优先按 LD 上色。
    """
    if not records:
        raise ValueError("records 为空")
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

    # ---------- 轨道与版面（按毫米分配，随画布大小自适应）
    w_mm, h_mm = float(figsize_mm[0]), float(figsize_mm[1])
    header_mm = 14.0 if (genes and show_reg and reg is not None) else 11.0
    footer_mm = 9.0
    gap_mm = 2.2
    rows_mm = []
    names = []
    names.append("assoc");  rows_mm.append(0.0)      # 占位，稍后按剩余高度算
    if genes:
        names.append("genes"); rows_mm.append(gene_rows * 4.6)
    if genes and show_reg and reg is not None:
        names.append("reg"); rows_mm.append(7.5)
    if show_encode and encode is not None:
        n_mk = len(encode_marks or ENC_COLORS)
        names.append("encode"); rows_mm.append(max(11.0, n_mk * 5.0))

    fixed = sum(rows_mm[1:]) + gap_mm * max(0, len(rows_mm) - 1)
    avail = h_mm - header_mm - footer_mm
    assoc_h = max(28.0, avail - fixed)
    rows_mm[0] = assoc_h
    total = sum(rows_mm)
    height_ratios = [r / total for r in rows_mm]
    bottom = footer_mm / h_mm
    top = 1.0 - header_mm / h_mm
    fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4))
    gs = fig.add_gridspec(len(names), 2, height_ratios=height_ratios,
                          width_ratios=[1, 0.030],
                          hspace=0.20, wspace=0.040,
                          left=0.142, right=0.885, top=top, bottom=bottom)
    axes = {}
    order = names
    for i, n in enumerate(order):
        axes[n] = fig.add_subplot(gs[i, 0], sharex=axes[order[0]] if i else None)
    cax = fig.add_subplot(gs[0, 1].subgridspec(10, 1)[0:7, 0])

    ax = axes["assoc"]

    # ---------- 关联
    ax.axhline(0, color="none")
    mode = color_by
    if mode == "auto":
        mode = "ld" if r2 is not None else ("pip" if credset else "none")
    if mode == "ld" and r2 is not None:
        r2 = np.asarray(r2, dtype=float)
        ok = ~np.isnan(r2)
        ax.scatter(pos, y, s=7.0, c=GREY, linewidths=0, alpha=.95, zorder=2)
        if ok.any():
            m = ax.scatter(pos[ok], y[ok], s=7.0, c=r2[ok], cmap=LD_CMAP, vmin=0, vmax=1,
                           linewidths=0, alpha=.95, zorder=3)
            cb = fig.colorbar(m, cax=cax)
            cb.set_label(ld_label, fontsize=5.8)
            cb.set_ticks([0, .2, .4, .6, .8, 1.0])
            cb.ax.tick_params(labelsize=5.2, length=1.5)
            cb.outline.set_linewidth(0.5)
    elif mode == "pip" and credset:
        pip = {}
        for r in credset:
            pip[(r.get("id") or "").strip()] = r.get("pip", float("nan"))
            if r.get("pos"):
                pip[r["pos"]] = r.get("pip", float("nan"))
        vals = np.array([pip.get(r["id"], pip.get(r["pos"], np.nan)) for r in recs], float)
        ok = ~np.isnan(vals)
        ax.scatter(pos, y, s=7.0, c=GREY, linewidths=0, alpha=.95, zorder=2)
        if ok.any():
            m = ax.scatter(pos[ok], y[ok], s=7.0, c=vals[ok], cmap="viridis", vmin=0, vmax=1,
                           linewidths=0, alpha=.95, zorder=3)
            cb = fig.colorbar(m, cax=cax)
            cb.set_label("credible-set PIP", fontsize=5.8)
            cb.ax.tick_params(labelsize=5.2, length=1.5)
            cb.outline.set_linewidth(0.5)
    else:
        ax.scatter(pos, y, s=7.0, c=GREY, linewidths=0, alpha=.95, zorder=2)

    # credible 圈
    if credset:
        cids = {(r.get("id") or "").strip() for r in credset}
        cpos = {r["pos"] for r in credset if r.get("pos")}
        ring = np.array([(r["id"] in cids) or (r["pos"] in cpos) for r in recs])
        if ring.any():
            ax.scatter(pos[ring], y[ring], s=13, facecolors="none", edgecolors=RING,
                       linewidths=.35, zorder=5)
    # index SNP
    if index_pos is not None:
        k = int(np.argmin(np.abs(pos - index_pos)))
        ax.scatter([index_pos], [y[k]], marker="D", s=15, c=IDX, edgecolors="white",
                   linewidths=.4, zorder=6)
    ax.set_xlim(lo, hi)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("\u2212log\u2081\u2080(P)", fontsize=6.8)
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
            axes[n].axvline(index_pos, color="#D7191C", lw=0.6, alpha=0.9, zorder=7,
                            solid_capstyle="butt")
    bottom = axes[order[-1]]
    bottom.set_xlabel(f"Chromosome {chrom} position (Mb, {build_label})", fontsize=6.8)
    bottom.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v/1e6:.2f}"))

    # ---------- 页眉与图例
    if title:
        ttl = title
    else:
        names = sorted({n for n in hl if n})
        ttl = str(index_snp) + (("  |  " + "-".join(names)) if names else "")
    y_title = 1.0 - 1.8 / h_mm
    y_sub = 1.0 - 4.6 / h_mm
    y_note = 1.0 - 6.9 / h_mm
    y_leg1 = 1.0 - 8.7 / h_mm
    y_leg2 = 1.0 - 11.4 / h_mm
    fig.text(0.142, y_title, ttl, fontsize=7.4, fontweight="bold", va="top")
    if subtitle:
        sub = subtitle
    elif len(recs):
        _k = int(np.nanargmin(pv))
        sub = f"{len(recs):,} variants   ·   min P {pv[_k]:.2g} ({recs[_k]['id']})"
    else:
        sub = "0 variants"
    fig.text(0.142, y_sub, sub, fontsize=5.9, va="top")
    if note:
        fig.text(0.142, y_note, note, fontsize=5.3, va="top", color="#a50f15")
    handles = []
    if credset:
        handles.append(Line2D([], [], marker="o", ls="", ms=3.6, mfc="none", mec="black",
                              mew=.35, label="credible-set member"))
    handles.append(Line2D([], [], marker="D", ls="", ms=3.2, mfc=IDX, mec="white",
                          mew=.45, label="index SNP"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.142, y_leg1),
               ncol=3, fontsize=5.4, handletextpad=.35, columnspacing=1.5,
               borderpad=0, borderaxespad=0)
    if reg_counts:
        rh = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="none") for _, c in REG_COLORS]
        fig.legend(handles=rh, labels=[REG_LABEL[t] for t, _ in REG_COLORS],
                   loc="upper left", bbox_to_anchor=(0.142, y_leg2), ncol=4, fontsize=5.4,
                   handlelength=1.1, handleheight=0.55, handletextpad=.35, columnspacing=1.2,
                   borderpad=0, borderaxespad=0)

    info = {"lo": lo, "hi": hi, "chrom": chrom, "index_pos": index_pos,
            "n_variants": len(recs), "n_genes": n_gene, "n_labelled": n_lab,
            "reg_counts": reg_counts, "color_mode": mode,
            "n_encode": sum(len(v) for v in enc_segments.values())}
    return fig, info
