"""style — 位点图的外观参数。

所有字段都有默认值，默认值等于"自动模式"（当前出图的样子）。
用法：

    plot_locus(..., style=Style(point_marker="*", index_color="#e31a1c"))
    plot_locus(..., style={"ld_cmap": ["#ffffff", "#08306b"], "point_size": 10})
    plot_locus(..., style=Style.preset("minimal"))

也可整段写进配置：

    style:
      point_marker: "s"
      ld_cmap: ["#d9d9d9", "#1f78b4", "#33a02c", "#fdb863", "#e31a1c"]
      style_preset: minimal        # 先套预设，再用上面的键覆盖
"""
from __future__ import annotations

import matplotlib

# ---------------------------------------------------------------- 预设色带
LD_STOPS = [(0.00, "#d9d9d9"), (0.20, "#d9d9d9"), (0.20, "#1f78b4"),
            (0.40, "#33a02c"), (0.60, "#fdb863"), (0.80, "#e31a1c"),
            (1.00, "#8c0d0d")]
REG_COLORS = [("promoter", "#e6550d"), ("enhancer", "#31a354"),
              ("CTCF_binding_site", "#3182bd"), ("open_chromatin_region", "#756bb1")]
ENC_COLORS = [("H3K27ac", "#D95F02"), ("H3K4me3", "#1B9E77"), ("H3K27me3", "#7570B3"),
              ("CTCF", "#4D4D4D"), ("DNase", "#E7298A")]

DEFAULTS = {
    # ---- 画布
    "width_mm": 90, "height_mm": 92, "dpi": 300,
    "facecolor": "white", "axes_facecolor": "white",
    # ---- 字体
    "font_family": "sans-serif",
    "font_sans": ("Arial", "Helvetica", "DejaVu Sans"),
    "font_size": 7, "title_size": 7.4, "subtitle_size": 5.9, "note_size": 5.3,
    "axis_label_size": 6.8, "tick_size": 6,
    "gene_label_size": 5.1, "gene_label_size_hl": 5.5,
    "legend_size": 5.4, "cbar_label_size": 5.8, "cbar_tick_size": 5.2,
    "title_weight": "bold", "title_color": "black",
    "subtitle_color": "black", "note_color": "#a50f15",
    # ---- 散点
    "point_marker": "o", "point_size": 7.0, "point_alpha": 0.95,
    "point_edge": "none", "point_edge_width": 0.0,
    "draw_background_points": True, "background_color": "#d9d9d9",
    "background_size": 7.0, "background_alpha": 0.95,
    # ---- LD / PIP 着色
    "ld_cmap": LD_STOPS, "ld_vmin": 0.0, "ld_vmax": 1.0,
    "ld_cbar": True, "ld_label": "LD r² with index SNP",
    "ld_ticks": (0, 0.2, 0.4, 0.6, 0.8, 1.0),
    "pip_cmap": "viridis", "pip_cbar": True, "pip_label": "credible-set PIP",
    # ---- index SNP
    "index_marker": "D", "index_size": 15, "index_color": "#6a3d9a",
    "index_edge": "white", "index_edge_width": 0.4,
    "index_label": False, "index_label_size": 5.4, "index_label_dy": 6,
    "index_label_color": "#6a3d9a",
    # ---- credible 圈
    "cred_marker": "o", "cred_size": 13, "cred_edge": "black",
    "cred_edge_width": 0.35, "cred_face": "none",
    # ---- 参考线 / 区域
    "vline": True, "vline_color": "#D7191C", "vline_width": 0.6,
    "vline_style": "-", "vline_alpha": 0.9,
    "region_span": True, "region_color": "black", "region_alpha": 0.03,
    # ---- 基因轨道
    "gene_color": "#4d4d4d", "gene_hl_color": "#b2182b",
    "gene_lw": 0.5, "gene_lw_hl": 0.75,
    "gene_exon_h": 0.144, "gene_body_h": 0.0,
    "gene_arrows": True, "gene_arrow_size": 2.2,
    # ---- 调控 / ENCODE 轨道
    "reg_colors": REG_COLORS, "reg_min_width_frac": 0.0016,
    "encode_colors": ENC_COLORS, "encode_datasets": 20,
    # ---- 布局（左侧留白、页眉页脚、行高）
    "left": 0.142, "right": 0.885, "right_no_cbar": 0.975, "bottom": 0.105,
    "header_mm": 14.0, "header_mm_small": 11.0, "footer_mm": 9.0, "gap_mm": 2.2,
    "gene_row_mm": 4.6, "reg_track_mm": 7.5, "encode_row_mm": 5.0,
    "assoc_min_mm": 28.0, "hspace": 0.20, "wspace": 0.040,
    "cbar_width_ratio": 0.030,
    # ---- 脊线 / 网格
    "spines": ("left", "right", "top", "bottom"),
    "spine_width": 0.6, "grid": False, "grid_color": "#e6e6e6", "grid_width": 0.4,
}

PRESETS = {
    # 默认：当前出版风格
    "default": {},
    # 极简：无背景层、无红线、无图例外的装饰
    "minimal": {"vline": False, "index_label": False, "region_span": False,
                "draw_background_points": True, "background_alpha": 0.7},
    # 大尺寸海报：整体放大
    "poster": {"width_mm": 180, "height_mm": 190, "font_size": 9,
               "title_size": 11, "subtitle_size": 8.4, "note_size": 7.6,
               "axis_label_size": 9.6, "tick_size": 8.4,
               "gene_label_size": 7.4, "gene_label_size_hl": 8.0,
               "legend_size": 7.8, "cbar_label_size": 8.2, "cbar_tick_size": 7.4,
               "point_size": 12, "background_size": 12,
               "index_size": 26, "cred_size": 22, "vline_width": 1.0,
               "gene_row_mm": 7.0, "encode_row_mm": 8.0, "reg_track_mm": 11.0,
               "header_mm": 20.0, "header_mm_small": 18.0, "footer_mm": 12.0},
    # 深色背景
    "dark": {"facecolor": "#111214", "axes_facecolor": "#111214",
             "gene_color": "#c8ccd0", "background_color": "#3a3f44",
             "vline_color": "#ff5a5f", "index_color": "#ffd166",
             "cred_edge": "#ffffff", "title_color": "#ffffff"},
}


class Style(dict):
    """带默认值的样式表；未知键也允许（透传给 matplotlib 的地方会用到）。"""

    def __init__(self, **over):
        base = dict(DEFAULTS)
        preset = over.pop("style_preset", None) or over.pop("preset", None)
        if preset:
            name = str(preset)
            if name not in PRESETS:
                raise KeyError(f"未知预设 {name!r}；可选 {sorted(PRESETS)}")
            base.update(PRESETS[name])
        base.update({k: v for k, v in over.items() if v is not None})
        super().__init__(base)

    # ---------------- 构造
    @classmethod
    def of(cls, obj=None):
        if obj is None:
            return cls()
        if isinstance(obj, Style):
            return obj
        if isinstance(obj, dict):
            return cls(**obj)
        raise TypeError(f"不能把 {type(obj)} 当作 Style")

    @classmethod
    def preset(cls, name, **over):
        return cls(style_preset=name, **over)

    # ---------------- 访问
    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError:
            raise AttributeError(k) from None

    def updated(self, **over):
        d = dict(self)
        d.update({k: v for k, v in over.items() if v is not None})
        return Style(**d)

    # ---------------- 派生对象
    def ld_colormap(self):
        """ld_cmap 可以是 matplotlib 色带名、[(pos,color),...] 或 [color,...]。"""
        cm = self["ld_cmap"]
        if isinstance(cm, str):
            return matplotlib.colormaps[cm] if hasattr(matplotlib, "colormaps") \
                else matplotlib.cm.get_cmap(cm)
        if isinstance(cm, (list, tuple)):
            if cm and isinstance(cm[0], (list, tuple)):
                return matplotlib.colors.LinearSegmentedColormap.from_list("custom", cm)
            return matplotlib.colors.LinearSegmentedColormap.from_list(
                "custom", list(cm), N=256)
        return cm

    def pip_colormap(self):
        cm = self["pip_cmap"]
        if isinstance(cm, str):
            return matplotlib.colormaps[cm] if hasattr(matplotlib, "colormaps") \
                else matplotlib.cm.get_cmap(cm)
        return self.ld_colormap() if cm is None else cm

    def reg_color_map(self):
        c = self["reg_colors"]
        return dict(c) if isinstance(c, (list, tuple)) else dict(c or {})

    def encode_color_map(self):
        c = self["encode_colors"]
        return dict(c) if isinstance(c, (list, tuple)) else dict(c or {})

    def apply_rcparams(self):
        """把字体类参数写进 matplotlib rcParams（每个进程调一次即可）。"""
        matplotlib.rcParams.update({
            "font.family": self["font_family"],
            "font.sans-serif": list(self["font_sans"]),
            "font.size": self["font_size"],
            "axes.labelsize": self["axis_label_size"],
            "axes.titlesize": self["title_size"],
            "xtick.labelsize": self["tick_size"],
            "ytick.labelsize": self["tick_size"],
            "pdf.fonttype": 42, "svg.fonttype": "none",
            "axes.linewidth": self["spine_width"],
            "xtick.major.width": 0.5, "ytick.major.width": 0.5,
            "xtick.major.size": 2.0, "ytick.major.size": 2.0,
        })
