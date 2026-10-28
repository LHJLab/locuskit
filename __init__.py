"""locuskit — 位点精细定位与调控证据作图工具包。

典型用法（面向对象 API）
----------------------
    from locuskit import (Chain, scan_windows, apply_liftover, LDPanel,
                          GeneIndex, RegFeatures, EncodePeaks, EvidenceTable,
                          plot_locus, detail_pages, style)

    style()
    recs, meta = scan_windows("PGC3.tsv.gz",
                              {"rs4702": ("5", 137_000_000, 138_000_000)},
                              columns={"chr": "CHROM", "pos": "POS", "a1": "A1",
                                       "a2": "A2", "p": "PVAL", "id": "ID"})
    ld = LDPanel("1000G.EUR.QC.bed")
    genes = GeneIndex.from_gtf("Homo_sapiens.GRCh37.87.gtf.gz", cache="cache/genes.pkl")
    fig, info = plot_locus(recs["rs4702"], chrom="5", index_snp="rs4702",
                           genes=genes.window("5", *win), ld=ld, r2=r2,
                           title="rs4702 | FURIN")
    tables = detail_pages(recs["rs4702"], title="rs4702")

命令行：`python -m locuskit.cli config.yaml`
"""
from .liftover import Chain
from .gwas import scan_windows, apply_liftover, read_header, resolve_columns
from .ld import LDPanel
from .annot import GeneIndex, RegFeatures, EncodePeaks, EvidenceTable
from .render import plot_locus, style, LD_CMAP, REG_COLORS, ENC_COLORS
from .style import Style, PRESETS
from .tables import detail_pages, DEFAULT_COLUMNS

__all__ = [
    "Chain", "scan_windows", "apply_liftover", "read_header", "resolve_columns",
    "LDPanel", "GeneIndex", "RegFeatures", "EncodePeaks", "EvidenceTable",
    "plot_locus", "style", "LD_CMAP", "REG_COLORS", "ENC_COLORS",
    "Style", "PRESETS",
    "detail_pages", "DEFAULT_COLUMNS",
]
__version__ = "0.1.0"
