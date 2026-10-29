#!/usr/bin/env python3
"""locuskit 命令行：读 YAML 配置，出图 + 附表。

    python -m locuskit.cli config.yaml
    python -m locuskit.cli config.yaml --only rs4702 --limit 3 --list

配置见 `locuskit/README.md`，示例见 `code/config_*.yaml`。
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
from pathlib import Path

import numpy as np
import yaml

import matplotlib.pyplot as plt

from . import (Chain, EncodePeaks, EvidenceTable, GeneIndex, LDPanel, RegFeatures,
               apply_liftover, detail_pages, plot_locus, scan_windows, style)
from .render import REG_COLORS

DEFAULTS = {
    "build_label": "GRCh37",
    "window": {"flank_kb": 150, "min_kb": 200, "max_kb": 1500},
    "figure": {"width_mm": 90, "height_mm": 92, "dpi": 300},
    "output": {"dir": "figures", "formats": "pdf", "preview": True, "gene_prefix": True},
    "genes": {"rows": 5, "name_style": "short", "biotypes": None, "highlight": []},
    "regulatory": {"types": ["promoter", "enhancer", "CTCF_binding_site",
                             "open_chromatin_region"]},
    "encode": {"show": False, "datasets": 20,
               "marks": ["H3K27ac", "H3K4me3", "H3K27me3", "CTCF", "DNase"]},
    "ld": {"match_by": "position"},
    "table": {"enabled": True, "source": "credset", "top_n": 50, "page": "A4L",
              "fontsize": 5.8, "ld_cut": 0.8},
    "cache_dir": "cache",
    "style": None,
    "fasta": None,
}


def deep_merge(base, over):
    out = dict(base)
    for k, v in (over or {}).items():
        out[k] = deep_merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def load_config(path):
    cfg = deep_merge(DEFAULTS, yaml.safe_load(open(path)) or {})

    def expand(sec):
        if isinstance(sec, dict):
            for k, v in sec.items():
                if k in ("path", "bed", "gtf", "gff", "peaks", "evidence", "chain", "dir") \
                        and isinstance(v, str):
                    sec[k] = os.path.expanduser(v)
                elif isinstance(v, dict):
                    expand(v)
    expand(cfg)
    return cfg


# ---------------------------------------------------------------- 输入装配
def read_table_dicts(path, columns, comment="##", key_col=None):
    """通用小表读取（loci / credset）。

    表头可以以 # 或 ## 开头（`##index_snp` 之类），列名匹配忽略前导 # 与大小写。
    """
    import gzip
    from .gwas import _looks_like_header, _norm_col
    op = gzip.open if str(path).endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        header = None
        for line in fh:
            if not line.strip():
                continue
            f = [c.strip() for c in line.rstrip("\n").rstrip("\r").split("\t")]
            if line.startswith(comment) and not _looks_like_header(f, min_hits=1):
                continue
            if header is None:
                header = f
                norm = {}
                for i, h in enumerate(header):
                    norm.setdefault(_norm_col(h), i)
                idx = {}
                for k, v in (columns or {}).items():
                    if v in header:
                        idx[k] = header.index(v)
                    elif v and _norm_col(v) in norm:
                        idx[k] = norm[_norm_col(v)]
                    else:
                        idx[k] = None
                if key_col and idx.get(key_col) is None:
                    raise SystemExit(f"{path}: 找不到列 {columns.get(key_col)}；实际列名={header}")
                continue
            if len(f) != len(header):
                continue
            out.append({k: (f[i] if i is not None else None) for k, i in idx.items()} | {"_f": f})
    return out


def build_loci(cfg, cred_rows):
    sec = cfg.get("loci") or {}
    loci = []
    if sec.get("path"):
        for r in read_table_dicts(sec["path"], sec.get("columns"), key_col="index_snp"):
            g = (r.get("genes") or "").replace(",", ";")
            st = int(float(r["start"])) if r.get("start") else None
            en = int(float(r["end"])) if r.get("end") else None
            wcol = r.get("window")
            if wcol and (st is None or en is None):
                import re as _re
                m = _re.search(r":([\d,]+)-([\d,]+)$", str(wcol))
                if m:
                    st = int(m.group(1).replace(",", ""))
                    en = int(m.group(2).replace(",", ""))
            loci.append({"index_snp": r.get("index_snp"),
                         "chr": str(r.get("chr") or "").replace("chr", ""),
                         "start": st, "end": en,
                         "genes": [x.strip() for x in g.split(";") if x.strip()]})
    if not loci and cred_rows:
        for idx, rs in cred_rows.items():
            loci.append({"index_snp": idx, "chr": rs[0]["chr"], "start": None, "end": None,
                         "genes": []})
    if not loci and (cfg.get("auto_loci") or {}).get("enabled"):
        loci = "AUTO"
    return loci


def auto_loci(windows_recs, p_thresh, clump_kb):
    out = []
    gap = int(clump_kb) * 1000
    for chrom, recs in windows_recs.items():
        hits = sorted((r for r in recs if r["p"] <= p_thresh), key=lambda r: r["p"])
        chosen = []
        for r in hits:
            if all(abs(r["pos"] - c["pos"]) > gap for c in chosen):
                chosen.append(r)
        out += [{"index_snp": r["id"], "chr": chrom, "start": None, "end": None,
                 "genes": [], "index_pos": r["pos"]} for r in chosen]
    return out


def window_for(locus, recs_all, cfg):
    if locus.get("start") and locus.get("end"):
        lead = locus.get("index_pos")
        if lead is None and recs_all:
            lead = next((r["pos"] for r in recs_all if r["id"] == locus["index_snp"]), None)
        return lead, int(locus["start"]), int(locus["end"])
    w = cfg["window"]
    lead = locus.get("index_pos")
    if lead is None and recs_all:
        lead = next((r["pos"] for r in recs_all if r["id"] == locus["index_snp"]), None)
    if lead is None and recs_all:
        lead = min(recs_all, key=lambda r: r["p"])["pos"]
    if lead is None:
        return None
    lo = locus.get("start") or lead - 100_000
    hi = locus.get("end") or lead + 100_000
    flank = int(w["flank_kb"]) * 1000
    lo, hi = min(lo, lead) - flank, max(hi, lead) + flank
    span = hi - lo
    if span < int(w["min_kb"]) * 1000:
        mid = (lo + hi) // 2
        lo, hi = mid - int(w["min_kb"]) * 500, mid + int(w["min_kb"]) * 500
    if hi - lo > int(w["max_kb"]) * 1000:
        lo, hi = lead - int(w["max_kb"]) * 500, lead + int(w["max_kb"]) * 500
    return lead, max(lo, 1), hi


# ---------------------------------------------------------------- 主流程
def main(argv=None):
    ap = argparse.ArgumentParser(description="locuskit 位点图")
    ap.add_argument("config")
    ap.add_argument("--only")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--outdir")
    ap.add_argument("--formats")
    ap.add_argument("--no-table", action="store_true")
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    if a.outdir:
        cfg["output"]["dir"] = a.outdir
    if a.formats:
        cfg["output"]["formats"] = a.formats
    if a.no_table:
        cfg["table"]["enabled"] = False
    from .style import Style
    Style.of(cfg.get("style")).apply_rcparams()
    cache = Path(cfg["cache_dir"])
    cache.mkdir(parents=True, exist_ok=True)
    outdir = Path(cfg["output"]["dir"])
    outdir.mkdir(parents=True, exist_ok=True)

    # --- 参考数据
    chain = Chain(cfg["liftover"]["chain"]) if (cfg.get("liftover") or {}).get("chain") else None
    ld = None
    if (cfg.get("ld") or {}).get("bed"):
        ld = LDPanel(cfg["ld"]["bed"])
        print(f"LD panel: {ld.samples} samples, {len(ld.chromosomes())} chroms", flush=True)
    genes_idx = None
    if (cfg.get("genes") or {}).get("gtf"):
        genes_idx = GeneIndex.from_gtf(cfg["genes"]["gtf"],
                                       cache=cache / "genes.pkl")
        print(f"genes: {len(genes_idx):,}", flush=True)
    reg = None
    if (cfg.get("regulatory") or {}).get("gff"):
        reg = RegFeatures.from_gff(cfg["regulatory"]["gff"],
                                   types=cfg["regulatory"].get("types"),
                                   cache=cache / "reg.pkl")
    enc_peaks = None
    ev = EvidenceTable({})
    if (cfg.get("encode") or {}).get("peaks"):
        enc_peaks = EncodePeaks(cfg["encode"]["peaks"], cfg["encode"]["marks"],
                                cache_dir=cache / "encode")
        ev = EvidenceTable.from_tsv(cfg["encode"].get("evidence") or "")

    # --- 精细定位
    cred_rows = {}
    cs = cfg.get("credset") or {}
    if cs.get("path"):
        for r in read_table_dicts(cs["path"], cs.get("columns"), key_col="index_snp"):
            idx = r.get("index_snp")
            try:
                pr = float(r.get("pip")) if r.get("pip") not in (None, "", "NA") else float("nan")
            except ValueError:
                pr = float("nan")
            cred_rows.setdefault(idx, []).append(
                {"id": r.get("id"), "pos": int(float(r["pos"])) if r.get("pos") else None,
                 "pip": pr, "chr": str(r.get("chr") or "").replace("chr", "")})

    loci = build_loci(cfg, cred_rows)

    # --- 窗口（在"源构建"里扫描）
    auto = loci == "AUTO"
    if auto:
        # 自动识别需要先扫全文件（不设窗口）
        recs, meta = scan_windows(cfg["gwas"]["path"], {"ALL": ("", 0, 1 << 62)},
                                  columns=cfg["gwas"].get("columns"),
                                  comment_prefix=cfg["gwas"].get("comment_prefix", "##"))
        allrecs = recs["ALL"]
        by_chr = {}
        for r in allrecs:
            by_chr.setdefault(r["chr"], []).append(r)
        loci = auto_loci(by_chr, float(cfg["auto_loci"].get("p_thresh", 5e-8)),
                         float(cfg["auto_loci"].get("clump_kb", 500)))
        print(f"auto loci: {len(loci)}", flush=True)
        windows = {}
        for l in loci:
            w = window_for(l, by_chr.get(l["chr"], []), cfg)
            if w:
                windows[l["index_snp"]] = (l["chr"], w[1], w[2])
    else:
        if not loci:
            sys.exit("没有位点：配置 loci.path 或 auto_loci.enabled")
        windows = {}
        for l in loci:
            w = window_for(l, None, cfg)
            if not w:
                continue
            l["index_pos"] = w[0]
            windows[l["index_snp"]] = (l["chr"], w[1], w[2])
        # 若有 liftover：显示坐标 → 源坐标（反推扫描窗口）
        if chain:
            rev = bool(cfg["liftover"].get("reverse", False))
            inv = {}
            for k, (c, lo, hi) in windows.items():
                a0 = chain.lift(c, lo, reverse=not rev)
                a1 = chain.lift(c, hi, reverse=not rev)
                inv[k] = (c, min(a0, a1) - 5000, max(a0, a1) + 5000) if (a0 and a1) else (c, 0, 1 << 62)
            windows_scan = inv
        else:
            windows_scan = windows
        recs, meta = scan_windows(cfg["gwas"]["path"], windows_scan,
                                  columns=cfg["gwas"].get("columns"),
                                  comment_prefix=cfg["gwas"].get("comment_prefix", "##"))

    print(f"scanned {meta['rows']:,} rows; windows={len(windows)}", flush=True)

    if a.only:
        loci = [l for l in loci if l["index_snp"] == a.only]
    if a.limit:
        loci = loci[: a.limit]
    if a.list:
        for l in loci:
            print(f"{l['chr']}\t{l['index_snp']}")
        return

    # --- 逐个位点
    ok = 0
    for i, l in enumerate(loci, 1):
        idx = l["index_snp"]
        if idx not in windows:
            continue
        chrom, lo, hi = windows[idx]
        recs_l = recs.get(idx, [])
        if chain:
            recs_l = apply_liftover(recs_l, chain, reverse=bool(cfg["liftover"].get("reverse", False)))
            recs_l = [r for r in recs_l if lo <= r["pos"] <= hi]
        if not recs_l:
            print(f"  [{i}] {idx}: 窗口内无变异，跳过", flush=True)
            continue
        ipos = next((r["pos"] for r in recs_l if r["id"] == idx), l.get("index_pos"))

        # LD
        r2 = None
        if ld is not None and ipos is not None:
            if str(cfg["ld"].get("match_by", "position")).lower() == "rsid":
                m = ld.lookup_rsid([r["id"] for r in recs_l] + [idx], chrom=chrom)
                if idx in m:
                    ref = ld.dosage(m[idx])
                    r2 = np.array([ld.r2(ref, ld.dosage(m[r["id"]])) if r["id"] in m else np.nan
                                   for r in recs_l])
            else:
                m = ld.lookup_pos(chrom, [r["pos"] for r in recs_l])
                if ipos in m:
                    ref = ld.dosage(m[ipos])
                    r2 = np.array([ld.r2(ref, ld.dosage(m[r["pos"]])) if r["pos"] in m else np.nan
                                   for r in recs_l])

        # 注释
        hl = set(l.get("genes") or []) | set(cfg["genes"].get("highlight") or [])
        gwin = genes_idx.window(chrom, lo, hi, biotypes=cfg["genes"].get("biotypes")) if genes_idx else None
        rwin = reg.window(chrom, lo, hi) if reg else None
        enc_seg = {}
        if enc_peaks and cfg["encode"].get("show"):
            for mk in cfg["encode"]["marks"]:
                enc_seg[mk] = enc_peaks.support(mk, chrom)
            enc_seg = {mk: [x for x in v if x[1] > lo and x[0] < hi] for mk, v in enc_seg.items()}

        fig, info = plot_locus(
            recs_l, chrom=chrom, index_snp=idx, index_pos=ipos,
            credset=cred_rows.get(idx), genes=gwin, reg=rwin,
            reg_types=[(t, dict(REG_COLORS)[t]) for t in cfg["regulatory"]["types"]
                       if t in dict(REG_COLORS)] if reg else None,
            encode=enc_seg or None, encode_marks=cfg["encode"]["marks"],
            show_encode=bool(enc_peaks and cfg["encode"].get("show")),
            r2=r2, highlight_genes=hl, gene_rows=cfg["genes"]["rows"],
            gene_name_style=cfg["genes"]["name_style"],
            build_label=cfg["build_label"],
            figsize_mm=(cfg["figure"]["width_mm"], cfg["figure"]["height_mm"]),
            dpi=cfg["figure"]["dpi"], style=cfg.get("style"))

        # 附表
        tfigs = []
        if cfg["table"]["enabled"]:
            src = cfg["table"]["source"]
            rows = []
            if src == "ld" and r2 is not None:
                cut = float(cfg["table"]["ld_cut"])
                rows = [dict(r, r2=float(v)) for r, v in zip(recs_l, r2) if not np.isnan(v) and v > cut]
                rows.sort(key=lambda r: (-r["r2"], r["p"]))
            elif src == "credset" and cred_rows.get(idx):
                rows = [dict(r) for r in recs_l if r["id"] in {c["id"] for c in cred_rows[idx]}]
                rows.sort(key=lambda r: r["p"])
            elif src == "window":
                rows = [dict(r) for r in recs_l]
            if cfg["table"].get("p_max"):
                rows = [r for r in rows if r["p"] <= float(cfg["table"]["p_max"])]
            if src == "top":
                rows = sorted((dict(r) for r in recs_l), key=lambda r: r["p"])[: int(cfg["table"]["top_n"])]
            getters = {}
            if ev.rows:
                for k in EvidenceTable.FIELDS:
                    getters[k[1]] = (lambda kk: (lambda row: (ev.get(row["id"]) or {}).get(kk)))(k[1])
                cols = [("id", "rsid", None), ("pos", "position", None), ("alleles", "ea/oa", None),
                        ("maf", "maf", 3), ("beta", "BETA", 4), ("se", "SE", 4), ("p", "P", None)]
                if src == "ld":
                    cols.append(("r2", "r²", 3))
                cols += [(c, c, None) for _, c in EvidenceTable.FIELDS] + [("info", "INFO", 3)]
            else:
                cols = None
            lines = [f"window: chr{chrom}:{lo:,}-{hi:,}   ·   {len(recs_l):,} variants"
                     + (f"   ·   {len(cred_rows.get(idx, []))} credible-set members" if cred_rows.get(idx) else ""),
                     f"gene models: {info['n_genes']} ({info['n_labelled']} labelled)"
                     + (f"   ·   regulatory elements: {sum(info['reg_counts'].values())}" if info["reg_counts"] else "")]
            tfigs = detail_pages(rows, columns=cols, extra_getters=getters or None,
                                 title=f"{idx}  |  {'-'.join(sorted(hl)) if hl else '—'}  —  "
                                       f"{'high-LD SNPs' if src == 'ld' else src}",
                                 lines=lines, page=cfg["table"]["page"],
                                 fontsize=cfg["table"]["fontsize"],
                                 source_note=f"source: {Path(cfg['gwas']['path']).name}")

        name = idx.replace(":", "_").replace("/", "_")
        if cfg["output"].get("gene_prefix") and (l.get("genes") or hl):
            name = "-".join(sorted(set(l.get("genes") or hl))) + "_" + name
        for fmt in str(cfg["output"]["formats"]).split(","):
            fmt = fmt.strip()
            if fmt == "pdf":
                from matplotlib.backends.backend_pdf import PdfPages
                with PdfPages(outdir / f"{name}.pdf") as pdf:
                    pdf.savefig(fig)
                    for tf in tfigs:
                        pdf.savefig(tf)
            else:
                fig.savefig(outdir / f"{name}.{fmt}")
        if cfg["output"].get("preview"):
            fig.savefig(outdir / f"{name}.png", dpi=int(cfg["figure"]["dpi"]))
        for tf in tfigs:
            plt.close(tf)
        plt.close(fig)
        ok += 1
        print(f"  [{i}/{len(loci)}] {name}"
              + (f"  (r²>{cfg['table']['ld_cut']}: {len(tfigs)} 页表)" if cfg["table"]["source"] == "ld" else ""),
              flush=True)
    print(f"完成 {ok}/{len(loci)} → {outdir}")


if __name__ == "__main__":
    main()
