"""gwas — 汇总统计的按窗口流式读取（列名全部由调用方映射）。

不依赖 pandas，内存只与"落在窗口内的变异数"成正比。
"""
from __future__ import annotations

import gzip
import math

# 内部字段 -> 猜列时的候选名（调用方给 columns 时优先使用调用方的）
CANDIDATES = {
    "chr": ["CHROM", "#CHROM", "chr", "CHR", "chromosome"],
    "pos": ["POS", "BP", "pos", "base_pair_location"],
    "id": ["ID", "SNP", "rsid", "MarkerName"],
    "p": ["PVAL", "P", "pval", "p_value", "P.value"],
    "a1": ["EA", "A1", "effect_allele", "ALT"],
    "a2": ["NEA", "A2", "other_allele", "REF"],
    "beta": ["BETA", "beta", "Effect"],
    "se": ["SE", "standard_error", "StdErr"],
    "maf": ["MAF", "maf", "EAF", "FRQ"],
    "fcas": ["FCAS", "FRQ_A", "effect_allele_frequency_cases"],
    "fcon": ["FCON", "FRQ_U", "effect_allele_frequency"],
    "info": ["IMPINFO", "INFO", "info", "R2"],
    "n": ["NEFF", "N", "n", "NEFFDIV2"],
}


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _num(x, default=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def read_header(path, comment_prefix="##"):
    """返回 (skiprows, header)。header 是第一条非注释行按 \\t 切分的列名。"""
    skip = 0
    with _open_text(path) as fh:
        for line in fh:
            if not line or line.startswith(comment_prefix):
                skip += 1
                continue
            return skip, [c.strip() for c in line.rstrip("\n").rstrip("\r").split("\t")]
    raise ValueError(f"{path}: 未找到表头")


def resolve_columns(header, columns=None, required=("chr", "pos", "p")):
    """把内部字段映射到列名；columns 里显式给了就用它，否则按候选名猜。"""
    columns = dict(columns or {})
    out = {}
    for key, cands in CANDIDATES.items():
        val = columns.get(key) or (columns.get(key) if key in columns else None)
        if val in (None, "", "null", "None"):
            val = next((c for c in cands if c in header), None)
        out[key] = val
    missing = [k for k in required if not out.get(k)]
    if missing:
        raise ValueError(f"无法确定列 {missing}；header={header}，可用 columns={{...}} 显式指定")
    return out


def scan_windows(path, windows, columns=None, comment_prefix="##",
                 p_min=0.0, p_max=1.0, or_from_beta=True):
    """单遍扫描，抽取落在窗口内的变异。

    windows: {key: (chrom, start, end)}，坐标为文件自身构建。
    返回 {key: [record, ...]}，record 含 pos/id/a1/a2/beta/se/p/maf/or/info/n/chrom。
    """
    skip, header = read_header(path, comment_prefix)
    col = resolve_columns(header, columns)
    idx = {k: (header.index(v) if v else None) for k, v in col.items()}
    by_chr = {}
    for key, (chrom, lo, hi) in windows.items():
        by_chr.setdefault(str(chrom).replace("chr", ""), []).append((lo, hi, key))
    out = {k: [] for k in windows}
    n_rows = 0
    with _open_text(path) as fh:
        for _ in range(skip):
            fh.readline()
        for line in fh:
            fields = line.rstrip("\n").rstrip("\r").split("\t")
            n_rows += 1
            if len(fields) != len(header):
                continue
            chrom = fields[idx["chr"]].replace("chr", "")
            hit = by_chr.get(chrom)
            if not hit:
                continue
            try:
                pos = int(float(fields[idx["pos"]]))
            except (TypeError, ValueError):
                continue
            keys = [k for (lo, hi, k) in hit if lo <= pos <= hi]
            if not keys:
                continue
            p = _num(fields[idx["p"]])
            if not (p_min < p <= p_max):
                continue
            a1 = fields[idx["a1"]].strip() if idx["a1"] is not None else ""
            a2 = fields[idx["a2"]].strip() if idx["a2"] is not None else ""
            beta = _num(fields[idx["beta"]]) if idx["beta"] is not None else float("nan")
            se = _num(fields[idx["se"]]) if idx["se"] is not None else float("nan")
            maf = _num(fields[idx["maf"]]) if idx["maf"] is not None else float("nan")
            if math.isnan(maf) and idx["fcas"] is not None and idx["fcon"] is not None:
                fr = (_num(fields[idx["fcas"]]) + _num(fields[idx["fcon"]])) / 2.0
                if not math.isnan(fr):
                    maf = min(fr, 1 - fr)
            rec = {
                "chr": chrom, "pos": pos,
                "id": (fields[idx["id"]].strip() if idx["id"] is not None else f"{chrom}:{pos}"),
                "a1": a1, "a2": a2,
                "beta": beta, "se": se, "p": p, "maf": maf,
                "info": _num(fields[idx["info"]]) if idx["info"] is not None else float("nan"),
                "n": _num(fields[idx["n"]]) if idx["n"] is not None else float("nan"),
                "or": math.exp(beta) if (or_from_beta and not math.isnan(beta)) else float("nan"),
            }
            for k in keys:
                out[k].append(dict(rec))
    for k in out:
        out[k].sort(key=lambda r: r["pos"])
    return out, {"rows": n_rows, "columns": col}


def apply_liftover(records, chain, reverse=False):
    """把 records 的 pos 换成另一套构建的坐标；无法映射的丢弃。返回新列表。"""
    kept = []
    for r in records:
        q = chain.lift(r["chr"], r["pos"], reverse=reverse)
        if q is None:
            continue
        r = dict(r)
        r["pos_src"] = r["pos"]
        r["pos"] = q
        kept.append(r)
    kept.sort(key=lambda r: r["pos"])
    return kept
