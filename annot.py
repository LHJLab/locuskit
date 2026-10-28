"""annot — 基因模型、Ensembl 调控元件、ENCODE peak 支持度与证据表。"""
from __future__ import annotations

import collections
import csv
import gzip
import pickle
import re
from pathlib import Path


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# =============================================================================
# 基因模型
# =============================================================================
class GeneIndex:
    """Ensembl GTF 基因模型索引，可按窗口查询。

    genes = GeneIndex.from_gtf("Homo_sapiens.GRCh38.116.gtf.gz", cache="cache/genes.pkl")
    genes.window("7", 8_300_000, 8_900_000)          # -> [gene dict, ...]
    genes.nearest("7", 8_640_813)                    # -> gene dict（含该点的优先）
    """

    def __init__(self, genes):
        self.genes = genes
        self._by_chr = collections.defaultdict(list)
        for g in genes:
            self._by_chr[g["chrom"]].append(g)
        for c in self._by_chr:
            self._by_chr[c].sort(key=lambda g: g["start"])

    @classmethod
    def from_gtf(cls, path, cache=None, chroms=None):
        path = Path(path)
        if cache and Path(cache).exists():
            return cls(pickle.loads(Path(cache).read_bytes()))
        pat = re.compile(r'(\S+)\s+"([^"]*)"')
        acc = {}
        with _open_text(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] not in ("gene", "exon"):
                    continue
                chrom = f[0].replace("chr", "")
                if chroms and chrom not in chroms:
                    continue
                attrs = dict(pat.findall(f[8]))
                gid = attrs.get("gene_id")
                if not gid:
                    continue
                g = acc.get(gid)
                if g is None:
                    g = acc[gid] = {"id": gid, "symbol": attrs.get("gene_name", ""),
                                    "biotype": attrs.get("gene_biotype", ""),
                                    "chrom": chrom, "strand": f[6],
                                    "start": int(f[3]), "end": int(f[4]), "exons": []}
                if f[2] == "exon":
                    g["exons"].append((int(f[3]), int(f[4])))
        genes = [g for g in acc.values() if g["exons"]]
        for g in genes:
            g["exons"].sort()
            g["start"] = min(g["start"], g["exons"][0][0])
            g["end"] = max(g["end"], g["exons"][-1][1])
        genes.sort(key=lambda g: (g["chrom"], g["start"]))
        if cache:
            Path(cache).parent.mkdir(parents=True, exist_ok=True)
            Path(cache).write_bytes(pickle.dumps(genes, protocol=4))
        return cls(genes)

    def window(self, chrom, lo, hi, biotypes=None, max_n=None):
        hit = [g for g in self._by_chr.get(str(chrom).replace("chr", ""), [])
               if g["end"] > lo and g["start"] < hi
               and (not biotypes or g["biotype"] in biotypes)]
        hit.sort(key=lambda g: -(min(g["end"], hi) - max(g["start"], lo)))
        return hit[:max_n] if max_n else hit

    def nearest(self, chrom, pos, biotypes=None):
        cand = [g for g in self._by_chr.get(str(chrom).replace("chr", ""), [])
                if not biotypes or g["biotype"] in biotypes]
        inside = [g for g in cand if g["start"] <= pos <= g["end"]]
        pool = inside or cand
        if not pool:
            return None
        return min(pool, key=lambda g: min(abs(g["start"] - pos), abs(g["end"] - pos)))

    def __len__(self):
        return len(self.genes)


# =============================================================================
# Ensembl 调控元件
# =============================================================================
class RegFeatures:
    """Ensembl Regulatory Build（promoter / enhancer / CTCF / open chromatin ...）。"""

    def __init__(self, feats):
        self.feats = feats                      # chrom -> [(start, end, type)]
        for c in self.feats:
            self.feats[c].sort()

    @classmethod
    def from_gff(cls, path, types=None, cache=None):
        if cache and Path(cache).exists():
            d = pickle.loads(Path(cache).read_bytes())
            return cls(d)
        feats = collections.defaultdict(list)
        with _open_text(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                if types and f[2] not in types:
                    continue
                try:
                    feats[f[0].replace("chr", "")].append((int(f[3]), int(f[4]), f[2]))
                except ValueError:
                    continue
        if cache:
            Path(cache).parent.mkdir(parents=True, exist_ok=True)
            Path(cache).write_bytes(pickle.dumps(dict(feats), protocol=4))
        return cls(dict(feats))

    def window(self, chrom, lo, hi):
        return [(s, e, t) for (s, e, t) in self.feats.get(str(chrom).replace("chr", ""), [])
                if e > lo and s < hi]


# =============================================================================
# ENCODE peak 支持度 + SNP 证据表
# =============================================================================
class EncodePeaks:
    """把一个 mark 的多个 bed.gz 合成"支持度区间"。

    peaks = EncodePeaks(".../01_peaks", marks=["H3K27ac", "H3K4me3", "CTCF", "H3K27me3", "DNase"],
                        cache_dir="cache/encode")
    peaks.support("H3K27ac", "chr7")      # -> [(start, end, n_datasets), ...]
    """

    def __init__(self, root, marks, cache_dir=None, dnase_subdir="hotspots"):
        self.root = Path(root)
        self.marks = list(marks)
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.dnase_subdir = dnase_subdir

    def files(self, mark):
        d = self.root / mark
        if mark.upper() == "DNASE" and (d / self.dnase_subdir).is_dir():
            d = d / self.dnase_subdir
        return sorted(d.glob("*.bed.gz"))

    def support(self, mark, chrom):
        c = str(chrom).replace("chr", "")
        if self.cache_dir:
            cf = self.cache_dir / f"peaks_{mark}_{c}.pkl"
            if cf.exists():
                return pickle.loads(cf.read_bytes())
        events = collections.defaultdict(int)
        for f in self.files(mark):
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    if not line.strip() or line[0] == "#":
                        continue
                    p = line.split("\t", 3)
                    if len(p) < 3 or p[0].replace("chr", "") != c:
                        continue
                    try:
                        s, e = int(p[1]), int(p[2])
                    except ValueError:
                        continue
                    events[s] += 1
                    events[e] -= 1
        rows, cur, prev = [], 0, None
        for pos in sorted(events):
            if prev is not None and cur > 0:
                rows.append((prev, pos, cur))
            cur += events[pos]
            prev = pos
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            (self.cache_dir / f"peaks_{mark}_{c}.pkl").write_bytes(
                pickle.dumps(rows, protocol=4))
        return rows


class EvidenceTable:
    """PGC 显著 SNP 的 dlPFC 证据表（snp_evidence_wide.tsv）。

    ev = EvidenceTable.from_tsv(".../snp_evidence_wide.tsv")
    ev.get("rs3824754")      # -> {'n_ev': 5, 'CTCF': 7, 'H3K27ac': 16, ...}
    """

    FIELDS = [("n_evidence_types", "n_ev"), ("CTCF_n", "CTCF"), ("H3K27ac_n", "H3K27ac"),
              ("H3K4me3_n", "H3K4me3"), ("H3K27me3_n", "H3K27me3"),
              ("DNase_n", "DNase"), ("cCRE_n", "cCRE")]

    def __init__(self, rows):
        self.rows = rows

    @classmethod
    def from_tsv(cls, path):
        rows = {}
        p = Path(path)
        if not p.exists():
            return cls(rows)
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    rows[r["rs"]] = {k: int(float(r.get(src) or 0)) for src, k in cls.FIELDS}
                except (TypeError, ValueError):
                    continue
        return cls(rows)

    def get(self, rsid):
        return self.rows.get(rsid)

    def __len__(self):
        return len(self.rows)
