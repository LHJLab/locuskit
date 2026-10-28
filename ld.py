"""ld — 直接读 PLINK .bed 计算 r²（不依赖 plink，自动读 .bim/.fam）。"""
from __future__ import annotations

import collections
import numpy as np


class LDPanel:
    """1000G 之类的参考面板。

    panel = LDPanel("1000G.EUR.QC.bed")
    panel.samples                     # 样本数
    panel.lookup_pos("7", [8_640_813, ...])   # {pos: bed 行号}
    panel.lookup_rsid(["rs123", ...])         # {rsid: bed 行号}
    panel.r2(ref_dosage, bed_index)
    """

    def __init__(self, bed_path):
        self.bed = str(bed_path)
        stem = self.bed[:-4] if self.bed.endswith(".bed") else self.bed
        self.bim = stem + ".bim"
        self.fam = stem + ".fam"
        with open(self.fam) as fh:
            self.samples = sum(1 for _ in fh)
        self.nbytes = (self.samples + 3) // 4
        pos = collections.defaultdict(list)
        idx = collections.defaultdict(list)
        with open(self.bim) as fh:
            for i, line in enumerate(fh):
                p = line.split("\t")
                if len(p) < 4:
                    continue
                pos[p[0].replace("chr", "")].append(int(p[3]))
                idx[p[0].replace("chr", "")].append(i)
        self._index = {c: (np.asarray(pos[c], dtype=np.int64),
                           np.asarray(idx[c], dtype=np.int64)) for c in pos}
        self._fh = None

    # ---------------- 索引 ----------------
    def chromosomes(self):
        return sorted(self._index, key=lambda c: (len(c), c))

    def lookup_pos(self, chrom, positions):
        """只保留面板里存在的位点，返回 {pos: bed 行号}。"""
        c = str(chrom).replace("chr", "")
        if c not in self._index:
            return {}
        bp, ix = self._index[c]
        out = {}
        for p in positions:
            j = int(np.searchsorted(bp, p))
            if j < len(bp) and bp[j] == p:
                out[int(p)] = int(ix[j])
        return out

    def lookup_rsid(self, ids, chrom=None):
        """{rsid: bed 行号}。给 chrom 时只扫该染色体，速度快很多。"""
        want = set(ids)
        out = {}
        with open(self.bim) as fh:
            for i, line in enumerate(fh):
                p = line.split("\t")
                if len(p) < 2:
                    continue
                if chrom is not None and p[0].replace("chr", "") != str(chrom).replace("chr", ""):
                    continue
                if p[1] in want:
                    out[p[1]] = i
        return out

    # ---------------- 基因型 ----------------
    def _open(self):
        if self._fh is None:
            self._fh = open(self.bed, "rb")
            if self._fh.read(3) != b"\x6c\x1b\x01":
                raise ValueError(f"{self.bed}: 不是 PLINK bed 文件")
        return self._fh

    def dosage(self, i):
        """第 i 个变异的 0/1/2 剂量；缺失为 nan。"""
        fh = self._open()
        fh.seek(3 + i * self.nbytes)
        raw = np.frombuffer(fh.read(self.nbytes), dtype=np.uint8)
        codes = np.empty(len(raw) * 4, dtype=np.uint8)
        for k in range(4):
            codes[k::4] = (raw >> (2 * k)) & 3
        codes = codes[: self.samples]
        dos = np.where(codes == 3, 2, np.where(codes == 2, 1, 0)).astype(float)
        dos[codes == 1] = np.nan
        return dos

    @staticmethod
    def r2(a, b, min_n=20):
        m = ~(np.isnan(a) | np.isnan(b))
        if m.sum() < min_n:
            return float("nan")
        x, y = a[m], b[m]
        if x.std() == 0 or y.std() == 0:
            return float("nan")
        return float(np.corrcoef(x, y)[0, 1] ** 2)

    def r2_to_index(self, index_bed, target_beds):
        """给定 index 的 bed 行号与一批目标行号，返回与之对齐的 r² 数组。"""
        ref = self.dosage(index_bed)
        return np.array([self.r2(ref, self.dosage(j)) if j is not None else np.nan
                         for j in target_beds], dtype=float)

    def close(self):
        if self._fh is not None:
            self._fh.close()
            self._fh = None
