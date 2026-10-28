"""liftover — 基于 UCSC chain 的坐标转换（GRCh37 <-> GRCh38）。

UCSC chain 文件的方向由文件头决定：`chain score tName ... qName ...`，
其中 t = target、q = query。本项目的 `GRCh37_to_GRCh38.chain.gz` 实际是
t=GRCh37、q=GRCh38（与文件名相反），`Chain.lift` 默认把 t 坐标映射到 q，
即 GRCh37 → GRCh38；需要反向时用 `reverse=True`。
"""
from __future__ import annotations

import bisect
import collections
import gzip


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


class Chain:
    """按染色体保存的 chain 块，支持单点与批量转换。"""

    def __init__(self, path, keep_strands="++"):
        self.path = str(path)
        self.blocks = collections.defaultdict(list)   # chrom -> [(t0, t1, q0)]
        self.target_sizes = {}
        self.query_sizes = {}
        self.skipped = 0
        cur = None
        with _open_text(path) as fh:
            for line in fh:
                if line.startswith("chain"):
                    p = line.split()
                    try:
                        chrom = p[2].replace("chr", "")
                        self.target_sizes[chrom] = int(p[3])
                        self.query_sizes[chrom] = int(p[8])
                        cur = dict(t=chrom, tstart=int(p[5]), qstart=int(p[10]),
                                   ts=p[4], qs=p[9])
                    except (ValueError, IndexError):
                        cur = None
                    continue
                v = line.split()
                if cur is None or len(v) < 3:
                    continue
                try:
                    size, dt, dq = int(v[0]), int(v[1]), int(v[2])
                except ValueError:
                    continue
                if cur["ts"] + cur["qs"] == keep_strands:
                    self.blocks[cur["t"]].append(
                        (cur["tstart"], cur["tstart"] + size, cur["qstart"]))
                else:
                    self.skipped += 1
                cur["tstart"] += size + dt
                cur["qstart"] += size + dq
        for c in self.blocks:
            self.blocks[c].sort()

    def lift(self, chrom, pos, reverse=False):
        """把 pos 从 t 坐标转到 q 坐标（reverse=True 则反过来）。不在块内返回 None。"""
        L = self.blocks.get(str(chrom).replace("chr", ""))
        if not L:
            return None
        if reverse:
            L = sorted((q0, q0 + (t1 - t0), t0) for t0, t1, q0 in L)
        i = bisect.bisect_right(L, (pos, 1 << 62, 1 << 62)) - 1
        if i >= 0 and L[i][0] <= pos < L[i][1]:
            return L[i][2] + (pos - L[i][0])
        return None

    def lift_many(self, positions, chrom, reverse=False):
        """批量转换，未映射为 None。"""
        return [self.lift(chrom, p, reverse) for p in positions]

    def __repr__(self):
        n = sum(len(v) for v in self.blocks.values())
        return f"<Chain {self.path!s} blocks={n:,} chroms={len(self.blocks)} skipped={self.skipped}>"
