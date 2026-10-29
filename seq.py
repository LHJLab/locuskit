#!/usr/bin/env python3
"""seq — 导出目标 SNP 两侧的参考序列与"效应等位"版本序列。

    # 单个 / 多个 SNP（批量）
    python -m locuskit.seq \\
        --gwas   data/003_SCZ_Nature2022/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv.gz \\
        --columns '{"chr":"CHROM","pos":"POS","id":"ID","a1":"A1","a2":"A2"}' \\
        --fasta  data/001_GRCh37_DNA/Homo_sapiens.GRCh37.dna.chromosomes.fa \\
        --snp    rs4702 rs3824754 rs1619283 \\
        --flank  500 \\
        --out    seqs.tsv  --fasta-out seqs.fa

    # 从文件批量读 rsID（每行一个，或第一列）
    python -m locuskit.seq ... --snp-file my_snps.txt --flank 200

    # 只想看某一段
    python -m locuskit.seq ... --snp rs4702 --flank 50 --mark brackets

输出
----
* TSV：每个 SNP 一行 —— 位置、参考碱基、GWAS 两个等位、效应等位、等位匹配情况、
  序列长度、SNP 在序列中的偏移、参考序列、效应等位序列
* FASTA：每个 SNP 一条记录，序列为**效应等位版本**（默认）；`--include-ref` 再加一条参考序列

标记方式 `--mark`
----------------
    upper     侧翼小写、SNP 位点大写（默认，适合做 motif/PWM 扫描）
    brackets  侧翼原样、SNP 位点用 [ ] 括起
    none      全部大写原样
"""
from __future__ import annotations

import argparse
import json
import os
import sys

from .gwas import _open_text, read_header, resolve_columns

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


# =============================================================================
# 参考基因组（按 .fai 随机访问）
# =============================================================================
class Fasta:
    """只读随机访问。需同目录同名 `.fai`（`samtools faidx` 或 `pyfaidx` 生成）。"""

    def __init__(self, path):
        self.path = str(path)
        self.fai = self.path + ".fai"
        if not os.path.exists(self.fai):
            raise SystemExit(f"缺少索引：{self.fai}\n先用 `samtools faidx {self.path}` 建索引")
        self.index = {}
        with open(self.fai) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 5:
                    self.index[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
        self._fh = None

    @property
    def chroms(self):
        return list(self.index)

    def _resolve(self, chrom):
        c = str(chrom).replace("chr", "")
        for cand in (c, "chr" + c, str(chrom)):
            if cand in self.index:
                return cand
        raise KeyError(f"参考序列里没有染色体 {chrom!r}（有 {self.chroms[:6]} …）")

    def _open(self):
        if self._fh is None:
            self._fh = open(self.path, "rb")
        return self._fh

    def fetch(self, chrom, start, end):
        """1-based、闭区间，返回大写序列。"""
        name = self._resolve(chrom)
        length, offset, lb, lw = self.index[name]
        start = max(1, int(start))
        end = min(int(end), length)
        if end < start:
            return ""

        def byte_of(pos):                     # 1-based 位置 -> 字节偏移
            return offset + (pos - 1) // lb * lw + (pos - 1) % lb

        fh = self._open()
        fh.seek(byte_of(start))
        raw = fh.read(byte_of(end) - byte_of(start) + 1)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode().upper()

    def close(self):
        if self._fh is not None:
            self._fh.close()
            self._fh = None


# =============================================================================
# 从汇总统计里取目标 SNP 的等位与位置
# =============================================================================
def lookup_variants(gwas_path, targets, columns=None, comment_prefix="##"):
    """单遍扫描，返回 {target: record}；record 含 chrom/pos/a1/a2。

    target 可以是 rsID，也可以是 "chr:pos"。
    """
    want = set(targets)
    want_pos = {}
    for t in targets:
        if ":" in str(t) and not str(t).startswith("rs"):
            c, p = str(t).split(":")[:2]
            want_pos[(c.replace("chr", ""), int(p.replace(",", "")))] = t

    skip, header = read_header(gwas_path, comment_prefix)
    col = resolve_columns(header, columns)
    idx = {k: (header.index(v) if v else None) for k, v in col.items()}
    out = {}
    seen = 0
    with _open_text(gwas_path) as fh:
        for _ in range(skip):
            fh.readline()
        for line in fh:
            f = line.rstrip("\n").rstrip("\r").split("\t")
            if len(f) != len(header):
                continue
            seen += 1
            rid = f[idx["id"]].strip() if idx.get("id") is not None else None
            chrom = f[idx["chr"]].replace("chr", "")
            try:
                pos = int(float(f[idx["pos"]]))
            except (TypeError, ValueError):
                continue
            hit = None
            if rid in want:
                hit = rid
            elif (chrom, pos) in want_pos:
                hit = want_pos[(chrom, pos)]
            if hit is None:
                continue
            out[hit] = {
                "chrom": chrom, "pos": pos, "id": rid or hit,
                "a1": f[idx["a1"]].strip() if idx.get("a1") is not None else "",
                "a2": f[idx["a2"]].strip() if idx.get("a2") is not None else "",
                "p": f[idx["p"]].strip() if idx.get("p") is not None else "",
            }
            if len(out) == len(want) + len(want_pos):
                break
    return out, seen


# =============================================================================
# 序列导出
# =============================================================================
def build_sequence(seq, pos_in_seq, allele, mark="upper"):
    """把 SNP 位点换成 allele 并按 mark 方式标注。"""
    i = pos_in_seq - 1
    left, right = seq[:i], seq[i + 1:]
    if mark == "upper":
        return left.lower() + allele.upper() + right.lower()
    if mark == "brackets":
        return f"{left}[{allele.upper()}]{right}"
    return left + allele.upper() + right


def export(gwas_path, targets, fasta_path, flank=500, columns=None,
           comment_prefix="##", build_label="GRCh37", mark="upper",
           include_ref=False, effect_col="a1"):
    """返回 (rows, fasta_records)。"""
    variants, n_rows = lookup_variants(gwas_path, targets, columns, comment_prefix)
    missing = [t for t in targets if t not in variants]
    fa = Fasta(fasta_path)
    rows, records = [], []
    try:
        for t in targets:
            v = variants.get(t)
            if not v:
                rows.append({"rsid": t, "status": "not_in_gwas"})
                continue
            lo, hi = v["pos"] - flank, v["pos"] + flank
            ref_seq_full = fa.fetch(v["chrom"], lo, hi)
            if len(ref_seq_full) < 2 * flank + 1:
                # 靠近染色体端点
                lo = max(1, v["pos"] - flank)
                hi = v["pos"] + flank
            ref_base = fa.fetch(v["chrom"], v["pos"], v["pos"])
            if not ref_base:
                rows.append({"rsid": t, "chrom": v["chrom"], "pos": v["pos"],
                             "status": "no_reference"})
                continue
            ref_base = ref_base[0]
            pos_in_seq = v["pos"] - lo + 1
            seq_len = len(ref_seq_full)

            a1, a2 = (v["a1"] or "").upper(), (v["a2"] or "").upper()
            # 哪一列是效应等位由调用方声明（SCZ 的 A1 常是参考等位，MDD 的 EA 是效应等位）
            if effect_col == "a2":
                effect, other = a2, a1
            else:
                effect, other = a1, a2
            # 定链向：哪个等位等于参考碱基（先正链，再反向互补）；
            # 再判断效应等位是参考等位还是替代等位。
            def _rel(a):
                if not a:
                    return None
                if a == ref_base:
                    return "forward"
                if revcomp(a) == ref_base:
                    return "reverse"
                return None

            r_eff, r_oth = _rel(effect), _rel(other)
            if not effect:
                match = "no_effect_allele"
            elif len(effect) != 1 or (other and len(other) != 1):
                match = "indel_or_multiallelic"
            elif r_oth == "forward":
                match = "forward/effect_is_alt"
            elif r_eff == "forward":
                match = "forward/effect_is_ref"
            elif r_oth == "reverse":
                match = "reverse/effect_is_alt"
            elif r_eff == "reverse":
                match = "reverse/effect_is_ref"
            else:
                match = "unmatched"           # 可疑：构建不符或多等位
            eff_base = effect[:1] if (effect and len(effect) == 1) else ref_base
            ref_marked = build_sequence(ref_seq_full, pos_in_seq, ref_base, mark)
            eff_marked = build_sequence(ref_seq_full, pos_in_seq, eff_base, mark)

            rows.append({
                "rsid": v["id"], "chrom": v["chrom"], "pos": v["pos"], "build": build_label,
                "ref_base": ref_base, "gwas_a1": a1, "gwas_a2": a2,
                "effect_allele": effect, "other_allele": other, "allele_match": match,
                "gwas_p": v["p"], "flank": flank, "seq_len": seq_len,
                "snp_offset": pos_in_seq, "ref_seq": ref_marked, "effect_seq": eff_marked,
                "status": "ok",
                "note": "" if len(effect) == 1 else "非单碱基等位，未做替换",
            })
            meta = (f"{v['id']} chr{v['chrom']}:{v['pos']:,} {build_label} "
                    f"ref={ref_base} a1={a1} a2={a2} effect={effect or '?'} "
                    f"flank={flank} match={match}")
            records.append((f"{v['id']}|effect|{meta}", eff_marked))
            if include_ref:
                records.append((f"{v['id']}|ref|{meta}", ref_marked))
    finally:
        fa.close()
    return rows, records, missing, n_rows


def write_outputs(rows, records, out_tsv, out_fasta):
    cols = ["rsid", "status", "chrom", "pos", "build", "ref_base", "gwas_a1", "gwas_a2",
            "effect_allele", "other_allele", "allele_match", "gwas_p",
            "flank", "seq_len", "snp_offset", "ref_seq", "effect_seq", "note"]
    if out_tsv:
        os.makedirs(os.path.dirname(os.path.abspath(out_tsv)) or ".", exist_ok=True)
        with open(out_tsv, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    if out_fasta and records:
        os.makedirs(os.path.dirname(os.path.abspath(out_fasta)) or ".", exist_ok=True)
        with open(out_fasta, "w") as fh:
            for name, seq in records:
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m locuskit.seq",
                                 description="导出目标 SNP 两侧参考序列与效应等位序列")
    ap.add_argument("--gwas", required=True, help="汇总统计（取位置与等位）")
    ap.add_argument("--columns", default=None,
                    help='列名映射 JSON，如 \'{"chr":"CHROM","pos":"POS","id":"ID","a1":"A1","a2":"A2"}\'')
    ap.add_argument("--comment-prefix", default="##")
    ap.add_argument("--fasta", required=True, help="参考基因组 FASTA（需同名 .fai）")
    ap.add_argument("--snp", nargs="*", default=[], help="rsID 或 chr:pos，可多个")
    ap.add_argument("--snp-file", default=None, help="每行一个 rsID（取第一列）")
    ap.add_argument("--flank", type=int, default=500, help="SNP 两侧各取多少 bp（默认 500）")
    ap.add_argument("--mark", choices=["upper", "brackets", "none"], default="upper")
    ap.add_argument("--build-label", default="GRCh37")
    ap.add_argument("--effect-col", choices=["a1", "a2"], default="a1",
                    help="哪一列是效应等位：MDD 的 EA→a1；SCZ 若把 A1(参考等位) 映射成 a1，则用 a2")
    ap.add_argument("--include-ref", action="store_true", help="FASTA 里同时输出参考序列")
    ap.add_argument("--out", default=None, help="输出 TSV")
    ap.add_argument("--fasta-out", default=None, help="输出 FASTA")
    a = ap.parse_args(argv)

    targets = list(a.snp)
    if a.snp_file:
        with open(a.snp_file) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    targets.append(line.split("\t")[0].split(",")[0].strip())
    targets = list(dict.fromkeys(targets))       # 去重保序
    if not targets:
        ap.error("至少给一个 --snp 或 --snp-file")
    if a.flank < 0:
        ap.error("--flank 不能为负")
    if not a.out and not a.fasta_out:
        a.out = "snp_sequences.tsv"

    cols = json.loads(a.columns) if a.columns else None
    print(f"[1/2] 扫描汇总统计取等位：{a.gwas}", flush=True)
    rows, records, missing, n_rows = export(
        a.gwas, targets, a.fasta, flank=a.flank, columns=cols,
        comment_prefix=a.comment_prefix, build_label=a.build_label,
        mark=a.mark, include_ref=a.include_ref, effect_col=a.effect_col)
    ok = sum(1 for r in rows if r.get("status") == "ok")
    print(f"      扫描 {n_rows:,} 行 → 命中 {ok}/{len(targets)} 个 SNP", flush=True)

    print("[2/2] 写出", flush=True)
    write_outputs(rows, records, a.out, a.fasta_out)
    if a.out:
        print(f"      TSV   {a.out}  ({len(rows)} 行)")
    if a.fasta_out:
        print(f"      FASTA {a.fasta_out}  ({len(records)} 条记录)")
    if missing:
        print(f"      未在汇总统计中找到：{', '.join(missing[:8])}"
              + (" …" if len(missing) > 8 else ""), flush=True)
    sus = [r["rsid"] for r in rows if r.get("allele_match") == "unmatched"]
    if sus:
        print(f"      注意：{len(sus)} 个 SNP 的效应等位既不是参考碱基也不是其互补"
              f"（可能是构建不符/多等位）：{', '.join(sus[:6])}", flush=True)
    bad = [r["rsid"] for r in rows if r.get("status") != "ok"]
    if bad:
        print(f"      序列导出失败：{', '.join(bad[:8])}", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
