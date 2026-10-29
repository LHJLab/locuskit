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
           include_genome=True, effect_col="a1", ref_col=None):
    """导出每个 SNP 的三条序列：基因组参考 / GWAS 参考等位 / GWAS 效应等位。

    返回 (rows, fasta_records, missing, n_rows_scanned)。
    三条序列都是**正链**（与 FASTA 同向）：若 GWAS 等位在反向链，会先取互补再写入。
    """
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
            lo, hi = max(1, v["pos"] - flank), v["pos"] + flank
            seq = fa.fetch(v["chrom"], lo, hi)
            ref_base = fa.fetch(v["chrom"], v["pos"], v["pos"])
            if not seq or not ref_base:
                rows.append({"rsid": t, "chrom": v["chrom"], "pos": v["pos"],
                             "status": "no_reference"})
                continue
            ref_base = ref_base[0]
            pos_in_seq = v["pos"] - lo + 1

            a1, a2 = (v["a1"] or "").upper(), (v["a2"] or "").upper()
            if effect_col == "a2":
                eff_rep, ref_rep = a2, a1
            else:
                eff_rep, ref_rep = a1, a2
            if ref_col == "a1":
                ref_rep = a1
            elif ref_col == "a2":
                ref_rep = a2

            def _rel(a):
                if not a:
                    return None
                if a == ref_base:
                    return "forward"
                return "reverse" if revcomp(a) == ref_base else None

            r_eff, r_ref = _rel(eff_rep), _rel(ref_rep)
            if not eff_rep:
                strand = "no_effect_allele"
            elif len(eff_rep) != 1 or (ref_rep and len(ref_rep) != 1):
                strand = "indel_or_multiallelic"
            elif r_eff == "forward" or r_ref == "forward":
                strand = "forward"
            elif r_eff == "reverse" or r_ref == "reverse":
                strand = "reverse"
            else:
                strand = "unmatched"

            # 统一到正链
            if strand == "reverse":
                eff_fwd = revcomp(eff_rep)[:1] if len(eff_rep) == 1 else ""
                ref_fwd = revcomp(ref_rep)[:1] if len(ref_rep) == 1 else ""
            else:
                eff_fwd = eff_rep[:1] if len(eff_rep) == 1 else ""
                ref_fwd = ref_rep[:1] if len(ref_rep) == 1 else ""

            seqs = {
                "genome": build_sequence(seq, pos_in_seq, ref_base, mark),
                "gwas_ref": build_sequence(seq, pos_in_seq, ref_fwd or ref_base, mark),
                "gwas_effect": build_sequence(seq, pos_in_seq, eff_fwd or ref_base, mark),
            }
            row = {
                "rsid": v["id"], "chrom": v["chrom"], "pos": v["pos"], "build": build_label,
                "genome_ref_base": ref_base,
                "gwas_a1": a1, "gwas_a2": a2,
                "gwas_ref_allele": ref_rep, "gwas_effect_allele": eff_rep,
                "ref_fwd": ref_fwd, "effect_fwd": eff_fwd,
                "strand_vs_genome": strand,
                "effect_is_ref": "yes" if (eff_fwd and eff_fwd == ref_base) else "no",
                "gwas_p": v["p"], "flank": flank, "seq_len": len(seq),
                "snp_offset": pos_in_seq,
                "genome_seq": seqs["genome"],
                "gwas_ref_seq": seqs["gwas_ref"],
                "gwas_effect_seq": seqs["gwas_effect"],
                "status": "ok",
                "note": "" if len(eff_rep) == 1 else "非单碱基等位，未替换",
            }
            rows.append(row)
            meta = (f"{v['id']} chr{v['chrom']}:{v['pos']:,} {build_label} "
                    f"genome={ref_base} a1={a1} a2={a2} "
                    f"ref={ref_rep or '?'} effect={eff_rep or '?'} "
                    f"strand={strand} flank={flank}")
            if include_genome:
                records.append((f"{v['id']}|genome|{meta}", seqs["genome"]))
            records.append((f"{v['id']}|gwas_ref|{meta}", seqs["gwas_ref"]))
            records.append((f"{v['id']}|gwas_effect|{meta}", seqs["gwas_effect"]))
    finally:
        fa.close()
    return rows, records, missing, n_rows


# =============================================================================
# 可视化
# =============================================================================
# 淡色系（ColorBrewer Paired 浅色），配深灰字，避免大片高饱和色
BASE_COLORS = {"A": "#b2df8a", "C": "#a6cee3", "G": "#fdbf6f", "T": "#fb9a99",
               "N": "#d9d9d9", "a": "#b2df8a", "c": "#a6cee3", "g": "#fdbf6f",
               "t": "#fb9a99", "n": "#d9d9d9"}
BASE_TEXT = "#333333"


def plot_sequence(row, view=30, out=None, dpi=300, width_mm=180, row_mm=7.0,
                  show_letters=None, font_size=6.0):
    """把三条序列画成碱基色块图（每个 SNP 一张）。返回 figure。

    row    : export() 产出的一行（需 status == ok）
    view   : SNP 两侧各显示多少 bp
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    if row.get("status") != "ok":
        raise ValueError(f"{row.get('rsid')}: 该行没有可用序列")

    def _plain(seq):
        """去掉 [ ] 标记（--mark brackets 时序列里带括号）"""
        return seq.replace("[", "").replace("]", "")

    def _center(seq, marked_offset0):
        """括号会让偏移错位，按剥掉标记后的位置重算。"""
        return len(_plain(seq[:marked_offset0]))

    raw = {"genome": row["genome_seq"],
           "gwas_ref": row["gwas_ref_seq"],
           "gwas_effect": row["gwas_effect_seq"]}
    groups = [("genome", "genome", row["genome_ref_base"]),
              ("GWAS ref", "gwas_ref", row["ref_fwd"] or row["genome_ref_base"]),
              ("GWAS effect", "gwas_effect", row["effect_fwd"] or row["genome_ref_base"])]
    seqs = [(label, _plain(raw[key]), base_letter) for label, key, base_letter in groups]
    n_full = max(len(x[1]) for x in seqs)
    center = _center(raw["genome"], int(row["snp_offset"]) - 1)
    lo = max(0, center - view)
    hi = min(n_full, center + view + 1)
    n = hi - lo
    shown = min(center - lo, hi - 1 - center)

    h_mm = 20 + row_mm * len(seqs)
    fig = plt.figure(figsize=(width_mm / 25.4, h_mm / 25.4))
    ax = fig.add_axes([0.13, 0.06, 0.83, 0.76])
    ax.set_xlim(lo - 0.5, hi - 0.5)
    ax.set_ylim(-0.6, len(seqs) - 0.05)
    ax.invert_yaxis()
    ax.set_yticks(range(len(seqs)))
    ax.set_yticklabels([x[0] for x in seqs], fontsize=font_size + 0.6)
    ax.set_xticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    letters = show_letters
    if letters is None:
        cell_pt = (width_mm / 25.4 * 0.84 * 72) / max(n, 1)
        letters = cell_pt >= 4.2

    for r, (_, seq, _b) in enumerate(seqs):
        seg = seq[lo:hi]
        for i, ch in enumerate(seg):
            up = ch.upper()
            base_col = BASE_COLORS.get(up, "#cccccc")
            is_snp = (lo + i) == center
            ax.add_patch(Rectangle((lo + i - 0.5, r - 0.42), 1.0, 0.84,
                                   facecolor=base_col, edgecolor="white",
                                   linewidth=0.5, zorder=2))
            if letters:
                ax.text(lo + i, r, up, ha="center", va="center",
                        fontsize=font_size, color=BASE_TEXT, zorder=3,
                        fontweight="bold" if is_snp else "normal")
        if not letters and len(seg):
            ax.text(lo - 1.0, r, "5\u2032", ha="right", va="center", fontsize=font_size)
    # SNP 竖线 + 方框
    for r in range(len(seqs)):
        ax.add_patch(Rectangle((center - 0.5, r - 0.44), 1.0, 0.88, facecolor="none",
                               edgecolor="#222222", linewidth=0.8, zorder=4))
    # 与基因组不同的列加标记
    g = seqs[0][1][lo:hi]
    for r in (1, 2):
        for i, ch in enumerate(seqs[r][1][lo:hi]):
            if ch.upper() != g[i].upper():
                ax.plot([lo + i], [len(seqs) - 0.17 + (r - 1) * 0.11],
                        marker="v", ms=2.4, color="#111111", zorder=5)

    ttl = (f"{row['rsid']}  |  chr{row['chrom']}:{int(row['pos']):,} ({row['build']})  |  "
           f"genome={row['genome_ref_base']}  GWAS ref={row.get('gwas_ref_allele') or '?'}"
           f"  effect={row.get('gwas_effect_allele') or '?'}")
    fig.text(0.02, 0.96, ttl, fontsize=font_size + 1.6, fontweight="bold", va="top")
    _p = row.get("gwas_p")
    sub = (f"P = {_p}   " if _p else "") + f"\u00b1{shown} bp"
    _strand = row.get("strand_vs_genome")
    if _strand and _strand != "forward":
        sub += f"   strand={_strand}"
    fig.text(0.02, 0.90, sub, fontsize=font_size, va="top", color="#555555")
    if out:
        os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
        fig.savefig(out, dpi=dpi)
        if str(out).lower().endswith((".png", ".jpg")):
            pass
        else:
            fig.savefig(os.path.splitext(out)[0] + ".png", dpi=dpi)
    return fig


def write_outputs(rows, records, out_tsv, out_fasta):
    cols = ["rsid", "status", "chrom", "pos", "build", "genome_ref_base",
            "gwas_a1", "gwas_a2", "gwas_ref_allele", "gwas_effect_allele",
            "ref_fwd", "effect_fwd", "strand_vs_genome", "effect_is_ref",
            "gwas_p", "flank", "seq_len", "snp_offset",
            "genome_seq", "gwas_ref_seq", "gwas_effect_seq", "note"]
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
    ap.add_argument("--no-genome", action="store_true",
                    help="FASTA 里不输出基因组参考序列（默认三条都输出）")
    ap.add_argument("--ref-col", choices=["a1", "a2"], default=None,
                    help="哪一列是 GWAS 的参考等位；缺省取 effect_col 的另一列")
    ap.add_argument("--figure", default=None, metavar="DIR",
                    help="每个 SNP 另存一张序列图到该目录（pdf + png）")
    ap.add_argument("--view", type=int, default=30, help="图上显示 SNP 两侧各多少 bp（默认 30）")
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
        mark=a.mark, include_genome=not a.no_genome,
        effect_col=a.effect_col, ref_col=a.ref_col)
    ok = sum(1 for r in rows if r.get("status") == "ok")
    print(f"      扫描 {n_rows:,} 行 → 命中 {ok}/{len(targets)} 个 SNP", flush=True)

    print("[2/2] 写出", flush=True)
    write_outputs(rows, records, a.out, a.fasta_out)
    if a.figure:
        ok_rows = [r for r in rows if r.get("status") == "ok"]
        if ok_rows:
            os.makedirs(a.figure, exist_ok=True)
            import matplotlib.pyplot as plt
            for r in ok_rows:
                p = os.path.join(a.figure, f"{r['rsid']}.pdf")
                try:
                    fig = plot_sequence(r, view=a.view, out=p)
                    plt.close(fig)
                except Exception as e:                        # noqa: BLE001
                    print(f"      出图失败 {r['rsid']}: {type(e).__name__}: {e}", flush=True)
            print(f"      图    {a.figure}/  ({len(ok_rows)} 个 SNP，±{a.view} bp)")
    if a.out:
        print(f"      TSV   {a.out}  ({len(rows)} 行)")
    if a.fasta_out:
        print(f"      FASTA {a.fasta_out}  ({len(records)} 条记录)")
    if missing:
        print(f"      未在汇总统计中找到：{', '.join(missing[:8])}"
              + (" …" if len(missing) > 8 else ""), flush=True)
    sus = [r["rsid"] for r in rows if r.get("strand_vs_genome") == "unmatched"]
    if sus:
        print(f"      注意：{len(sus)} 个 SNP 的效应等位既不是参考碱基也不是其互补"
              f"（构建不符或多等位，请核对）：{', '.join(sus[:6])}", flush=True)
    bad = [r["rsid"] for r in rows if r.get("status") != "ok"]
    if bad:
        print(f"      序列导出失败：{', '.join(bad[:8])}", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
