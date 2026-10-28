# locuskit

位点精细定位与调控证据作图。三个入口：**网页 GUI**、**命令行**、**Python 库**。

依赖：`numpy`、`matplotlib`、`PyYAML`（GUI 另需 `Flask`）。**不需要** pandas / plink。

> 本文件先讲"要准备什么数据"，再讲"怎么用"，最后是配置字段与限制。
>
> ⚠️ **最容易踩的坑是构建（GRCh37 vs GRCh38）** —— 见 [第二节](#二构建grch37-还是-grch38)。

---

## 一、需要的数据

### 1.1 一览

| # | 数据 | 用途 | 必需性 | 本项目中的文件 |
|---|---|---|---|---|
| 1 | **GWAS 汇总统计** | 关联曲线（P/BETA/SE/等位） | **必需** | `data/003_SCZ_Nature2022/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv.gz`<br>`data/010_MDD_Cell2025/Cell2025/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv` |
| 2 | **位点定义**（位点表 *或* `auto_loci`） | 画哪几个位点、窗口多大 | **必需**（二选一） | `results/locus_figures.data_integrity.tsv`；或用 `auto_loci` 从汇总统计自动识别 |
| 3 | LD 参考面板 | r² 着色、`r²>阈值` 附表 | 可选 | `data/009_LD_REF/1000G_EUR_Phase3_plink_merge/1000G.EUR.QC.{bed,bim,fam}` |
| 4 | 基因注释 GTF | 基因轨道 | 可选 | GRCh37：`data/001_GRCh37_DNA/Homo_sapiens.GRCh37.87.gtf.gz`<br>GRCh38：`data/002_GRCh38_DNA/Homo_sapiens.GRCh38.116.gtf.gz` |
| 5 | 精细定位 credible set | 空心圈、PIP 着色、附表来源 | 可选 | `results/table2_SCZ106_loci.credible_set_variants.tsv` |
| 6 | Ensembl 调控元件 | 调控轨道 | 可选 | `data/001_GRCh37_DNA/Homo_sapiens.GRCh37.regulatory_features.v114.gff3.gz`<br>`data/002_GRCh38_DNA/Homo_sapiens.GRCh38.regulatory_features.v112.gff3.gz` |
| 7 | ENCODE dlPFC peak | 五条表观证据轨道 | 可选 | `data/014_ENCODE_snp_scan/snp_annotation/ENCODE_database/dlPFC_panel/processed_data/01_peaks/` |
| 8 | ENCODE SNP 证据表 | 附表里的 `n_ev / CTCF / H3K27ac / …` 列 | 可选 | `…/dlPFC_panel/snp_evidence/snp_evidence_wide.tsv` |
| 9 | liftover chain | 汇总统计与出图坐标不同构建时 | 可选 | `data/001_GRCh37_DNA/GRCh37_to_GRCh38.chain.gz` |

**最小可跑**：只要 1 + 2（2 也可以用 `auto_loci` 代替）。此时出纯灰点图 + 标题。

### 1.2 各类数据长什么样

以下都是本项目真实文件的前几行。**代码块里多个空格只是为了对齐，实际是制表符分隔**；
过长的列用 `…` 截断。

**(1) GWAS 汇总统计** —— 制表符分隔，可 `.gz`；开头可有 `##` 注释行。

```
# SCZ (PGC3 core)，注意列顺序：POS 在第 3 列，A1 是参考等位
CHROM  ID          POS        A1  A2  FCAS   FCON   IMPINFO  BETA      SE      PVAL    NGT  DIRE  NCAS  NCON  NEFF
8      rs62513865  101592213  C   T   0.932  0.928  0.954    0.01340   0.0171  0.4319  0    -?-+…  0     0     0

# MDD (PGC-MDD 2025)，注意列名带 #、POS 在第 2 列、EA 是效应等位
#CHROM  POS     ID           EA  NEA  BETA       SE      PVAL    FCAS   FCON   IMPINFO  NEFF       NCAS    NCON    HETI  HETDF  HETPVAL
1       752566  rs3094315    G   A    0.00280    0.0047  0.5535  0.177  0.179  0.843    927713.52  341197  989052  39.1  12     0.073
```

两份文件的列名和顺序完全不同 —— 所以要在 `gwas.columns` 里逐列映射，代码不猜。

**表头的井号数量不限**：`CHROM` / `#CHROM` / `##CHROM` / `##Chr` 都能识别。
判定规则是"这一行是否至少有两个已知列名"，所以 `##fileFormat=…` 这类真注释行不会被误当表头；
列名匹配时忽略**前导 `#`** 与大小写。位点表、credible set 表同理（`##index_snp` 也认）。

**(2) 位点表** —— 每行一个位点，给窗口边界。

`code/loci_enc_scz.tsv`（最简形式）：

```
index_snp   chrom  start       end
rs4146428   10     102441964   102841964
rs3824754   10     102654593   103054593
```

`results/locus_figures.data_integrity.tsv`（本项目已有的表，窗口是合并列）：

```
index_snp       chrom  window_grch37      window_kb  genes   n_credible  …
8:4180090_T_A   8      8:4028791-4902494  874        CSMD1   9           …
```

两种都支持：分开的 `start`/`end`，或合并的 `window`（`chr:start-end`）。

**(3) LD 面板** —— PLINK 三件套，同名同目录。

```
1000G.EUR.QC.bim      # 第 4 列必须是位置；第 2 列是 rsID（跨构建时按它匹配）
1   rs575272151  0  11008  G  C
1   rs544419019  0  11012  G  C

1000G.EUR.QC.fam      # 每行一个样本，只用来数样本数
HG00096 HG00096 0 0 0 -9
HG00097 HG00097 0 0 0 -9
```

**(4) 基因注释 GTF** —— Ensembl GTF，`gene` 与 `exon` 行都要。

```
1  ensembl_havana  gene  11869  14412  .  +  .  gene_id "ENSG00000223972"; gene_version "4"; gene_name "DDX11L1"; gene_biotype "pseudogene";
1  ensembl_havana  exon  11869  12227  .  +  .  gene_id "ENSG00000223972"; transcript_id "ENST00000456328";
```

**(5) credible set** —— 每行一个可信变异；PIP 可有可无。

`results/table2_SCZ106_loci.credible_set_variants.tsv`：

```
index_snp       expected_causals_k  extended_gwas  rsid         chromosome  position  GWAS_effect_allele  other_allele  maf    or      se        pval        finemap_posterior_probability  gene_symbol  …
8:4180090_T_A   1.45                YES            rs10108772   8           4178791   G                   A             0.243  0.9458  0.0091    9.287e-10   0.0164                         CSMD1        …
```

最小要求只有 `index_snp` + rsID 列；`position` 和 PIP 列缺了也能跑（缺 PIP 就不上 PIP 色）。

**(6) Ensembl 调控元件 GFF3**

```
##gff-version 3
1  Ensembl  promoter  10936  11436  .  .  .  ID=ENSR1_958;extended_start=10001;extended_end=11436;gene_id=ENSG00000290825;…
```
第 3 列是元件类型（`promoter` / `enhancer` / `CTCF_binding_site` / `open_chromatin_region`），
第 4、5 列是起止。

**(7) ENCODE peak** —— 每个 mark 一个目录，内含若干 `.bed.gz`（BED6+ 格式，前 3 列即 chr/start/end）。

```
chr2  46680  48272  Peak_20281  141  .  6.25571  14.10579  11.77212  219
chr2  51282  51664  Peak_84932  51   .  3.35656  5.19958   3.42724   153
```
每个文件是一个数据集（供体）；同一 mark 的多个文件叠加后的区间数 = 支持度（0–20）。

**(8) SNP 证据表** —— 宽表，每行一个 SNP，`*_n` 列是命中的数据集数。

```
rs         chr   pos      n_evidence_types  CTCF_n  CTCF_src     H3K27ac_n  H3K27ac_src  H3K4me3_n  H3K4me3_src  H3K27me3_n  H3K27me3_src  DNase_n  DNase_src    cCRE_n  cCRE_src                        cCRE_class
rs4592207  chr1  2438059  3                 0                    0                      0                      1           ENCSR391SNO    1        ENCSR849WGE  1       cCRE_registry_v3_highconf.bed  CTCF-bound,dELS
```

**(9) liftover chain** —— UCSC chain。

```
chain 1 1 249250621 + 10000 267719 1 248956422 + 10000 297968 2
167417 50000 80249
```
字段含义：`chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id`。
上例的 `tSize=249250621` 是 **GRCh37** 的 1 号染色体长度，`qSize=248956422` 是 GRCh38 —— 
说明本项目的 `GRCh37_to_GRCh38.chain.gz` 实际是「t=GRCh37，q=GRCh38」，
`Chain.lift()` 默认把 t 映射到 q，正好是 GRCh37 → GRCh38；反向用 `reverse=True`。

---

## 二、构建：GRCh37 还是 GRCh38

这是最容易出错的地方。**一张图里有四个坐标系，必须各自声明清楚**：

| 角色 | 由谁决定 | 是否必须等于 `build_label` |
|---|---|---|
| GWAS 汇总统计 | 文件本身（看它的 `##genomeReference`） | ❌ 可以不同，用 `liftover.chain` 桥接 |
| 位点表 `start` / `end` | 你的配置 | ✅ **必须**是 `build_label` 的坐标 |
| 基因 GTF / 调控 GFF | 你的配置 | ✅ **必须**与 `build_label` 一致 |
| LD 参考面板 | 文件本身 | ❌ 可以不同，但必须 `match_by: rsid` |

**工具不会自动猜构建** —— 必须由你声明。`build_label` 只影响坐标轴标签与注释轨道的解读，
不会改变任何输入文件的坐标。

### 2.1 四种组合

| GWAS | 图上要显示 | 需要的配置 |
|---|---|---|
| GRCh37 | GRCh37 | 全默认：`build_label: GRCh37`、无 `liftover`、`ld.match_by: position` |
| **GRCh37** | **GRCh38** | `build_label: GRCh38` + `liftover.chain` + **GRCh38 的 GTF/GFF** + `ld.match_by: rsid` |
| GRCh38 | GRCh38 | 无 `liftover`；若 LD 面板还是 GRCh37 → `ld.match_by: rsid` |
| GRCh38 | GRCh37 | `build_label: GRCh37` + `liftover.chain` + `liftover.reverse: true` + GTF/GFF 用 37 |

现成例子：`config_scz_v2.yaml`（第 1 种）、`config_nxph1_grch38.yaml` 与
`config_enc_scz.yaml`（第 2 种，PGC 的 37 数据 + 38 注释 + rsID 匹配 LD）。

### 2.2 位点表也要换构建

`loci` 的坐标必须已经是 `build_label`。如果手上只有另一套构建的位点表，先转一下：

```python
from locuskit import Chain
c = Chain("data/001_GRCh37_DNA/GRCh37_to_GRCh38.chain.gz")   # t=37 → q=38

rows = [l for l in loci if l.get("start") and l.get("end")]
for l in rows:
    l["start"] = c.lift(l["chr"], l["start"])
    l["end"]   = c.lift(l["chr"], l["end"])
```

（`lift` 返回 `None` 表示该点不在 chain 块内，这类位点要么换端点，要么整段跳过。）

### 2.3 怎么判断一个文件是哪版

1. 看头部：PGC3 SCZ 是 `##genomeReference=GRCh37`，PGC-MDD 是 `##genomeReference="GRCh37"`；
2. 拿一个已知 SNP 对坐标：`rs3094315` 在 GRCh37 是 `chr1:752,566`，在 GRCh38 是 `chr1:817,186`；
3. 差值几十 kb 的成对位置（如 `chr7:8,680,443` vs `chr7:8,640,813`）基本可判定为 37↔38；
4. **最硬的验证**：取该位置 ±30 bp 的序列，分别到两套参考 FASTA 里比。
   例如 `chr1:752,566` 附近的 `ATGAAAGAGCT…` 在 GRCh38 只出现在 `chr1:817,186`，
   换一个候选位置序列就完全不同 —— 一次比对即可确认方向与偏移。

本项目实例：NXPH1 的 lead SNP `rs10238024`
— GRCh37 `7:8,680,443` → GRCh38 `7:8,640,813`（偏移 −39,630），
`config_nxph1_grch38.yaml` 的 loci 表写的是 GRCh38 坐标，`gwas` 仍是 GRCh37，靠 `liftover` 衔接。

### 2.4 配错了会怎样

| 症状 | 原因 |
|---|---|
| 窗口内 0 个变异、图空白 | 位点表用了另一套构建的坐标 |
| 峰不在 index SNP 附近 | GWAS 与位点表构建不一致，且没配 `liftover` |
| LD 几乎全灰，只有 index SNP 有色 | 面板与数据不同构建，却用了 `match_by: position` |
| 基因名对不上 / 轨道空 | GTF 构建与 `build_label` 不一致 |
| 所有位置整体偏 ~40 kb | 就是 37↔38 的差异，换注释并补 `liftover` |

---

## 三、怎么用

三种入口等价：PDF 首页是图，其后是附表；另存 PNG / SVG。

### 3.1 网页 GUI

```bash
./plot gui                  # 等价于 python3 -m locuskit.gui
# → http://127.0.0.1:8088
```

选配置 → 填 index SNP / 窗口 → 勾选轨道 → 选附表 → 出图。
右侧实时预览 PNG，下方给 PDF / PNG / SVG 下载，底部是运行日志。
参考数据与扫描结果有缓存，同一位点第二次约 1–2 秒。

### 3.2 命令行

```bash
python3 -m locuskit.cli config.yaml                   # 配置里的全部位点
python3 -m locuskit.cli config.yaml --only rs4702     # 单个
python3 -m locuskit.cli config.yaml --limit 3         # 前 3 个
python3 -m locuskit.cli config.yaml --list            # 只列位点
python3 -m locuskit.cli config.yaml --no-table        # 不出附表
python3 -m locuskit.cli config.yaml --outdir out --formats pdf,svg
```

现成配置（`code/` 下，也是 GUI 下拉里的四项）：

| 配置 | 数据 | 构建 | 特点 |
|---|---|---|---|
| `config_scz_v2.yaml` | PGC3 SCZ core | GRCh37 | LD + credible 圈 + 基因 + 调控 |
| `config_enc_scz.yaml` | SCZ european | GRCh38 | liftover + rsID 匹配 LD + **ENCODE 五轨** + `r²>0.8` 附表 |
| `config_mdd_v2.yaml` | PGC-MDD 2025 EUR | GRCh37 | 无 credible set，`auto_loci` 自动识别 |
| `config_nxph1_grch38.yaml` | MDD EUR · NXPH1 | GRCh38 | liftover + 显著 SNP 附表（P<5e-8） |

### 3.3 桌面可执行文件

不想装 Python 时，可以打包成单文件可执行程序（`locuskit-gui` / `locuskit-gui.exe`）：

```bash
bash packaging/build.sh      # Linux / macOS  → dist/locuskit-gui
packaging\build.bat          # Windows        → dist\locuskit-gui.exe
```

把产物放到含 `code/` 与 `data/` 的目录后双击，浏览器会自动打开
`http://127.0.0.1:8088/`。**数据不打进 exe**（通常几十 GB）。
构建细节、跨平台限制与踩坑记录见 `packaging/README.md`。

### 3.4 Python 库

```python
from locuskit import (scan_windows, apply_liftover, Chain, LDPanel,
                      GeneIndex, RegFeatures, plot_locus, detail_pages, style)

style()                                                   # matplotlib 样式

# 读窗口（只留窗口内变异，内存与窗口大小成正比）
recs, meta = scan_windows(
    "PGC3.tsv.gz",
    {"rs4702": ("5", 136_900_000, 137_900_000)},
    columns={"chr": "CHROM", "pos": "POS", "id": "ID", "p": "PVAL",
             "a1": "A1", "a2": "A2", "beta": "BETA", "se": "SE"})
r = recs["rs4702"]

r = apply_liftover(r, Chain("GRCh37_to_GRCh38.chain.gz"))   # 可选：换构建

ld = LDPanel("1000G.EUR.QC.bed")
genes = GeneIndex.from_gtf("GRCh38.116.gtf.gz", cache="cache/genes.pkl")
reg = RegFeatures.from_gff("regulatory_features.v112.gff3.gz", cache="cache/reg.pkl")

m = ld.lookup_rsid([x["id"] for x in r] + ["rs4702"])       # r²
r2 = ld.r2_to_index(m["rs4702"], [m.get(x["id"]) for x in r])

lo, hi = r[0]["pos"], r[-1]["pos"]
fig, info = plot_locus(r, chrom="5", index_snp="rs4702",
                       genes=genes.window("5", lo, hi), reg=reg.window("5", lo, hi),
                       r2=r2, highlight_genes={"FURIN"},
                       build_label="GRCh38", figsize_mm=(90, 92))
fig.savefig("FURIN_rs4702.png", dpi=300)

pages = detail_pages([dict(x, r2=float(v)) for x, v in zip(r, r2) if v > 0.8],
                     title="rs4702 — high-LD SNPs")
```

---

## 四、写配置

### 4.1 最小配置

```yaml
gwas:
  path: data/003_SCZ_Nature2022/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv.gz
  columns: {chr: CHROM, pos: POS, id: ID, p: PVAL}

loci:
  path: my_loci.tsv
  columns: {index_snp: index_snp, chr: chrom, start: start, end: end}

output: {dir: out}
```

### 4.2 完整配置

```yaml
build_label: GRCh38              # 图上显示的构建
cache_dir: cache                 # 基因/peak 索引缓存

gwas:
  path: ...
  comment_prefix: "##"
  columns: {chr: CHROM, pos: POS, id: ID, p: PVAL, a1: A1, a2: A2,
            beta: BETA, se: SE, info: IMPINFO, n: NEFF, fcas: FCAS, fcon: FCON}

loci:                            # 三选一：loci / credset / auto_loci
  path: my_loci.tsv
  columns: {index_snp: index_snp, chr: chrom, start: start, end: end, genes: genes}

credset:                         # 可选
  path: table2_credible.tsv
  columns: {index_snp: index_snp, id: rsid, chr: chromosome, pos: position,
            pip: finemap_posterior_probability}

liftover:                        # 汇总统计与 build_label 不同构建时
  chain: data/001_GRCh37_DNA/GRCh37_to_GRCh38.chain.gz
  reverse: false                 # chain 方向与预期相反时置 true

ld:
  bed: 1000G.EUR.QC.bed
  match_by: rsid                 # position（同构建）| rsid（跨构建）

genes: {gtf: Homo_sapiens.GRCh38.116.gtf.gz, rows: 5,
        biotypes: null,          # 例 [protein_coding] 只看编码基因
        name_style: short}       # 无 symbol 的基因显示 …229970

regulatory:
  gff: Homo_sapiens.GRCh38.regulatory_features.v112.gff3.gz
  types: [promoter, enhancer, CTCF_binding_site, open_chromatin_region]

encode:                          # dlPFC 证据轨道 + 附表证据列
  show: true
  peaks: .../dlPFC_panel/processed_data/01_peaks
  evidence: .../dlPFC_panel/snp_evidence/snp_evidence_wide.tsv
  marks: [H3K27ac, H3K4me3, H3K27me3, CTCF, DNase]
  datasets: 20

window: {flank_kb: 150, min_kb: 200, max_kb: 1500}
figure: {width_mm: 90, height_mm: 92, dpi: 300}

output:
  dir: figures
  formats: pdf                   # pdf / pdf,svg / png …
  preview: true                  # 额外存 PNG
  gene_prefix: true              # 文件名用「基因名_位点」

table:
  enabled: true
  source: ld                     # credset | top | window | ld | none
  ld_cut: 0.8                    # source=ld 时的 r² 阈值
  top_n: 50                      # source=top 时取前 N
  p_max: null                    # 例 5e-8 只保留显著变异
  page: A4L                      # A4L | A4 | A3L
  fontsize: 5.8
```

### 4.3 缺省行为

| 段 | 不写会怎样 |
|---|---|
| `gwas` | 必须有；其中只有 `chr/pos/p` 是硬要求 |
| `loci` | 没写就要求 `auto_loci.enabled: true` |
| `credset` | 无空心圈、无 PIP 着色 |
| `liftover` | 不做坐标变换（汇总统计与图上视为同一构建） |
| `ld` | 纯灰点，无 r² 色标、无 `r²>阈值` 附表 |
| `genes` | 无基因轨道 |
| `regulatory` | 无调控轨道 |
| `encode` | 无证据轨道、附表无证据列 |
| `table` | `enabled: false` 等价于不出附表 |

**任一可选项缺失，对应轨道/列自动消失，不报错。**

---

## 五、常见操作

| 想做的事 | 怎么做 |
|---|---|
| 换数据集 | 改 `gwas.path` + `gwas.columns` |
| 换构建 | 见第二节；要点：`build_label` + `liftover.chain` + 换 GTF/GFF + `ld.match_by: rsid` |
| 调窗口 | 改 `window.flank_kb`（或 GUI 的"窗口 kb"）；有 `loci.start/end` 时以它为准 |
| 只列 r²>0.9 | `table: {source: ld, ld_cut: 0.9}` |
| 只列显著 SNP | `table: {source: window, p_max: 5e-8}` |
| 不要附表 | `table.enabled: false` 或 `--no-table` |
| 只看编码基因 | `genes.biotypes: [protein_coding]` |
| 双栏大图 | `figure: {width_mm: 180, height_mm: 190}`（版面按毫米自适应） |
| 没有位点表 | `auto_loci: {enabled: true, p_thresh: 5e-8, clump_kb: 500}` |

---

## 六、代码结构

```
locuskit/
├── liftover.py   Chain            GRCh37 ↔ GRCh38（UCSC chain，双向）
├── gwas.py       scan_windows()   汇总统计按窗口流式读取
├── ld.py         LDPanel          PLINK .bed 直算 r²（位置 / rsID 匹配）
├── annot.py      GeneIndex        Ensembl GTF 基因模型（pickle 缓存）
│                 RegFeatures      Ensembl 调控元件
│                 EncodePeaks      多数据集 peak 支持度（缓存）
│                 EvidenceTable    SNP → 各 mark 命中数
├── render.py     plot_locus()     1–4 条轨道的位点图（版面按毫米自适应）
├── tables.py     detail_pages()   附表（列可自定义，缺字段自动省略）
├── cli.py                         YAML 批量出图
└── gui.py                         Flask 单页界面
```

## 七、限制

* `ld.match_by: position` 要求面板与汇总统计同构建；跨构建请用 `rsid`。
* 用 EUR 面板给东亚数据上色属近似：面板里没有的变异保持灰点。
* `EncodePeaks` 首次读该 mark 全部 `bed.gz`（每类 20 个），之后走缓存。
* 只做"以某个 index SNP 为参照的 r²"，不含 LD 聚类或条件分析。
* 附表 `position` 列是图上所用构建的坐标；原始坐标保留在记录的 `pos_src`。
