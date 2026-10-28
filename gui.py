#!/usr/bin/env python3
"""locuskit.gui — 本地网页 GUI（Flask 单页，无前端构建）。

    python -m locuskit.gui                 # 默认 127.0.0.1:8088
    python -m locuskit.gui --port 9000 --host 0.0.0.0

功能：选配置 → 填位点/窗口 → 勾选轨道 → 出图 → 页面预览 + 下载 PDF/PNG/SVG。
参考数据（GTF / LD 面板 / peak）按配置缓存；每个位点的汇总统计按窗口缓存，
所以第二次出同一个位点是秒级。
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import pickle
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from . import (Chain, EncodePeaks, EvidenceTable, GeneIndex, LDPanel, RegFeatures,
               apply_liftover, detail_pages, plot_locus, scan_windows, style)
from .cli import build_loci, load_config, read_table_dicts, window_for


def yaml_load(path):
    import yaml
    return yaml.safe_load(open(path))
from .render import REG_COLORS

try:
    from flask import Flask, jsonify, request, send_from_directory
except ImportError:                                  # pragma: no cover
    sys.exit("需要 Flask：pip install flask")

ROOT = Path(__file__).resolve().parent.parent
APP = Flask(__name__)
CACHE = ROOT / "cache" / "gui"
OUT = ROOT / "figures_gui"


def set_root(path):
    """切换工作目录（打包成可执行文件后，工作目录 = exe 所在目录）。"""
    global ROOT, CACHE, OUT
    ROOT = Path(path).resolve()
    CACHE = ROOT / "cache" / "gui"
    OUT = ROOT / "figures_gui"
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    return ROOT

_STATE = {"cfg_path": None, "cfg": None, "objs": None, "t0": 0.0}


# ---------------------------------------------------------------- 参考数据
def _load_objects(cfg_path):
    if _STATE["cfg_path"] == cfg_path and _STATE["objs"]:
        return _STATE["cfg"], _STATE["objs"]
    cfg = load_config(cfg_path)
    cache = Path(cfg["cache_dir"])
    cache.mkdir(parents=True, exist_ok=True)
    objs = {}
    if (cfg.get("ld") or {}).get("bed"):
        objs["ld"] = LDPanel(cfg["ld"]["bed"])
    if (cfg.get("genes") or {}).get("gtf"):
        objs["genes"] = GeneIndex.from_gtf(cfg["genes"]["gtf"], cache=cache / "genes.pkl")
    if (cfg.get("regulatory") or {}).get("gff"):
        objs["reg"] = RegFeatures.from_gff(cfg["regulatory"]["gff"],
                                           types=cfg["regulatory"].get("types"),
                                           cache=cache / "reg.pkl")
    if (cfg.get("encode") or {}).get("peaks"):
        objs["peaks"] = EncodePeaks(cfg["encode"]["peaks"], cfg["encode"]["marks"],
                                    cache_dir=cache / "encode")
        objs["ev"] = EvidenceTable.from_tsv(cfg["encode"].get("evidence") or "")
    if (cfg.get("liftover") or {}).get("chain"):
        objs["chain"] = Chain(cfg["liftover"]["chain"])
    _STATE.update({"cfg_path": cfg_path, "cfg": cfg, "objs": objs, "t0": time.time()})
    return cfg, objs


def _credset(cfg):
    cs = cfg.get("credset") or {}
    if not cs.get("path"):
        return {}
    out = {}
    for r in read_table_dicts(cs["path"], cs.get("columns"), key_col="index_snp"):
        try:
            pip = float(r.get("pip")) if r.get("pip") not in (None, "", "NA") else float("nan")
        except ValueError:
            pip = float("nan")
        out.setdefault(r.get("index_snp"), []).append(
            {"id": r.get("id"), "pos": int(float(r["pos"])) if r.get("pos") else None,
             "pip": pip, "chr": str(r.get("chr") or "").replace("chr", "")})
    return out


def _scan_cached(cfg, windows, force=False):
    """按 (汇总统计, 窗口) 缓存扫描结果。"""
    key = hashlib.md5(json.dumps([cfg["gwas"]["path"],
                                  {k: list(v) for k, v in windows.items()}],
                                 sort_keys=True).encode()).hexdigest()[:16]
    cf = CACHE / f"scan_{key}.pkl"
    if cf.exists() and not force:
        return pickle.loads(cf.read_bytes())
    recs, meta = scan_windows(cfg["gwas"]["path"], windows,
                              columns=cfg["gwas"].get("columns"),
                              comment_prefix=cfg["gwas"].get("comment_prefix", "##"))
    CACHE.mkdir(parents=True, exist_ok=True)
    cf.write_bytes(pickle.dumps((recs, meta), protocol=4))
    return recs, meta


# ---------------------------------------------------------------- 出图
def run(params):
    log = []
    t0 = time.time()
    cfg_path = params["config"]
    cfg, objs = _load_objects(cfg_path)
    cfg = json.loads(json.dumps(cfg))          # 拷贝，允许本次覆盖
    log.append(f"配置: {cfg_path}")

    # 覆盖项
    if params.get("window_kb"):
        kb = int(params["window_kb"])
        cfg["window"] = {"flank_kb": kb // 2, "min_kb": max(200, kb), "max_kb": 3000}
    if params.get("width_mm"):
        cfg["figure"]["width_mm"] = int(params["width_mm"])
    if params.get("height_mm"):
        cfg["figure"]["height_mm"] = int(params["height_mm"])
    style_in = None
    _sp = (params.get("style_preset") or "").strip()
    _sj = (params.get("style_json") or "").strip()
    if _sp or _sj:
        style_in = {}
        if _sp and _sp != "default":
            style_in["style_preset"] = _sp
        if _sj:
            try:
                extra = json.loads(_sj)
                if not isinstance(extra, dict):
                    raise ValueError("style_json 必须是对象 {…}")
                style_in.update(extra)
            except (ValueError, json.JSONDecodeError) as e:
                raise SystemExit(f"style_json 解析失败: {e}")
        log.append(f"样式: 预设={_sp or 'default'}"
                   + (f"，覆盖 {len(json.loads(_sj))} 项" if _sj else ""))

    tracks = set(params.get("tracks") or [])
    cfg["encode"]["show"] = "encode" in tracks
    show_reg = "reg" in tracks
    show_genes = "genes" in tracks

    # 位点
    idx = (params.get("snp") or "").strip()
    locus = None
    if idx:
        chrom = (params.get("chrom") or "").strip().replace("chr", "")
        start = params.get("start")
        end = params.get("end")
        if start and end:
            locus = {"index_snp": idx, "chr": chrom, "start": int(start), "end": int(end),
                     "genes": []}
        else:
            # 从配置文件里找同名位点（拿到区间与基因）
            cred = _credset(cfg)
            for l in build_loci(cfg, cred) if not isinstance(build_loci(cfg, cred), str) else []:
                if l["index_snp"] == idx:
                    locus = l
                    break
            if locus is None and (cfg.get("auto_loci") or {}).get("enabled") and chrom:
                # 从汇总统计自动识别，找出该 SNP 的位置
                allrecs, _m = _scan_cached(cfg, {"ALL": (chrom, 0, 1 << 62)},
                                           force=bool(params.get("rescan")))
                hit = next((r for r in allrecs.get("ALL", []) if r["id"] == idx), None)
                if hit is None:
                    raise SystemExit(f"{idx} 不在 chr{chrom} 的汇总统计里")
                locus = {"index_snp": idx, "chr": chrom, "start": None, "end": None,
                         "genes": [], "index_pos": hit["pos"]}
            if locus is None:
                raise SystemExit(f"配置里没有位点 {idx}；请填 start/end，或给带位点表的配置")
    else:
        raise SystemExit("请填 index SNP")

    cfg_cred = _credset(cfg)
    wins = {}
    lead, lo, hi = window_for(locus, None, cfg)
    if lo is None:
        raise SystemExit("无法确定窗口")
    wins[idx] = (locus["chr"], lo, hi)
    scan_wins = wins
    chain = objs.get("chain")
    rev = bool(cfg.get("liftover", {}).get("reverse", False)) if chain else False
    if chain:
        a0 = chain.lift(locus["chr"], lo, reverse=not rev)
        a1 = chain.lift(locus["chr"], hi, reverse=not rev)
        if a0 and a1:
            scan_wins = {idx: (locus["chr"], min(a0, a1) - 5000, max(a0, a1) + 5000)}
    recs, meta = _scan_cached(cfg, scan_wins, force=bool(params.get("rescan")))
    recs_l = recs.get(idx, [])
    log.append(f"扫描 {meta['rows']:,} 行 → 窗口内 {len(recs_l):,} 个变异")
    if chain:
        recs_l = apply_liftover(recs_l, chain, reverse=rev)
        recs_l = [r for r in recs_l if lo <= r["pos"] <= hi]
        log.append(f"liftover 后 {len(recs_l):,} 个变异")
    if not recs_l:
        raise SystemExit("窗口内没有变异")

    ipos = next((r["pos"] for r in recs_l if r["id"] == idx), None)

    # LD
    r2 = None
    ld = objs.get("ld")
    if ld is not None and "ld" in tracks and ipos is not None:
        if str(cfg["ld"].get("match_by", "position")).lower() == "rsid":
            m = ld.lookup_rsid([r["id"] for r in recs_l] + [idx], chrom=locus["chr"])
            if idx in m:
                ref = ld.dosage(m[idx])
                r2 = np.array([ld.r2(ref, ld.dosage(m[r["id"]])) if r["id"] in m else np.nan
                               for r in recs_l])
        else:
            m = ld.lookup_pos(locus["chr"], [r["pos"] for r in recs_l])
            if ipos in m:
                ref = ld.dosage(m[ipos])
                r2 = np.array([ld.r2(ref, ld.dosage(m[p])) if p in m else np.nan
                               for p in (r["pos"] for r in recs_l)])
        if r2 is not None:
            log.append(f"LD: {int(np.isfinite(r2).sum())}/{len(r2)} 个变异取到 r²")

    # 注释
    hl = set(locus.get("genes") or [])
    gwin = objs["genes"].window(locus["chr"], lo, hi,
                                biotypes=cfg["genes"].get("biotypes")) if (show_genes and "genes" in objs) else None
    rwin = objs["reg"].window(locus["chr"], lo, hi) if (show_reg and "reg" in objs) else None
    enc = None
    if cfg["encode"]["show"] and "peaks" in objs:
        enc = {}
        for mk in cfg["encode"]["marks"]:
            enc[mk] = [x for x in objs["peaks"].support(mk, locus["chr"])
                       if x[1] > lo and x[0] < hi]

    fig, info = plot_locus(
        recs_l, chrom=locus["chr"], index_snp=idx, index_pos=ipos,
        credset=cfg_cred.get(idx), genes=gwin, reg=rwin,
        reg_types=[(t, dict(REG_COLORS)[t]) for t in cfg["regulatory"]["types"]],
        encode=enc, encode_marks=cfg["encode"]["marks"], show_encode=cfg["encode"]["show"],
        r2=r2, highlight_genes=hl,
        gene_rows=cfg["genes"]["rows"], gene_name_style=cfg["genes"]["name_style"],
        build_label=cfg["build_label"],
        figsize_mm=(cfg["figure"]["width_mm"], cfg["figure"]["height_mm"]),
        style=style_in)

    # 附表
    tfigs = []
    src = params.get("table_source") or "none"
    if src != "none":
        rows = []
        if src == "ld" and r2 is not None:
            cut = float(params.get("ld_cut") or 0.8)
            rows = [dict(r, r2=float(v)) for r, v in zip(recs_l, r2)
                    if np.isfinite(v) and v > cut]
            rows.sort(key=lambda r: (-r["r2"], r["p"]))
        elif src == "credset":
            ids = {c["id"] for c in cfg_cred.get(idx, [])}
            rows = sorted((dict(r) for r in recs_l if r["id"] in ids), key=lambda r: r["p"])
        else:
            rows = sorted((dict(r) for r in recs_l), key=lambda r: r["p"])
        if params.get("p_max"):
            rows = [r for r in rows if r["p"] <= float(params["p_max"])]
        getters = None
        cols = None
        if "ev" in objs and objs["ev"].rows:
            getters = {c: (lambda kk: (lambda row: (objs["ev"].get(row["id"]) or {}).get(kk)))(c)
                       for _, c in EvidenceTable.FIELDS}
            cols = [("id", "rsid", None), ("pos", "position", None), ("alleles", "ea/oa", None),
                    ("maf", "maf", 3), ("beta", "BETA", 4), ("se", "SE", 4), ("p", "P", None)]
            if src == "ld":
                cols.append(("r2", "r²", 3))
            cols += [(c, c, None) for _, c in EvidenceTable.FIELDS] + [("info", "INFO", 3)]
        tfigs = detail_pages(rows, columns=cols, extra_getters=getters,
                             title=f"{idx}  —  {src}", page=cfg["table"]["page"],
                             lines=[f"window: chr{locus['chr']}:{lo:,}-{hi:,}   ·   "
                                    f"{len(recs_l):,} variants"]
                                   + ([f"{len(rows)} rows"] if rows else []),
                             source_note=f"source: {Path(cfg['gwas']['path']).name}")
        log.append(f"附表: {src} → {len(rows)} 行 / {len(tfigs)} 页")

    # 保存
    OUT.mkdir(parents=True, exist_ok=True)
    stem = idx.replace(":", "_").replace("/", "_")
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    pdf = OUT / f"{stem}.pdf"
    with PdfPages(pdf) as f:
        f.savefig(fig)
        for tf in tfigs:
            f.savefig(tf)
    png = OUT / f"{stem}.png"
    fig.savefig(png, dpi=int(cfg["figure"]["dpi"]))
    fig.savefig(OUT / f"{stem}.svg")
    for tf in tfigs:
        plt.close(tf)
    plt.close(fig)
    log.append(f"完成，用时 {time.time() - t0:.1f}s")
    return {"ok": True, "png": f"/out/{png.name}", "pdf": f"/out/{pdf.name}",
            "svg": f"/out/{stem}.svg", "log": log,
            "info": {"window": f"chr{locus['chr']}:{lo:,}-{hi:,}",
                     "variants": len(recs_l), "genes": info["n_genes"],
                     "labelled": info["n_labelled"],
                     "pages": 1 + len(tfigs)}}


# ---------------------------------------------------------------- 页面
PAGE = """<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>locuskit GUI</title>
<style>
 body{font:13px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif;margin:0;background:#f6f7f9;color:#222}
 .wrap{display:flex;gap:16px;padding:16px;align-items:flex-start}
 .panel{background:#fff;border:1px solid #e3e5e8;border-radius:10px;padding:14px 16px;box-shadow:0 1px 2px rgba(0,0,0,.04)}
 .left{width:330px;flex:none} .right{flex:1;min-width:420px}
 h1{font-size:15px;margin:0 0 12px} h2{font-size:12px;margin:14px 0 6px;color:#555;text-transform:uppercase;letter-spacing:.04em}
 label{display:block;font-size:12px;color:#444;margin:8px 0 3px}
 input,select{width:100%;box-sizing:border-box;padding:6px 8px;border:1px solid #d5d8dc;border-radius:6px;font-size:12px;background:#fff}
 .row{display:flex;gap:8px} .row>*{flex:1}
 .chk{display:flex;flex-wrap:wrap;gap:6px 12px;font-size:12px;color:#333}
 .chk label{display:flex;align-items:center;gap:5px;margin:0;width:auto}
 button{width:100%;padding:9px;border:0;border-radius:8px;background:#2f6f4f;color:#fff;font-size:13px;font-weight:600;cursor:pointer;margin-top:14px}
 button:disabled{background:#9bb3a6;cursor:default}
 img{max-width:100%;border:1px solid #e3e5e8;border-radius:8px;background:#fff}
 .links a{display:inline-block;margin-right:12px;font-size:12px;color:#2f6f4f;text-decoration:none;border-bottom:1px dashed #9ec3ae}
 pre{background:#1e2327;color:#d7e0d9;padding:10px;border-radius:8px;font-size:11px;line-height:1.45;max-height:190px;overflow:auto;white-space:pre-wrap}
 .meta{font-size:12px;color:#555;margin:8px 0}
 .hint{font-size:11px;color:#888;margin-top:4px}
</style></head><body><div class="wrap">
<div class="panel left">
 <h1>locuskit · 位点图</h1>
 <label>配置文件</label>
 <select id="config">__CONFIGS__</select>
 <h2>位点</h2>
 <label>index SNP</label><input id="snp" placeholder="rs72986630">
 <div class="row"><div><label>染色体</label><input id="chrom" placeholder="19"></div>
 <div><label>窗口 (kb)</label><input id="window_kb" value="600"></div></div>
 <div class="row"><div><label>start（可选）</label><input id="start"></div>
 <div><label>end（可选）</label><input id="end"></div></div>
 <div class="hint">不填 start/end 时用配置文件里的区间或 SNP ± 窗口</div>
 <h2>轨道</h2>
 <div class="chk">
  <label><input type="checkbox" id="t_ld" checked>LD r²</label>
  <label><input type="checkbox" id="t_genes" checked>基因</label>
  <label><input type="checkbox" id="t_reg" checked>调控元件</label>
  <label><input type="checkbox" id="t_encode">ENCODE 证据</label>
 </div>
 <h2>附表</h2>
 <div class="row"><div><label>来源</label><select id="table_source">
  <option value="none">不要</option><option value="ld">LD &gt; 阈值</option>
  <option value="credset">credible set</option><option value="window">窗口全部</option>
 </select></div><div><label>LD 阈值</label><input id="ld_cut" value="0.8"></div></div>
 <div class="row"><div><label>P 上限（可选）</label><input id="p_max" placeholder="5e-8"></div></div>
 <h2>样式</h2>
 <div class="row"><div><label>预设</label><select id="style_preset">
  <option value="default">default</option><option value="minimal">minimal</option>
  <option value="poster">poster</option><option value="dark">dark</option>
 </select></div></div>
 <label>自定义（JSON，覆盖预设）</label>
 <textarea id="style_json" rows="3" placeholder='{"point_marker":"s","ld_cmap":["#ffffff","#08306b"],"vline_style":"--"}'
   style="width:100%;box-sizing:border-box;padding:6px 8px;border:1px solid #d5d8dc;border-radius:6px;font:11px/1.4 monospace"></textarea>
 <h2>画布 (mm)</h2>
 <div class="row"><div><label>宽</label><input id="width_mm" value="90"></div>
 <div><label>高</label><input id="height_mm" value="92"></div></div>
 <label style="margin-top:10px"><input type="checkbox" id="rescan"> 忽略缓存重扫</label>
 <button id="go">出图</button>
 <div class="hint" id="status"></div>
</div>
<div class="panel right">
 <h2>预览</h2>
 <div class="meta" id="meta">—</div>
 <img id="prev" style="display:none">
 <div class="links" id="links"></div>
 <h2>日志</h2>
 <pre id="log">等待出图…</pre>
</div></div>
<script>
const $ = id => document.getElementById(id);
const tracks = () => ['ld','genes','reg','encode'].filter(t => $('t_'+t).checked);
$('go').onclick = async () => {
  const body = {
    config: $('config').value, snp: $('snp').value.trim(), chrom: $('chrom').value.trim(),
    window_kb: $('window_kb').value, start: $('start').value, end: $('end').value,
    tracks: tracks(), table_source: $('table_source').value,
    ld_cut: $('ld_cut').value, p_max: $('p_max').value,
    width_mm: $('width_mm').value, height_mm: $('height_mm').value,
    style_preset: $('style_preset').value, style_json: $('style_json').value,
    rescan: $('rescan').checked
  };
  $('go').disabled = true; $('status').textContent = '运行中…'; $('log').textContent = '';
  try {
    const r = await fetch('/plot', {method:'POST', headers:{'Content-Type':'application/json'},
                                   body: JSON.stringify(body)});
    const d = await r.json();
    $('log').textContent = (d.log || d.error || '').join ? (d.log||[]).join('\\n') : (d.error||'');
    if (!d.ok) { $('status').textContent = '失败'; return; }
    $('prev').src = d.png + '?t=' + Date.now(); $('prev').style.display = 'block';
    $('links').innerHTML = `<a href="${d.pdf}" target="_blank">下载 PDF（含附表）</a>
      <a href="${d.png}" target="_blank">PNG</a><a href="${d.svg}" target="_blank">SVG</a>`;
    $('meta').textContent = `${d.info.window} · ${d.info.variants} 变异 · ${d.info.genes} 基因`
      + ` (${d.info.labelled} 标注) · ${d.info.pages} 页`;
    $('status').textContent = '完成';
  } catch (e) { $('log').textContent = String(e); $('status').textContent = '失败'; }
  finally { $('go').disabled = false; }
};
</script></body></html>"""


@APP.route("/")
def index():
    cfgs = []
    code_dir = ROOT / "code"
    for p in sorted(code_dir.glob("config*.yaml")) if code_dir.is_dir() else []:
        try:
            if "cache_dir" in (yaml_load(p) or {}):
                cfgs.append(str(p.relative_to(ROOT)))
        except Exception:                      # noqa: BLE001
            continue
    opts = "".join(f'<option value="{c}">{c}</option>' for c in cfgs)
    return PAGE.replace("__CONFIGS__", opts)


@APP.route("/plot", methods=["POST"])
def plot_endpoint():
    params = request.get_json(force=True)
    try:
        return jsonify(run(params))
    except SystemExit as e:
        return jsonify({"ok": False, "log": [str(e)], "error": str(e)})
    except Exception as e:                            # noqa: BLE001
        return jsonify({"ok": False, "log": traceback.format_exc().splitlines(),
                        "error": f"{type(e).__name__}: {e}"})


@APP.route("/out/<path:name>")
def out_file(name):
    return send_from_directory(OUT, name)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8088)
    a = ap.parse_args(argv)
    style()
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"locuskit GUI → http://{a.host}:{a.port}", flush=True)
    APP.run(host=a.host, port=a.port, threaded=True, debug=False)


if __name__ == "__main__":
    main()
