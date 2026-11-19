#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figure S19: 三基因 KEGG GSEA（经典 running-ES 图），NES/padj 取自管线表 26.GSEA"""
import os, csv, io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tifffile as tif
from PIL import Image

PKG = r"F:\图复现代码与原始数据包_2026-11-19"
MTX = os.path.join(PKG, "D_共享输入数据", "merge.normalize.txt")
GMT = r"C:\Users\Administrator\AppData\Local\hermes\cache\scratch\fgsea_pkg\2_管线GSEA完整表_KEGG\c2.cp.kegg.Hs.symbols.gmt"
TAB = r"C:\Users\Administrator\AppData\Local\hermes\cache\scratch\fgsea_pkg\2_管线GSEA完整表_KEGG"
DST = os.path.join(PKG, "图S19_GSEA_20261119")
os.makedirs(DST, exist_ok=True)
DPI = 300

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.7,
                     "axes.labelsize": 7.5, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "xtick.direction": "out", "ytick.direction": "out",
                     "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42})

# ------------------------------------------------------------------ data
def load_matrix():
    with open(MTX, encoding="utf-8", errors="ignore") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        names = hdr[1:]
        keep = [i for i, n in enumerate(names) if n.endswith("_Treat")]
        genes, rows = [], []
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue
            try:
                v = np.array([float(x) if x not in ("", "NA") else 0.0 for x in p[1:]])
            except ValueError:
                continue
            genes.append(p[0]); rows.append(v[keep])
    X = np.vstack(rows)
    # avereps: 同名基因取均值
    from collections import defaultdict
    idx = defaultdict(list)
    for i, g in enumerate(genes):
        idx[g].append(i)
    uniq = sorted(idx)
    Y = np.vstack([X[idx[g]].mean(axis=0) for g in uniq])
    mask = Y.mean(axis=1) > 0
    return uniq_arr(uniq, mask), Y[mask], [names[i] for i in keep]


def uniq_arr(names, mask):
    return [n for n, m in zip(names, mask) if m]


def load_gmt():
    sets = {}
    for ln in open(GMT, encoding="utf-8", errors="ignore"):
        p = ln.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        sets[p[0]] = [g for g in p[2:] if g]
    return sets


def ranking(gene_idx, data, genes):
    v = data[gene_idx]
    med = np.median(v)
    lo = v < med; hi = ~lo
    diff = data[:, hi].mean(axis=1) - data[:, lo].mean(axis=1)
    order = np.argsort(-diff)
    return diff[order], [genes[i] for i in order]


def es_curve(stats_sorted, genes_sorted, members):
    N = len(stats_sorted)
    hits = np.isin(genes_sorted, list(members))
    Nh = hits.sum()
    if Nh == 0 or Nh == N:
        return None
    w = np.abs(stats_sorted) * hits
    nrm = w.sum()
    if nrm == 0:
        w = hits.astype(float); nrm = w.sum()
    run = np.cumsum(w) / nrm - np.cumsum(~hits) / (N - Nh)
    run = np.concatenate([[0.0], run])
    return run, hits


def table_row(gene, term):
    rows = list(csv.DictReader(open(os.path.join(TAB, f"{gene}.result.txt"),
                                  encoding="utf-8", errors="ignore"), delimiter="\t"))
    for r in rows:
        if r["ID"] == term:
            return float(r["enrichmentScore"]), float(r["NES"]), float(r["p.adjust"]), int(r["rank"])
    return None


PANELS = [("GUCA2A", "KEGG_RETINOL_METABOLISM", "GUCA2A | retinol metabolism"),
          ("CA7", "KEGG_RETINOL_METABOLISM", "CA7 | retinol metabolism"),
          ("SST", "KEGG_CELL_CYCLE", "SST | cell cycle"),
          ("CA7", "KEGG_ECM_RECEPTOR_INTERACTION", "CA7 | ECM-receptor interaction")]


def draw(ax_top, ax_tick, ax_rank, stats, genes, members, color, title, es_tab, nes, padj):
    out = es_curve(stats, genes, members)
    run, hits = out
    x = np.arange(len(run))
    ax_top.plot(x, run, color=color, lw=1.0)
    ax_top.axhline(0, color="#888888", lw=0.5, ls=":")
    imax = int(np.argmax(np.abs(run)))
    ax_top.plot([imax, imax], [0, run[imax]], color="#555555", lw=0.6, ls="--")
    ax_top.set_ylabel("Enrichment score")
    ax_top.set_xticks([]); ax_top.set_xlim(0, len(run) - 1)
    ax_top.set_title(title, fontsize=7.4, pad=3)
    ax_top.text(0.03, 0.05, f"NES = {nes:+.2f}\nP$_{{adj}}$ = {padj:.1e}",
                transform=ax_top.transAxes, fontsize=6.4, va="bottom", ha="left",
                linespacing=1.3, color="#222222")

    idx = np.where(hits)[0]
    ax_tick.vlines(idx, 0, 1, color=color, lw=0.5)
    ax_tick.set_yticks([]); ax_tick.set_xticks([]); ax_tick.set_xlim(0, len(run) - 1)
    ax_tick.set_ylim(0, 1)

    ax_rank.fill_between(np.arange(len(stats)), 0, stats, color="#BFBFBF", lw=0)
    ax_rank.set_xlim(0, len(stats) - 1)
    ax_rank.set_ylabel("Ranked metric")
    ax_rank.set_xlabel("Rank in ordered gene list", fontsize=6.8)
    ax_rank.tick_params(axis="x", labelsize=5.5)


def build():
    genes, data, samples = load_matrix()
    sets = load_gmt()
    print("matrix:", data.shape, "sets:", len(sets))
    fig = plt.figure(figsize=(7.09, 6.2))
    gs = GridSpec(7, 2, figure=fig, height_ratios=[1.35, 0.22, 0.95, 0.30, 1.35, 0.22, 0.95],
                  hspace=0.055, wspace=0.28,
                  left=0.105, right=0.975, top=0.935, bottom=0.085)
    for k, (gene, term, title) in enumerate(PANELS):
        col = k % 2
        base = 0 if k < 2 else 4
        ax_top = fig.add_subplot(gs[base, col])
        ax_tick = fig.add_subplot(gs[base + 1, col], sharex=ax_top)
        ax_rank = fig.add_subplot(gs[base + 2, col], sharex=ax_top)
        gi = genes.index(gene)
        stats, gorder = ranking(gi, data, genes)
        row = table_row(gene, term)
        es_tab, nes, padj, rk = row
        out = es_curve(stats, gorder, sets[term])
        es_mine = float(out[0][np.argmax(np.abs(out[0]))])
        print(f"{gene:8s}{term:44s} ES_table={es_tab:+.5f} ES_mine={es_mine:+.5f} NES={nes:+.3f}")
        color = "#C0392B" if nes > 0 else "#2563A8"
        draw(ax_top, ax_tick, ax_rank, stats, gorder, sets[term], color, title, es_tab, nes, padj)
        fig.text(0.012 if col == 0 else 0.520, 0.992 - (k // 2) * 0.495,
                 ("A" if k == 0 else "B" if k == 1 else "C" if k == 2 else "D"),
                 fontsize=10, fontweight="bold", va="top")
    return fig


def save(fig, stem):
    png = os.path.join(DST, stem + ".png")
    fig.savefig(png, dpi=DPI, bbox_inches="tight", pad_inches=0.045, facecolor="white")
    fig.savefig(os.path.join(DST, stem + ".pdf"), dpi=DPI, bbox_inches="tight",
                pad_inches=0.045, facecolor="white")
    arr = np.array(Image.open(png).convert("RGB"))
    tif.imwrite(os.path.join(DST, stem + ".tif"), arr, compression="zlib",
                photometric="rgb", resolution=(DPI, DPI))
    print("saved", stem, arr.shape)


if __name__ == "__main__":
    f = build()
    save(f, "FigureS19_GSEA_3genes_KEGG")
    print("done")
