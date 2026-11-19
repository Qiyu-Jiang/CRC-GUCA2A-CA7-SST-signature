#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""图8 重绘（A 绝对信息流 top20 / B Top12 差异，方向标签修正）
   + S15 相对信息流（全通路，归一化份额）
   + S17 Top12（1.52 宽高比）
   + 7B / S13 计数轴修正
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tifffile as tif
from PIL import Image

PKG = r"F:\图复现代码与原始数据包_2026-11-19"
CSV8 = os.path.join(PKG, "A_正文图脚本与数据", "Figure_8AB_CellChat", "pathway_diff_EpiSigHigh_vs_Low.csv")
MUTCNA = os.path.join(PKG, "A_正文图脚本与数据", "Figure_7AB_TCGA_KM突变CNA", "tcga_mut_cna.csv")
DST = os.path.join(PKG, "图8_Fig7B_修正_20261119")
os.makedirs(DST, exist_ok=True)
DPI = 300

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7,
    "axes.linewidth": 0.8, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "xtick.direction": "out", "ytick.direction": "out",
    "legend.frameon": False, "legend.fontsize": 6.5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
C_HI, C_LO = "#D96C5F", "#4FB0AE"


def save(fig, stem):
    png = os.path.join(DST, stem + ".png")
    fig.savefig(png, dpi=DPI, bbox_inches="tight", pad_inches=0.02, facecolor="white")

    fig.savefig(os.path.join(DST, stem + ".pdf"), dpi=DPI, bbox_inches="tight", pad_inches=0.02, facecolor="white")
    arr = np.array(Image.open(png).convert("RGB"))
    tif.imwrite(os.path.join(DST, stem + ".tif"), arr, compression="zlib",
                photometric="rgb", resolution=(DPI, DPI))
    print("saved", stem, arr.shape)


def load8():
    rows = list(csv.DictReader(open(CSV8, encoding="utf-8-sig")))
    d = [(r["pathway"], float(r["high_total"]), float(r["low_total"]),
          float(r["total_diff"]), float(r["epiHigh_out"]), float(r["epiLow_out"]))
         for r in rows]
    return d


# ---------------------------------------------------------------- 8A absolute
def panel_A(ax, d, topn=20):
    d = sorted(d, key=lambda r: max(r[1], r[2]), reverse=True)[:topn]
    labels = [r[0] for r in d]
    hi = np.array([r[1] for r in d]); lo = np.array([r[2] for r in d])
    y = np.arange(len(d))[::-1]
    h = 0.36
    ax.barh(y + h / 2, hi, h, color=C_HI, edgecolor="none", label="EpiSigHigh")
    ax.barh(y - h / 2, lo, h, color=C_LO, edgecolor="none", label="EpiSigLow")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Information flow")
    ax.set_ylim(-0.7, len(d) - 0.3)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="lower right", handlelength=1.2, borderaxespad=0.2)
    ax.set_title("Total information flow (top 20 pathways)", fontsize=8, pad=4)


# ---------------------------------------------------------------- 8B differences
def panel_B(ax, d, topn=12):
    d = sorted(d, key=lambda r: abs(r[3]), reverse=True)[:topn]
    d = sorted(d, key=lambda r: abs(r[3]))
    labels = [r[0] for r in d]; dif = np.array([r[3] for r in d])
    y = np.arange(len(d))
    ax.barh(y, dif, 0.62, color="#3B6FB8", edgecolor="none")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7)
    ax.tick_params(axis="y", length=0)
    lim = max(abs(dif.min()), abs(dif.max()))
    ax.set_xlim(-lim * 1.22, 0.4)
    ax.set_xlabel("Communication difference\n(EpiSigHigh − EpiSigLow)")
    for yi, v in zip(y, dif):
        ax.text(v - lim * 0.03, yi, f"{v:.2f}", va="center", ha="right",
                fontsize=6.2, color="#333333")
    ax.set_ylim(-0.6, len(d) - 0.4)


# ---------------------------------------------------------------- S15 relative
def panel_rel(ax, d, topn=24):
    hi = np.array([r[1] for r in d]); lo = np.array([r[2] for r in d])
    sh_hi = 100 * hi / hi.sum(); sh_lo = 100 * lo / lo.sum()
    idx = np.argsort(-np.maximum(sh_hi, sh_lo))[:topn]
    labels = [d[i][0] for i in idx]
    y = np.arange(topn)[::-1]; h = 0.36
    ax.barh(y + h / 2, sh_hi[idx], h, color=C_HI, edgecolor="none", label="EpiSigHigh")
    ax.barh(y - h / 2, sh_lo[idx], h, color=C_LO, edgecolor="none", label="EpiSigLow")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Relative information flow (% of group total)")
    ax.legend(loc="lower right", handlelength=1.2, borderaxespad=0.2)
    ax.set_ylim(-0.7, topn - 0.3)
    ax.set_title(f"Relative information flow (top {topn} pathways of {len(d)})",
                 fontsize=8, pad=4)


# ---------------------------------------------------------------- 7B / S13
def mutcna(ax, letter=False, caption="TCGA-COADREAD mutation/CNA (n = 257)"):
    rows = list(csv.DictReader(open(MUTCNA, encoding="utf-8-sig"), delimiter="\t"))
    genes = [r["gene"] for r in rows]
    mut = np.array([int(r["mutationCount"] or 0) for r in rows], dtype=float)
    amp = np.array([float(r["CNA_amp"]) if r["CNA_amp"] not in ("", None) else 0 for r in rows])
    dele = np.array([float(r["CNA_del"]) if r["CNA_del"] not in ("", None) else 0 for r in rows])
    x = np.arange(len(genes)); w = 0.24
    ax.bar(x - w, mut, w, color="#8C8C8C", edgecolor="black", lw=0.5,
           label="Non-synonymous mutation")
    ax.bar(x, amp, w, color="#C0392B", edgecolor="black", lw=0.5,
           label="CNA amplification (GISTIC2)")
    ax.bar(x + w, dele, w, color="#2563A8", edgecolor="black", lw=0.5,
           label="CNA deep deletion")
    for xi, v in zip(x - w, mut):
        if v > 0:
            ax.text(xi, v + 0.05, f"{int(v)}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(genes)
    ax.set_ylabel("Number of samples")
    ax.set_ylim(0, 1.45)
    ax.set_yticks([0, 1])
    ax.set_yticks([0, 0.5, 1.0], minor=True)
    ax.set_title(caption, fontsize=8, pad=4)
    ax.legend(loc="upper right", handlelength=1.1, borderaxespad=0.3, fontsize=6.2)
    ax.text(0.5, -0.215, "No copy-number events; one SST missense SNV (R88I)",
            transform=ax.transAxes, ha="center", va="top", fontsize=6.4, color="#444444")


# ================================================================ build
def fig8(d):
    fig = plt.figure(figsize=(7.09, 4.35))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.18, 1.0], wspace=0.42,
                  left=0.155, right=0.985, top=0.925, bottom=0.175)
    axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1])
    panel_A(axA, d); panel_B(axB, d)
    fig.text(0.012, 0.965, "A", fontsize=10, fontweight="bold")
    fig.text(0.545, 0.965, "B", fontsize=10, fontweight="bold")
    return fig


def figS15(d):
    fig = plt.figure(figsize=(7.09, 4.9))
    ax = fig.add_subplot(111)
    panel_rel(ax, d)
    fig.subplots_adjust(left=0.185, right=0.975, top=0.93, bottom=0.115)
    return fig


def figS17(d):
    """全部 61 条通路的差异（棒棒糖图），与 8B 的 Top-12 互补"""
    d = sorted(d, key=lambda r: r[3])
    fig = plt.figure(figsize=(7.09, 6.6))
    ax = fig.add_subplot(111)
    y = np.arange(len(d)); v = np.array([r[3] for r in d])
    ax.hlines(y, 0, v, color="#9DB4D4", lw=0.8)
    ax.plot(v, y, "o", ms=2.6, color="#3B6FB8", mec="none")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in d], fontsize=5.0)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.8, len(d) - 0.2)
    lim = max(abs(v.min()), abs(v.max()))
    ax.set_xlim(-lim * 1.12, lim * 0.45)
    ax.set_xlabel("Communication difference (EpiSigHigh − EpiSigLow)", fontsize=8)
    ax.set_title(f"Communication differences across all {len(d)} pathways (positive = higher in EpiSigHigh)",
                 fontsize=7.2, pad=3)
    fig.subplots_adjust(left=0.235, right=0.975, top=0.945, bottom=0.085)
    return fig


def fig7B():
    fig = plt.figure(figsize=(5.6, 4.05))
    ax = fig.add_subplot(111)
    mutcna(ax)
    fig.subplots_adjust(left=0.145, right=0.97, top=0.92, bottom=0.245)
    p = ax.get_position()
    fig.text(p.x0 - 0.105, p.y1 + 0.045, "B", fontsize=10, fontweight="bold")
    return fig


def figS13():
    fig = plt.figure(figsize=(5.6, 4.05))
    ax = fig.add_subplot(111)
    mutcna(ax)
    fig.subplots_adjust(left=0.145, right=0.97, top=0.92, bottom=0.245)
    return fig


if __name__ == "__main__":
    d = load8()
    print("pathways:", len(d), "| sum hi=%.3f lo=%.3f" % (sum(r[1] for r in d), sum(r[2] for r in d)))
    for fn, stem in ((lambda: fig8(d), "Figure8_v3_A-B"),
                     (lambda: figS15(d), "FigureS15_relative_information_flow"),
                     (lambda: figS17(d), "FigureS17_Top12_pathway_diff_v2"),
                     (fig7B, "Figure7B_mutation_CNA_v2"),
                     (figS13, "FigureS13_mutation_CNA_n257_v2")):
        f = fn(); save(f, stem); plt.close(f)
    print("ALL DONE")
