#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Figure 7 重建：A KM（n=242，HR 参照组标注）/ B 突变-CNA（计数轴）/ C UALCAN 甲基化"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import tifffile as tif
from PIL import Image

OUT = r"C:\Users\Administrator\AppData\Local\hermes\cache\scratch\crc_redraw\out"
BOXF = r"F:\图复现代码与原始数据包_2026-11-19\A_正文图脚本与数据\Figure_7C_UALCAN甲基化\ualcan_methylation_boxstats.csv"
DST = r"F:\图复现代码与原始数据包_2026-11-19\图8_Fig7B_修正_20261119"
DPI = 300
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.8,
                     "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "xtick.direction": "out", "ytick.direction": "out",
                     "legend.frameon": False, "legend.fontsize": 6.4,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
C_HI, C_LO = "#C0392B", "#2563A8"
GENES = ["GUCA2A", "CA7", "SST"]
PSTAT = {"GUCA2A": ("p = 0.140", "HR 0.54 (Low vs High)"),
         "CA7": ("p = 0.066", "HR 0.47 (Low vs High)"),
         "SST": ("p = 0.752", "HR 0.87 (Low vs High)")}


def save(fig, stem):
    png = os.path.join(DST, stem + ".png")
    fig.savefig(png, dpi=DPI, bbox_inches="tight", pad_inches=0.02, facecolor="white")

    fig.savefig(os.path.join(DST, stem + ".pdf"), dpi=DPI, bbox_inches="tight", pad_inches=0.02, facecolor="white")
    arr = np.array(Image.open(png).convert("RGB"))
    tif.imwrite(os.path.join(DST, stem + ".tif"), arr, compression="zlib",
                photometric="rgb", resolution=(DPI, DPI))
    print("saved", stem, arr.shape)


def km(ax, g, legend=False):
    rows = [r for r in csv.DictReader(open(os.path.join(OUT, "km_steps.csv"))) if r["gene"] == g]
    for grp, col in (("grp=High", C_HI), ("grp=Low", C_LO)):
        st = [r for r in rows if r["strata"] == grp]
        t = [0.0] + [float(r["time"]) for r in st]
        s = [1.0] + [float(r["surv"]) for r in st]
        ax.step(t, s, where="post", color=col, lw=1.0,
                label=("High (n=121)" if grp.endswith("High") else "Low (n=121)"))
    ax.set_xlim(0, 75); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Overall survival (months)")
    ax.set_ylabel("Probability")
    ax.set_title(f"{g} (median split)", fontsize=7.6, pad=3)
    p, hr = PSTAT[g]
    ax.text(0.96, 0.06, f"log-rank {p}\n{hr}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.2, linespacing=1.3, color="#333333")
    if legend:
        ax.legend(loc="upper right", handlelength=1.4, borderaxespad=0.25)


def mutcna(ax, letter=True):
    rows = list(csv.DictReader(open(
        r"F:\图复现代码与原始数据包_2026-11-19\A_正文图脚本与数据\Figure_7AB_TCGA_KM突变CNA\tcga_mut_cna.csv".replace("tcga_mut", "tcga_mut"), encoding="utf-8-sig"), delimiter="\t"))
    genes = [r["gene"] for r in rows]
    mut = np.array([float(r["mutationCount"] or 0) for r in rows])
    amp = np.array([float(r["CNA_amp"]) if r["CNA_amp"] else 0 for r in rows])
    de = np.array([float(r["CNA_del"]) if r["CNA_del"] else 0 for r in rows])
    x = np.arange(len(genes)); w = 0.24
    ax.bar(x - w, mut, w, color="#8C8C8C", edgecolor="black", lw=0.5, label="Non-synonymous mutation")
    ax.bar(x, amp, w, color="#C0392B", edgecolor="black", lw=0.5, label="CNA amplification (GISTIC2)")
    ax.bar(x + w, de, w, color="#2563A8", edgecolor="black", lw=0.5, label="CNA deep deletion")
    ax.text(x[2] - w, mut[2] + 0.04, "1", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(genes)
    ax.set_ylabel("Number of samples"); ax.set_ylim(0, 1.45); ax.set_yticks([0, 1])
    ax.set_title("TCGA-COADREAD mutation/CNA (n = 257)", fontsize=7.6, pad=3)
    ax.legend(loc="upper right", handlelength=1.1, borderaxespad=0.3, fontsize=6.0)
    if letter:
        ax.text(-0.20, 1.09, "B", transform=ax.transAxes, fontsize=10, fontweight="bold")


def meth(ax, g, legend=False):
    rows = list(csv.DictReader(open(BOXF, encoding="utf-8-sig")))
    stats, labels, cols = [], [], []
    for grp, col, lab in (("Normal (n=37)", "#8C8C8C", "Normal (n=37)"),
                          ("Primary tumor (n=313)", "#E1895F", "Tumour (n=313)")):
        d = {r["stat"]: float(r["value"]) for r in rows if r["gene"] == g and r["group"] == grp}
        stats.append({"med": d["median"], "q1": d["q1"], "q3": d["q3"],
                      "whislo": d["low"], "whishi": d["high"], "fliers": []})
        labels.append(lab); cols.append(col)
    bp = ax.bxp(stats, positions=[1, 2], widths=0.5, patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.9)
    for med in bp["medians"]:
        med.set_color("black"); med.set_linewidth(1.0)
    for w in bp["whiskers"]:
        w.set_linewidth(0.8)
    ax.set_xticks([1, 2]); ax.set_xticklabels(labels, fontsize=6.4)
    ax.set_ylabel("Promoter methylation (β value)")
    med_n = [s["med"] for s in stats]
    ax.set_title(f"{g}  (Normal {med_n[0]:.3f} | Tumour {med_n[1]:.3f})", fontsize=7.2, pad=3)
    ax.set_xlim(0.4, 2.6)


def fig7():
    fig = plt.figure(figsize=(7.09, 6.45))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[1.25, 0.80, 1.05],
                  hspace=0.50, wspace=0.30,
                  left=0.085, right=0.985, top=0.945, bottom=0.075)
    axA1 = fig.add_subplot(gs[0, 0]); axA2 = fig.add_subplot(gs[0, 1]); axA3 = fig.add_subplot(gs[0, 2])
    for ax, g, leg in ((axA1, "GUCA2A", True), (axA2, "CA7", False), (axA3, "SST", False)):
        km(ax, g, legend=leg)
    axB = fig.add_subplot(gs[1, :]); mutcna(axB)
    axC1 = fig.add_subplot(gs[2, 0]); axC2 = fig.add_subplot(gs[2, 1]); axC3 = fig.add_subplot(gs[2, 2])
    for ax, g in ((axC1, "GUCA2A"), (axC2, "CA7"), (axC3, "SST")):
        meth(ax, g)
    fig.text(0.012, 0.995, "A", fontsize=10, fontweight="bold", va="top")
    fig.text(0.012, 0.315, "C", fontsize=10, fontweight="bold", va="top")
    return fig


if __name__ == "__main__":
    f = fig7(); save(f, "Figure7_v3_A-B-C")
    print("done")
