#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""CRC 稿件 图4(A/B/C) 与 图5(A/B/C) 出版级重绘
输入: crc_redraw/out/*.csv (由 regen.R 生成)
输出: F:/图复现代码与原始数据包_2026-11-19/图4_图5_重绘_20261119/
      整图 + 分面板, 每张 PNG(300dpi) + TIFF(LZW, 300dpi)
"""
import os, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle
import tifffile as tif
from PIL import Image

OUT_CSV = r"C:\Users\Administrator\AppData\Local\hermes\cache\scratch\crc_redraw\out"
DST = r"F:\图复现代码与原始数据包_2026-11-19\图4_图5_重绘_20261119"
os.makedirs(DST, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 2.6,
    "ytick.major.size": 2.6,
    "legend.frameon": False,
    "legend.fontsize": 6.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
    "axes.spines.top": True,
    "axes.spines.right": True,
})
C_RED, C_BLUE, C_GREY = "#C0392B", "#2563A8", "#7F8C8D"
C_VENN = ("#D95F5F", "#4A7EBB", "#E8A33D")
DPI = 300


def read_csv(name):
    with open(os.path.join(OUT_CSV, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fl(rows, key):
    return np.array([float(r[key]) for r in rows])


# ---------------------------------------------------------------- panel A (4A)
def draw_4A(ax):
    d = read_csv("fig4A_cvfit.csv")
    mn = {r["metric"]: float(r["value"]) for r in read_csv("fig4A_min.csv")}
    x, cvm, cvsd, nz = fl(d, "loglambda"), fl(d, "cvm"), fl(d, "cvsd"), fl(d, "nzero")

    ax.errorbar(x, cvm, yerr=cvsd, fmt="none", ecolor="#9E9E9E",
                elinewidth=0.7, capsize=1.8, capthick=0.7, zorder=2)
    ax.plot(x, cvm, "-", color=C_RED, lw=0.8, zorder=2)
    ax.plot(x, cvm, "o", ms=2.6, mfc=C_RED, mec=C_RED, mew=0.3, zorder=3)

    for key in ("lambda.min", "lambda.1se"):
        ax.axvline(math.log(mn[key]), color="black", lw=0.6, ls=(0, (2.5, 2.5)), zorder=1)

    ax.set_xlabel("log(λ)")
    ax.set_ylabel("Binomial deviance")
    ax.set_xlim(x.min() - 0.15, x.max() + 0.15)   # 左＝小 λ（多变量），右＝大 λ
    lo = float(np.min(cvm - cvsd)); hi = float(np.max(cvm + cvsd))
    ax.set_ylim(lo - 0.03, hi + 0.03)
    ax.set_xticks([-8, -6, -4, -2])

    # 顶部非零系数个数（在 ~15 个等距位置上标注，避免重叠）
    idx = np.linspace(0, len(x) - 1, 15).round().astype(int)
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    top.set_xticks(x[idx])
    top.set_xticklabels([f"{int(nz[i])}" for i in idx], fontsize=5.2)
    top.tick_params(axis="x", length=0, pad=1.2)
    top.spines["top"].set_linewidth(0.8)
    for s in ("bottom", "left", "right"):
        top.spines[s].set_visible(False)
    top.set_xlabel("")
    ax.text(0.0, 1.115, "Non-zero coefficients", transform=ax.transAxes,
            fontsize=5.6, ha="left", va="bottom")


# ---------------------------------------------------------------- panel B (4B)
def draw_4B(ax):
    R = 0.74
    cA, cB, cC = (-0.40, 0.26), (0.40, 0.26), (0.00, -0.16)
    for (cx, cy), col in zip((cA, cB, cC), C_VENN):
        ax.add_patch(Circle((cx, cy), R, facecolor=col, edgecolor="none", alpha=0.42, zorder=1))
        ax.add_patch(Circle((cx, cy), R, facecolor="none", edgecolor=col, lw=0.7, zorder=2))

    lab = dict(size=7.6)
    ax.text(cA[0] - 0.30, cA[1] + R + 0.10, "LASSO", **lab, color=C_VENN[0],
            ha="center", va="bottom", fontweight="bold")
    ax.text(cB[0] + 0.30, cB[1] + R + 0.10, "RF", **lab, color=C_VENN[1],
            ha="center", va="bottom", fontweight="bold")
    ax.text(0.00, cC[1] - R - 0.10, "SVM-RFE", **lab, color=C_VENN[2],
            ha="center", va="top", fontweight="bold")

    s = dict(fontsize=7.4, ha="center", va="center", zorder=4)
    ax.text(-0.74, 0.52, "0", **s)
    ax.text(0.74, 0.52, "0", **s)
    ax.text(0.00, -0.70, "5\n(45.5%)", fontsize=6.4, ha="center", va="center", zorder=4)
    ax.text(0.00, 0.78, "0", **s)
    ax.text(-0.44, -0.26, "3\n(27.3%)", fontsize=6.4, ha="center", va="center", zorder=4)
    ax.text(0.44, -0.26, "0", **s)
    ax.text(0.00, 0.06, "3\n(27.3%)", fontsize=6.4, ha="center", va="center",
            zorder=4, fontweight="bold")

    ax.text(0.00, -1.30, "Intersection: GUCA2A / CA7 / SST", fontsize=6.6,
            ha="center", va="top", color="#333333")
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.55, 1.22)
    ax.set_aspect("equal"); ax.axis("off")


# ---------------------------------------------------------------- panel C (4C)
def draw_4C(ax_l, ax_r):
    d = read_csv("fig4C_freq.csv")
    genes = [r["gene"] for r in d]
    lf, rf = fl(d, "LASSO_freq"), fl(d, "RF_freq")
    xi = np.arange(len(genes))

    for ax, vals, col, ttl in ((ax_l, lf, "#3B6FB8", "LASSO"),
                               (ax_r, rf, "#B04A3A", "Random forest")):
        ax.bar(xi, vals, width=0.68, color=col, edgecolor="black", linewidth=0.5, zorder=2)
        ax.set_ylim(0, 1.06)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticks(xi)
        ax.set_xticklabels(genes, rotation=90, fontsize=5.6)
        ax.set_ylabel("Selection frequency")
        ax.set_title(f"{ttl}  (100 × 80% subsamples)", fontsize=7.2, fontweight="bold", pad=3)
        ax.set_xlim(-0.6, len(genes) - 0.4)
        ax.grid(axis="y", lw=0.35, color="#DDDDDD", zorder=0, alpha=0.8)


# ---------------------------------------------------------------- ROC helper
def draw_roc(ax, series, xlab=True, ylab=True, legend=False):
    """series: list of (rows, color, linestyle, label)"""
    for rows, col, lsty, lbl in series:
        fpr, tpr = fl(rows, "fpr"), fl(rows, "tpr")
        o = np.argsort(fpr, kind="mergesort")        # pROC 输出为 fpr 递减，需按 fpr 升序
        fpr, tpr = fpr[o], tpr[o]
        ax.plot(np.r_[0, fpr, 1], np.r_[0, tpr, 1], color=col, lw=1.1, ls=lsty,
                drawstyle="steps-post", label=lbl, zorder=3)
    ax.plot([0, 1], [0, 1], color="#AAAAAA", lw=0.6, ls=(0, (3, 3)), zorder=1)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_aspect("equal")
    if xlab: ax.set_xlabel("1 − Specificity")
    if ylab: ax.set_ylabel("Sensitivity")
    if legend:
        ax.legend(loc="lower right", handlelength=1.6, borderaxespad=0.3,
                  labelspacing=0.25, fontsize=6.2)


# ---------------------------------------------------------------- save helpers
def save(fig, stem):
    png = os.path.join(DST, stem + ".png")
    fig.savefig(png, dpi=DPI, bbox_inches="tight", pad_inches=0.045, facecolor="white")

    fig.savefig(os.path.join(DST, stem + ".pdf"), dpi=DPI, bbox_inches="tight", pad_inches=0.045, facecolor="white")
    arr = np.array(Image.open(png).convert("RGB"))
    tif.imwrite(os.path.join(DST, stem + ".tif"), arr, compression="zlib",
                photometric="rgb", resolution=(DPI, DPI))
    print("saved", stem, arr.shape)
    return png


def letter(fig, ax, ch, dx=-0.075, dy=0.022, size=10):
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y0 + bb.height + dy, ch, fontsize=size,
             fontweight="bold", ha="left", va="top")


# ================================================================ 图4 整图
def fig4_composite():
    fig = plt.figure(figsize=(7.09, 4.30))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 0.80], wspace=0.40,
                  hspace=0.62, left=0.085, right=0.985, top=0.905, bottom=0.175)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    draw_4A(axA); draw_4B(axB)
    sub = gs[1, :].subgridspec(1, 2, wspace=0.34)
    axC1 = fig.add_subplot(sub[0, 0]); axC2 = fig.add_subplot(sub[0, 1])
    draw_4C(axC1, axC2)
    fig.text(0.010, 0.985, "A", fontsize=10, fontweight="bold", va="top")
    fig.text(0.520, 0.985, "B", fontsize=10, fontweight="bold", va="top")
    fig.text(0.010, 0.425, "C", fontsize=10, fontweight="bold", va="top")
    return fig


# ================================================================ 图5 整图
def fig5_composite():
    fig = plt.figure(figsize=(6.60, 2.95))
    gs = GridSpec(1, 3, figure=fig, wspace=0.22,
                  left=0.078, right=0.988, top=0.90, bottom=0.285)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[0, 2])

    draw_roc(axA, [(read_csv("fig5A_roc.csv"), C_RED, "-", "")])
    axA.text(0.96, 0.20, "AUC 0.987\n95% CI 0.965–1.000", fontsize=6.4,
             ha="right", va="bottom", color=C_RED, linespacing=1.35)
    axA.text(0.5, -0.235, "In-sample fit (descriptive)", transform=axA.transAxes,
             fontsize=6.6, ha="center", va="top", color="#444444")

    draw_roc(axB, [(read_csv("fig5B_roc.csv"), C_BLUE, "-", "")])
    axB.text(0.96, 0.20, "AUC 0.907\nAccuracy 0.900", fontsize=6.4,
             ha="right", va="bottom", color=C_BLUE, linespacing=1.35)
    axB.text(0.5, -0.235, "10-fold cross-validation", transform=axB.transAxes,
             fontsize=6.6, ha="center", va="top", color="#444444")

    draw_roc(axC, [(read_csv("fig5C_ANN_z.csv"), C_RED, "-", None),
                   (read_csv("fig5C_GLM.csv"), C_BLUE, "-", None),
                   (read_csv("fig5C_ANN_raw.csv"), C_GREY, (0, (4, 2)), None)], legend=False)
    axC.legend([plt.Line2D([], [], color=C_RED, lw=1.1),
                plt.Line2D([], [], color=C_BLUE, lw=1.1),
                plt.Line2D([], [], color=C_GREY, lw=1.0, ls=(0, (4, 2)))],
               [f"ANN (z-score), AUC 0.820", f"Logistic (z-score), AUC 0.841",
                f"ANN (raw), AUC 0.772"],
               loc="lower right", fontsize=6.0, handlelength=1.5,
               borderaxespad=0.35, labelspacing=0.28)
    axC.text(0.5, -0.235, "External: GSE110224 (n = 34)", transform=axC.transAxes,
             fontsize=6.6, ha="center", va="top", color="#444444")

    letter(fig, axA, "A", dx=-0.048, dy=0.055)
    letter(fig, axB, "B", dx=-0.048, dy=0.055)
    letter(fig, axC, "C", dx=-0.048, dy=0.055)
    return fig


# ================================================================ 分面板
def single(fn, size=(2.55, 2.45), stem="x"):
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    fn(fig, ax)
    return fig


def f4a(fig, ax):
    draw_4A(ax)
    ax.text(-0.15, 1.16, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")


def f4b(fig, ax):
    draw_4B(ax)


def f4c(fig, ax=None):
    ax1 = fig.add_subplot(121); ax2 = fig.add_subplot(122)
    draw_4C(ax1, ax2)


def f5(fig, ax, csvname, color, label, note):
    draw_roc(ax, [(read_csv(csvname), color, "-", "")])
    ax.text(0.96, 0.18, label, fontsize=6.4, ha="right", va="bottom",
            color=color, linespacing=1.35)
    ax.text(0.5, -0.24, note, transform=ax.transAxes, fontsize=6.6,
            ha="center", va="top", color="#444444")


if __name__ == "__main__":
    fig = fig4_composite(); save(fig, "Figure4_v2_A-B-C"); plt.close(fig)
    fig = fig5_composite(); save(fig, "Figure5_v2_A-B-C"); plt.close(fig)

    fig = single(f4a, (2.55, 2.55)); save(fig, "Figure4A_LASSO_CV"); plt.close(fig)
    fig = single(f4b, (2.55, 2.75)); save(fig, "Figure4B_Venn"); plt.close(fig)
    fig = plt.figure(figsize=(6.6, 2.35)); f4c(fig); save(fig, "Figure4C_Selection_Frequency"); plt.close(fig)
    fig = single(lambda f, a: f5(f, a, "fig5A_roc.csv", C_RED, "AUC 0.987\n95% CI 0.965–1.000",
                                 "In-sample fit (descriptive)"), (2.50, 2.80))
    save(fig, "Figure5A_ANN_Train_ROC"); plt.close(fig)
    fig = single(lambda f, a: f5(f, a, "fig5B_roc.csv", C_BLUE, "AUC 0.907\nAccuracy 0.900",
                                 "10-fold cross-validation"), (2.50, 2.80))
    save(fig, "Figure5B_ANN_CV_ROC"); plt.close(fig)
    print("ALL DONE")
