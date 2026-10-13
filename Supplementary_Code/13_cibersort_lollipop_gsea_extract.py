# -*- coding: utf-8 -*-
"""Extract CIBERSORT group diffs, lollipop correlations, per-gene GSEA -> JSON for manuscript."""
import csv, json, os, re
import statistics

BASE = r"F:\232机器学习联合神经网络\CRC"
out = {}

# ---------- 1) CIBERSORT fractions: group means + Wilcoxon ----------
cs_path = os.path.join(BASE, "28.barplot", "CIBERSORT-Results.txt")
lines = open(cs_path, encoding="utf-8", errors="replace").read().splitlines()
data_lines = [l for l in lines if l.strip() and not l.startswith("//")]
rdr = list(csv.reader(data_lines, delimiter="\t"))
header = rdr[0]
cols = header[1:]
samples = [r0[0] for r0 in rdr[1:]]
vals = {c: [] for c in cols}
grp = {"Control": {c: [] for c in cols}, "Treat": {c: [] for c in cols}}
for r0 in rdr[1:]:
    g = "Control" if r0[0].endswith("Control") else ("Treat" if r0[0].endswith("Treat") else None)
    for j, c in enumerate(cols):
        try:
            v = float(r0[j + 1])
        except Exception:
            continue
        vals[c].append(v)
        if g: grp[g][c].append(v)
# stats for cell-type columns only (exclude P-value / correlation columns)
cell_cols = [c for c in cols if c not in ("P-value", "Pearson.CC", "P.value", "Pearson.CC.")]
try:
    from scipy.stats import mannwhitneyu
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False
cs_out = []
for c in cell_cols:
    a, b = grp["Control"][c], grp["Treat"][c]
    if not a or not b: continue
    if HAVE_SCIPY:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
    else:
        p = float("nan")
    cs_out.append({"cell": c,
                   "mean_control": round(sum(a)/len(a), 4),
                   "mean_tumor": round(sum(b)/len(b), 4),
                   "p": p})
if HAVE_SCIPY:
    # BH
    m = len(cs_out)
    ranked = sorted(cs_out, key=lambda x: x["p"])
    prev = 1.0
    for i in range(m - 1, -1, -1):
        ranked[i]["q"] = min(prev, ranked[i]["p"] * m / (i + 1))
        prev = ranked[i]["q"]
    for x in cs_out: x["q"] = round(x["q"], 4)
    for x in cs_out: x["p"] = round(x["p"], 5)
cs_out.sort(key=lambda x: -abs(x["mean_tumor"] - x["mean_control"]))
out["cibersort"] = cs_out
out["cibersort_n"] = {"control": len(grp["Control"][cell_cols[0]]) if cell_cols else 0,
                      "tumor": len(grp["Treat"][cell_cols[0]]) if cell_cols else 0}
print("CIBERSORT cell cols:", len(cell_cols), "| n:", out["cibersort_n"])

# ---------- 2) Lollipop: per-gene immune correlations ----------
lol = {}
for g in ("GUCA2A", "CA7", "SST"):
    pth = os.path.join(BASE, "29.Lollipop", f"corr.{g}.result.txt")
    rr = list(csv.reader(open(pth, encoding="utf-8", errors="replace"), delimiter="\t"))
    hdr = rr[0]
    ii_c = next((i for i, h in enumerate(hdr) if "cell" in h.lower() or "Immune" in h), 0)
    ii_r = next((i for i, h in enumerate(hdr) if "rho" in h.lower() or "cor" in h.lower() or "estimate" in h.lower()), None)
    ii_p = next((i for i, h in enumerate(hdr) if h.strip().lower() in ("pvalue", "p.value", "p") or "pvalue" in h.lower()), None)
    lol[g] = []
    for r0 in rr[1:]:
        try:
            rho = float(r0[ii_r]) if ii_r is not None else float("nan")
            pv = float(r0[ii_p]) if ii_p is not None else float("nan")
        except Exception:
            continue
        lol[g].append({"cell": r0[ii_c], "rho": round(rho, 3), "p": pv})
    lol[g] = sorted(lol[g], key=lambda x: -abs(x["rho"]))
out["lollipop"] = lol
print("lollipop genes:", list(lol.keys()), "| cells per gene:", len(lol.get("CA7", [])))

# ---------- 3) per-gene GSEA ----------
gsea = {}
for g in ("GUCA2A", "CA7", "SST"):
    pth = os.path.join(BASE, "26.GSEA", f"{g}.result.txt")
    rr = list(csv.reader(open(pth, encoding="utf-8", errors="replace"), delimiter="\t"))
    hdr = rr[0]
    iN, iP, iQ, iAdj = hdr.index("NES"), hdr.index("pvalue"), hdr.index("qvalue"), hdr.index("p.adjust")
    sets_ = []
    for r0 in rr[1:]:
        try:
            nes = float(r0[iN]); padj = float(r0[iAdj])
        except Exception:
            continue
        sets_.append({"set": r0[1], "NES": nes, "pvalue": r0[iP], "p.adjust": padj})
    sig = [s for s in sets_ if s["p.adjust"] < 0.05]
    pos = sorted([s for s in sig if s["NES"] > 0], key=lambda x: -x["NES"])
    neg = sorted([s for s in sig if s["NES"] < 0], key=lambda x: x["NES"])
    gsea[g] = {"n_sig": len(sig), "pos": pos[:6], "neg": neg[:6]}
    print(f"GSEA {g}: sig={len(sig)} pos={len(pos)} neg={len(neg)}")
out["gsea"] = gsea

json.dump(out, open(r"D:\Autoclaw\.openclaw\tmp\results_extra.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("saved results_extra.json")
