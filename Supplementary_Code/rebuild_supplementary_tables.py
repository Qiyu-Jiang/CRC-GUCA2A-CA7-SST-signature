#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rebuild the journal-format supplementary document from the authoritative CSVs.

The submitted Supplementary_Materials_COMPLETE.docx originally carried abbreviated
copies of nine per-sample tables ("Showing first N of M rows/columns") plus a
mis-ordered S25 block and a caption/table mismatch around S6/S6b. This script
replaces those tables with the full data held in Supplementary_Tables/, so that
the docx and the CSV files cannot drift apart.

Usage
-----
    python rebuild_supplementary_tables.py --docx path/to/Supplementary_Materials_COMPLETE.docx \
                                           --csv-dir Supplementary_Tables

Requirements: python-docx (>= 1.1).

Notes
-----
* S9 and S11 are transposed (samples as rows) so that 150 sample columns do not
  have to be squeezed into a portrait Word page; the untransposed matrices remain
  available as the original CSV files.
* S38 is appended if it is missing (HPA protein-level evidence).
* Run once per document; the script is not idempotent for S38 (it will add a
  second copy if one already exists).
"""
import argparse
import csv
import os
import re
import sys

import docx
from docx.shared import Pt


# caption prefix -> (csv file, font size, transpose)
MAPPING = {
    "Table S3.": ("Table_S3_downsampling_stability_feature_selection.csv", 7, False),
    "Table S5.": ("Table_S5_seeded_feature_lists_full_training.csv", 7, False),
    "Table S6.": ("Table_S6_ANN_performance_summary.csv", 9, False),
    "Table S6b.": ("Table_S6b_LOCO_cohort_blocked_predictions.csv", 7, False),
    "Table S7.": ("Table_S7_ANN_external_sample_level_predictions.csv", 8, False),
    "Table S9.": ("Table_S9_ssGSEA_scores_Senbabaoglu20.csv", 6, True),
    "Table S11.": ("Table_S11_ssGSEA_scores_Charoentong28.csv", 6, True),
    "Table S12.": ("Table_S12_ssGSEA_correlations_Charoentong28.csv", 7, False),
    "Table S18.": ("Table_S18_UALCAN_methylation_boxstats.csv", 8, False),
    "Table S25.": ("Table_S25_ten_seed_consensus_full.csv", 6, False),
    "Table S26.": ("Table_S26_ten_seed_consensus_downsample.csv", 7, False),
    "Table S27.": ("Table_S27_nested_cv_predictions.csv", 7, False),
    "Table S32.": ("Table_S32_CIBERSORT_per_sample_scores.csv", 6, False),
}


def read_csv(path):
    raw = open(path, encoding="utf-8-sig", errors="ignore").read().splitlines()
    delim = "\t" if ("\t" in raw[0] and "," not in raw[0]) else ","
    return list(csv.reader(raw, delimiter=delim))


def transpose(rows):
    out = [["Sample"] + [r[0] for r in rows[1:]]]
    for j in range(1, len(rows[0])):
        out.append([rows[0][j]] + [r[j] if j < len(r) else "" for r in rows[1:]])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docx", required=True)
    ap.add_argument("--csv-dir", default="Supplementary_Tables")
    ap.add_argument("--no-s38", action="store_true", help="skip adding Table S38")
    args = ap.parse_args()

    doc = docx.Document(args.docx)
    body = doc.element.body

    def find_caption(prefix):
        for p in doc.paragraphs:
            if p.text.strip().startswith(prefix):
                return p
        return None

    def table_for(cap_par):
        """Next table before the following 'Table S' caption (None if absent)."""
        nxt = cap_par._p.getnext()
        while nxt is not None:
            if nxt.tag.endswith("}tbl"):
                return nxt
            if nxt.tag.endswith("}p") and "".join(nxt.itertext()).strip().startswith("Table S"):
                return None
            nxt = nxt.getnext()
        return None

    def build(rows, pt):
        ncol = max(len(r) for r in rows)
        tb = doc.add_table(rows=len(rows), cols=ncol)
        try:
            tb.style = "Table Grid"
        except Exception:
            pass
        for i, r in enumerate(rows):
            for j in range(ncol):
                cell = tb.rows[i].cells[j]
                cell.text = r[j] if j < len(r) else ""
                for par in cell.paragraphs:
                    for run in par.runs:
                        run.font.size = Pt(pt)
                        run.font.bold = (i == 0)
        return tb

    for cap_prefix, (fname, pt, trans) in MAPPING.items():
        par = find_caption(cap_prefix)
        if par is None:
            print("caption not found:", cap_prefix, file=sys.stderr)
            continue
        rows = read_csv(os.path.join(args.csv_dir, fname))
        if trans:
            rows = transpose(rows)
        old = table_for(par)
        new = build(rows, pt)
        par._p.addnext(new._tbl)
        if old is not None:
            body.remove(old)
        if trans and "samples in rows" not in par.text:
            par.runs[0].text = (par.text.rstrip() + " Samples are given in rows and cell types in "
                                "columns; the untransposed matrix is provided as " + fname + ".")
        print(f"{cap_prefix:12s} {fname:55s} {len(rows)}x{max(len(r) for r in rows)}")

    # strap table that used to sit after S29 (old S25 data, no caption)
    s29 = find_caption("Table S29.")
    if s29 is not None:
        t29 = table_for(s29)
        nxt = t29.getnext() if t29 is not None else None
        while nxt is not None and not nxt.tag.endswith("}tbl"):
            nxt = nxt.getnext()
        if nxt is not None:
            body.remove(nxt)
            print("removed stray table after S29")

    # stale "Showing first N of M" notes
    removed = 0
    for p in list(doc.paragraphs):
        if p.text.strip().startswith(("Showing first", "Abbreviated view")):
            p._p.getparent().remove(p._p)
            removed += 1
    print("removed abbreviation notes:", removed)

    if not args.no_s38 and find_caption("Table S38.") is None:
        cap = doc.add_paragraph()
        run = cap.add_run("Table S38. Protein-level evidence for GUCA2A, CA7 and SST from the Human "
                          "Protein Atlas (antibody identifiers, reliability scores, annotated normal-colon "
                          "cell types and colorectal cancer IHC counts).")
        run.bold = True
        run.font.size = Pt(9)
        rows = read_csv(os.path.join(args.csv_dir, "Table_S38_HPA_protein_level_evidence.csv"))
        build(rows, 7)
        print("added S38")

    doc.save(args.docx)
    print("saved", args.docx)


if __name__ == "__main__":
    main()
