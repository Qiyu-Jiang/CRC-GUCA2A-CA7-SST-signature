# Supplementary tables: provenance notes

- `Table_S23_nested_reselection_cv_heldout_metrics.csv`,
  `Table_S24_nested_reselection_cv_features_per_fold.csv` and
  `Table_S27_nested_cv_predictions.csv` are the **historical record** of the original
  single-seed nested-CV run. They are kept unchanged and are superseded by
  `Supplementary_Code/nested_stable.R` (see `Supplementary_Code/ANALYSIS_CORRECTION.md`).
- `Table_S28_fixed_gene_repeated_cv_summary.csv` and
  `Table_S29_fixed_gene_repeated_cv_per_partition.csv` are unaffected by the correction.
- `Table_S33 = Table_S9` and `Table_S33b = Table_S11` (identical files, kept for number continuity).

## TCGA cross-platform validation (2026-11)

- `Table_S42_clinical_associations.csv`: clinical associations of GUCA2A/CA7/SST in TCGA-COADREAD, **tumor samples only** (n = 380; 15 tests, Benjamini-Hochberg q values). Histologic grade was not available in the Xena clinical matrix.
- `Table_S43_TCGA_crossplatform_AUC.csv`: AUC registry for the frozen three-gene models in TCGA-COADREAD (380 tumors + 51 adjacent normals = 431 samples; DeLong 95% CI).
- `Supplementary_Figures/Figure_S15.pdf|png`: four-panel ROC/AUC summary for the same cohort.
