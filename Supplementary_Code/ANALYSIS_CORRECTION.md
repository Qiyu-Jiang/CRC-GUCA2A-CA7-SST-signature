# Correction notice: nested re-selection cross-validation (2026-11)

## What was wrong

`P0_03_nested_cv.R` ran the nested re-selection analysis once (one seed) and fitted an
`nnet` classifier on the training fold. Two defects made the held-out predictions
unusable:

1. **No re-standardisation of the held-out matrix.** The training matrix was z-scored
   once, before the fold split, and the same scaling was applied to held-out samples
   whose gene set differed from fold to fold.
2. **Degenerate single-seed fit.** `nnet` with default settings converged to a
   near-constant output in 7 of 10 folds: the stored probabilities in
   `Table_S27_nested_cv_predictions.csv` take few distinct values per fold, and the
   per-fold AUC is exactly 0 or 1 in those folds.

Pooling those predictions gives **AUC 0.528**, which reads as "re-selection destroys
the signature" but is an artefact of the fit, not a property of the selection
procedure. `verify_s27.R` reproduces the 0.528 fold by fold from the stored table.

## What replaces it

`nested_stable.R` re-estimates the whole chain inside every training fold:

- limma tumour-vs-normal (|log2FC| > 1, BH P < 0.05)
- LASSO (`cv.glmnet`, `lambda.min`)
- z-scoring fitted on the training fold only
- logistic regression fitted on the training fold

Design: 5 partitions x 10 folds = **50 paired folds**, seed 20261114. Both the fixed
three-gene model and the re-selected model are scored on the same held-out fold, so
their difference is paired.

## Results (50 folds)

| Quantity | Value |
|---|---|
| Fixed three-gene AUC | 0.980 (SD 0.040) |
| Nested re-selection AUC | 0.999 (SD 0.009) |
| Paired difference, Delta = fixed - nested | **-0.019 (95% CI -0.030 to -0.009)** |
| GUCA2A selected | 50/50 folds |
| CA7 selected | 0/50 folds |
| SST selected | 0/50 folds |

A negative Delta means the nested estimate is marginally **higher**: re-selection does
not degrade performance. The manuscript reports this as
"paired DeltaAUC = -0.02; the nested estimate was marginally higher, a difference
without practical significance" (Table 5, Methods 2.7).

## Which files are authoritative

| File | Status |
|---|---|
| `nested_stable.R`, `nested_stable_50folds.csv`, `nested_stable_summary.txt` | **authoritative** |
| `P0_03_nested_cv.R` | historical; header flagged as superseded |
| `Table_S23`, `Table_S24`, `Table_S27` | historical record, retained unchanged |
| `Table_S28`, `Table_S29` | fixed-gene repeated CV, unaffected |

## Reproducing

```bash
# 1. place the three GEO series matrices in GEO/
#    GSE22598, GSE33113, GSE41328 (.txt.gz)
Rscript Supplementary_Code/build_matrix.R      # -> GEO/train_matrix.rds
Rscript Supplementary_Code/nested_stable.R     # -> nested_stable_50folds.csv
Rscript Supplementary_Code/verify_s27.R        # reproduces the 0.528 artefact
```

Required R packages: limma, glmnet, pROC, hgu133plus2.db (Bioconductor).
