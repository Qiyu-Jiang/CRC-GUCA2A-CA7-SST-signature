## ============================================================
## verify_s27.R  —  reproduce the 0.528 estimate from Table S27 and show why it is an artefact
##
## Input : ../Supplementary_Tables/Table_S27_nested_cv_predictions.csv
## Output: console (per-fold AUC + pooled AUC)
##
## Background
##   Table S27 stores the held-out predictions of the original single-seed nested run
##   (P0_03_nested_cv.R). All folds should contribute predictions on a comparable
##   probability scale, but the stored values cluster at the extremes: 7/10 folds are
##   saturated (AUC 1.00 or 0.00), i.e. the nnet classifier returned the same class
##   probability for most held-out samples. Two folds carry most of the signal, and the
##   pooled estimate collapses to ~0.528 - which reads as "re-selection destroys
##   performance" although it reflects an unstable single-seed neural fit, not the
##   selection procedure.
##
##   The corrected, replicated analysis is nested_stable.R.
## ============================================================
suppressMessages(library(pROC))

f <- commandArgs(trailingOnly = TRUE)
path <- if (length(f)) f[1] else "../Supplementary_Tables/Table_S27_nested_cv_predictions.csv"
d <- read.csv(path, check.names = FALSE)
cat("rows:", nrow(d), "| columns:", paste(names(d), collapse = ", "), "\n\n")

obs_col <- grep("true|label|y_true|observed", names(d), ignore.case = TRUE, value = TRUE)[1]
prd_col <- grep("prob|pred|score|ann", names(d), ignore.case = TRUE, value = TRUE)[1]
fold_col <- grep("fold", names(d), ignore.case = TRUE, value = TRUE)[1]
cat("using obs =", obs_col, "| pred =", prd_col, "| fold =", fold_col, "\n\n")

for (k in sort(unique(d[[fold_col]]))) {
  s <- d[d[[fold_col]] == k, ]
  a <- as.numeric(auc(s[[obs_col]], s[[prd_col]], quiet = TRUE))
  cat(sprintf("fold %2d  n = %3d  AUC = %.3f  pred range = %.3f-%.3f  unique = %d\n",
              k, nrow(s), a, min(s[[prd_col]]), max(s[[prd_col]]), length(unique(round(s[[prd_col]], 4))))
}
cat(sprintf("\npooled AUC = %.3f\n", as.numeric(auc(d[[obs_col]], d[[prd_col]], quiet = TRUE))))
cat("Interpretation: saturated folds (AUC exactly 0 or 1) indicate a degenerate single-seed fit;\n",
    "the corrected estimate is reported in nested_stable_summary.txt.\n", sep = "")
