## ============================================================
## P0-4 Fixed-gene nested CV (addresses review M2 properly)
## Features are FIXED (GUCA2A/CA7/SST) - no selection inside folds, so this is
## a standard (leak-free) CV of the 3-gene model. 5 fold-partition seeds x 10 folds.
## Fold-level z-scoring: parameters computed within each training fold.
## ANN identical architecture to the original (3-5-2 -> hidden=5).
## Outputs -> outputs/phase0/04_fixed_nested_cv/
##   fixed_nested_by_seed.csv   per-seed pooled AUC/acc/sens/spec/PR-AUC
##   fixed_nested_summary.csv   mean +/- SD across 5 seeds
##   fixed_predictions_seed123.csv  sample-level predictions (seed 123)
## ============================================================
suppressMessages({library(neuralnet); library(pROC)})
root <- "./"
outdir <- file.path(root, "outputs/phase0", "04_fixed_nested_cv")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("== P0-4 fixed-gene nested CV start ==\n")

readExpr <- function(path){
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "", comment.char = "")
  m <- as.matrix(rt[, -1]); storage.mode(m) <- "numeric"
  rownames(m) <- as.character(rt[[1]]); m
}
genes3 <- c("GUCA2A", "CA7", "SST")
trainRaw <- t(readExpr(file.path(root, "06.sva", "merge.normalize.txt"))[genes3, ])
grp <- sub(".*_(Control|Treat)$", "\\1", rownames(trainRaw))
y <- ifelse(grp == "Treat", 1L, 0L)
cat("train:", nrow(trainRaw), "samples; y:", table(y), "\n")

pr_auc <- function(yy, s){
  o <- order(s, decreasing = TRUE); yy <- yy[o]
  tp <- cumsum(yy == 1); fp <- cumsum(yy != 1)
  prec <- tp / (tp + fp); rec <- tp / sum(yy == 1)
  x <- c(0, rec); v <- c(1, prec)
  sum(diff(x) * (head(v, -1) + tail(v, -1)) / 2)
}
one_cv <- function(seed, model = c("ann", "glm")){
  model <- match.arg(model)
  set.seed(seed)
  i1 <- sample(which(y == 1)); i0 <- sample(which(y == 0))
  f1 <- split(i1, cut(seq_along(i1), 10)); f0 <- split(i0, cut(seq_along(i0), 10))
  pred <- rep(NA_real_, length(y))
  for (k in 1:10) {
    te <- c(f1[[k]], f0[[k]]); tr <- setdiff(seq_along(y), te)
    mu <- colMeans(trainRaw[tr, , drop = FALSE])
    sdv <- apply(trainRaw[tr, , drop = FALSE], 2, sd)
    Ztr <- sweep(sweep(trainRaw[tr, , drop = FALSE], 2, mu, "-"), 2, sdv, "/")
    Zte <- sweep(sweep(trainRaw[te, , drop = FALSE], 2, mu, "-"), 2, sdv, "/")
    if (model == "ann") {
      set.seed(600 + seed %% 1000 * 10 + k)
      dd <- cbind(as.data.frame(Ztr, check.names = FALSE),
                  Control = 1 - y[tr], Treat = y[tr])
      fit <- neuralnet(Control + Treat ~ GUCA2A + CA7 + SST, dd,
                       hidden = 5, stepmax = 1e6)
      nn <- compute(fit, as.data.frame(Zte, check.names = FALSE))$net.result
      if (is.null(dim(nn))) nn <- matrix(nn, ncol = 1)
      if (ncol(nn) == 1) nn <- cbind(1 - nn, nn)
      colnames(nn) <- c("Control", "Treat")
      pred[te] <- nn[, "Treat"]
    } else {
      df.tr <- data.frame(Ztr, check.names = FALSE); df.tr$y <- factor(y[tr])
      fit <- glm(y ~ ., family = binomial, data = df.tr)
      pred[te] <- as.numeric(predict(fit, newdata = data.frame(Zte, check.names = FALSE),
                                     type = "response"))
    }
  }
  roc.o <- roc(y, pred); ci <- as.numeric(ci.auc(roc.o, method = "delong"))
  thr <- coords(roc.o, "best", transpose = FALSE)
  assign(if (model == "ann") "one_cv_pred" else "one_cv_pred_glm", pred, envir = parent.frame())
  data.frame(seed = seed, auc = as.numeric(roc.o$auc), ci_low = ci[1], ci_high = ci[3],
             accuracy = mean((pred >= thr$threshold) == (y == 1)),
             sensitivity = thr$sensitivity, specificity = thr$specificity,
             pr_auc = pr_auc(y, pred))
}
res_ann <- data.frame(); res_glm <- data.frame()
pred_cols <- data.frame(row.names = rownames(trainRaw))
pred_cols$y <- y
for (s in c(123, 456, 789, 2024, 314)) {
  r1 <- one_cv(s, "ann");  res_ann <- rbind(res_ann, r1)
  r2 <- one_cv(s, "glm");  res_glm <- rbind(res_glm, r2)
  if (s == 123) pred_cols$ANN_seed123 <- one_cv_pred
  if (s == 123) pred_cols$GLM_seed123 <- one_cv_pred_glm
  cat("seed", s, "ANN AUC", round(r1$auc, 3), "| GLM AUC", round(r2$auc, 3), "\n")
}
write.csv(pred_cols, file.path(outdir, "fixed_predictions_seed123.csv"), row.names = TRUE)
write.csv(res_ann, file.path(outdir, "fixed_nested_by_seed.csv"), row.names = FALSE)
write.csv(res_glm, file.path(outdir, "fixed_nested_glm_by_seed.csv"), row.names = FALSE)

summ <- data.frame(
  model = c("ANN fixed-gene nested CV (5 partitions)", "Logistic fixed-gene nested CV (5 partitions)"),
  auc_mean = c(mean(res_ann$auc), mean(res_glm$auc)),
  auc_sd = c(sd(res_ann$auc), sd(res_glm$auc)),
  auc_min = c(min(res_ann$auc), min(res_glm$auc)),
  auc_max = c(max(res_ann$auc), max(res_glm$auc)),
  sens_mean = c(mean(res_ann$sensitivity), mean(res_glm$sensitivity)),
  spec_mean = c(mean(res_ann$specificity), mean(res_glm$specificity)),
  pr_auc_mean = c(mean(res_ann$pr_auc), mean(res_glm$pr_auc)))
write.csv(summ, file.path(outdir, "fixed_nested_summary.csv"), row.names = FALSE)
cat("== fixed-gene nested summary ==\n"); print(summ, digits = 3)
cat("== P0-4 done ==\n")
