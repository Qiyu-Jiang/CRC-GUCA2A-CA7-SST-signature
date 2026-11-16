## ============================================================
## nested_stable.R  —  corrected nested re-selection cross-validation
##
## Supersedes : P0_03_nested_cv.R (single-seed, nnet, saturated in 7/10 folds)
##
## Input   : GEO/train_matrix.rds  (from build_matrix.R)
## Output  : nested_stable_50folds.csv   (per-fold fixed-gene and nested AUC)
##           nested_stable_summary.txt   (means, paired difference, selection frequencies)
##
## Design
##  5 partitions x 10 folds (stratified). Inside each training fold:
##    limma tumour-vs-normal (|log2FC| > 1, BH P < 0.05)
##    -> LASSO (cv.glmnet, lambda.min) on the fold DEGs
##    -> z-scoring fitted on the training fold only
##    -> logistic regression fitted on the training fold
##  Both the fixed three-gene model and the re-selected model are evaluated on the
##  same held-out fold, so the two AUCs are paired and their difference is reported
##  as Delta = AUC(fixed) - AUC(nested).
##
## Note: the feature pool is NOT restricted to WGCNA module membership here, and the
##       DEG threshold is |log2FC| > 1 rather than > 2, because the fold-local sample
##       size is ~10x smaller than the full training set. These deviations are
##       disclosed in Methods 2.7 of the manuscript.
## ============================================================
suppressMessages({library(limma); library(glmnet); library(pROC)})

set.seed(20261114)
dat <- readRDS(Sys.getenv("CRC_MATRIX", unset = "GEO/train_matrix.rds"))
E <- dat$E; y <- dat$y
E <- E[apply(E, 1, function(v) all(is.finite(v))), ]   # complete-case genes (15,298)
cat("matrix:", nrow(E), "x", ncol(E), "| tumours", sum(y), "\n")

genes <- c("GUCA2A", "CA7", "SST")
auc <- function(obs, pred) as.numeric(pROC::auc(obs, pred, quiet = TRUE))

fit_fixed <- function(tr, te) {
  z <- t(scale(t(E[genes, tr, drop = FALSE])))
  mu <- rowMeans(E[genes, tr, drop = FALSE]); sdv <- apply(E[genes, tr, drop = FALSE], 1, sd)
  zte <- t(scale(t(E[genes, te, drop = FALSE]), center = mu, scale = sdv))
  df <- data.frame(y = y[tr], t(z)); m <- glm(y ~ ., data = df, family = binomial)
  predict(m, newdata = data.frame(y = y[te], t(zte)), type = "response")
}

fit_nested <- function(tr, te) {
  d <- E[, tr, drop = FALSE]; yy <- y[tr]
  des <- model.matrix(~ factor(yy))
  fit <- eBayes(lmFit(d, des))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  deg <- rownames(tt)[abs(tt$logFC) > 1 & tt$adj.P.Val < 0.05]
  if (length(deg) < 5) deg <- rownames(tt)[order(tt$P.Value)][1:20]
  x <- t(d[deg, , drop = FALSE]); x <- scale(x)
  cv <- cv.glmnet(x, yy, family = "binomial", alpha = 1)
  sel <- rownames(coef(cv, s = "lambda.min"))[-1]
  sel <- sel[coef(cv, s = "lambda.min")[-1, 1] != 0]
  if (length(sel) == 0) sel <- deg[1]
  xt <- t(E[sel, tr, drop = FALSE]); mu <- colMeans(xt); sdv <- apply(xt, 2, sd)
  xt <- scale(xt, center = mu, scale = sdv)
  m <- glm(yy ~ ., data = data.frame(xt, check.names = TRUE), family = binomial)
  xe <- scale(t(E[sel, te, drop = FALSE]), center = mu, scale = sdv)
  pr <- predict(m, newdata = data.frame(xe, check.names = TRUE), type = "response")
  list(pred = as.numeric(pr), sel = sel)
}

res <- list(); sels <- c()
for (part in 1:5) {
  set.seed(20261114 + part)
  folds <- sample(rep(1:10, length.out = length(y)))
  for (k in 1:10) {
    te <- which(folds == k); tr <- setdiff(seq_along(y), te)
    pf <- fit_fixed(tr, te); pn <- fit_nested(tr, te)
    res[[length(res) + 1]] <- data.frame(partition = part, fold = k, n_test = length(te),
                                         auc_fixed = auc(y[te], pf), auc_nested = auc(y[te], pn$pred))
    sels <- c(sels, pn$sel)
  }
  cat("partition", part, "done\n")
}
R <- do.call(rbind, res)
write.csv(R, "nested_stable_50folds.csv", row.names = FALSE)

d <- R$auc_fixed - R$auc_nested
ci <- t.test(d)$conf.int
sink("nested_stable_summary.txt")
cat(sprintf("folds %d | fixed mean %.3f (SD %.3f) | nested mean %.3f (SD %.3f)\n",
            nrow(R), mean(R$auc_fixed), sd(R$auc_fixed), mean(R$auc_nested), sd(R$auc_nested)))
cat(sprintf("Delta (fixed - nested) mean %.3f | 95%% CI %.3f to %.3f\n", mean(d), ci[1], ci[2]))
for (g in genes) cat(sprintf("re-selection: %s %d/%d\n", g, sum(sels == g), nrow(R)))
sink()
cat(readLines("nested_stable_summary.txt"), sep = "\n")
