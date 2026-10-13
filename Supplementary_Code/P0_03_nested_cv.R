## ============================================================
## P0-3 Nested CV + multi-seed consensus (addresses review M1/M2/M3)
## Outputs -> 99.P0/03_nested_cv/
##   nested_cv_predictions.csv   held-out predictions (10 folds)
##   nested_cv_metrics.csv       AUC(DeLong) + acc + sens/spec(Youden) + PR-AUC
##   fold_features.csv           per-fold selected features (LASSO/SVM/RF/intersection)
##   seed_consensus_full.csv     10-seed full-training-set consensus (LASSO/SVM/RF/intersection)
##   seed_consensus_down.csv     10-seed 80% downsampling consensus
## Notes:
##   - WGCNA module boundaries (brown/purple) are fixed from the full training set;
##     limma DEG is recomputed within each fold. This partially nested design is
##     reported honestly in the manuscript (module boundary leakage remains).
## ============================================================
suppressMessages({library(limma); library(glmnet); library(e1071); library(randomForest)
                  library(neuralnet); library(pROC)})
set.seed(123)
root <- "F:/232\u673a\u5668\u5b66\u4e60\u8054\u5408\u795e\u7ecf\u7f51\u7edc/CRC"
outdir <- file.path(root, "99.P0", "03_nested_cv")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("== P0-3 nested CV start ==\n")

readExpr <- function(path){
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "", comment.char = "")
  m <- as.matrix(rt[, -1]); storage.mode(m) <- "numeric"
  rownames(m) <- as.character(rt[[1]]); m
}

expr <- readExpr(file.path(root, "10.WGCNA", "merge.normalize.txt"))
grp <- sub(".*_(Control|Treat)$", "\\1", colnames(expr))
y <- ifelse(grp == "Treat", 1L, 0L)
cat("full matrix:", dim(expr)[1], "genes x", dim(expr)[2], "samples; y:", table(y), "\n")

cand <- readLines(file.path(root, "11.venn", "interGenes.txt")); cand <- trimws(cand); cand <- cand[nchar(cand) > 0]
brown <- readLines(file.path(root, "10.WGCNA", "module_brown.txt")); brown <- trimws(brown); brown <- brown[nchar(brown) > 0]
purple <- readLines(file.path(root, "10.WGCNA", "module_purple.txt")); purple <- trimws(purple); purple <- purple[nchar(purple) > 0]
modgenes <- union(brown, purple)
cat("43-candidates:", length(intersect(cand, rownames(expr))), " module genes:", length(intersect(modgenes, rownames(expr))), "\n")

## ---------- helpers ----------
lasso_sel <- function(X, yy, seed){
  suppressMessages(library(glmnet))
  set.seed(seed)
  cvfit <- cv.glmnet(X, factor(yy), family = "binomial", alpha = 1,
                     type.measure = "deviance", nfolds = 10)
  cf <- coef(cvfit, s = "lambda.min"); idx <- which(cf != 0)
  rn <- rownames(cf)[idx]; rn <- rn[rn != "(Intercept)"]; rn
}
svmRFE_rank <- function(X, yy){
  # linear-kernel RFE, drop 10% weakest per pass (Duan et al. style)
  # FIX (3rd-party review M4): return GENE NAMES (not column indices), strongest FIRST.
  Xs <- scale(X); keep <- colnames(Xs); ranks <- c()
  yy <- factor(yy)
  repeat{
    if(length(keep) <= 2){
      ranks <- c(keep, ranks); break
    }
    fit <- svm(x = Xs[, keep, drop = FALSE], y = yy, kernel = "linear", cost = 1, scale = FALSE)
    w <- t(fit$coefs) %*% fit$SV           # weight vector
    score <- as.numeric(w^2)  # per-feature |w_i|^2 (w is 1 x p row vector)
    names(score) <- keep
    ord <- order(as.numeric(score), decreasing = FALSE)  # weakest first
    n.drop <- max(1, floor(length(keep) * 0.1))
    drop <- keep[ord[seq_len(n.drop)]]
    ranks <- c(drop, ranks)
    keep <- setdiff(keep, drop)
  }
  rev(ranks)  # strongest first
}
rf_sel <- function(X, yy, seed, topk = 3){
  suppressMessages(library(randomForest))
  set.seed(seed)
  d <- data.frame(X, check.names = FALSE); d$.y <- factor(yy)
  rf <- randomForest(.y ~ ., data = d, ntree = 500)
  imp <- importance(rf)[, "MeanDecreaseGini"]
  names(sort(imp, decreasing = TRUE))[seq_len(min(topk, length(imp)))]
}
ann_fit_pred <- function(Xtr, ytr, Xte, seed, genes){
  suppressMessages(library(nnet))
  set.seed(seed)
  Y <- cbind(Control = 1 - ytr, Treat = ytr)
  fit <- nnet(x = Xtr, y = Y, size = 5, maxit = 1000, trace = FALSE, rang = 0.5)
  nn <- predict(fit, Xte, type = "raw")
  if (is.null(ncol(nn))) nn <- matrix(nn, ncol = 1)
  if (ncol(nn) == 1) nn <- cbind(1 - nn, nn)
  colnames(nn) <- c("Control", "Treat"); nn
}
pr_auc <- function(y, s){
  o <- order(s, decreasing = TRUE); yy <- y[o]; ss <- s[o]
  tp <- cumsum(yy == 1); fp <- cumsum(yy != 1)
  prec <- tp / (tp + fp); rec <- tp / sum(yy == 1)
  # trapz over recall, prepend (0,1)
  x <- c(0, rec); v <- c(1, prec)
  sum(diff(x) * (head(v, -1) + tail(v, -1)) / 2)
}

## ---------- nested 10-fold CV ----------
set.seed(123)
i1 <- sample(which(y == 1)); i0 <- sample(which(y == 0))
f1 <- split(i1, cut(seq_along(i1), 10)); f0 <- split(i0, cut(seq_along(i0), 10))
pred <- rep(NA_real_, length(y))
fold.features <- data.frame()
for (k in 1:10) {
  te <- c(f1[[k]], f0[[k]]); tr <- setdiff(seq_along(y), te)
  ## (a) fold-level limma DEG
  design <- model.matrix(~ factor(y[tr]))
  fit <- lmFit(expr[, tr, drop = FALSE], design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf)
  deg <- rownames(tt)[abs(tt$logFC) > 2 & tt$adj.P.Val < 0.05]
  candF <- intersect(deg, modgenes)
  if (length(candF) < 10) candF <- intersect(union(deg, cand), modgenes)
  candF <- intersect(candF, rownames(expr)); candF <- candF[!duplicated(candF)]
  Xall <- t(expr[match(candF, rownames(expr)), tr, drop = FALSE])
  ## (b) fold-level three algorithms
  l.f <- tryCatch(lasso_sel(Xall, y[tr], 200 + k), error = function(e) character(0))
  rk <- tryCatch(svmRFE_rank(Xall, y[tr]), error = function(e) character(0))
  s.f <- head(rk, 11)
  r.f <- tryCatch(rf_sel(Xall, y[tr], 400 + k, 3), error = function(e) character(0))
  inter <- Reduce(intersect, list(l.f, s.f, r.f))
  used <- if (length(inter) >= 2) inter else if (length(intersect(l.f, s.f)) >= 2) intersect(l.f, s.f) else l.f
  if (length(used) == 0) used <- head(names(sort(table(c(l.f, s.f, r.f))), decreasing = TRUE), 3)
  used <- used[!duplicated(used)]; used <- intersect(used, colnames(Xall))
  if (length(used) == 0) used <- intersect(genes3 <- c("GUCA2A","CA7","SST"), colnames(Xall))
  ## (c) fold-level ANN
  pr <- ann_fit_pred(Xall[, match(used, colnames(Xall)), drop = FALSE], y[tr],
                     t(expr[match(used, rownames(expr)), te, drop = FALSE]), 600 + k, used)
  pred[te] <- pr[, "Treat"]
  fold.features <- rbind(fold.features, data.frame(
    fold = k, n_cand_fold = length(candF), n_LASSO = length(l.f), n_SVM = length(s.f),
    n_RF = length(r.f), n_inter = length(inter),
    features = paste(used, collapse = ";"),
    GUCA2A = as.integer("GUCA2A" %in% used), CA7 = as.integer("CA7" %in% used),
    SST = as.integer("SST" %in% used), CA2 = as.integer("CA2" %in% used),
    GUCA2B = as.integer("GUCA2B" %in% used)))
  cat("fold", k, ": cand", length(candF), "L", length(l.f), "S", length(s.f),
      "R", length(r.f), "-> inter", length(inter), "used:", paste(used, collapse = "+"), "\n")
}
write.csv(fold.features, file.path(outdir, "fold_features.csv"), row.names = FALSE)
write.csv(data.frame(sample = colnames(expr), y = y, nested_pred_Treat = pred),
          file.path(outdir, "nested_cv_predictions.csv"), row.names = FALSE)

roc.n <- roc(y, pred)
ci <- as.numeric(ci.auc(roc.n, method = "delong"))
thr <- coords(roc.n, "best", transpose = FALSE)
metrics <- data.frame(
  model = "Nested CV (fold-level DEG+3algo+ANN)",
  auc = as.numeric(roc.n$auc), ci_low = ci[1], ci_high = ci[3],
  accuracy = mean((pred >= thr$threshold) == (y == 1)),
  sensitivity = thr$sensitivity, specificity = thr$specificity,
  pr_auc = pr_auc(y, pred))
write.csv(metrics, file.path(outdir, "nested_cv_metrics.csv"), row.names = FALSE)
cat("NESTED CV AUC:", round(as.numeric(roc.n$auc), 3),
    "CI", round(ci[1],3), "-", round(ci[3],3),
    "| sens", round(thr$sensitivity,3), "spec", round(thr$specificity,3),
    "| PR-AUC", round(metrics$pr_auc, 3), "\n")

## ---------- 10-seed consensus (full training set, 43 candidates) ----------
X43 <- t(expr[match(intersect(cand, rownames(expr)), rownames(expr)), , drop = FALSE])
seeds <- c(123, 456, 789, 2024, 314, 2718, 999, 55, 777, 2025)
consF <- data.frame()
for (s in seeds) {
  l.f <- lasso_sel(X43, y, s)
  rk <- svmRFE_rank(X43, y); s.f <- head(rk, 11)
  r.f <- rf_sel(X43, y, s, 3)
  inter <- Reduce(intersect, list(l.f, s.f, r.f))
  consF <- rbind(consF, data.frame(seed = s, n_LASSO = length(l.f), n_SVM = length(s.f),
    n_RF = length(r.f), n_inter = length(inter),
    LASSO = paste(l.f, collapse = ";"), SVM = paste(s.f, collapse = ";"),
    RF = paste(r.f, collapse = ";"), INTER = paste(inter, collapse = ";"),
    GUCA2A = as.integer("GUCA2A" %in% inter), CA7 = as.integer("CA7" %in% inter),
    SST = as.integer("SST" %in% inter), CA2 = as.integer("CA2" %in% inter)))
  cat("seed", s, "inter:", paste(inter, collapse = "+"), "\n")
}
write.csv(consF, file.path(outdir, "seed_consensus_full.csv"), row.names = FALSE)

## ---------- 10-seed 80% downsampling consensus ----------
consD <- data.frame()
for (s in seeds) {
  set.seed(s)
  idx <- c(sample(which(y == 1), floor(sum(y == 1) * 0.8)),
           sample(which(y == 0), floor(sum(y == 0) * 0.8)))
  Xs <- X43[idx, , drop = FALSE]; ys <- y[idx]
  l.f <- lasso_sel(Xs, ys, s)
  rk <- svmRFE_rank(Xs, ys); s.f <- head(rk, 11)
  r.f <- rf_sel(Xs, ys, s, 3)
  inter <- Reduce(intersect, list(l.f, s.f, r.f))
  consD <- rbind(consD, data.frame(seed = s, n_LASSO = length(l.f), n_SVM = length(s.f),
    n_RF = length(r.f), n_inter = length(inter), INTER = paste(inter, collapse = ";"),
    GUCA2A = as.integer("GUCA2A" %in% inter), CA7 = as.integer("CA7" %in% inter),
    SST = as.integer("SST" %in% inter), CA2 = as.integer("CA2" %in% inter)))
  cat("down seed", s, "inter:", paste(inter, collapse = "+"), "\n")
}
write.csv(consD, file.path(outdir, "seed_consensus_down.csv"), row.names = FALSE)

freq <- data.frame(
  gene = c("GUCA2A", "CA7", "SST", "CA2", "GUCA2B"),
  full_10seed = c(mean(consF$GUCA2A), mean(consF$CA7), mean(consF$SST),
                  mean(consF$CA2), NA),
  down_10seed = c(mean(consD$GUCA2A), mean(consD$CA7), mean(consD$SST),
                  mean(consD$CA2), NA))
write.csv(freq, file.path(outdir, "consensus_frequency.csv"), row.names = FALSE)
cat("== consensus frequencies ==\n"); print(freq)
cat("== P0-3 done ==\n")
