## ============================================================
## P0-2  敏感性分析（修正版）
##   A. 合并后 PCA：按数据集着色 + 按分组标形（preNorm vs ComBat 后）
##   B. Leave-one-dataset-out：每次用 2 个数据集重做 ComBat+limma，
##      核查 GUCA2A/CA7/SST 的 logFC/adj.P
##   C. GSE33113 肿瘤抽样降权(30例) × 3 个种子：ComBat+limma 核查三基因，
##      并在 ComBat 后矩阵上重跑 LASSO/SVM-RFE/RF 特征共识
## 输出: outputs/phase0/02.sensitivity/
## ============================================================
suppressMessages({library(limma); library(sva); library(glmnet);
                  library(randomForest); library(e1071)})
root   <- "./"
outdir <- file.path(root, "outputs/phase0", "02.sensitivity")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
genes3 <- c("GUCA2A", "CA7", "SST")
hub <- readLines(file.path(root, "16.hubGene", "hubGenes.txt"))
tp <- function(cn) sub(".*_(Control|Treat)$", "\\1", cn)
cat("== P0-2 start ==\n")

readExpr <- function(path, average = FALSE){
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "", comment.char = "")
  m <- as.matrix(rt[, -1]); storage.mode(m) <- "numeric"
  rownames(m) <- as.character(rt[[1]])
  if (average) limma::avereps(m) else m
}
combatExpr <- function(ex, batch){
  batchf <- factor(batch)
  if (nlevels(batchf) > 1) ComBat(ex, batchf, par.prior = TRUE) else ex
}
degTable <- function(ex, label){
  grp <- tp(colnames(ex))
  design <- model.matrix(~0 + factor(grp)); colnames(design) <- c("Control", "Treat")
  fit <- eBayes(contrasts.fit(lmFit(ex, design),
             makeContrasts(Treat - Control, levels = design)))
  tt <- topTable(fit, number = Inf, adjust = "fdr", sort.by = "none")
  data.frame(analysis = label, gene = genes3,
             logFC = tt[genes3, "logFC"], adjP = tt[genes3, "adj.P.Val"],
             direction = ifelse(tt[genes3, "logFC"] > 0, "up_in_tumor", "down_in_tumor"),
             nDEG = sum(abs(tt$logFC) > 2 & tt$adj.P.Val < 0.05, na.rm = TRUE),
             row.names = NULL)
}

ds22598 <- readExpr(file.path(root, "05.normalize", "GSE22598.normalize.txt"), average = TRUE)
ds33113 <- readExpr(file.path(root, "05.normalize", "GSE33113.normalize.txt"), average = TRUE)
ds41328 <- readExpr(file.path(root, "05.normalize", "GSE41328.normalize.txt"), average = TRUE)

## ---------- A. PCA ----------
pcaPanel <- function(pre, post, file){
  mk <- function(m){
    md <- apply(m, 1, mad); top <- names(sort(md, decreasing = TRUE))[1:2000]
    prcomp(t(m[top, ]))
  }
  p1 <- mk(pre); p2 <- mk(post)
  ds <- sub("^(GSE\\d+)_.*", "\\1", colnames(pre))
  gp <- tp(colnames(pre))
  pdf(file, width = 11, height = 5.2)
  par(mfrow = c(1, 2))
  for (pc in list(p1, p2)) {
    v <- summary(pc)$importance[2, 1:2] * 100
    plot(pc$x[, 1], pc$x[, 2], pch = ifelse(gp == "Treat", 17, 19), cex = 1.1,
         col = ifelse(ds == "GSE22598", "#c0392b",
                      ifelse(ds == "GSE33113", "#2563a8", "#0e7c66")),
         xlab = sprintf("PC1 (%.1f%%)", v[1]), ylab = sprintf("PC2 (%.1f%%)", v[2]))
    legend("bottomleft", pch = 19, col = c("#c0392b", "#2563a8", "#0e7c66"),
           legend = c("GSE22598", "GSE33113", "GSE41328"), bty = "n", cex = 0.9)
    legend("topright", pch = c(19, 17), col = "grey30",
           legend = c("Control", "Treat"), bty = "n", cex = 0.9)
  }
  dev.off()
}
pre  <- readExpr(file.path(root, "06.sva", "merge.preNorm.txt"))
post <- readExpr(file.path(root, "06.sva", "merge.normalize.txt"))
pcaPanel(pre, post, file.path(outdir, "PCA_by_dataset.pdf"))
cat("A. PCA done\n")

## ---------- B. LODO ----------
lodo <- rbind(
  degTable(combatExpr(cbind(ds22598, ds41328),
                      c(rep("A", ncol(ds22598)), rep("B", ncol(ds41328)))),
           "leaveOut_GSE33113"),
  degTable(combatExpr(cbind(ds33113, ds41328),
                      c(rep("A", ncol(ds33113)), rep("B", ncol(ds41328)))),
           "leaveOut_GSE22598"),
  degTable(combatExpr(cbind(ds22598, ds33113),
                      c(rep("A", ncol(ds22598)), rep("B", ncol(ds33113)))),
           "leaveOut_GSE41328"))
write.csv(lodo, file.path(outdir, "LODO_3gene.csv"), row.names = FALSE)
print(lodo); cat("B. LODO done\n")

## ---------- C. 降权抽样 + 特征重选 ----------
source(file.path(root, "18.SVM", "neuralML18.msvmRFE.R"), encoding = "UTF-8")
featSel <- function(ex, grp, seed, tag){
  hb <- intersect(hub, rownames(ex))
  X <- t(ex[hb, , drop = FALSE])
  g <- factor(grp, levels = c("Control", "Treat"))
  set.seed(seed + 1)
  cv <- cv.glmnet(X, g, family = "binomial", alpha = 1,
                  type.measure = "deviance", nfolds = 10)
  cf <- coef(cv, s = cv$lambda.min)
  l <- setdiff(rownames(cf)[as.numeric(cf) != 0], "(Intercept)")
  set.seed(seed + 2)
  d <- as.data.frame(X); d$group <- g
  nfold <- 10
  folds <- rep(1:nfold, len = nrow(d))[sample(nrow(d))]
  folds <- lapply(1:nfold, function(x) which(folds == x))
  rl <- lapply(folds, svmRFE.wrap, d, k = 10, halve.above = 50)
  tf <- WriteFeatures(rl, d, save = FALSE)
  nm <- min(30, ncol(d) - 1)
  fs <- lapply(1:nm, FeatSweep.wrap, rl, d)
  er <- sapply(fs, function(x) ifelse(is.null(x), NA, x$error))
  sv <- as.character(tf[1:which.min(er), 1])
  set.seed(seed + 3)
  rf <- randomForest(g ~ ., data = d, ntree = 500)
  opt <- which.min(rf$err.rate[, 1])
  set.seed(seed + 4)
  rf2 <- randomForest(g ~ ., data = d, ntree = opt)
  imp <- importance(rf2)
  imp <- imp[order(imp[, "MeanDecreaseGini"], decreasing = TRUE), , drop = FALSE]
  rfG <- rownames(imp)[imp[, "MeanDecreaseGini"] > 5]
  inter <- sort(Reduce(intersect, list(l, sv, rfG)))
  data.frame(analysis = tag, seed = seed,
             nLASSO = length(l), nSVM = length(sv), nRF = length(rfG),
             LASSO = paste(sort(l), collapse = ";"),
             SVM = paste(sort(sv), collapse = ";"),
             RF = paste(sort(rfG), collapse = ";"),
             Intersect = paste(inter, collapse = ";"),
             GUCA2A_in = "GUCA2A" %in% inter, CA7_in = "CA7" %in% inter,
             SST_in = "SST" %in% inter,
             stringsAsFactors = FALSE)
}
selList <- list(); degList <- list()
k <- 0
for (sd in c(123, 456, 789)) {
  set.seed(sd)
  tum <- colnames(ds33113)[tp(colnames(ds33113)) == "Treat"]
  keep <- c(colnames(ds33113)[tp(colnames(ds33113)) == "Control"], sample(tum, 30))
  d3 <- cbind(ds22598, ds33113[, keep], ds41328)
  batch <- c(rep("A", ncol(ds22598)), rep("B", length(keep)), rep("C", ncol(ds41328)))
  exc <- combatExpr(d3, batch)
  g3 <- tp(colnames(d3))
  cat("round seed", sd, " nControl:", sum(g3 == "Control"),
      " nTreat:", sum(g3 == "Treat"), "\n")
  degList[[length(degList) + 1]] <- degTable(exc, paste0("downsample_seed", sd))
  selList[[length(selList) + 1]] <- featSel(exc, g3, sd, "downsample")
}
sall <- do.call(rbind, selList)
dall <- do.call(rbind, degList)
write.csv(sall, file.path(outdir, "downsample_selection.csv"), row.names = FALSE)
write.csv(dall, file.path(outdir, "downsample_3gene_deg.csv"), row.names = FALSE)
print(sall[, c("analysis", "seed", "nLASSO", "nSVM", "nRF",
               "Intersect", "GUCA2A_in", "CA7_in", "SST_in")])
print(dall)
cat("== P0-2 done ==\n")
