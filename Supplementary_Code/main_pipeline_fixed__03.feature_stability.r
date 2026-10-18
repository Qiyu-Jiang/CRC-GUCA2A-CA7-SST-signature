## ============================================================
## P0-3  特征筛选可复现性核查：锁定 set.seed 后重跑 LASSO/SVM-RFE/RF
##       在全量训练集（150 样本）上 × 3 个种子，与原始输出对比
## 输出: outputs/phase0/03.feature_stability/
## ============================================================
suppressMessages({library(glmnet); library(randomForest); library(e1071)})
root   <- "./"
outdir <- file.path(root, "outputs/phase0", "03.feature_stability")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
hub <- readLines(file.path(root, "16.hubGene", "hubGenes.txt"))
source(file.path(root, "18.SVM", "neuralML18.msvmRFE.R"), encoding = "UTF-8")
cat("== P0-3 start ==\n")

rt <- read.table(file.path(root, "06.sva", "merge.normalize.txt"),
                 header = TRUE, sep = "\t", check.names = FALSE,
                 quote = "", comment.char = "", row.names = 1)
ex <- as.matrix(rt); storage.mode(ex) <- "numeric"
hb <- intersect(hub, rownames(ex))
grp <- sub(".*_(Control|Treat)$", "\\1", colnames(ex))
X <- t(ex[hb, ])
g <- factor(grp, levels = c("Control", "Treat"))
cat("hub genes available:", length(hb), " samples:", nrow(X), "\n")

orig <- data.frame(LASSO = "GUCA2A;GUCA2B;CHGA;CA2;CA7;SST",
                   SVM = "CA7;SST;CA2;CLCA4;GUCA2B;CA1;GUCA2A;CA4;CHGA;ZG16;SLC26A3",
                   RF = "CA7;GUCA2A;SST",
                   Intersect = "GUCA2A;CA7;SST", stringsAsFactors = FALSE)

out <- list()
for (sd in c(123, 456, 789)) {
  set.seed(sd + 1)
  cv <- cv.glmnet(X, g, family = "binomial", alpha = 1,
                  type.measure = "deviance", nfolds = 10)
  cf <- coef(cv, s = cv$lambda.min)
  l <- setdiff(rownames(cf)[as.numeric(cf) != 0], "(Intercept)")
  set.seed(sd + 2)
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
  set.seed(sd + 3)
  rf <- randomForest(g ~ ., data = d, ntree = 500)
  opt <- which.min(rf$err.rate[, 1])
  set.seed(sd + 4)
  rf2 <- randomForest(g ~ ., data = d, ntree = opt)
  imp <- importance(rf2)
  imp <- imp[order(imp[, "MeanDecreaseGini"], decreasing = TRUE), , drop = FALSE]
  rfG <- rownames(imp)[imp[, "MeanDecreaseGini"] > 5]
  inter <- sort(Reduce(intersect, list(l, sv, rfG)))
  out[[length(out) + 1]] <- data.frame(
    seed = sd,
    LASSO = paste(sort(l), collapse = ";"), nLASSO = length(l),
    SVM = paste(sort(sv), collapse = ";"), nSVM = length(sv),
    RF = paste(sort(rfG), collapse = ";"), nRF = length(rfG),
    Intersect = paste(inter, collapse = ";"),
    match_original_LASSO = identical(sort(l), sort(strsplit(orig$LASSO, ";")[[1]])),
    match_original_SVM = identical(sort(sv), sort(strsplit(orig$SVM, ";")[[1]])),
    match_original_RF = identical(sort(rfG), sort(strsplit(orig$RF, ";")[[1]])),
    match_original_intersect = identical(inter, sort(strsplit(orig$Intersect, ";")[[1]])),
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, out)
write.csv(res, file.path(outdir, "feature_lists_seeded.csv"), row.names = FALSE)
print(res[, c("seed", "nLASSO", "nSVM", "nRF", "Intersect",
              "match_original_LASSO", "match_original_SVM",
              "match_original_RF", "match_original_intersect")])
cat("== P0-3 done ==\n")
