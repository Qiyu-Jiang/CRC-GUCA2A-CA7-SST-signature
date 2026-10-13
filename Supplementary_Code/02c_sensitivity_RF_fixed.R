## P0-2c RF 泄漏修复版（降权抽样）：重建抽样+ComBat，重算 RF/LASSO，复用有效 SVM 结果
## 用法: Rscript P0_02c.rf_fix.R <seed>
suppressMessages({library(limma); library(sva); library(glmnet); library(randomForest)})
args <- commandArgs(trailingOnly = TRUE); sd <- as.integer(args[1])
root   <- "F:/232机器学习联合神经网络/CRC"
outdir <- file.path(root, "99.P0", "02.sensitivity")
genes3 <- c("GUCA2A", "CA7", "SST")
hub <- readLines(file.path(root, "16.hubGene", "hubGenes.txt"))
tp <- function(cn) sub(".*_(Control|Treat)$", "\\1", cn)
cat("== P0-2c seed", sd, "start ==\n")

readExpr <- function(path){
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "", comment.char = "")
  m <- as.matrix(rt[, -1]); storage.mode(m) <- "numeric"
  rownames(m) <- as.character(rt[[1]]); limma::avereps(m)
}
ds22598 <- readExpr(file.path(root, "05.normalize", "GSE22598", "GSE22598.normalize.txt"))
ds33113 <- readExpr(file.path(root, "05.normalize", "GSE33113", "GSE33113.normalize.txt"))
ds41328 <- readExpr(file.path(root, "05.normalize", "GSE41328", "GSE41328.normalize.txt"))

set.seed(sd)
tum <- colnames(ds33113)[tp(colnames(ds33113)) == "Treat"]
keep <- c(colnames(ds33113)[tp(colnames(ds33113)) == "Control"], sample(tum, 30))
d3 <- cbind(ds22598, ds33113[, keep], ds41328)
batch <- c(rep("A", ncol(ds22598)), rep("B", length(keep)), rep("C", ncol(ds41328)))
g3 <- tp(colnames(d3))
exc <- ComBat(d3, factor(batch), par.prior = TRUE)
hb <- intersect(hub, rownames(exc))
X <- t(exc[hb, , drop = FALSE])
g <- factor(g3, levels = c("Control", "Treat"))

## LASSO（复核，同种子可复现）
set.seed(sd + 1)
cv <- cv.glmnet(X, g, family = "binomial", alpha = 1,
                type.measure = "deviance", nfolds = 10)
cf <- coef(cv, s = cv$lambda.min)
l <- setdiff(rownames(cf)[as.numeric(cf) != 0], "(Intercept)")
cat("LASSO done:", length(l), "\n")
## RF（修复泄漏）
set.seed(sd + 3)
rf <- randomForest(x = X, y = g, ntree = 500)
opt <- which.min(rf$err.rate[, "OOB"])
set.seed(sd + 4)
rf2 <- randomForest(x = X, y = g, ntree = opt)
imp <- importance(rf2)
imp <- imp[order(imp[, "MeanDecreaseGini"], decreasing = TRUE), , drop = FALSE]
rfG <- rownames(imp)[imp[, "MeanDecreaseGini"] > 5]
cat("RF done:", length(rfG), "\n")
## SVM（沿用同种子有效结果）
saved <- read.csv(file.path(outdir, sprintf("downsample_selection_seed%d.csv", sd)))
sv <- strsplit(as.character(saved$SVM), ";")[[1]]
inter <- sort(Reduce(intersect, list(l, sv, rfG)))
sel <- data.frame(seed = sd,
                  LASSO = paste(sort(l), collapse = ";"), nLASSO = length(l),
                  SVM = paste(sort(sv), collapse = ";"), nSVM = length(sv),
                  RF = paste(sort(rfG), collapse = ";"), nRF = length(rfG),
                  Intersect = paste(inter, collapse = ";"),
                  GUCA2A_in = "GUCA2A" %in% inter, CA7_in = "CA7" %in% inter,
                  SST_in = "SST" %in% inter)
write.csv(sel, file.path(outdir, sprintf("downsample_selection_fixed_seed%d.csv", sd)),
          row.names = FALSE)
print(sel)
cat("== P0-2c seed", sd, "done ==\n")
