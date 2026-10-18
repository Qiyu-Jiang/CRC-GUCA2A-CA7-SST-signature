## ============================================================
## P0-1  ANN 内部 10 折交叉验证 + GSE110224 外部验证 + logistic 对比
## 输出: outputs/phase0/01.ANN_CV_external/
##   ROC_insample_repro.pdf   原模型(raw)训练集内 AUC 复现
##   ROC_CV.pdf               10 折交叉验证 ROC（模型级）
##   ROC_external.pdf         外部验证 ROC：ANN(z-score) / ANN(raw直推) / logistic
##   results_summary.csv      全部 AUC(95%CI DeLong) 与准确率
##   predictions_external.csv 外部集样本级预测值
## 说明: 训练集 = 06.sva/merge.normalize.txt (150 样本, 3 特征基因)
##       外部集 = 24.testDiff/GSE110224.normalize.txt (34 样本)
##       z-score 方案: 训练集与外部集各自内部标准化后建模/预测
## ============================================================
suppressMessages({library(neuralnet); library(pROC)})
set.seed(123)
root   <- "./"
outdir <- file.path(root, "outputs/phase0", "01.ANN_CV_external")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
genes3 <- c("GUCA2A", "CA7", "SST")
cat("== P0-1 start ==\n")

readExpr <- function(path){
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "", comment.char = "")
  m <- as.matrix(rt[, -1]); storage.mode(m) <- "numeric"
  rownames(m) <- as.character(rt[[1]]); m
}

trainRaw <- t(readExpr(file.path(root, "06.sva", "merge.normalize.txt"))[genes3, ])
grp <- sub(".*_(Control|Treat)$", "\\1", rownames(trainRaw))
y   <- ifelse(grp == "Treat", 1L, 0L)
trainZ <- scale(trainRaw)
cat("train dim:", dim(trainRaw), " Control:", sum(y==0), " Treat:", sum(y==1), "\n")

ext <- readExpr(file.path(root, "24.testDiff", "GSE110224.normalize.txt"))
testRaw <- t(ext[genes3, ])
grp2 <- sub(".*_(Control|Treat)$", "\\1", colnames(ext))
y2   <- ifelse(grp2 == "Treat", 1L, 0L)
testZ <- scale(testRaw)
cat("external dim:", dim(testRaw), " Control:", sum(y2==0), " Treat:", sum(y2==1), "\n")

mkX <- function(X){ d <- as.data.frame(X); colnames(d) <- genes3; d }
netfit <- function(X, yy, seed){
  set.seed(seed)
  neuralnet(Control + Treat ~ GUCA2A + CA7 + SST,
            cbind(mkX(X), Control = 1 - yy, Treat = yy),
            hidden = 5, stepmax = 1e6)
}
netpred <- function(fit, X){
  nn <- compute(fit, mkX(X))$net.result
  colnames(nn) <- c("Control", "Treat"); nn
}
aucRow <- function(name, rocObj, acc){
  ci <- as.numeric(ci.auc(rocObj, method = "delong"))
  data.frame(model = name, auc = as.numeric(rocObj$auc),
             ci_low = ci[1], ci_high = ci[3], accuracy = acc)
}
results <- data.frame()

## ---------- (a) 复现原模型（raw scale, 训练集内） ----------
fit.raw <- netfit(trainRaw, y, 123)
pr <- netpred(fit.raw, trainRaw)
r <- roc(y, pr[, "Treat"])
acc <- mean(max.col(pr) == (y + 1))
results <- rbind(results, aucRow("ANN raw in-sample (repro)", r, acc))
pdf(file.path(outdir, "ROC_insample_repro.pdf"), width = 5.5, height = 5)
plot(r, print.auc = TRUE, col = "red", legacy.axes = TRUE,
     main = "ANN in-sample (raw inputs, repro)")
dev.off()
cat("raw in-sample AUC:", round(as.numeric(r$auc), 3), "(原文为 0.987)\n")

## ---------- (b) 10 折交叉验证（z-score, 分层抽样） ----------
set.seed(123)
i1 <- sample(which(y == 1)); i0 <- sample(which(y == 0))
f1 <- split(i1, cut(seq_along(i1), 10)); f0 <- split(i0, cut(seq_along(i0), 10))
pm <- matrix(NA_real_, nrow = length(y), ncol = 2,
             dimnames = list(NULL, c("Control", "Treat")))
for (i in 1:10) {
  te <- c(f1[[i]], f0[[i]]); tr <- setdiff(seq_along(y), te)
  fit <- netfit(trainZ[tr, , drop = FALSE], y[tr], 100 + i)
  pr  <- netpred(fit, trainZ[te, , drop = FALSE])
  pm[te, ] <- pr
}
roc.cv <- roc(y, pm[, "Treat"])
acc.cv <- mean(max.col(pm) == (y + 1))
results <- rbind(results, aucRow("ANN 10-fold CV (z-score)", roc.cv, acc.cv))
pdf(file.path(outdir, "ROC_CV.pdf"), width = 5.5, height = 5)
plot(roc.cv, print.auc = TRUE, col = "#2563a8", legacy.axes = TRUE,
     main = "ANN 10-fold cross-validation")
mtext(sprintf("Accuracy = %.3f (10-fold, held-out)", acc.cv), side = 1, line = 2.5, cex = 0.8)
dev.off()
cat("CV AUC:", round(as.numeric(roc.cv$auc), 3), " accuracy:", round(acc.cv, 3), "\n")

## ---------- (c) 外部验证 GSE110224 ----------
fit.z <- netfit(trainZ, y, 123)
pe.ann  <- netpred(fit.z, testZ)
pe.raw  <- netpred(fit.raw, testRaw)
df.tr <- data.frame(trainZ, check.names = FALSE); df.tr$y <- factor(y)
glmfit <- glm(y ~ ., family = binomial, data = df.tr)
pr.glm <- as.numeric(predict(glmfit, newdata = data.frame(testZ, check.names = FALSE),
                             type = "response"))
roc.ann  <- roc(y2, pe.ann[, "Treat"]);  acc.ann  <- mean(max.col(pe.ann) == (y2 + 1))
roc.raw2 <- roc(y2, pe.raw[, "Treat"]);  acc.raw2 <- mean(max.col(pe.raw) == (y2 + 1))
roc.glm  <- roc(y2, pr.glm);             acc.glm  <- mean((pr.glm >= 0.5) == (y2 == 1))
results <- rbind(results,
                 aucRow("ANN external GSE110224 (z-score)", roc.ann, acc.ann),
                 aucRow("ANN external GSE110224 (raw direct)", roc.raw2, acc.raw2),
                 aucRow("Logistic external GSE110224 (z-score)", roc.glm, acc.glm))
for (g in genes3) {
  r1 <- roc(y,  trainRaw[, g]); results <- rbind(results, aucRow(paste0(g, " train AUC"), r1, NA))
  r2 <- roc(y2, testRaw[, g]);  results <- rbind(results, aucRow(paste0(g, " external AUC"), r2, NA))
}
pdf(file.path(outdir, "ROC_external.pdf"), width = 6, height = 5.5)
plot(roc.ann, col = "#c0392b", lwd = 2, legacy.axes = TRUE, main = "External validation: GSE110224 (n=34)")
plot(roc.glm, col = "#2563a8", lwd = 2, add = TRUE)
plot(roc.raw2, col = "#7f8c8d", lwd = 2, lty = 2, add = TRUE)
legend("bottomright", bty = "n", lwd = 2, col = c("#c0392b", "#2563a8", "#7f8c8d"), lty = c(1, 1, 2),
       legend = c(sprintf("ANN (z-score) AUC=%.3f", as.numeric(roc.ann$auc)),
                  sprintf("Logistic (z-score) AUC=%.3f", as.numeric(roc.glm$auc)),
                  sprintf("ANN (raw direct) AUC=%.3f", as.numeric(roc.raw2$auc))))
dev.off()

write.csv(results, file.path(outdir, "results_summary.csv"), row.names = FALSE)
write.csv(data.frame(sample = colnames(ext), group = grp2,
                     ANN_z = pe.ann[, "Treat"], ANN_raw = pe.raw[, "Treat"], GLM = pr.glm),
          file.path(outdir, "predictions_external.csv"), row.names = FALSE)
cat("== results ==\n"); print(results, digits = 3)
cat("== P0-1 done ==\n")
