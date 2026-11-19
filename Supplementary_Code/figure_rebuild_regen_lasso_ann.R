## ============================================================
## CRC 图4 / 图5 重绘：重新生成绘图所需的原始数值
## 输入: D_共享输入数据/merge.normalize.txt, GSE110224.normalize.txt,
##       A_正文图脚本与数据/Figure_4A_LASSO/hubGenes.txt,
##       A_正文图脚本与数据/Figure_4C_特征选择频率/selection_frequency_输出参考.csv
## 输出: out/*.csv (供 Python 渲染)
## 说明: 完全复刻 17.lasso.R 与 P0_01.ANN_CV_external.R 的建模设置
## ============================================================
suppressMessages({library(glmnet); library(neuralnet); library(pROC)})
set.seed(123)

pkg    <- "F:/图复现代码与原始数据包_2026-11-19"
datdir <- file.path(pkg, "D_共享输入数据")
out    <- "C:/Users/Administrator/AppData/Local/hermes/cache/scratch/crc_redraw/out"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

## ---------- 4A: LASSO cv ----------
hub <- readLines(file.path(pkg, "A_正文图脚本与数据/Figure_4A_LASSO/hubGenes.txt"))
hub <- trimws(hub); hub <- hub[nzchar(hub)]

rt <- read.table(file.path(datdir, "merge.normalize.txt"), header = TRUE, sep = "\t",
                 check.names = FALSE, quote = "", comment.char = "")
m  <- as.matrix(rt[, -1]); storage.mode(m) <- "numeric"
rownames(m) <- as.character(rt[[1]])
gene <- rownames(m)
keep <- !duplicated(gene)
m <- m[keep, , drop = FALSE]; rownames(m) <- gene[keep]

sub <- m[hub, , drop = FALSE]
x <- t(sub)
y <- gsub("(.*)_(.*)_(.*)", "\\3", rownames(x))          # 与脚本一致
cat("4A dim:", dim(x), "| classes:", paste(names(table(y)), table(y), collapse = " "), "\n")

fit   <- glmnet(x, y, family = "binomial", alpha = 1)
cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1, type.measure = "deviance", nfolds = 10)

cvd <- data.frame(loglambda = log(cvfit$lambda), lambda = cvfit$lambda,
                  cvm = cvfit$cvm, cvsd = cvfit$cvsd, nzero = cvfit$nzero)
write.csv(cvd, file.path(out, "fig4A_cvfit.csv"), row.names = FALSE)
write.csv(data.frame(metric = c("lambda.min", "lambda.1se", "nzero.min", "nzero.1se"),
                     value  = c(cvfit$lambda.min, cvfit$lambda.1se,
                                cvfit$nzero[which(cvfit$lambda == cvfit$lambda.min)],
                                cvfit$nzero[which(cvfit$lambda == cvfit$lambda.1se)])),
          file.path(out, "fig4A_min.csv"), row.names = FALSE)

co <- as.matrix(coef(fit, s = cvfit$lambda.min))
sel <- rownames(co)[as.numeric(co) != 0]
sel <- setdiff(sel, "(Intercept)")
write.table(data.frame(gene = sel, coef = as.numeric(co[sel, 1])),
            file.path(out, "fig4A_selected.csv"), sep = ",", quote = FALSE, row.names = FALSE)
cat("4A selected:", paste(sel, collapse = "/"), "\n")
cat("4A nzero head:", paste(head(cvfit$nzero, 20), collapse = " "), "\n")
cat("4A lambda.min:", cvfit$lambda.min, " 1se:", cvfit$lambda.1se, "\n")

## ---------- 4C: selection frequency (直接读已算好的参考表) ----------
sf <- read.csv(file.path(pkg, "A_正文图脚本与数据/Figure_4C_特征选择频率/selection_frequency_输出参考.csv"),
               stringsAsFactors = FALSE)
write.csv(sf, file.path(out, "fig4C_freq.csv"), row.names = FALSE)
cat("4C rows:", nrow(sf), "\n")

## ---------- 5A/5B/5C: ANN ----------
genes3 <- c("GUCA2A", "CA7", "SST")
readExpr <- function(path){
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "", comment.char = "")
  mm <- as.matrix(rt[, -1]); storage.mode(mm) <- "numeric"
  rownames(mm) <- as.character(rt[[1]]); mm
}
trainRaw <- t(readExpr(file.path(datdir, "merge.normalize.txt"))[genes3, ])
grp <- sub(".*_(Control|Treat)$", "\\1", rownames(trainRaw))
yv  <- ifelse(grp == "Treat", 1L, 0L)
trainZ <- scale(trainRaw)

ext <- readExpr(file.path(datdir, "GSE110224.normalize.txt"))
testRaw <- t(ext[genes3, ])
grp2 <- sub(".*_(Control|Treat)$", "\\1", colnames(ext))
y2   <- ifelse(grp2 == "Treat", 1L, 0L)
testZ <- scale(testRaw)
cat("5 train dim:", dim(trainRaw), "ext dim:", dim(testRaw), "\n")

mkX <- function(X){ d <- as.data.frame(X); colnames(d) <- genes3; d }
netfit <- function(X, yy, seed){
  set.seed(seed)
  neuralnet(Control + Treat ~ GUCA2A + CA7 + SST,
            cbind(mkX(X), Control = 1 - yy, Treat = yy),
            hidden = 5, stepmax = 1e6)
}
netpred <- function(f, X){
  nn <- compute(f, mkX(X))$net.result
  colnames(nn) <- c("Control", "Treat"); nn
}

## (a) in-sample, raw
fit.raw <- netfit(trainRaw, yv, 123)
pr <- netpred(fit.raw, trainRaw)
r_in <- roc(yv, pr[, "Treat"])
write.csv(data.frame(fpr = 1 - r_in$specificities, tpr = r_in$sensitivities),
          file.path(out, "fig5A_roc.csv"), row.names = FALSE)
cat("5A AUC:", as.numeric(r_in$auc), " acc:", mean(max.col(pr) == (yv + 1)), "\n")

## (b) 10-fold CV, z-score, 分层
set.seed(123)
i1 <- sample(which(yv == 1)); i0 <- sample(which(yv == 0))
f1 <- split(i1, cut(seq_along(i1), 10)); f0 <- split(i0, cut(seq_along(i0), 10))
pm <- matrix(NA_real_, nrow = length(yv), ncol = 2, dimnames = list(NULL, c("Control", "Treat")))
for (i in 1:10) {
  te <- c(f1[[i]], f0[[i]]); tr <- setdiff(seq_along(yv), te)
  fiti <- netfit(trainZ[tr, , drop = FALSE], yv[tr], 100 + i)
  pm[te, ] <- netpred(fiti, trainZ[te, , drop = FALSE])
}
roc.cv <- roc(yv, pm[, "Treat"])
write.csv(data.frame(fpr = 1 - roc.cv$specificities, tpr = roc.cv$sensitivities),
          file.path(out, "fig5B_roc.csv"), row.names = FALSE)
cat("5B AUC:", as.numeric(roc.cv$auc), " acc:", mean(max.col(pm) == (yv + 1)), "\n")
write.csv(data.frame(sample = rownames(trainRaw), group = grp, p_treat = pm[, "Treat"]),
          file.path(out, "fig5B_predictions.csv"), row.names = FALSE)

## (c) external: 直接从已保存的预测值重算（与脚本同源）
pe <- read.csv(file.path(pkg, "A_正文图脚本与数据/Figure_5ABC_ANN_ROC/predictions_external_输出参考.csv"),
               stringsAsFactors = FALSE)
lab <- ifelse(pe$group == "Treat", 1L, 0L)
for (nm in c("ANN_z", "ANN_raw", "GLM")) {
  rr <- roc(lab, pe[[nm]])
  write.csv(data.frame(fpr = 1 - rr$specificities, tpr = rr$sensitivities),
            file.path(out, paste0("fig5C_", nm, ".csv")), row.names = FALSE)
  cat("5C", nm, "AUC:", round(as.numeric(rr$auc), 3), "\n")
}
write.csv(pe, file.path(out, "fig5C_predictions.csv"), row.names = FALSE)
cat("== done ==\n")
