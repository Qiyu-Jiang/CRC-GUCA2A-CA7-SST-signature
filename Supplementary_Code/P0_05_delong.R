# -*- coding: utf-8 -*-
# P0-5: DeLong tests + external PR-AUC (3rd-party M5/P0-7)
suppressMessages({library(pROC)})
root <- "F:/232\u673a\u5668\u5b66\u4e60\u8054\u5408\u795e\u7ecf\u7f51\u7edc/CRC"
outdir <- file.path(root, "99.P0", "05_delong_tests")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("== P0-5 DeLong tests ==\n")

readExpr <- function(path){
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "", comment.char = "")
  m <- as.matrix(rt[, -1]); storage.mode(m) <- "numeric"
  rownames(m) <- as.character(rt[[1]]); m
}
genes3 <- c("GUCA2A", "CA7", "SST")
ext <- readExpr(file.path(root, "24.testDiff", "GSE110224.normalize.txt"))
ext3 <- ext[genes3, ]
grp2 <- sub(".*_(Control|Treat)$", "\\1", colnames(ext))
y2 <- ifelse(grp2 == "Treat", 1L, 0L)

# per-sample predictions from P0_01
pred <- read.csv(file.path(root, "99.P0", "01.ANN_CV_external", "predictions_external.csv"))
# align by sample order (predictions_external uses colnames(ext) order)
stopifnot(nrow(pred) == length(y2))

roc.ann  <- roc(y2, pred$ANN_z)
roc.glm  <- roc(y2, pred$GLM)
roc.raw  <- roc(y2, pred$ANN_raw)
roc.guca <- roc(y2, ext3["GUCA2A", ])
roc.ca7  <- roc(y2, ext3["CA7", ])
roc.sst  <- roc(y2, ext3["SST", ])

pr_auc <- function(yy, s){
  o <- order(s, decreasing = TRUE); yy <- yy[o]
  tp <- cumsum(yy == 1); fp <- cumsum(yy != 1)
  prec <- tp / (tp + fp); rec <- tp / sum(yy == 1)
  x <- c(0, rec); v <- c(1, prec)
  sum(diff(x) * (head(v, -1) + tail(v, -1)) / 2)
}

# DeLong paired tests
tests <- data.frame()
t1 <- roc.test(roc.ann, roc.glm, method = "delong", paired = TRUE)
tests <- rbind(tests, data.frame(comparison = "ANN vs logistic", D = t1$statistic, p = t1$p.value))
t2 <- roc.test(roc.ann, roc.sst, method = "delong", paired = TRUE)
tests <- rbind(tests, data.frame(comparison = "ANN vs SST", D = t2$statistic, p = t2$p.value))
t3 <- roc.test(roc.ann, roc.guca, method = "delong", paired = TRUE)
tests <- rbind(tests, data.frame(comparison = "ANN vs GUCA2A", D = t3$statistic, p = t3$p.value))
t4 <- roc.test(roc.ann, roc.ca7, method = "delong", paired = TRUE)
tests <- rbind(tests, data.frame(comparison = "ANN vs CA7", D = t4$statistic, p = t4$p.value))
t5 <- roc.test(roc.glm, roc.sst, method = "delong", paired = TRUE)
tests <- rbind(tests, data.frame(comparison = "logistic vs SST", D = t5$statistic, p = t5$p.value))
t6 <- roc.test(roc.glm, roc.guca, method = "delong", paired = TRUE)
tests <- rbind(tests, data.frame(comparison = "logistic vs GUCA2A", D = t6$statistic, p = t6$p.value))
t7 <- roc.test(roc.glm, roc.ca7, method = "delong", paired = TRUE)
tests <- rbind(tests, data.frame(comparison = "logistic vs CA7", D = t7$statistic, p = t7$p.value))

write.csv(tests, file.path(outdir, "delong_external.csv"), row.names = FALSE)
cat("\n== DeLong tests (external, paired) ==\n"); print(tests, digits = 4)

# External PR-AUC
pr <- data.frame(
  model = c("ANN z-scored", "ANN raw", "Logistic", "GUCA2A", "CA7", "SST"),
  pr_auc = c(pr_auc(y2, pred$ANN_z), pr_auc(y2, pred$ANN_raw), pr_auc(y2, pred$GLM),
             pr_auc(y2, ext3["GUCA2A", ]), pr_auc(y2, ext3["CA7", ]), pr_auc(y2, ext3["SST", ])),
  auc = c(as.numeric(roc.ann$auc), as.numeric(roc.raw$auc), as.numeric(roc.glm$auc),
          as.numeric(roc.guca$auc), as.numeric(roc.ca7$auc), as.numeric(roc.sst$auc)))
write.csv(pr, file.path(outdir, "external_PR-AUC.csv"), row.names = FALSE)
cat("\n== External PR-AUC ==\n"); print(pr, digits = 4)

# AUC values for reference
cat("\n== External AUC values ==\n")
cat("  ANN-z:", round(as.numeric(roc.ann$auc), 3), "\n")
cat("  GLM:  ", round(as.numeric(roc.glm$auc), 3), "\n")
cat("  SST:  ", round(as.numeric(roc.sst$auc), 3), "\n")
cat("  GUCA2A:", round(as.numeric(roc.guca$auc), 3), "\n")
cat("  CA7:  ", round(as.numeric(roc.ca7$auc), 3), "\n")
cat("== P0-5 done ==\n")
