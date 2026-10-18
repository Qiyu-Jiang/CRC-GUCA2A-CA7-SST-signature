## ============================================================
## P1-5  特征选择频率稳定性：LASSO + RF × 100 次（80% 分层子样本）
## 输出: outputs/phase1/05.selection_freq/{selection_frequency.csv, selection_freq_barplot.pdf}
## ============================================================
suppressMessages({library(glmnet); library(randomForest)})
root   <- "./"
outdir <- file.path(root, "outputs/phase1", "05.selection_freq")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
hub <- readLines(file.path(root, "16.hubGene", "hubGenes.txt"))
genes3 <- c("GUCA2A", "CA7", "SST")
cat("== P1-5 start ==\n")

rt <- read.table(file.path(root, "06.sva", "merge.normalize.txt"), header = TRUE,
                 sep = "\t", check.names = FALSE, quote = "", comment.char = "", row.names = 1)
ex <- as.matrix(rt); storage.mode(ex) <- "numeric"
hb <- intersect(hub, rownames(ex))
grp <- sub(".*_(Control|Treat)$", "\\1", colnames(ex))
X <- t(ex[hb, ])
g <- factor(grp, levels = c("Control", "Treat"))
N_REP <- 100; FRAC <- 0.8

freqL <- setNames(rep(0, length(hb)), hb)
freqR <- setNames(rep(0, length(hb)), hb)
rankR <- setNames(rep(0, length(hb)), hb)  # RF 平均重要性排名
for (i in seq_len(N_REP)) {
  set.seed(20000 + i)
  i1 <- sample(which(g == "Treat"), round(FRAC * sum(g == "Treat")))
  i0 <- sample(which(g == "Control"), round(FRAC * sum(g == "Control")))
  idx <- c(i1, i0)
  Xi <- X[idx, , drop = FALSE]; gi <- g[idx]
  ## LASSO
  set.seed(30000 + i)
  cv <- cv.glmnet(Xi, gi, family = "binomial", alpha = 1, nfolds = 10)
  cf <- coef(cv, s = cv$lambda.min)
  l <- setdiff(rownames(cf)[as.numeric(cf) != 0], "(Intercept)")
  freqL[l] <- freqL[l] + 1
  ## RF
  set.seed(31000 + i)
  rf <- randomForest(x = Xi, y = gi, ntree = 500)
  imp <- importance(rf)[, "MeanDecreaseGini"]
  freqR[names(imp)[imp > 5]] <- freqR[names(imp)[imp > 5]] + 1
  rk <- rank(-imp); rankR <- rankR + rk
  if (i %% 20 == 0) cat("rep", i, "\n")
}
out <- data.frame(gene = hb,
                  LASSO_freq = freqL / N_REP,
                  RF_freq = freqR / N_REP,
                  RF_mean_rank = round(rankR / N_REP, 1))
out <- out[order(-out$LASSO_freq, -out$RF_freq), ]
write.csv(out, file.path(outdir, "selection_frequency.csv"), row.names = FALSE)
print(out)

pdf(file.path(outdir, "selection_freq_barplot.pdf"), width = 9, height = 5.5)
par(mfrow = c(1, 2), mar = c(9, 4.2, 3, 1))
barplot(sort(out$LASSO_freq, decreasing = TRUE), las = 2, cex.names = 0.85,
        ylab = "Selection frequency (LASSO, 100 reps)",
        main = sprintf("LASSO selection frequency\n(%d%% subsample x %d)", FRAC * 100, N_REP),
        col = ifelse(names(sort(out$LASSO_freq, decreasing = TRUE)) %in% genes3, "#c0392b", "grey70"))
barplot(sort(out$RF_freq, decreasing = TRUE), las = 2, cex.names = 0.85,
        ylab = "Selection frequency (RF Gini>5, 100 reps)",
        main = "Random forest selection frequency",
        col = ifelse(names(sort(out$RF_freq, decreasing = TRUE)) %in% genes3, "#c0392b", "grey70"))
dev.off()
cat("== P1-5 done ==\n")
