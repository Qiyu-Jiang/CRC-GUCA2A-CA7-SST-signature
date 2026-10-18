# P0_07.postreview_sensitivity.R  (2026-09-28)
## (a) BH correction for immune-checkpoint correlations (per cohort)
## (b) Frozen-parameter z-score external validation (deployment-oriented)
root  <- "./"
out99 <- file.path(root, "outputs/phase0", "07.postreview_sensitivity")
dir.create(out99, showWarnings = FALSE, recursive = TRUE)
supp  <- file.path(root, "DELIVERY_2026-09-15", "3_GitHub公开资料与投稿上传包",
                   "SUBMISSION_SUPP", "Supplementary_Tables")

## ---------- (a) BH per cohort ----------
bh1 <- function(d, cohort){
  d <- d[order(d$pvalue), ]
  d$q_BH <- p.adjust(d$pvalue, method = "BH")
  d$cohort <- cohort
  d
}
tr <- read.csv(file.path(supp, "Table_S13_checkpoint_correlations_training_cohort.csv"), check.names = FALSE)
tc <- read.csv(file.path(supp, "Table_S14_checkpoint_correlations_TCGA_cohort.csv"),        check.names = FALSE)
tr <- bh1(tr, "training"); tc <- bh1(tc, "TCGA")
bh <- rbind(tr, tc)
bh <- bh[, c("cohort","gene","checkpoint","rho","pvalue","q_BH")]
write.csv(bh, file.path(supp, "Table_S36_checkpoint_correlations_BH.csv"), row.names = FALSE)

cat("=== training: tests n =", nrow(tr), "===\n")
print(tr[tr$q_BH < 0.05, c("gene","checkpoint","rho","pvalue","q_BH")])
cat("=== TCGA: tests n =", nrow(tc), "===\n")
print(tc[tc$q_BH < 0.05, c("gene","checkpoint","rho","pvalue","q_BH")])

m <- merge(tr[, c("gene","checkpoint","rho","pvalue","q_BH")],
           tc[, c("gene","checkpoint","rho","pvalue","q_BH")],
           by = c("gene","checkpoint"), suffixes = c("_tr","_tc"))
surv <- m[m$q_BH_tr < 0.05 & m$q_BH_tc < 0.05, ]
cat("=== survivors in BOTH cohorts (q<0.05 both) ===\n")
print(surv)
nom  <- m[m$pvalue_tr < 0.05 & m$pvalue_tc < 0.05, ]
cat("=== nominal in BOTH (raw P<0.05 both) ===\n")
print(nom[, c("gene","checkpoint","rho_tr","p_tr","q_BH_tr","rho_tc","p_tc","q_BH_tc")])

## ---------- (b) frozen z-score external validation ----------
## Exact replication of P0_01 model calls, only the scaling changes:
##   original: trainZ = scale(trainRaw); testZ = scale(testRaw)  (cohort-internal)
##   frozen:   extZ_f = (extRaw - mu_train) / sd_train           (train-set parameters)
genes3 <- c("GUCA2A", "CA7", "SST")
readExpr <- function(path){
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "", comment.char = "", stringsAsFactors = FALSE)
  m <- as.matrix(rt[, -1]); storage.mode(m) <- "numeric"
  rownames(m) <- as.character(rt[[1]]); m
}
trainRaw <- t(readExpr(file.path(root, "06.sva", "merge.normalize.txt"))[genes3, ])
ext      <- readExpr(file.path(root, "24.testDiff", "GSE110224.normalize.txt"))
testRaw  <- t(ext[genes3, ])
grp  <- sub(".*_(Control|Treat)$", "\\1", rownames(trainRaw))
y    <- ifelse(grp  == "Treat", 1L, 0L)
grp2 <- sub(".*_(Control|Treat)$", "\\1", colnames(ext))
y2   <- ifelse(grp2 == "Treat", 1L, 0L)
cat("train dim:", dim(trainRaw), " Control:", sum(y==0), " Treat:", sum(y==1), "\n")
cat("external dim:", dim(testRaw), " Control:", sum(y2==0), " Treat:", sum(y2==1), "\n")

mu  <- colMeans(trainRaw); sdv <- apply(trainRaw, 2, sd)
fz  <- function(X) sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
trainZ_f <- fz(trainRaw); extZ_f <- fz(testRaw)

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

library(neuralnet)
library(pROC)
nn_f <- netfit(trainZ_f, y, 123)
pn <- netpred(nn_f, extZ_f)
auc.ann_f <- as.numeric(pROC::roc(y2, pn[, "Treat"], quiet = TRUE)$auc)
ci.ann_f  <- as.numeric(pROC::ci.auc(y2, pn[, "Treat"], quiet = TRUE))

glm_f <- glm(y ~ ., data = cbind(mkX(trainZ_f), y = y), family = binomial)
pg <- as.numeric(predict(glm_f, newdata = data.frame(mkX(extZ_f), check.names = FALSE), type = "response"))
auc.glm_f <- as.numeric(pROC::roc(y2, pg, quiet = TRUE)$auc)
ci.glm_f  <- as.numeric(pROC::ci.auc(y2, pg, quiet = TRUE))

res <- data.frame(
  model = c("ANN external (frozen train-set z-params)",
            "Logistic external (frozen train-set z-params)"),
  auc = c(auc.ann_f, auc.glm_f),
  ci_low = c(ci.ann_f[1], ci.glm_f[1]),
  ci_high = c(ci.ann_f[3], ci.glm_f[3]))
write.csv(res, file.path(supp, "Table_S37_frozen_zscore_external.csv"), row.names = FALSE)
write.csv(res, file.path(out99, "frozen_zscore_external.csv"),        row.names = FALSE)
cat("=== frozen z-score external AUC ===\n"); print(res)
cat("reference: cohort-internal z ANN 0.820 / logistic 0.841; raw direct 0.772\n")
