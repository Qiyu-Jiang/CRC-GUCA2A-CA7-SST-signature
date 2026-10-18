suppressMessages({library(neuralnet); library(pROC)})
root <- "./"
outdir <- file.path(root, "outputs/phase0", "06_patient_blocked_cv")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("== P0-8 v2: cohort-blocked CV ==\n")

readExpr <- function(path){
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "", comment.char = "")
  m <- as.matrix(rt[, -1]); storage.mode(m) <- "numeric"
  rownames(m) <- as.character(rt[[1]]); m
}
genes3 <- c("GUCA2A", "CA7", "SST")
trainRaw <- t(readExpr(file.path(root, "06.sva", "merge.normalize.txt"))[genes3, ])
samples <- rownames(trainRaw)
cohort <- sub("_GSM.*", "", samples)
y <- ifelse(sub(".*_(Control|Treat)$", "\\1", samples) == "Treat", 1L, 0L)
cohorts <- unique(cohort)
cat("cohorts:", cohorts, "\n")
for (cc in cohorts) cat("  ", cc, ":", sum(cohort == cc), "samples (", sum(y[cohort==cc]==1), "tumour,", sum(y[cohort==cc]==0), "normal )\n")

# Leave-one-cohort-out CV: each cohort is held out once
pred <- rep(NA_real_, length(y))
for (cc in cohorts) {
  te <- which(cohort == cc); tr <- which(cohort != cc)
  mu <- colMeans(trainRaw[tr, , drop = FALSE])
  sdv <- apply(trainRaw[tr, , drop = FALSE], 2, sd)
  Ztr <- sweep(sweep(trainRaw[tr, , drop = FALSE], 2, mu, "-"), 2, sdv, "/")
  Zte <- sweep(sweep(trainRaw[te, , drop = FALSE], 2, mu, "-"), 2, sdv, "/")
  set.seed(600 + which(cohorts == cc))
  dd <- cbind(as.data.frame(Ztr, check.names = FALSE), Control = 1 - y[tr], Treat = y[tr])
  fit <- neuralnet(Control + Treat ~ GUCA2A + CA7 + SST, dd, hidden = 5, stepmax = 1e6)
  nn <- compute(fit, as.data.frame(Zte, check.names = FALSE))$net.result
  if (is.null(dim(nn))) nn <- matrix(nn, ncol = 1)
  if (ncol(nn) == 1) nn <- cbind(1 - nn, nn)
  colnames(nn) <- c("Control", "Treat")
  pred[te] <- nn[, "Treat"]
  cat("held out:", cc, "n =", length(te), "\n")
}

roc.p <- roc(y, pred)
ci <- as.numeric(ci.auc(roc.p, method = "delong"))
result <- data.frame(
  model = "Leave-one-cohort-out CV (patient-blocked approximation)",
  auc = as.numeric(roc.p$auc), ci_low = ci[1], ci_high = ci[3],
  n = length(y))
write.csv(result, file.path(outdir, "cohort_blocked_cv.csv"), row.names = FALSE)
write.csv(data.frame(sample = samples, cohort = cohort, y = y, pred = pred),
          file.path(outdir, "cohort_blocked_predictions.csv"), row.names = FALSE)
cat("\n== Cohort-blocked CV (LOCO) result ==\n")
cat("AUC:", round(as.numeric(roc.p$auc), 3), " CI:", round(ci[1], 3), "-", round(ci[3], 3), "\n")
cat("== P0-8 v2 done ==\n")
