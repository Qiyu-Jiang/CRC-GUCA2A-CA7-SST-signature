## P0-3c RF 泄漏修复版（全量训练集）：重算 RF（无 group 泄漏）+ 复核 LASSO + 复用有效 SVM 结果
## 用法: Rscript P0_03c.rf_fix.R <seed>
suppressMessages({library(glmnet); library(randomForest)})
args <- commandArgs(trailingOnly = TRUE); sd <- as.integer(args[1])
root   <- "./"
outdir <- file.path(root, "outputs/phase0", "03.feature_stability")
genes3 <- c("GUCA2A", "CA7", "SST")
hub <- readLines(file.path(root, "16.hubGene", "hubGenes.txt"))
cat("== P0-3c seed", sd, "start ==\n")

rt <- read.table(file.path(root, "06.sva", "merge.normalize.txt"),
                 header = TRUE, sep = "\t", check.names = FALSE,
                 quote = "", comment.char = "", row.names = 1)
ex <- as.matrix(rt); storage.mode(ex) <- "numeric"
hb <- intersect(hub, rownames(ex))
grp <- sub(".*_(Control|Treat)$", "\\1", colnames(ex))
X <- t(ex[hb, ])
g <- factor(grp, levels = c("Control", "Treat"))
orig <- list(LASSO = strsplit("GUCA2A;GUCA2B;CHGA;CA2;CA7;SST", ";")[[1]],
             SVM = strsplit("CA7;SST;CA2;CLCA4;GUCA2B;CA1;GUCA2A;CA4;CHGA;ZG16;SLC26A3", ";")[[1]],
             RF = strsplit("CA7;GUCA2A;SST", ";")[[1]],
             Intersect = strsplit("GUCA2A;CA7;SST", ";")[[1]])

## LASSO（复核）
set.seed(sd + 1)
cv <- cv.glmnet(X, g, family = "binomial", alpha = 1,
                type.measure = "deviance", nfolds = 10)
cf <- coef(cv, s = cv$lambda.min)
l <- setdiff(rownames(cf)[as.numeric(cf) != 0], "(Intercept)")
## RF（修复：x/y 接口，绝无 group 列泄漏）
set.seed(sd + 3)
rf <- randomForest(x = X, y = g, ntree = 500)
opt <- which.min(rf$err.rate[, "OOB"])
set.seed(sd + 4)
rf2 <- randomForest(x = X, y = g, ntree = opt)
imp <- importance(rf2)
imp <- imp[order(imp[, "MeanDecreaseGini"], decreasing = TRUE), , drop = FALSE]
rfG <- rownames(imp)[imp[, "MeanDecreaseGini"] > 5]
## SVM（沿用同种子已跑出的有效结果）
saved <- read.csv(file.path(outdir, sprintf("feature_lists_seeded_seed%d.csv", sd)))
sv <- strsplit(as.character(saved$SVM), ";")[[1]]
inter <- sort(Reduce(intersect, list(l, sv, rfG)))
res <- data.frame(seed = sd,
                  LASSO = paste(sort(l), collapse = ";"), nLASSO = length(l),
                  SVM = paste(sort(sv), collapse = ";"), nSVM = length(sv),
                  RF = paste(sort(rfG), collapse = ";"), nRF = length(rfG),
                  Intersect = paste(inter, collapse = ";"),
                  GUCA2A_in = "GUCA2A" %in% inter, CA7_in = "CA7" %in% inter,
                  SST_in = "SST" %in% inter,
                  lasso_setequal_orig = setequal(l, orig$LASSO),
                  rf_setequal_orig = setequal(rfG, orig$RF),
                  intersect_setequal_orig = setequal(inter, orig$Intersect),
                  stringsAsFactors = FALSE)
write.csv(res, file.path(outdir, sprintf("feature_lists_fixed_seed%d.csv", sd)),
          row.names = FALSE)
print(res)
cat("== P0-3c seed", sd, "done ==\n")
