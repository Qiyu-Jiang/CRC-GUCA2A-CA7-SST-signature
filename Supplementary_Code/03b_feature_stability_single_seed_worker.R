## P0-3b 单种子版：全量训练集(150样本)上锁定种子重跑 LASSO/SVM-RFE/RF
## 用法: Rscript P0_03b.stability_seed.R <seed>
suppressMessages({library(glmnet); library(randomForest); library(e1071)})
args <- commandArgs(trailingOnly = TRUE); sd <- as.integer(args[1])
root   <- "F:/232机器学习联合神经网络/CRC"
outdir <- file.path(root, "99.P0", "03.feature_stability")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
hub <- readLines(file.path(root, "16.hubGene", "hubGenes.txt"))
source(file.path(root, "18.SVM", "neuralML18.msvmRFE.R"), encoding = "UTF-8")
cat("== P0-3b seed", sd, "start ==\n")

rt <- read.table(file.path(root, "06.sva", "merge.normalize.txt"),
                 header = TRUE, sep = "\t", check.names = FALSE,
                 quote = "", comment.char = "", row.names = 1)
ex <- as.matrix(rt); storage.mode(ex) <- "numeric"
hb <- intersect(hub, rownames(ex))
grp <- sub(".*_(Control|Treat)$", "\\1", colnames(ex))
X <- t(ex[hb, ])
g <- factor(grp, levels = c("Control", "Treat"))
cat("hub:", length(hb), " samples:", nrow(X), "\n")

orig <- list(LASSO = strsplit("GUCA2A;GUCA2B;CHGA;CA2;CA7;SST", ";")[[1]],
             SVM = strsplit("CA7;SST;CA2;CLCA4;GUCA2B;CA1;GUCA2A;CA4;CHGA;ZG16;SLC26A3", ";")[[1]],
             RF = strsplit("CA7;GUCA2A;SST", ";")[[1]],
             Intersect = strsplit("GUCA2A;CA7;SST", ";")[[1]])

set.seed(sd + 1)
cv <- cv.glmnet(X, g, family = "binomial", alpha = 1,
                type.measure = "deviance", nfolds = 10)
cf <- coef(cv, s = cv$lambda.min)
l <- setdiff(rownames(cf)[as.numeric(cf) != 0], "(Intercept)")
cat("LASSO done:", length(l), "genes\n")
set.seed(sd + 2)
d <- data.frame(group = g, X, check.names = FALSE)
nfold <- 10
folds <- rep(1:nfold, len = nrow(d))[sample(nrow(d))]
folds <- lapply(1:nfold, function(x) which(folds == x))
rl <- lapply(folds, svmRFE.wrap, d, k = 10, halve.above = 50)
tf <- WriteFeatures(rl, d, save = FALSE)
nm <- min(30, ncol(d) - 1)
fs <- lapply(1:nm, FeatSweep.wrap, rl, d)
er <- sapply(fs, function(x) ifelse(is.null(x), NA, x$error))
sv <- as.character(tf[1:which.min(er), 1])
cat("SVM-RFE done:", length(sv), "genes\n")
set.seed(sd + 3)
rf <- randomForest(g ~ ., data = d, ntree = 500)
opt <- which.min(rf$err.rate[, 1])
set.seed(sd + 4)
rf2 <- randomForest(g ~ ., data = d, ntree = opt)
imp <- importance(rf2)
imp <- imp[order(imp[, "MeanDecreaseGini"], decreasing = TRUE), , drop = FALSE]
rfG <- rownames(imp)[imp[, "MeanDecreaseGini"] > 5]
cat("RF done:", length(rfG), "genes\n")
inter <- sort(Reduce(intersect, list(l, sv, rfG)))
res <- data.frame(seed = sd,
                  LASSO = paste(sort(l), collapse = ";"), nLASSO = length(l),
                  SVM = paste(sort(sv), collapse = ";"), nSVM = length(sv),
                  RF = paste(sort(rfG), collapse = ";"), nRF = length(rfG),
                  Intersect = paste(inter, collapse = ";"),
                  match_original_LASSO = identical(sort(l), orig$LASSO),
                  match_original_SVM = identical(sort(sv), orig$SVM),
                  match_original_RF = identical(sort(rfG), orig$RF),
                  match_original_intersect = identical(inter, orig$Intersect))
write.csv(res, file.path(outdir, sprintf("feature_lists_seeded_seed%d.csv", sd)),
          row.names = FALSE)
print(res[, c("seed", "nLASSO", "nSVM", "nRF", "Intersect",
              "match_original_LASSO", "match_original_SVM",
              "match_original_RF", "match_original_intersect")])
cat("== P0-3b seed", sd, "done ==\n")
