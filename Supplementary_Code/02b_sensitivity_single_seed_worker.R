## P0-2b 单种子版：GSE33113 肿瘤抽样降权 + ComBat + limma + 特征重选
## 用法: Rscript P0_02b.sensitivity_seed.R <seed>
suppressMessages({library(limma); library(sva); library(glmnet);
                  library(randomForest); library(e1071)})
args <- commandArgs(trailingOnly = TRUE); sd <- as.integer(args[1])
root   <- "F:/232机器学习联合神经网络/CRC"
outdir <- file.path(root, "99.P0", "02.sensitivity")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
genes3 <- c("GUCA2A", "CA7", "SST")
hub <- readLines(file.path(root, "16.hubGene", "hubGenes.txt"))
tp <- function(cn) sub(".*_(Control|Treat)$", "\\1", cn)
cat("== P0-2b seed", sd, "start ==\n")

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
cat("nControl:", sum(g3 == "Control"), " nTreat:", sum(g3 == "Treat"), "\n")

exc <- ComBat(d3, factor(batch), par.prior = TRUE)
grp <- g3
design <- model.matrix(~0 + factor(grp)); colnames(design) <- c("Control", "Treat")
fit <- eBayes(contrasts.fit(lmFit(exc, design),
           makeContrasts(Treat - Control, levels = design)))
tt <- topTable(fit, number = Inf, adjust = "fdr", sort.by = "none")
deg <- data.frame(seed = sd, gene = genes3,
                  logFC = tt[genes3, "logFC"], adjP = tt[genes3, "adj.P.Val"],
                  direction = ifelse(tt[genes3, "logFC"] > 0, "up_in_tumor", "down_in_tumor"),
                  nDEG = sum(abs(tt$logFC) > 2 & tt$adj.P.Val < 0.05, na.rm = TRUE))
write.csv(deg, file.path(outdir, sprintf("downsample_3gene_deg_seed%d.csv", sd)),
          row.names = FALSE)
print(deg)

source(file.path(root, "18.SVM", "neuralML18.msvmRFE.R"), encoding = "UTF-8")
hb <- intersect(hub, rownames(exc))
X <- t(exc[hb, , drop = FALSE])
g <- factor(g3, levels = c("Control", "Treat"))
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
sel <- data.frame(seed = sd,
                  LASSO = paste(sort(l), collapse = ";"), nLASSO = length(l),
                  SVM = paste(sort(sv), collapse = ";"), nSVM = length(sv),
                  RF = paste(sort(rfG), collapse = ";"), nRF = length(rfG),
                  Intersect = paste(inter, collapse = ";"),
                  GUCA2A_in = "GUCA2A" %in% inter, CA7_in = "CA7" %in% inter,
                  SST_in = "SST" %in% inter)
write.csv(sel, file.path(outdir, sprintf("downsample_selection_seed%d.csv", sd)),
          row.names = FALSE)
print(sel[, c("seed", "nLASSO", "nSVM", "nRF", "Intersect",
              "GUCA2A_in", "CA7_in", "SST_in")])
cat("== P0-2b seed", sd, "done ==\n")
