## ============================================================
## P1-3  免疫检查点/共刺激分子相关性（训练集 150 样本 + TCGA 肿瘤队列）
## 输出: 99.P1/03.checkpoint/{checkpoint_corr_train.pdf/csv, checkpoint_corr_tcga.pdf/csv}
## ============================================================
suppressMessages({library(pheatmap)})
root   <- "F:/232机器学习联合神经网络/CRC"
outdir <- file.path(root, "99.P1", "03.checkpoint")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
genes3 <- c("GUCA2A", "CA7", "SST")
ck <- c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "BTLA", "IDO1",
        "CD276", "VSIR", "SIGLEC15", "TNFRSF9", "ICOS", "CD27", "CD28")
cat("== P1-3 start ==\n")

corrBlock <- function(ex, tag){
  hb <- intersect(ck, rownames(ex))
  cat(tag, "checkpoint genes available:", length(hb), "\n")
  mat <- matrix(NA_real_, nrow = length(genes3), ncol = length(hb),
                dimnames = list(genes3, hb))
  pmat <- mat
  for (g in genes3) for (c0 in hb) {
    tt <- suppressWarnings(cor.test(ex[g, ], ex[c0, ], method = "spearman"))
    mat[g, c0] <- as.numeric(tt$estimate); pmat[g, c0] <- tt$p.value
  }
  stars <- ifelse(pmat < 0.001, "***", ifelse(pmat < 0.01, "**", ifelse(pmat < 0.05, "*", "")))
  out.csv <- data.frame(gene = rep(genes3, each = length(hb)),
                        checkpoint = rep(hb, times = 3),
                        rho = as.vector(mat), pvalue = as.vector(pmat))
  write.csv(out.csv, file.path(outdir, sprintf("checkpoint_corr_%s.csv", tag)), row.names = FALSE)
  pdf(file.path(outdir, sprintf("checkpoint_corr_%s.pdf", tag)), width = 11.5, height = 3.2)
  pheatmap(mat, cluster_rows = FALSE, display_numbers = stars,
           color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(50),
           main = sprintf("Spearman: feature genes vs immune checkpoints (%s)", tag),
           fontsize_number = 9, border_color = NA)
  dev.off()
  out.csv
}

## 训练集
rt <- read.table(file.path(root, "06.sva", "merge.normalize.txt"), header = TRUE,
                 sep = "\t", check.names = FALSE, quote = "", comment.char = "", row.names = 1)
ex <- as.matrix(rt); storage.mode(ex) <- "numeric"
corrBlock(ex, "train_merge")

## TCGA 肿瘤队列（cBioPortal RSEM）
f <- file.path(root, "99.P1", "tcga_expr_wide.csv")
if (file.exists(f)) {
  w <- read.csv(f, sep = "\t", check.names = FALSE, row.names = 1)
  ex2 <- as.matrix(w); storage.mode(ex2) <- "numeric"
  corrBlock(ex2, "tcga_tumor")
} else cat("TCGA expression not found, skip\n")
cat("== P1-3 done ==\n")
