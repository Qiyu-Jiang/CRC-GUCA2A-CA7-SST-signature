## ============================================================
## P1-2  ssGSEA 免疫细胞浸润（Senbabaoglu 2016 细胞签名，inRNAseq=YES）
##       + 三基因 × 细胞类型 Spearman 相关
## 输入: F:/.../CRC/06.sva/merge.normalize.txt + immune_cell_type_list.csv
## 输出: 99.P1/02.ssgsea/{ssGSEA_scores.csv, ssGSEA_boxplot.pdf, ssGSEA_corr_heatmap.pdf, ssGSEA_corr.csv}
## ============================================================
suppressMessages({library(GSVA); library(GSEABase); library(limma); library(pheatmap)})
root   <- "F:/232机器学习联合神经网络/CRC"
outdir <- file.path(root, "99.P1", "02.ssgsea")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
genes3 <- c("GUCA2A", "CA7", "SST")
cat("== P1-2 ssGSEA start ==\n")

sets.df <- read.delim(file.path(root, "99.P1", "immune_cell_type_list.tsv"), stringsAsFactors = FALSE)
sets.df <- subset(sets.df, inRNAseq == "YES")
sets.df$Symbol <- sub(" .*", "", sets.df$Symbol)
cellSets <- lapply(split(sets.df$Symbol, sets.df$Cell_type), function(x) unique(toupper(x)))
cellSets <- cellSets[sapply(cellSets, length) >= 5]
cat("cell types:", length(cellSets), "\n")

rt <- read.table(file.path(root, "06.sva", "merge.normalize.txt"), header = TRUE,
                 sep = "\t", check.names = FALSE, quote = "", comment.char = "", row.names = 1)
ex <- as.matrix(rt); storage.mode(ex) <- "numeric"
grp <- sub(".*_(Control|Treat)$", "\\1", colnames(ex))

## 新版 GSVA (>=1.50) 用参数对象；旧版回退
es <- tryCatch({
  par <- ssgseaParam(as.matrix(ex), cellSets)
  gsva(par)
}, error = function(e) {
  gsva(as.matrix(ex), cellSets, method = "ssgsea", verbose = FALSE)
})
es <- es[!grepl("^Version", rownames(es)), , drop = FALSE]
write.csv(es, file.path(outdir, "ssGSEA_scores.csv"))
cat("ssGSEA matrix:", dim(es), "\n")

## 箱图：Control vs Treat
pdf(file.path(outdir, "ssGSEA_boxplot.pdf"), width = 13, height = 5.5)
par(mar = c(9, 4.2, 3, 1))
boxplot(t(es) ~ grp, las = 2, col = c("#5C88DAFF", "#F05C3BFF"), cex.axis = 0.8,
        ylab = "ssGSEA enrichment", main = "ssGSEA immune cell types (Control vs Treat)")
dev.off()

## 三基因 × 细胞类型 Spearman
mat <- matrix(NA_real_, nrow = 3, ncol = nrow(es), dimnames = list(genes3, rownames(es)))
pmat <- mat
for (g in genes3) for (ct in rownames(es)) {
  ct.test <- suppressWarnings(cor.test(ex[g, ], es[ct, ], method = "spearman"))
  mat[g, ct] <- as.numeric(ct.test$estimate); pmat[g, ct] <- ct.test$p.value
}
stars <- ifelse(pmat < 0.001, "***", ifelse(pmat < 0.01, "**", ifelse(pmat < 0.05, "*", "")))
write.csv(cbind(t(mat), t(pmat)), file.path(outdir, "ssGSEA_corr.csv"))
pdf(file.path(outdir, "ssGSEA_corr_heatmap.pdf"), width = 12, height = 3.4)
pheatmap(mat, cluster_rows = FALSE, display_numbers = stars,
         color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(50),
         main = "Spearman: feature genes vs immune cell types (ssGSEA)",
         fontsize_number = 9, border_color = NA)
dev.off()
cat("== P1-2 done ==\n")
