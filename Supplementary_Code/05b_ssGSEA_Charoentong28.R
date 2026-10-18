## P1-9 Charoentong 28 免疫细胞签名 ssGSEA（对齐模板的原始基因集）
suppressMessages({library(GSVA); library(GSEABase); library(pheatmap)})
root   <- "./"
outdir <- file.path(root, "outputs/phase1", "02.ssgsea")
genes3 <- c("GUCA2A", "CA7", "SST")
cat("== P1-9 start ==\n")

## 读 GMT
sets <- list()
for (ln in readLines(file.path(outdir, "charoentong28.gmt"), encoding = "UTF-8")) {
  p <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  sets[[p[1]]] <- p[3:length(p)]
}
cat("cell types:", length(sets), "\n")

rt <- read.table(file.path(root, "06.sva", "merge.normalize.txt"), header = TRUE,
                 sep = "\t", check.names = FALSE, quote = "", comment.char = "", row.names = 1)
ex <- as.matrix(rt); storage.mode(ex) <- "numeric"
grp <- sub(".*_(Control|Treat)$", "\\1", colnames(ex))

es <- tryCatch({
  par <- ssgseaParam(as.matrix(ex), sets)
  gsva(par)
}, error = function(e) gsva(as.matrix(ex), sets, method = "ssgsea", verbose = FALSE))
es <- es[!grepl("^Version", rownames(es)), , drop = FALSE]
write.csv(es, file.path(outdir, "charoentong28_scores.csv"))
cat("ssGSEA matrix:", dim(es), "\n")

## 箱图
pdf(file.path(outdir, "charoentong28_boxplot.pdf"), width = 14, height = 5.5)
par(mar = c(11, 4.2, 3, 1))
boxplot(t(es) ~ grp, las = 2, col = c("#5C88DAFF", "#F05C3BFF"), cex.axis = 0.75,
        ylab = "ssGSEA enrichment", main = "Charoentong 28 immune cell types (Control vs Treat)")
dev.off()

## 三基因 × 28 细胞类型 Spearman
mat <- matrix(NA_real_, nrow = 3, ncol = nrow(es), dimnames = list(genes3, rownames(es)))
pmat <- mat
for (g in genes3) for (ct in rownames(es)) {
  tt <- suppressWarnings(cor.test(ex[g, ], es[ct, ], method = "spearman"))
  mat[g, ct] <- as.numeric(tt$estimate); pmat[g, ct] <- tt$p.value
}
stars <- ifelse(pmat < 0.001, "***", ifelse(pmat < 0.01, "**", ifelse(pmat < 0.05, "*", "")))
write.csv(cbind(t(mat), t(pmat)), file.path(outdir, "charoentong28_corr.csv"))
png(file.path(outdir, "charoentong28_corr_heatmap.png"), width = 2400, height = 900, res = 230)
pheatmap(mat, cluster_rows = FALSE, display_numbers = stars,
         color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(50),
         main = "Spearman: feature genes vs Charoentong 28 immune cell types (ssGSEA)",
         fontsize_number = 8, border_color = NA)
dev.off()
cat("== P1-9 done ==\n")
