## 用管线同款 clusterProfiler::GSEA 复核 NES
suppressMessages({library(clusterProfiler); library(limma)})
mtx <- "F:/图复现代码与原始数据包_2026-11-19/D_共享输入数据/merge.normalize.txt"
gmt <- "C:/Users/Administrator/AppData/Local/hermes/cache/scratch/fgsea_pkg/2_管线GSEA完整表_KEGG/c2.cp.kegg.Hs.symbols.gmt"

rt <- read.table(mtx, header = TRUE, sep = "\t", check.names = FALSE,
                 quote = "", comment.char = "", row.names = 1)
data <- as.matrix(rt); storage.mode(data) <- "numeric"
data <- avereps(data)
data <- data[rowMeans(data) > 0, ]
Type <- gsub("(.*)_(.*)_(.*)", "\\3", colnames(data))
data <- data[, Type == "Treat", drop = FALSE]
cat("tumour-only:", dim(data), "\n")

ln <- readLines(gmt)
parts <- strsplit(ln, "\t", fixed = TRUE)
t2g <- do.call(rbind, lapply(parts, function(x) {
  gs <- x[3:length(x)]; gs <- gs[gs != ""]
  if (length(gs)) data.frame(term = x[1], gene = gs, stringsAsFactors = FALSE)
}))
cat("terms:", length(unique(t2g$term)), " pairs:", nrow(t2g), "\n")

KEY <- c("KEGG_RETINOL_METABOLISM", "KEGG_DRUG_METABOLISM_CYTOCHROME_P450",
         "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450",
         "KEGG_CELL_CYCLE", "KEGG_ECM_RECEPTOR_INTERACTION")

for (gene in c("GUCA2A", "CA7", "SST")) {
  v <- data[gene, ]
  meanL <- rowMeans(data[, v <  median(v), drop = FALSE])
  meanH <- rowMeans(data[, v >= median(v), drop = FALSE])
  meanL[meanL < 1e-05] <- 1e-05; meanH[meanH < 1e-05] <- 1e-05
  stats <- sort(meanH - meanL, decreasing = TRUE)

  eg <- GSEA(stats, TERM2GENE = t2g, pvalueCutoff = 1, pAdjustMethod = "BH",
             eps = 0, verbose = FALSE, seed = FALSE)
  df <- as.data.frame(eg)
  cat("\n=====", gene, "rows:", nrow(df), "p<0.05:", sum(df$pvalue < 0.05), "\n")
  for (k in KEY) {
    r <- df[df$ID == k, ]
    if (nrow(r)) cat(sprintf("  %-46s NES=%+6.3f p=%.3g padj=%.3g\n", k, r$NES, r$pvalue, r$p.adjust))
  }
  top <- df[order(df$pvalue), ][1:3, ]
  cat("  top3:", paste(sprintf("%s(%.2f)", top$ID, top$NES), collapse = " | "), "\n")
}
cat("\n== done ==\n")
