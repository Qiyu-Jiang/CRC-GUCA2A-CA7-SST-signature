## P1-8 UALCAN 启动子甲基化五数概括 -> 出版级分组箱图（bxp 版）
root <- "F:/232机器学习联合神经网络/CRC/99.P1"
outdir <- file.path(root, "06.ualcan")
d <- read.csv(file.path(outdir, "ualcan_methylation_boxstats.csv"), stringsAsFactors = FALSE)
genes <- c("GUCA2A", "CA7", "SST")
cat("== P1-8 start ==\n")

meds <- data.frame()
pdf(file.path(outdir, "UALCAN_methylation_3genes.pdf"), width = 12, height = 4.2)
par(mfrow = c(1, 3), mar = c(7, 4.5, 3.5, 1))
for (g in genes) {
  dg <- d[d$gene == g, ]
  mk <- function(grp) {
    x <- dg[grepl(grp, dg$group, fixed = TRUE), ]
    s <- x$value; names(s) <- x$stat
    unname(c(s["low"], s["q1"], s["median"], s["q3"], s["high"]))
  }
  normal <- mk("Normal"); tumor <- mk("Primary tumor")
  nLab <- gsub("[^n=0-9]", "", grep("Normal", dg$group, value = TRUE)[1])
  tLab <- gsub("[^n=0-9]", "", grep("tumor", dg$group, value = TRUE)[1])
  sm <- rbind(normal, tumor)   # 2 x 5
  meds <- rbind(meds, data.frame(gene = g,
                                 normal_median = normal[3], tumor_median = tumor[3]))
  bx <- bxp(list(stats = t(sm), n = c(nLab, tLab)),
            names = c(paste0("Normal\n(", nLab, ")"), paste0("Tumor\n(", tLab, ")")),
            col = c("grey80", "#F05C3BFF"), border = c("grey30", "#8B2E14"),
            ylab = "Promoter methylation (Beta value)",
            main = g, medcol = "black", medlwd = 2.5,
            boxfill = c("grey80", "#F05C3BFF"), outpch = NA,
            ylim = c(0, max(sm) * 1.12), cex.axis = 0.85, las = 1)
  mtext(sprintf("median: Normal=%.3f | Tumor=%.3f", normal[3], tumor[3]),
        side = 3, line = 0.4, cex = 0.75)
}
dev.off()
write.csv(meds, file.path(outdir, "UALCAN_methylation_medians.csv"), row.names = FALSE)
print(meds)

## PNG 输出（每基因一张）
for (g in genes) {
  png(file.path(outdir, paste0("UALCAN_methylation_", g, ".png")), width = 1600, height = 1300, res = 220)
  dg <- d[d$gene == g, ]
  mk <- function(grp) {
    x <- dg[grepl(grp, dg$group, fixed = TRUE), ]
    s <- x$value; names(s) <- x$stat
    unname(c(s["low"], s["q1"], s["median"], s["q3"], s["high"]))
  }
  normal <- mk("Normal"); tumor <- mk("Primary tumor")
  nLab <- gsub("[^n=0-9]", "", grep("Normal", dg$group, value = TRUE)[1])
  tLab <- gsub("[^n=0-9]", "", grep("tumor", dg$group, value = TRUE)[1])
  sm <- rbind(normal, tumor)
  par(mar = c(6.5, 4.5, 3.5, 1))
  bxp(list(stats = t(sm), n = c(nLab, tLab)),
      names = c(paste0("Normal\n(", nLab, ")"), paste0("Tumor\n(", tLab, ")")),
      col = c("grey80", "#F05C3BFF"), border = c("grey30", "#8B2E14"),
      ylab = "Promoter methylation (Beta value)",
      main = paste0(g, " — promoter methylation (UALCAN, TCGA COAD)"),
      medcol = "black", medlwd = 2.5, outpch = NA,
      ylim = c(0, max(sm) * 1.12), cex.axis = 0.9, las = 1)
  mtext(sprintf("median: Normal=%.3f | Tumor=%.3f", normal[3], tumor[3]),
        side = 3, line = 0.4, cex = 0.8)
  dev.off()
}
cat("== P1-8 done ==\n")
