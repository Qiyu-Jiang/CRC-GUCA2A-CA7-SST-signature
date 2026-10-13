## ============================================================
## P1-4  TCGA-COADREAD (cBioPortal) 三基因：表达分布 + KM 生存 + 突变/CNA 汇总
## 输出: 99.P1/04.tcga_cbio/{expr_boxplot.pdf, KM_*.pdf x3, KM_summary.csv, mut_cna_bar.pdf}
## ============================================================
suppressMessages({library(survival)})
root   <- "F:/232机器学习联合神经网络/CRC"
outdir <- file.path(root, "99.P1", "04.tcga_cbio")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
genes3 <- c("GUCA2A", "CA7", "SST")
cat("== P1-4 start ==\n")

w <- read.csv(file.path(root, "99.P1", "tcga_expr_wide.csv"), sep = "\t", check.names = FALSE, row.names = 1)
ex <- as.matrix(w); storage.mode(ex) <- "numeric"
s2p <- read.csv(file.path(root, "99.P1", "tcga_sample_patient.tsv"), sep = "\t")
os  <- read.csv(file.path(root, "99.P1", "tcga_os.csv"), sep = "\t")
os$status <- ifelse(grepl("DECEASED", os$OS_STATUS), 1, 0)
os$months <- suppressWarnings(as.numeric(os$OS_MONTHS))
os <- os[!is.na(os$months) & !is.na(os$status), ]
cat("patients with OS:", nrow(os), "\n")

## 样本级基因表达 → 患者级（一个患者多个样本取均值, 只保留 primary tumor 01）
s2p <- s2p[grepl("^TCGA-", s2p$sampleId), ]
s2p$sampleType <- substr(s2p$sampleId, 14, 15)
s2p <- s2p[s2p$sampleType == "01", ]
cat("primary tumor samples:", nrow(s2p), "\n")

res <- data.frame(gene = genes3, n = NA, pval = NA, hr = NA, high_med = NA, low_med = NA, stringsAsFactors = FALSE)
pdf(file.path(outdir, "KM_3genes.pdf"), width = 12, height = 4)
par(mfrow = c(1, 3))
for (g in genes3) {
  sub <- s2p[s2p$sampleId %in% colnames(ex), c("sampleId", "patientId")]
  sub$expr <- as.numeric(ex[g, sub$sampleId])
  sub$expr[is.na(sub$expr)] <- NA
  agg <- aggregate(expr ~ patientId, data = sub, FUN = mean, na.rm = TRUE)
  m <- merge(agg, os, by = "patientId")
  m <- m[!is.na(m$expr), ]
  med <- median(m$expr)
  m$grp <- ifelse(m$expr >= med, "High", "Low")
  res$n[res$gene == g] <- nrow(m)
  res$high_med[res$gene == g] <- sum(m$grp == "High")
  res$low_med[res$gene == g] <- sum(m$grp == "Low")
  fit <- survfit(Surv(months, status) ~ grp, data = m)
  lr <- survdiff(Surv(months, status) ~ grp, data = m)
  pv <- 1 - pchisq(lr$chisq, 1)
  cox <- tryCatch(exp(coef(coxph(Surv(months, status) ~ grp, data = m)))[1], error = function(e) NA)
  res$pval[res$gene == g] <- pv; res$hr[res$gene == g] <- cox
  plot(fit, col = c("#c0392b", "#2563a8"), lwd = 2, xlab = "Overall survival (months)",
       ylab = "Probability", main = sprintf("%s (median split)", g))
  legend("topright", bty = "n", col = c("#c0392b", "#2563a8"), lwd = 2, lty = 1,
         legend = c(sprintf("High (n=%d)", sum(m$grp == "High")), sprintf("Low (n=%d)", sum(m$grp == "Low"))))
  mtext(sprintf("log-rank p=%.3g, HR=%.2f", pv, cox), side = 3, cex = 0.8)
}
dev.off()
write.csv(res, file.path(outdir, "KM_summary.csv"), row.names = FALSE)
print(res)

## 表达分布（肿瘤队列内 3 基因 RSEM）
pdf(file.path(outdir, "expr_boxplot.pdf"), width = 6, height = 5)
e3 <- ex[genes3, , drop = FALSE]
boxplot(t(e3), col = c("#F05C3BFF", "#5C88DAFF", "#5CB85CFF"), log = "y",
        ylab = "RSEM (log scale)", main = "TCGA-COADREAD tumor expression")
dev.off()

## 突变/CNA
mc <- read.csv(file.path(root, "99.P1", "tcga_mut_cna.csv"), sep = "\t", check.names = FALSE)
pdf(file.path(outdir, "mut_cna_bar.pdf"), width = 8, height = 4.5)
barplot(rbind(mc$mutatedSamples, mc$CNA_amp, mc$CNA_del),
        beside = TRUE, names.arg = mc$gene, col = c("#7f8c8d", "#c0392b", "#2563a8"),
        ylab = "Number of samples", main = "TCGA-COADREAD (n=~637)")
legend("topright", bty = "n", fill = c("#7f8c8d", "#c0392b", "#2563a8"),
       legend = c("Non-synonymous mutation", "CNA amplification (GISTIC2)", "CNA deep deletion"))
dev.off()
print(mc)
cat("== P1-4 done ==\n")
