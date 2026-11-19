## 复现 TCGA KM（median split），导出绘图数据
suppressMessages({library(survival)})
setwd("F:/图复现代码与原始数据包_2026-11-19/A_正文图脚本与数据/Figure_7AB_TCGA_KM突变CNA")
out <- "C:/Users/Administrator/AppData/Local/hermes/cache/scratch/crc_redraw/out"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

expr <- read.delim("tcga_expr_wide.csv", check.names=FALSE, stringsAsFactors=FALSE)
osd  <- read.csv("tcga_os.csv", sep="\t", stringsAsFactors=FALSE)
osd$event <- ifelse(grepl("DECEASED", osd$OS_STATUS), 1, 0)

res <- list()
for (g in c("GUCA2A","CA7","SST")) {
  v <- as.numeric(expr[expr$gene==g, -1]); ids <- colnames(expr)[-1]
  pat <- sub("(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}).*","\\1", ids)
  df <- data.frame(pat=pat, expr=v, stringsAsFactors=FALSE)
  df <- merge(df, osd[,c("patientId","OS_MONTHS","event")], by.x="pat", by.y="patientId")
  df <- df[!is.na(df$expr) & !is.na(df$OS_MONTHS), ]
  grp <- ifelse(df$expr > median(df$expr), "High", "Low")
  df$grp <- factor(grp, levels=c("Low","High"))
  sd <- survdiff(Surv(OS_MONTHS, event) ~ grp, data=df)
  p  <- 1 - pchisq(sd$chisq, length(sd$n)-1)
  hr <- exp(coef(coxph(Surv(OS_MONTHS, event) ~ grp, data=df)))[[1]]
  cat(sprintf("%s: n=%d (%s/%s) p=%.6f HR=%.4f\n", g, nrow(df),
              sum(df$grp=="High"), sum(df$grp=="Low"), p, hr))
  # KM step data per group
  km <- survfit(Surv(OS_MONTHS, event) ~ grp, data=df)
  st <- data.frame(time=km$time, surv=km$surv, n.risk=km$n.risk, n.event=km$n.event,
                   strata=rep(names(km$strata), km$strata))
  st$gene <- g
  res[[g]] <- st
  write.csv(df[,c("pat","OS_MONTHS","event","expr","grp")],
            file.path(out, paste0("km_", g, "_patients.csv")), row.names=FALSE)
}
all <- do.call(rbind, res)
write.csv(all, file.path(out, "km_steps.csv"), row.names=FALSE)
cat("KM n rows:", nrow(all), "\n")
