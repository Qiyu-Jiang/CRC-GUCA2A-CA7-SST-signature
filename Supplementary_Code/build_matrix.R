## ============================================================
## build_matrix.R  —  build the three-cohort training matrix from GEO series matrices
##
## Inputs  : GEO/  containing GSE22598.txt.gz, GSE33113.txt.gz, GSE41328.txt.gz
##           (series matrix files; download from
##            https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22598/matrix/GSE22598_series_matrix.txt.gz etc.)
## Outputs : GEO/train_matrix.rds   (21356 genes x 150 samples; 117 tumour / 33 normal)
##
## Notes   : requires hgu133plus2.db (Bioconductor) for probe -> symbol mapping.
##           GUCA2A is carried by probe 207003_at, which the annotation package does not
##           map to a symbol in some builds; it is therefore assigned explicitly from the
##           platform probe table so that all three signature genes are present.
## ============================================================
suppressMessages({library(hgu133plus2.db); library(AnnotationDbi)})

data_dir <- Sys.getenv("CRC_GEO_DIR", unset = "GEO")
coh <- c("GSE22598", "GSE33113", "GSE41328")

read_series <- function(gse) {
  ln <- readLines(gzfile(file.path(data_dir, paste0(gse, ".txt.gz"))), warn = FALSE)
  ttl <- gsub('"', "", strsplit(grep("^!Sample_title", ln, value = TRUE), "\t")[[1]][-1])
  b <- grep("^!series_matrix_table_begin", ln); e <- grep("^!series_matrix_table_end", ln)
  tab <- read.delim(text = paste(ln[(b + 1):(e - 1)], collapse = "\n"), check.names = FALSE)
  m <- as.matrix(tab[, -1]); rownames(m) <- tab[[1]]; mode(m) <- "numeric"; colnames(m) <- ttl
  m
}

classify <- function(t) {
  low <- tolower(t); g <- rep(NA_character_, length(t))
  g[grepl("cell line|ht29|hct15", low)] <- NA
  g[grepl("normal|mucosa|non-?tumou?r|adjacent", low) & is.na(g)] <- "Normal"
  g[grepl("cancer|tumou?r|adenocarcinoma|^col[0-9]", low) & is.na(g)] <- "Tumor"
  g
}

E <- NULL; lab <- NULL; cohvec <- NULL
for (g in coh) {
  m <- read_series(g)
  grp <- classify(colnames(m)); keep <- !is.na(grp)
  guca_probe <- if ("207003_at" %in% rownames(m)) as.numeric(m["207003_at", keep]) else NULL
  sym <- mapIds(hgu133plus2.db, keys = rownames(m), keytype = "PROBEID", column = "SYMBOL",
                multiVals = "first")
  ok <- !is.na(sym); m <- m[ok, keep, drop = FALSE]; sym <- sym[ok]
  ord <- order(rowMeans(m), decreasing = TRUE); m <- m[ord, ]; sym <- sym[ord]
  dup <- duplicated(sym); m <- m[!dup, , drop = FALSE]; rownames(m) <- sym[!dup]
  if (max(m, na.rm = TRUE) > 100) {
    m <- log2(m + 1)
    if (!is.null(guca_probe)) guca_probe <- log2(guca_probe + 1)
  }
  if (!"GUCA2A" %in% rownames(m) && !is.null(guca_probe)) m <- rbind(m, GUCA2A = guca_probe)
  E <- cbind(E, m); lab <- c(lab, as.integer(grp[keep] == "Tumor")); cohvec <- c(cohvec, rep(g, ncol(m)))
  cat(g, ":", ncol(m), "samples,", sum(grp[keep] == "Tumor"), "tumour\n")
}
saveRDS(list(E = E, y = lab, cohort = cohvec), file.path(data_dir, "train_matrix.rds"))
cat("saved train_matrix.rds:", nrow(E), "genes x", ncol(E), "samples\n")
