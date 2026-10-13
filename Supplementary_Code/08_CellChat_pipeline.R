## ============================================================
## P1-6  CellChat 细胞通讯（GSE289314 结癌 TME，59,632 细胞）
## 设计：上皮细胞按 GUCA2A+CA7+SST 联合评分分为 Epi_SigHigh / Epi_SigLow，
##       其余 9 种细胞为共享背景；每组下采样 ≤2000 细胞；
##       跑 high / low 两个 CellChat 对象并 merge 比较（rankNet + diff + bubble）
## 输出: 99.P1/07.cellchat/
## ============================================================
suppressMessages({library(Seurat); library(CellChat); library(future); library(ggplot2)})
root   <- "F:/232机器学习联合神经网络/CRC"
outdir <- file.path(root, "99.P1", "07.cellchat")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sig <- c("GUCA2A", "CA7", "SST")
plan("multisession", workers = 4)
options(future.globals.maxSize = 4 * 1024^3)
cat("== P1-6 start ==\n")

env <- new.env(parent = emptyenv())
load(file.path(root, "32.Seurat", "Seurat.Rdata"), envir = env)
so <- get("pbmc", envir = env)
DefaultAssay(so) <- "RNA"
cat("loaded:", ncol(so), "cells\n")

## 确保有 normalized data
dl <- tryCatch(nrow(LayerData(so, assay = "RNA", layer = "data")), error = function(e) 0)
if (is.null(dl) || dl == 0) so <- NormalizeData(so, verbose = FALSE)

## 上皮细胞联合评分 + 中位数分组
cat("Type values:", paste(unique(as.character(so$Type)), collapse = " | "), "\n"); epi_cells <- colnames(so)[as.character(Idents(so)) == "Epithelial cells"]; cat("epi cells found:", length(epi_cells), "\n"); epi <- subset(so, cells = epi_cells)
epi <- AddModuleScore(epi, features = list(sig), name = "EpiSig", seed = 123)
med <- median(epi$EpiSig1)
epi$CellGroup <- ifelse(epi$EpiSig1 >= med, "EpiSigHigh", "EpiSigLow")
cat("epithelial high/low:", sum(epi$CellGroup == "EpiSigHigh"), "/", sum(epi$CellGroup == "EpiSigLow"), "\n")

so$CellGroup <- as.character(Idents(so))
so$CellGroup[match(colnames(epi), colnames(so))] <- epi$CellGroup
so$CellGroup <- gsub("\\+", "", so$CellGroup)   # 避免正则元字符
print(table(so$CellGroup))

## 每组下采样 ≤2000
set.seed(123)
cells.use <- unlist(lapply(split(colnames(so), so$CellGroup),
                           function(x) if (length(x) > 2000) sample(x, 2000) else x))
sub <- subset(so, cells = cells.use)
cat("downsampled object:", ncol(sub), "cells,", length(unique(sub$CellGroup)), "groups\n")

runCC <- function(keepGroup, tag){
  keep <- colnames(sub)[sub$CellGroup %in% keepGroup]
  o <- subset(sub, cells = keep)
  o$CellGroup <- factor(o$CellGroup)
  cc <- createCellChat(object = o, group.by = "CellGroup", assay = "RNA")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  saveRDS(cc, file.path(outdir, paste0("cellchat_", tag, ".rds")))
  cat(tag, "done: pathways", nrow(cc@netP$prob), "\n")
  cc
}
others <- setdiff(unique(sub$CellGroup), c("EpiSigHigh", "EpiSigLow"))
ccH <- runCC(c("EpiSigHigh", others), "high")
ccL <- runCC(c("EpiSigLow",  others), "low")

object.list <- list(EpiSigHigh = ccH, EpiSigLow = ccL)
merged <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(merged, file.path(outdir, "cellchat_merged.rds"))

pdf(file.path(outdir, "cellchat_comparison.pdf"), width = 11, height = 7)
gg1 <- rankNet(merged, mode = "comparison", stacked = FALSE, do.stat = TRUE, return.data = FALSE)
print(gg1)
gg2 <- netVisual_diff(merged, weight.scale = TRUE)
print(gg2)
dev.off()

pdf(file.path(outdir, "cellchat_bubble.pdf"), width = 12, height = 8)
evidx <- which(names(object.list) == "EpiSigHigh")
gg3 <- netVisual_bubble(merged, sources.use = grep("EpiSig", levels(merged@idents[[1]])),
                        targets.use = grep("EpiSig", levels(merged@idents[[1]]), invert = TRUE),
                        comparison = c(1, 2), angle.x = 45)
print(gg3)
dev.off()

## 关键通路汇总表
pw <- rbind(data.frame(group = "EpiSigHigh", pathway = rownames(ccH@netP$prob),
                       prob = apply(ccH@netP$prob, 1, max)),
            data.frame(group = "EpiSigLow", pathway = rownames(ccL@netP$prob),
                       prob = apply(ccL@netP$prob, 1, max)))
write.csv(pw, file.path(outdir, "pathway_prob_summary.csv"), row.names = FALSE)
cat("== P1-6 done ==\n")
