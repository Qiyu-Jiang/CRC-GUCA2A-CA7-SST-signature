## P1-6b CellChat 比较（修复版）：读入已算好的 high/low 对象，cell.prefix=TRUE 合并出图
suppressMessages({library(CellChat); library(ggplot2)})
root   <- "F:/232机器学习联合神经网络/CRC"
outdir <- file.path(root, "99.P1", "07.cellchat")
cat("== P1-6b start ==\n")

ccH <- readRDS(file.path(outdir, "cellchat_high.rds"))
ccL <- readRDS(file.path(outdir, "cellchat_low.rds"))
object.list <- list(EpiSigHigh = ccH, EpiSigLow = ccL)

merged <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
saveRDS(merged, file.path(outdir, "cellchat_merged.rds"))
cat("merged saved\n")

## rankNet：两组整体信号强度对比
pdf(file.path(outdir, "cellchat_rankNet.pdf"), width = 11, height = 7)
gg1 <- rankNet(merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)
print(gg1)
gg1b <- rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
print(gg1b)
dev.off()
cat("rankNet plots done\n")

## netVisual_diff：通讯数量/强度差异
pdf(file.path(outdir, "cellchat_diff.pdf"), width = 9.5, height = 8)
par(mfrow = c(1, 2))
gg2 <- tryCatch(netVisual_diff(merged, weight.scale = TRUE), error = function(e) { cat("netVisual_diff 不可用:", conditionMessage(e), "\n"); NULL })
gg3 <- tryCatch(netVisual_heatmap(merged, weight.scale = TRUE, color.heatmap = "Reds"), error = function(e) { cat("netVisual_heatmap 不可用:", conditionMessage(e), "\n"); NULL })
dev.off()
cat("diff/heatmap done\n")

## bubble：EpiSigHigh / EpiSigLow 作为 source 的比较
pdf(file.path(outdir, "cellchat_bubble.pdf"), width = 13, height = 9)
src <- grep("EpiSig", levels(merged@idents[[1]]), value = TRUE)
tgt <- setdiff(levels(merged@idents[[1]]), src)
gg4 <- netVisual_bubble(merged, sources.use = src, targets.use = tgt,
                        comparison = c(1, 2), angle.x = 45)
print(gg4)
dev.off()
cat("bubble done\n")

## 通路概率汇总表
pw <- rbind(
  data.frame(group = "EpiSigHigh", pathway = rownames(ccH@netP$prob),
             prob = apply(ccH@netP$prob, 1, max)),
  data.frame(group = "EpiSigLow", pathway = rownames(ccL@netP$prob),
             prob = apply(ccL@netP$prob, 1, max)))
write.csv(pw, file.path(outdir, "pathway_prob_summary.csv"), row.names = FALSE)
print(pw[order(-pw$prob), ][1:20, ])
cat("== P1-6b done ==\n")
