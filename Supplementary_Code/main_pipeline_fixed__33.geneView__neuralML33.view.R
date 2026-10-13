## [FIXED P0-3] 修正说明: 仅路径规范化
## 原始文件: F:/232机器学习联合神经网络/CRC/33.geneView/neuralML33.view.R (原编码 gbk)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("clustree")
#install.packages("Seurat")
#install.packages("harmony")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("monocle")

#install.packages("devtools")
#devtools::install_github('immunogenomics/presto')

#install.packages("assertthat")
#install.packages("SCpubr")



#######################01.数据前期处理和矫正#######################
#引用包
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(clustree)
library(harmony)
library(assertthat)
library(SCpubr)

#设置工作目录
workDir="C:\\biowolf\\neuralML\\33.geneView"
setwd(workDir)
load("Seurat.Rdata")

#读取交集特征基因的文件
geneRT=read.table("interFeatureGenes.txt", header=F, sep="\t", check.names=F)
hubGenes=as.vector(geneRT[,1])
hubGenes=hubGenes[hubGenes %in% rownames(pbmc)]

#绘制散点图
pdf(file="Scatter.pdf", width=10, height=8)
FeaturePlot(object = pbmc, features = hubGenes, cols = c("green", "red"), ncol=3)
dev.off()

#绘制小提琴图
pdf(file="Violin.pdf", width=15, height=10)
VlnPlot(object = pbmc, features = hubGenes, ncol=4)
dev.off()

#绘制气泡图
pdf(file="Bubble.pdf", width=8, height=5)
DotPlot(object = pbmc, features = hubGenes) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
dev.off()

# 对细胞进行循环,观察特征基因在哪些细胞中是具有差异的
pbmc$group2=pbmc@active.ident
for(cellName in unique(pbmc@active.ident)){
	target_cells <- subset(pbmc, subset = group2 == cellName)
	pdf(file=paste0("cell_",cellName,".Violin.pdf"), width=9, height=6)
	p1=VlnPlot(object = target_cells, features = hubGenes, group.by = "Type", ncol=4)
	print(p1)
	dev.off()
}


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


