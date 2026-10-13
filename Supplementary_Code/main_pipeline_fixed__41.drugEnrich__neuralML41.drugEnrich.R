## [FIXED P0-3] 修正说明: setwd -> F:/.../41.drugEnrich (x1)
## 原始文件: F:/232机器学习联合神经网络/CRC/41.drugEnrich/neuralML41.drugEnrich.R (原编码 gbk)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("tidytable")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("tidytable")
library("ggplot2")

pvalueFilter=0.05      #p值过滤条件
adjPvalFilter=0.05     #矫正后p值的过滤条件

#定义图形的颜色
colorSel="p.adjust"
if(adjPvalFilter>0.05){
	colorSel="pvalue"
}

drugFile="DSigDB_All_detailed.txt"      #药物和基因的关系文件
hubFile="interFeatureGenes.txt"          #交集特征基因的文件
setwd("F:/232机器学习联合神经网络/CRC/41.drugEnrich")      #设置工作目录

#读取交集特征基因的列表文件
rt=read.table(hubFile, header=F, sep="\t", check.names=F)
#提取基因的名称
genes=unique(as.vector(rt[,1]))

#读取药物和基因的关系文件
drugRT=read.table(drugFile, header=T, sep="\t", check.names=F, quote="", comment.char="")
drugRT=drugRT[,1:2]

#药物富集分析
kk=enricher(genes,
	pvalueCutoff=1, qvalueCutoff=1,
	minGSSize = 10, maxGSSize = 500,
	TERM2GENE=drugRT)

#输出显著富集的结果
DRUG=as.data.frame(kk)
DRUG=DRUG[(DRUG$pvalue<pvalueFilter & DRUG$p.adjust<adjPvalFilter),]
write.table(DRUG[,-c(3,4)], file="DRUG.enrich.xls", sep="\t", quote=F, row.names = F)

#设置展示药物的数目
showNum=30
if(nrow(DRUG)<showNum){
	showNum=nrow(DRUG)
}

#柱状图
pdf(file="barplot.pdf", width=7, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=7, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


