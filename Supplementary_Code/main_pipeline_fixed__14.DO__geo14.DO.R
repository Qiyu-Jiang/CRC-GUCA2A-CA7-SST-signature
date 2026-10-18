## [FIXED P0-3] 修正说明: setwd -> (project)/14.DO (x1)
## 原始文件: .//14.DO/geo14.DO.R (原编码 gbk)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("GSEABase")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("GSEABase")
library("DOSE")

pvalueFilter=0.05      #p值过滤条件
qvalueFilter=1          #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
	
setwd(".//14.DO")      #设置工作目录
rt=read.table("interGenes.txt", header=F, sep="\t", check.names=F)    #读取交集基因的文件

#基因名字转换为基因id
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#DO富集分析(疾病)
kk=enrichDO(gene=gene, pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
DO=as.data.frame(kk)
DO=DO[(DO$pvalue<pvalueFilter & DO$qvalue<qvalueFilter),]
#保存富集分析的结果
write.table(DO, file="DO.txt", sep="\t", quote=F, row.names = F)

#定义显示疾病数目
showNum=30
if(nrow(DO)<showNum){
	showNum=nrow(DO)
}

#柱状图
pdf(file="barplot.pdf", width = 7.5, height = 7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width = 7.5, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


