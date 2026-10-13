## [FIXED P0-3] 修正说明: setwd -> F:/.../26.GSEA (x1)
## 原始文件: F:/232机器学习联合神经网络/CRC/26.GSEA/neuralML26.GSEA.R (原编码 gbk)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

expFile="merge.normalize.txt"           #表达数据文件
geneFile="interFeatureGenes.txt"       #基因列表文件
gmtFile="c2.cp.kegg.Hs.symbols.gmt"    #基因集文件
setwd("F:/232机器学习联合神经网络/CRC/26.GSEA")      #设置工作目录

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#去除对照组的样品
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,Type=="Treat",drop=F]

#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#对基因进行循环
for(gene in as.vector(geneRT[,1])){
	#根据目标基因表达量对样品进行分组，得到高低表达组的logFC
	dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #低表达组的数据
	dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #高表达组的数据
	meanL=rowMeans(dataL)
	meanH=rowMeans(dataH)
	meanL[meanL<0.00001]=0.00001
	meanH[meanH<0.00001]=0.00001
	logFC=meanH-meanL
	#根据logFC对基因进行排序
	logFC=sort(logFC, decreasing=T)
	genes=names(logFC)

	#读取基因集文件
	gmt=read.gmt(gmtFile)

	#GSEA富集分析
	kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
	kkTab=as.data.frame(kk)
	#根据pvalue<0.05对数据进行过滤,得到显著富集的结果
	kkTab=kkTab[kkTab$pvalue<0.05,]
	write.table(kkTab,file=paste0(gene,".result.txt"),sep="\t",quote=F,row.names = F)

	#绘制高表达组富集的图形
	termNum=5     #设置展示通路的数目
	kkUp=kkTab[kkTab$NES>0,]
	termNum=ifelse(nrow(kkUp)>=termNum, termNum, nrow(kkUp))
	showTerm=row.names(kkUp)[1:termNum]     #展示通路的名称
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in high ", gene, " group"))
	pdf(file=paste0(gene,".highExp.pdf"), width=6.5, height=5.5)
	print(gseaplot)
	dev.off()

	#绘制低表达组富集的图形
	termNum=5     #设置展示通路的数目
	kkDown=kkTab[kkTab$NES<0,]
	termNum=ifelse(nrow(kkDown)>=termNum, termNum, nrow(kkDown))
	showTerm=row.names(kkDown)[1:termNum]     #展示通路的名称
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in low ", gene, " group"))
	pdf(file=paste0(gene,".lowExp.pdf"), width=6.5, height=5.5)
	print(gseaplot)
	dev.off()
}


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


