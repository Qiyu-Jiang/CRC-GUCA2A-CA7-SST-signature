## [FIXED P0-3] 修正说明: setwd -> (project)/24.testDiff (x1); expFile GSE36980 -> GSE110224
## 原始文件: .//24.testDiff/neuralML24.testDiff.R (原编码 gbk)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("pROC")


#引用包
library(limma)
library(reshape2)
library(ggpubr)
library(pROC)

expFile = "GSE110224.normalize.txt"      #表达数据文件(需要修改)
geneFile="interFeatureGenes.txt"      #交集特征基因的文件
setwd(".//24.testDiff")     #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件, 提取交集特征基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
data=t(data[sameGene,])

#获取样品的分组信息
Type=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type)

#将数据转换为箱线图的输入文件
data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color = "Type",
	     xlab="", ylab="Gene expression", legend.title="Type",
	     palette = c("blue", "red"), add="jitter", width=0.8)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")
#输出图形
pdf(file="boxplot.pdf", width=6, height=4.5)
print(p1)
dev.off()


#定义图形的颜色
bioCol=c("#F05C3BFF","#5C88DAFF","#5CB85CFF", "#EEA236FF", "#9632B8FF", "#17BECFFF", "#BCBD22FF")
bioCol=bioCol[1:nrow(geneRT)]
#对交集特征基因进行循环，绘制ROC曲线
aucText=c()
rocSigGenes=c()
k=0
y=ifelse(rt$Type=="Control", 0, 1)
for(x in colnames(rt)[1:(ncol(rt)-1)]){
	k=k+1
	#绘制ROC曲线
	roc1=roc(y, as.numeric(rt[,x]))     #得到ROC曲线的参数
	if(k==1){
		pdf(file="ROC.genes.pdf", width=5.5, height=5)
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}else{
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3, add=TRUE)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}
}
#绘制图例, 展示ROC曲线下的面积
legend("bottomright", aucText, lwd=3, bty="n", col=bioCol[1:(ncol(rt)-1)])
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


