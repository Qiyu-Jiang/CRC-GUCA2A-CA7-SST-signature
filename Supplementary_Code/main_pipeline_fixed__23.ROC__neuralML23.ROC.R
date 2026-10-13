## [FIXED P0-3] 修正说明: setwd -> F:/.../23.ROC (x1)
## 原始文件: F:/232机器学习联合神经网络/CRC/23.ROC/neuralML23.ROC.R (原编码 gbk)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("pROC")


#引用包
library(pROC)

expFile="merge.normalize.txt"          #表达数据文件
geneFile="interFeatureGenes.txt"      #交集特征基因的文件
setwd("F:/232机器学习联合神经网络/CRC/23.ROC")    #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

#提取样品的分组信息(对照组和实验组)
y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(rt))
y=ifelse(y=="Control", 0, 1)

#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#定义图形的颜色
bioCol=c("#F05C3BFF","#5C88DAFF","#5CB85CFF", "#EEA236FF", "#9632B8FF", "#17BECFFF", "#BCBD22FF")
bioCol=bioCol[1:nrow(geneRT)]

#对基因进行循环，绘制ROC曲线
aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
	k=k+1
	#绘制ROC曲线
	roc1=roc(y, as.numeric(rt[x,]))     #得到ROC曲线的参数
	if(k==1){
		pdf(file="ROC.genes.pdf", width=5.5, height=5)
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}else{
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE, lwd=3)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText, lwd=3, bty="n", cex=1, col=bioCol[1:(ncol(rt)-1)])
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


