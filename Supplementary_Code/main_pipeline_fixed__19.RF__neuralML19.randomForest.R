## [FIXED P0-3] 修正说明: setwd -> F:/.../19.RF (x1); enabled set.seed(123) (x1)
## 原始文件: F:/232机器学习联合神经网络/CRC/19.RF/neuralML19.randomForest.R (原编码 utf-8-sig)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("randomForest")
#install.packages("ggplot2")


#引用包
set.seed(123)
library(limma)
library(randomForest)
library(ggplot2)

expFile="merge.normalize.txt"      #表达数据文件
geneFile="hubGenes.txt"             #核心基因的文件
setwd("F:/232机器学习联合神经网络/CRC/19.RF")      #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件,提取核心基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),,drop=F]
data=t(data)
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))

#随机森林树
rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=5, height=5)
plot(rf, main="Random forest", lwd=2)
dev.off()

#找出误差最小的点
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

#获取基因重要性的评分
importance=importance(x=rf2)

#挑选疾病的特征基因
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>5])       #根据重要性评分大于5对基因进行过滤
write.table(rfGenes, file="RF.gene.txt", sep="\t", quote=F, col.names=F, row.names=F)

##绘制基因重要性的气泡图
importance=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),,drop=F]
rt=cbind(row.names(importance), as.data.frame(importance))
colnames(rt)=c("ID", "Importance")
if(nrow(rt)>30){rt=rt[1:30,]}
rt$ID=factor(rt$ID, levels=rev(rt$ID))
rt$Importance=as.numeric(rt$Importance)
p=ggplot(rt, aes(Importance, ID)) + 
    geom_segment(aes(x=0, xend=Importance, y=ID, yend=ID), color="grey", cex=1.5) + 
    geom_point(aes(color=Importance), cex=3.5)+
    scale_colour_gradient(low="#5050FFFF", high="#BB0021FF") + 
    labs(x="Importance", y="")+
    theme(axis.text.x=element_text(color="black", size=10),
          axis.text.y=element_text(color="black", size=10))+
    theme_bw()

#输出图形
pdf(file="geneImportance.pdf", width=5, height=4.5)
print(p)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


