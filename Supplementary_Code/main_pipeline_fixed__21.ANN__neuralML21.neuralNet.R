## [FIXED P0-3] 修正说明: setwd -> F:/.../21.ANN (x1)
## 原始文件: F:/232机器学习联合神经网络/CRC/21.ANN/neuralML21.neuralNet.R (原编码 utf-8-sig)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("neuralnet")
#install.packages("NeuralNetTools")
#install.packages("pROC")


#引用包
library(limma)
library(neuralnet)
library(NeuralNetTools)
library(pROC)
set.seed(123)

expFile="merge.normalize.txt"           #表达数据文件
geneFile="interFeatureGenes.txt"       #交集特征基因的文件
setwd("F:/232机器学习联合神经网络/CRC/21.ANN")      #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件,提取交集特征基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),,drop=F]
data=as.data.frame(t(data))

#获取样品分组信息(对照组和实验组)
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))
data$Control=ifelse(group=="Control", 1, 0)
data$Treat=ifelse(group=="Treat", 1, 0)

#构建人工神经网络的模型
fit=neuralnet(Control+Treat~., data, hidden=5, stepmax=1000000)
fit$result.matrix
fit$weight
plot(fit)
dev.copy2pdf(file="neuralnet1.pdf", width=9, height=7, out.type="pdf")

pdf(file="neuralnet2.pdf", width=9, height=7)
plotnet(fit)
dev.off()

#利用模型预测结果
net.predict=compute(fit, data)$net.result
net.prediction=c("Control", "Treat")[apply(net.predict, 1, which.max)]
predict.table=table(group, net.prediction)
predict.table
ControlAccuracy=predict.table[1,1]/(predict.table[1,1]+predict.table[1,2])
TreatAccuracy=predict.table[2,2]/(predict.table[2,1]+predict.table[2,2])
paste0("Control accuracy: ", sprintf("%.3f", ControlAccuracy))
paste0("Treat accuracy: ", sprintf("%.3f", TreatAccuracy))

#准备ROC曲线的数据
colnames(net.predict)=c("Control", "Treat")
y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(net.predict))
y=ifelse(y=="Control", 0, 1)

#绘制ROC曲线
roc1=roc(y, as.numeric(net.predict[,2]))
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T)
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


