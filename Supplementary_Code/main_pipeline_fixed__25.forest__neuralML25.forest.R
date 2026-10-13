## [FIXED P0-3] 修正说明: setwd -> F:/.../25.forest (x1)
## 原始文件: F:/232机器学习联合神经网络/CRC/25.forest/neuralML25.forest.R (原编码 gbk)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("glmnet")


#引用包
library(glmnet)
expFile="merge.normalize.txt"         #表达数据文件
geneFile="interFeatureGenes.txt"     #基因列表文件
setwd("F:/232机器学习联合神经网络/CRC/25.forest")    #设置工作目录

#读取表达数据文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

#读取基因列表文件,提取交集特征基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#获取样品分组信息(对照组和实验组)
data=t(data)
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))
rt=as.data.frame(data)
rt$Type=ifelse(group=="Control", 0, 1)

############定义森林图函数############
bioForest=function(logisticFile=null, forestFile=null, forestCol=null){
	#读取输入文件
	rt <- read.table(logisticFile, header=T, sep="\t", check.names=F, row.names=1)
	gene <- rownames(rt)
	or <- sprintf("%.3f",rt$"OR")
	orLow  <- sprintf("%.3f",rt$"OR.95L")
	orHigh <- sprintf("%.3f",rt$"OR.95H")
	Hazard.ratio <- paste0(or,"(",orLow,"-",orHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#输出图形
	pdf(file=forestFile, width=6.5, height=4.5)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#绘制森林图左边的基因信息
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Odds ratio',cex=text.cex,font=2,adj=1)
		
	#绘制右边的森林图
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(orLow),as.numeric(orHigh),1.1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Odds ratio")
	arrows(as.numeric(orLow),n:1,as.numeric(orHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=3)
	abline(v=1, col="black", lty=2, lwd=2)
	boxcolor = ifelse(as.numeric(or) > 1, forestCol, forestCol)
	points(as.numeric(or), n:1, pch = 15, col = boxcolor, cex=2)
	axis(1)
	dev.off()
}

##############单因素逻辑回归分析##############
outTab=data.frame()
for(i in colnames(rt)[1:(ncol(rt)-1)]){
	rt2=rt[,c(i, "Type")]
	colnames(rt2)=c("Gene", "Type")
	fit=glm(Type ~ Gene, family="binomial", data=rt2)
	summ2=summary(fit)
	conf=confint(fit, level=0.95)
	geneTab=cbind(Gene=i,
	               OR=exp(summ2$coefficients[,"Estimate"])[2],
	               OR.95L=exp(conf[,1])[2],
	               OR.95H=exp(conf[,2])[2],
	               pvalue=summ2$coefficients[,"Pr(>|z|)"][2])
	outTab=rbind(outTab, geneTab)
}
#输出单因素逻辑回归分析结果文件
write.table(outTab, file="uni.logistic.txt", sep="\t", quote=F, row.names=F)
#绘制单因素的森林图
bioForest(logisticFile="uni.logistic.txt", forestFile="uni.forest.pdf", forestCol="blue")


##############多因素逻辑回归分析##############
fit=glm(Type ~ ., family="binomial", data=rt)
summ2=summary(fit)
conf=confint(fit, level=0.95)
multiTab=cbind(
	         OR=exp(summ2$coefficients[,"Estimate"])[-1],
	         OR.95L=exp(conf[,1])[-1],
	         OR.95H=exp(conf[,2])[-1],
	         pvalue=summ2$coefficients[,"Pr(>|z|)"][-1])
#输出多因素逻辑回归分析的结果
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab, file="multi.logistic.txt", sep="\t", quote=F, row.names=F)
#绘制多因素的森林图
bioForest(logisticFile="multi.logistic.txt", forestFile="multi.forest.pdf", forestCol="red")


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


