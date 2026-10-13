## [FIXED P0-3] 修正说明: setwd -> F:/.../17.lasso (x1); enabled set.seed(123) (x1)
## 原始文件: F:/232机器学习联合神经网络/CRC/17.lasso/neuralML17.lasso.R (原编码 gbk)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("glmnet")
#install.packages("reshape2")
#install.packages("ggplot2")
#install.packages("ggrepel")


#引用包
set.seed(123)
library(glmnet)
library(reshape2)
library(ggplot2)
library(ggrepel)

expFile="merge.normalize.txt"      #表达数据文件
geneFile="hubGenes.txt"             #核心基因的文件
setwd("F:/232机器学习联合神经网络/CRC/17.lasso")      #设置工作目录

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
rt=t(data)

#构建Lasso回归模型
x=as.matrix(rt)
y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)

#绘制交叉验证的图形
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf", width=6, height=5.5)
plot(cvfit,las=1)
dev.off()

#输出lasso回归找到的特征基因
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=as.data.frame(lassoGene[-1])
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)

###############绘制lasso回归的图形###############
#定义图形颜色
colors <- c('#F39B7F', '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85',
            '#3B4992', '#008B45', '#631879', '#008280', '#5F559B', '#A20056', '#BC3C29', '#0072B5', '#E18727', '#AD494A',
            '#20854E', '#7876B1', '#6F99AD', '#FFDC91', '#EE4C97', '#00468B', '#42B540', '#0099B4', '#925E9F', '#FDAF91', '#D6616B',
            '#AD002A', '#ADB6B6', '#374E55', '#DF8F44', '#00A1D5', '#B24745', '#79AF97', '#6A6599', '#80796B', '#EFC000', '#E7969C',
            '#868686', '#CD534C', '#7AA6DC', '#003C67', '#8F7700', '#3B3B3B', '#A73030', '#4A6990', '#FF9900', '#FFCC00', '#756BB1')

#准备图形的数据
melt.data = melt(data.frame(x=log(fit$lambda), t(apply(fit$beta,2,c))),id="x")
textdata = lapply(as.character(unique(melt.data$variable)),function(x){
  genedatai = melt.data[melt.data$variable==x,,drop=F];
  if(sum(genedatai$value)<0){
    w = which.min(genedatai$value)
  }else{
    w = which.max(genedatai$value)
  }
  return(data.frame(x=genedatai$x[w],variable=x,value=genedatai$value[w]))
})
textdata = do.call('rbind',textdata)

#绘制图形
p = ggplot(melt.data,aes(x=x,y=value,group=variable,color=variable))+
  geom_vline(xintercept = log(cvfit$lambda.min),lty = 2, lwd = 1,colour = "red")+
  geom_vline(xintercept = log(cvfit$lambda.1se),lty = 2, lwd = 1,colour = "red")+
  geom_line(size=1.2)+     #线条宽度
  theme_bw()+     #背景和网格
  xlab('Log lambda')+      #x轴的名称
  ylab('Coefficients')+    #y轴的名称
  theme(legend.position = "none")+    #图例的位置
  geom_point(aes(x=-10,y=value),data=textdata,size=2)+
  geom_text_repel(aes(x=x,y=value,label=variable),data=textdata)+     #标注基因的名称
  scale_color_manual(values=colors)      #颜色的设置
 
#输出lasso回归图形
pdf(file="lasso.pdf", width=6, height=5)
print(p)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


