## [FIXED P0-3] 修正说明: setwd -> (project)/27.CIBERSORT (x1)
## 原始文件: .//27.CIBERSORT/neuralML27.run.R (原编码 gbk)
#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")


inputFile="merge.normalize.txt"      #表达数据文件
setwd(".//27.CIBERSORT")      #设置工作目录
source("neuralML27.CIBERSORT.R")       #引用包

#免疫细胞浸润分析
outTab=CIBERSORT("ref.txt", inputFile, perm=1000)

#对免疫浸润结果过滤，并且保存免疫细胞浸润结果
outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


