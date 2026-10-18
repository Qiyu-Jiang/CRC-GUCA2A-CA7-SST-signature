## [FIXED P0-3] 修正说明: setwd -> (project)/38.exposure (x1)
## 原始文件: .//38.exposure/neuralML38.URL.R (原编码 gbk)
refFile="IEU.ref.txt"        #数据库的参考文件
geneFile="interFeatureGenes.txt"     #交集特征基因的文件
setwd(".//38.exposure")     #设置工作目录

#读取数据库的参考文件
refRT=read.csv(refFile, header=T, sep="\t", check.names=F)

#读取交集特征基因的文件
geneRT=read.csv(geneFile, header=F, sep="\t", check.names=F)

#提取交集特征基因的数据
outTab=refRT[refRT$Symbol %in% as.vector(geneRT[,1]),]

#在表格中添加暴露数据的下载链接
url=paste0("https://gwas.mrcieu.ac.uk/files/", outTab$id, "/", outTab$id, ".vcf.gz")
outTab=cbind(outTab, url)
write.csv(outTab, file="URL.csv", row.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


