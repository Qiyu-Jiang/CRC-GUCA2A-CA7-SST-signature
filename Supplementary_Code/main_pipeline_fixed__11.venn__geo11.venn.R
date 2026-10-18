## [FIXED P0-3] 修正说明: setwd -> (project)/11.venn (x1)
## 原始文件: .//11.venn/geo11.venn.R (原编码 gbk)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("ggvenn")


#引用包
library(ggvenn)

diffFile="diff.txt"                  #差异分析的结果文件
moduleFile="module_green.txt"      #模块基因的文件
setwd(".//11.venn")     #设置工作目录
geneList=list()

#读取差异分析的结果文件
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])               #提取差异基因的名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)                #对差异基因取unique
geneList[["DEG"]]=uniqGene

#读取模块基因的文件
rt=read.table(moduleFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])                #提取模块基因的名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)                 #对模块基因取unique
geneList[["WGCNA"]]=uniqGene

#绘制venn图
pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList,show_percentage = T,
	stroke_color = "white", stroke_size = 0.5,
	fill_color = c("#E41A1C","#1E90FF"),
	set_name_color =c("#E41A1C","#1E90FF"),
	set_name_size=6, text_size=4.5)
dev.off()

#输出交集基因
interGenes=Reduce(intersect, geneList)
write.table(file="interGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


