## [FIXED P0-3] 修正说明: setwd -> F:/.../20.venn (x1)
## 原始文件: F:/232机器学习联合神经网络/CRC/20.venn/neuralML20.venn.R (原编码 gbk)
######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("ggvenn")


#引用包
library(ggvenn)

setwd("F:/232机器学习联合神经网络/CRC/20.venn")    #设置工作目录
files=list.files(pattern="*.txt$")     #获取目录下txt结尾的文件
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(inputFile in files){
	if(inputFile=="interFeatureGenes.txt"){next}
    rt=read.table(inputFile, header=F, sep="\t", check.names=F)      #读取输入文件
    geneNames=unlist(strsplit(as.vector(rt[,1]), " "))       #提取基因名称
    geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
    uniqGene=unique(geneNames)                 #基因取unique
    header=unlist(strsplit(inputFile,"\\.|\\-"))
    geneList[[header[1]]]=uniqGene
}
#绘制venn图
pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList, show_percentage = T,
  stroke_color = "white", stroke_size = 0.5,
  fill_color = c("#E41A1C","#1E90FF","#FF8C00"),
  set_name_color =c("#E41A1C","#1E90FF","#FF8C00"),
  set_name_size=6, text_size=4.5)
dev.off()

#输出三种方法的交集特征基因
interGenes=Reduce(intersect, geneList)
write.table(file="interFeatureGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


