## [FIXED P0-3] 修正说明: 仅路径规范化
## 原始文件: F:/232机器学习联合神经网络/CRC/16.hubGene/neuralML16.venn.R (原编码 gbk)
#install.packages("ggvenn")


library(ggvenn)      #引用包
scoreFile="score.csv"     #输入文件
setwd("C:\\Users\\Administrator\\Desktop\\neuralML\\16.hubGene")     #设置工作目录

#读取输入文件
rt=read.csv(scoreFile, header=T, sep=",", check.names=F, row.names=1)
colnames(rt)=c(colnames(rt)[2:ncol(rt)], "N")
scoreType=c("Degree", "EPC", "MCC", "MNC")

#根据打分对基因进行排序
geneList=c()
for(i in scoreType){
	data=rt[order(rt[,i],decreasing=T),]
	hubGene=row.names(data)[1:20]
	geneList[[i]]=hubGene
}

#绘制venn图
pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList,show_percentage = T,
	stroke_color = "white", stroke_size = 0.5,
	fill_color = c("#E41A1C","#1E90FF","#4DAF4A","#984EA3"),
	set_name_color = c("#E41A1C","#1E90FF","#4DAF4A","#984EA3"),
	set_name_size=6, text_size=4.5)
dev.off()

#得到四种方法的交集基因, 作为网络的核心基因
interGenes=Reduce(intersect, geneList)
write.table(interGenes, file="hubGenes.txt", sep="\t", quote=F, row.names=F, col.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


