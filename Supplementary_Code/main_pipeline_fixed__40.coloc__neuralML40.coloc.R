## [FIXED P0-3] 修正说明: setwd -> F:/.../40.coloc (x1)
## 原始文件: F:/232机器学习联合神经网络/CRC/40.coloc/neuralML40.coloc.R (原编码 utf-8-sig)
#install.packages("dplyr")
#install.packages("data.table")
#install.packages("coloc")

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")

#install.packages("devtools")
#devtools::install_github("mrcieu/ieugwasr",force = TRUE)
#devtools::install_github("mrcieu/gwasglue", force = TRUE)
#devtools::install_github("mrcieu/gwasvcf", force = TRUE)

#library(devtools)
#devtools::install_github("boxiangliu/locuscomparer")


#引用包
library(dplyr)
library(data.table)
library(coloc)
library(VariantAnnotation)
library(gwasglue)
library(locuscomparer)

#输入文件
drugFile="IEU.ref.txt"      #IEU数据库的参考文件
outcomeFile="finngen_R12_AD_EO_EXMORE.gz"     #结局数据文件
finnAnnFile="finngen_R12_manifest.tsv"         #芬兰数据库的注释表格文件

#设置工作目录
setwd("F:/232机器学习联合神经网络/CRC/40.coloc")

#读取IEU数据库的参考文件
drugRT=read.table(drugFile, header=T, sep="\t", check.names=F)
row.names(drugRT)=drugRT[,"id"]

#获取目录下所有vcf.gz结尾的文件(基因的数据)
allVcf=list.files(pattern="*.vcf.gz$")

#读取结局数据文件
data0=data.table::fread(outcomeFile, header=T, sep="\t", check.names=F)
finn_info <- fread(finnAnnFile, data.table=F)
trait_row <- finn_info[grepl(outcomeFile, finn_info$path_https),]
data0$ncase.outcome <- trait_row$num_cases
data0$ncontrol.outcome <- trait_row$num_controls
data0$samplesize.outcome <- trait_row$num_cases + trait_row$num_controls
	
#提取结局数据
data1 <- data0 %>% dplyr::select("rsids","#chrom","pos","alt","ref",
	                                 "af_alt","beta","sebeta","pval","samplesize.outcome","ncase.outcome")
colnames(data1) <- c('SNP','chrom',"pos",'effect_allele','other_allele',
	                     "eaf","beta","se","P","samplesize","number_cases")
data2 <- as.data.frame(data1)
data2$varbeta <- data2$se^2
data2$MAF <- ifelse(data2$eaf<0.5,data2$eaf,1-data2$eaf)
data3 <- subset(data2, !duplicated(SNP))
data3$s <- data3$number_cases/data3$samplesize
data3$z = data3$beta/data3$se
GWASdataAll <- data3 %>% na.omit()


#对vcf文件循环(基因的数据)
outTab=data.frame()
for(eqtlFile in allVcf){
	possibleError=tryCatch({
	#获取基因位置信息
	eqtlID=gsub(".vcf.gz", "", eqtlFile)
	geneChr=drugRT[eqtlID, "Chr"]
	geneStart=drugRT[eqtlID, "Start"]
	geneEnd=drugRT[eqtlID, "End"]
	geneName=drugRT[eqtlID, "Symbol"]
	
	#读取基因的数据文件
	vcfRT <- readVcf(eqtlFile)
	data1=gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="exposure")
	
	#挑选暴露数据
	data2 <- data1 %>% dplyr::select("SNP","chr.exposure","pos.exposure","effect_allele.exposure",
	                                 "other_allele.exposure", "eaf.exposure","beta.exposure",
	                                 "se.exposure", "pval.exposure", "samplesize.exposure")
	colnames(data2) <- c("SNP","chrom","Pos","effect_allele","other_allele",
	                     "MAF","beta","se","P","samplesize")
	#整理共定位所需数据
	data3 <- as.data.frame(data2)
	data3$varbeta <- data3$se^2
	data3$z = data3$beta/data3$se
	data3 <- subset(data3, !duplicated(SNP))
	geneData <- data3 %>% filter(chrom==geneChr, Pos>geneStart-100000, Pos<geneEnd+100000) %>% na.omit()
	lead <- geneData %>% dplyr::arrange(P)
	leadPos <- lead$Pos[1]
	QTLdata <-  geneData %>% filter(Pos>leadPos-50000,Pos<leadPos+50000) %>% na.omit()

	#共定位分析	
	sameSNP <- intersect(QTLdata$SNP, GWASdataAll$SNP)
	QTLdata <- QTLdata[QTLdata$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
	GWASdata <- GWASdataAll[GWASdataAll$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
	coloc_data <- list(dataset1=list(snp=QTLdata$SNP,beta=QTLdata$beta,varbeta=QTLdata$varbeta,
	                                 N=QTLdata$samplesize,MAF=QTLdata$MAF,z = QTLdata$z,
	                                 pvalues=QTLdata$P,type="quant"), 
	                   dataset2=list(snp=GWASdata$SNP,beta=GWASdata$beta,varbeta=GWASdata$varbeta,
	                                 N=GWASdata$samplesize,MAF=GWASdata$MAF,z = GWASdata$z,
	                                 pvalues=GWASdata$P,type="cc"))
	result <- coloc.abf(dataset1=coloc_data$dataset1, dataset2=coloc_data$dataset2)
	#保存共定位分析的结果
	result$results %>% dplyr::arrange(desc(SNP.PP.H4))
	outTab=rbind(outTab, data.frame(ID=eqtlID,
								   Symbol=geneName,
		                                PP.H0=result$summary[2],
		                                PP.H1=result$summary[3],
		                                PP.H2=result$summary[4],
		                                PP.H3=result$summary[5],
		                                PP.H4=result$summary[6]))
		
	#绘制图形
	GWAS_fn <- GWASdata[,c('SNP','P')] %>% dplyr::rename(rsid = SNP, pval = P)
	eQTL_fn <- QTLdata[,c('SNP','P')] %>% dplyr::rename(rsid = SNP, pval = P)
	pdf(file=paste0(eqtlID, "." , geneName, ".pdf"), width=8, height=6)
	print(locuscompare(in_fn1 = GWAS_fn, 
		                   in_fn2 = eQTL_fn,
		                   snp = NULL,
		                   title1 = 'GWAS', 
		                   title2 = 'eQTL'))
	dev.off()
	},error=function(e) e)
    if(inherits(possibleError, "error")){next}
}
#输出共定位分析的结果
write.csv(outTab, file="coloc.result.csv", row.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056


