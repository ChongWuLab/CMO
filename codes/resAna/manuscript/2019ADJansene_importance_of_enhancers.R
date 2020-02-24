setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/manuscript/")

source("ggvenn.R")
source("geom_venn.R")

setwd("/home/chong/Dropbox/FSU_research/R21 Grant Application/MWAS/resAna")
setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/AD_res/")
library(VennDiagram)

cpgRes = readRDS("res3_cross_CpG_AD_Jansene.rds")
cpgRes = cpgRes[!is.na(cpgRes[,"cross_CMO"]),]
cpgRes = cpgRes[cpgRes[,"cross_CMO"]>=0,]
CMO = cpgRes[cpgRes[,"cross_CMO"]<0.05/dim(cpgRes)[1],]
CMO.gene = CMO[,2]

write.table(CMO[CMO[,"CHR"]==19,2],"CMO.gene.list.txt",sep=",",quote=F,row.names=F,col.names=F)


gene = cpgRes[cpgRes[,"gene_CMO"]>=0,]
gene = gene[gene[,"gene_CMO"]<0.05/dim(gene)[1],]
gene = gene[,2]
gene.list = CMO.gene[!CMO.gene%in%gene ]

tmp = cpgRes[cpgRes[,2] %in%gene.list, ]

tmp[,c(1:4,29:31)]
