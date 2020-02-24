setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/manuscript/")

source("ggvenn.R")
source("geom_venn.R")

setwd("/home/chong/Dropbox/FSU_research/R21 Grant Application/MWAS/resAna")
setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/AD_res/")
library(VennDiagram)
cisMWAS = readRDS("summary_best_cis_MWAS_AD_Jansene.rds")


cpgRes = readRDS("res3_cross_CpG_AD_Jansene.rds")
cpgRes = cpgRes[!is.na(cpgRes[,"cross_CMO"]),]
cpgRes = cpgRes[cpgRes[,"cross_CMO"]>=0,]
cpgRes = cpgRes[cpgRes[,"cross_CMO"]<0.05/dim(cpgRes)[1],]
CMO.gene = cpgRes[,2]
length(CMO.gene)

tmp = cisMWAS[cisMWAS[,"Gene"] %in% CMO.gene,]
tmp[,"most_sig_500Kb"] = as.numeric(tmp[,"most_sig_500Kb"])
tmp = tmp[tmp[,"most_sig_500Kb"]>5e-8,]
CMO.gene = tmp[,"Gene"]
length(CMO.gene)

tab = cpgRes[cpgRes[,2] %in% CMO.gene,c(1,2,4,5,29)]
library(xtable)
xtable(tab,digits = -2)
