setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/manuscript/")

source("ggvenn.R")
source("geom_venn.R")

setwd("/home/chong/Dropbox/FSU_research/R21 Grant Application/MWAS/resAna")
setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/AD_res/")
library(VennDiagram)

cisTWAS = readRDS("summary_best_cis_AD_Jansene.rds")
cisMWAS = readRDS("summary_best_cis_MWAS_AD_Jansene.rds")

crosRes = readRDS("res_cross_tissue_AD_Jansene_1000G_SPrediXcan.rds")
single = readRDS("res_signle_tissue_AD_Jansene_1000G_SPrediXcan.rds")

dim(single)

crosRes = crosRes[!is.na(crosRes[, "MultiXcan"]),]
multixcan = crosRes[crosRes[, "MultiXcan"] < 0.05/dim(crosRes)[1],]

#multixcan = multixcan[multixcan[, "best_p-value"] < 1e-4,]
dim(multixcan)
multixcan.gene = multixcan[, "Gene"]
multixcan.gene = as.character(multixcan.gene)

tmp = cisTWAS[cisTWAS[,"Gene"] %in% multixcan.gene,]
tmp[,"most_sig_500Kb"] = as.numeric(tmp[,"most_sig_500Kb"])
tmp = tmp[tmp[,"most_sig_500Kb"]>5e-8,]
multixcan.gene = tmp[,"Gene"]

single = single[!is.na(single[,"p-value"]),]
twas = single[single[, "p-value"] < 0.05/dim(single)[1],]

twas.gene = unique(twas[, "Gene_name"])
twas.gene = as.character(twas.gene)

tmp = cisTWAS[cisTWAS[,"Gene"] %in% twas.gene,]
tmp[,"most_sig_500Kb"] = as.numeric(tmp[,"most_sig_500Kb"])
tmp = tmp[tmp[,"most_sig_500Kb"]>5e-8,]
twas.gene = tmp[,"Gene"]
length(twas.gene)


cpgRes = readRDS("res3_cross_CpG_AD_Jansene.rds")
cpgRes = cpgRes[!is.na(cpgRes[,"cross_CMO"]),]
cpgRes = cpgRes[cpgRes[,"cross_CMO"]>=0,]
cpgRes = cpgRes[cpgRes[,"cross_CMO"]<0.05/dim(cpgRes)[1],]
CMO.gene = cpgRes[,2]

tmp = cisMWAS[cisMWAS[,"Gene"] %in% CMO.gene,]
tmp[,"most_sig_500Kb"] = as.numeric(tmp[,"most_sig_500Kb"])
tmp = tmp[tmp[,"most_sig_500Kb"]>5e-8,]
CMO.gene = tmp[,"Gene"]
length(CMO.gene)

CpG = readRDS("res3_signle_CpG_AD_Jansene.rds")
CpG = CpG[CpG[,"Feature"] %in% c("GeneBody","Promoter"),]
n.CpG = length(CpG[,"CpG"])
CpG = CpG[!is.na(CpG[,"SUM"]),]
CpG = CpG[CpG[,"SUM"]<0.05/n.CpG,]
MWAS.gene = unique(CpG[,"geneID"])

tmp = cisMWAS[cisMWAS[,"Gene"] %in% MWAS.gene,]
tmp[,"most_sig_500Kb"] = as.numeric(tmp[,"most_sig_500Kb"])
tmp = tmp[tmp[,"most_sig_500Kb"]>5e-8,]
MWAS.gene = tmp[,"Gene"]

library(ggvenn)
library(ggplot2)
library("ggsci")
library("ggplot2")
library("gridExtra")

# use list as input
a <- list(`MWAS` = MWAS.gene,
`MultiXcan` = multixcan.gene,
`TWAS` = twas.gene,
`CMO` = CMO.gene)
ggvenn(a,set_name_size = 7,text_size = 7)+scale_fill_npg() + scale_color_npg()

outd = "/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/"
fig.name = "Sfig2_AD_Jansene_venn_all_novel"
ggsave(paste(outd,fig.name, ".pdf", sep = ""))

