setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/manuscript/")

source("ggvenn.R")
source("geom_venn.R")

setwd("/home/chong/Dropbox/FSU_research/R21 Grant Application/MWAS/resAna")
setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/AD_res/")
library(VennDiagram)

multiXcan = readRDS("res_cross_tissue_AD_Jansene_1000G_SPrediXcan.rds")
multiXcan = multiXcan[,c(1:3,15)]
multiXcan = multiXcan[!is.na(multiXcan[,4]),]

gene.list = multiXcan[,3]

TWAS = readRDS("res_signle_tissue_AD_Jansene_1000G_SPrediXcan.rds")

gene.list = gene.list[gene.list %in% TWAS[,"Gene_name"]]

CTO = readRDS("res3_cross_CpG_AD_Jansene.rds")
CTO = CTO[,c(1:3,29)]
CTO = CTO[!is.na(CTO[,4]),]

gene.list = gene.list[gene.list %in% CTO[,"geneID"]]
CpG = readRDS("res3_signle_CpG_AD_Jansene.rds")

gene.list = gene.list[gene.list %in% CpG[,"geneID"]]
length(gene.list)


crosRes = readRDS("res_cross_tissue_AD_Jansene_1000G_SPrediXcan.rds")
single = readRDS("res_signle_tissue_AD_Jansene_1000G_SPrediXcan.rds")

dim(single)

crosRes = crosRes[!is.na(crosRes[, "MultiXcan"]),]
multiXcan = multiXcan[multiXcan[,"Gene"] %in%gene.list,]

multixcan = crosRes[crosRes[, "MultiXcan"] < 0.05/dim(crosRes)[1],]

#multixcan = multixcan[multixcan[, "best_p-value"] < 1e-4,]
dim(multixcan)
multixcan.gene = multixcan[, "Gene"]
multixcan.gene = as.character(multixcan.gene)

single = single[!is.na(single[,"p-value"]),]
twas = single[single[, "p-value"] < 0.05/dim(single)[1],]

twas = twas[twas[,"Gene_name"]%in% gene.list,]


twas.gene = unique(twas[, "Gene_name"])
twas.gene = as.character(twas.gene)
length(twas.gene)


cpgRes = readRDS("res3_cross_CpG_AD_Jansene.rds")
cpgRes = cpgRes[!is.na(cpgRes[,"cross_CMO"]),]
cpgRes = cpgRes[cpgRes[,"geneID"]%in% gene.list,]

cpgRes = cpgRes[cpgRes[,"cross_CMO"]>=0,]
cpgRes = cpgRes[cpgRes[,"cross_CMO"]<0.05/dim(cpgRes)[1],]
CMO.gene = cpgRes[,2]

length(CMO.gene)

CpG = readRDS("res3_signle_CpG_AD_Jansene.rds")
CpG = CpG[CpG[,"Feature"] %in% c("GeneBody","Promoter"),]
CpG = CpG[CpG[,"geneID"]%in% gene.list,]

n.CpG = length(CpG[,"CpG"])
CpG = CpG[!is.na(CpG[,"SUM"]),]
CpG = CpG[CpG[,"SUM"]<0.05/n.CpG,]
MWAS.gene = unique(CpG[,"geneID"])

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
fig.name = "Sfig3_AD_Jansene_venn_common"
ggsave(paste(outd,fig.name, ".pdf", sep = ""))



