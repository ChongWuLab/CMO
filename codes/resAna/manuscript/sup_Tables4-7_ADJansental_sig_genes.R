
setwd("/home/chong/Dropbox/FSU_research/R21 Grant Application/MWAS/resAna")
setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/AD_res/")
library(VennDiagram)
library(xtable)

# get GWAS Catalog information
## replication by GWAS Catalog
res = as.data.frame(matrix(NA, 1, 5))
library(data.table)

gene.inf = fread("gwas_catalog_v1.0-associations_e98_r2019-12-16.tsv")
gene.inf = as.data.frame(gene.inf)

tmp.name = gene.inf[, "MAPPED_GENE"]
tmp.name = as.character(tmp.name)
tmp.name = unlist(strsplit(tmp.name, " - "))
tmp.name = unlist(strsplit(tmp.name, ","))
PopSize = length(unique(tmp.name))


gene.inf = gene.inf[grepl("Alzheimer", gene.inf[, 8]),]


for (i in 1:dim(gene.inf)[1]) {
    gene.inf.tmp = gene.inf[i,]
    tmp.name = gene.inf.tmp["MAPPED_GENE"]
    tmp.name = as.character(tmp.name)
    tmp.name = unlist(strsplit(tmp.name, " - "))
    tmp.name = unlist(strsplit(tmp.name, ","))
    if (is.null(tmp.name)) next
    tmp.name = unlist(strsplit(tmp.name, ";"))
    tmp.name = trimws(tmp.name)
    tmp.name1 = tmp.name
    
    tmp.name = gene.inf.tmp[14]
    tmp.name = as.character(tmp.name)
    tmp.name = unlist(strsplit(tmp.name, " - "))
    tmp.name = unlist(strsplit(tmp.name, ", "))
    tmp.name = unlist(strsplit(tmp.name, ","))
    if (is.null(tmp.name)) next
    tmp.name = unlist(strsplit(tmp.name, ";"))
    tmp.name = trimws(tmp.name)
    tmp.name = c(tmp.name, tmp.name1)
    tmp.name = unique(tmp.name)
    
    res.tmp = as.data.frame(matrix(NA, length(tmp.name), 5))
    res.tmp[, 2:5] = gene.inf.tmp[c(28, 8, 6, 7)] #c("P.VALUE","DISEASE.TRAIT","LINK",Study
    res.tmp[, 1] = tmp.name
    res = rbind(res, res.tmp)
}

res = res[-1,]

gene.list = unique(res[, 1])

cpgRes = readRDS("res3_cross_CpG_AD_Jansene.rds")
cpgRes = cpgRes[!is.na(cpgRes[,"cross_CMO"]),]
CMO = cpgRes[cpgRes[,"cross_CMO"]>=0,]


CMO = CMO[CMO[,"cross_CMO"]<0.05/dim(CMO)[1],]

#build table
out = cbind(CMO[,1:6],CMO[,"cross_CMO"])

gwas.catalog = res[res[,1] %in% out[,2],]

final.output = as.data.frame(matrix(-1,dim(out)[1],8))
final.output[,8] = "remNA"
final.output[,1:7] = out
gene.name = out[,2]
tmp.index = NULL
for(index in 1:dim(out)[1]) {
    if(gene.name[index] %in% gene.list) {
        cat("index:" ,index,"  ")
        tmp = res[res[,1] %in% gene.name[index],]
        tmp.index = rbind(tmp.index,tmp)
        tmp = tmp[1,]
        final.output[index,8] = tmp[4]
    }
}

write.csv(final.output,"AD_Jansene_sig_genes.csv")


CpG = readRDS("res3_signle_CpG_AD_Jansene.rds")
CpG = CpG[CpG[,"Feature"] %in% c("GeneBody","Promoter"),]
n.CpG = length(CpG[,"CpG"])
CpG = CpG[!is.na(CpG[,"SUM"]),]
CpG = CpG[CpG[,"SUM"]<0.05/n.CpG,]

CpG = CpG[order(CpG[,"geneID"],CpG[,"SUM"]),]
CpG = CpG[!duplicated(CpG[,"geneID"]),]
CpG = CpG[,c("geneID","SUM")]

tmp = cpgRes[cpgRes[,2] %in% CpG[,1],]
tmp = tmp[,1:6]

rownames(tmp) = tmp[,2]
rownames(CpG) = CpG[,1]
CpG = CpG[rownames(tmp),]

out = cbind(tmp,CpG[,2])


final.output = as.data.frame(matrix(-1,dim(out)[1],8))
final.output[,8] = "remNA"
final.output[,1:7] = out
gene.name = out[,2]
tmp.index = NULL
for(index in 1:dim(out)[1]) {
    if(gene.name[index] %in% gene.list) {
        cat("index:" ,index,"  ")
        tmp = res[res[,1] %in% gene.name[index],]
        tmp.index = rbind(tmp.index,tmp)
        tmp = tmp[1,]
        final.output[index,8] = tmp[4]
    }
}

write.csv(final.output,"AD_Jansene_sig_genes_MWAS.csv")

#cpgRes = cpgRes[cpgRes[,"cross_CMO"]<0.05/dim(cpgRes)[1],]
##########
# read the data
####
crosRes = readRDS("res_cross_tissue_AD_Jansene_1000G_SPrediXcan.rds")
single = readRDS("res_signle_tissue_AD_Jansene_1000G_SPrediXcan.rds")

single = single[single[,"p-value"]<0.05/dim(single)[1],]

out = single[,c(1,2:6,11)]


final.output = as.data.frame(matrix(-1,dim(out)[1],8))
final.output[,8] = "remNA"
final.output[,1:7] = out
gene.name = out[,4]
tmp.index = NULL
for(index in 1:dim(out)[1]) {
    if(gene.name[index] %in% gene.list) {
        cat("index:" ,index,"  ")
        tmp = res[res[,1] %in% gene.name[index],]
        tmp.index = rbind(tmp.index,tmp)
        tmp = tmp[1,]
        final.output[index,8] = tmp[4]
    }
}

write.csv(final.output,"AD_Jansene_sig_genes_TWAS.csv")

dim(single)

crosRes = crosRes[!is.na(crosRes[, "MultiXcan"]),]
multixcan = crosRes[crosRes[, "MultiXcan"] < 0.05/dim(crosRes)[1],]

tmp = multixcan[,c(1,3,4,5,6)]

out = cbind(tmp,multixcan[,"MultiXcan"])


final.output = as.data.frame(matrix(-1,dim(out)[1],8))
final.output[,8] = "remNA"
final.output[,1:6] = out
gene.name = out[,2]
tmp.index = NULL
for(index in 1:dim(out)[1]) {
    if(gene.name[index] %in% gene.list) {
        cat("index:" ,index,"  ")
        tmp = res[res[,1] %in% gene.name[index],]
        tmp.index = rbind(tmp.index,tmp)
        tmp = tmp[1,]
        final.output[index,8] = tmp[4]
    }
}

write.csv(final.output,"AD_Jansene_sig_genes_MultiXcan.csv")
