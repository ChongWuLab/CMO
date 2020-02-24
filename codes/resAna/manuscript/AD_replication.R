setwd("/home/chong/Dropbox/FSU_research/R21 Grant Application/MWAS/resAna")
setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/AD_res")
library(VennDiagram)

cpgRes = readRDS("res3_cross_CpG_IGAP1.rds")
cpgRes = cpgRes[!is.na(cpgRes[,"cross_CMO"]),]
cpgRes = cpgRes[cpgRes[,"cross_CMO"]>=0,]
cpgRes = cpgRes[cpgRes[,"cross_CMO"]<0.05/dim(cpgRes)[1],]
CMO.gene = cpgRes[,2]


list.len = length(CMO.gene)


gwax = readRDS("res3_cross_CpG_AD_GWAX.rds")
gwax = gwax[!is.na(gwax[,"cross_CMO"]),]
gwax = gwax[gwax[,"cross_CMO"]>=0,]

PopSize = dim(gwax)[1]

list2 = sum(gwax[, "cross_CMO"] < 0.05, na.rm = T)
list = sum(gwax[, "cross_CMO"] < 0.05 / PopSize, na.rm = T)


gwax = gwax[gwax[, "geneID"] %in% CMO.gene,]
dim(gwax)
overlap2 = sum(gwax[, "cross_CMO"] < 0.05, na.rm = T)
overlap = sum(gwax[, "cross_CMO"] < 0.05 / PopSize, na.rm = T)

phyper(overlap - 1, list.len, PopSize - list.len, list, lower.tail = FALSE)
overlap
list

phyper(overlap2 - 1, list.len, PopSize - list.len, list2, lower.tail = FALSE)
overlap2
list2


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


CMO.gene[CMO.gene %in% gene.list]

phyper(sum(CMO.gene %in% gene.list) - 1, list.len, PopSize - list.len, length(gene.list), lower.tail = FALSE)

CMO.gene[!CMO.gene %in% gene.list]
