setwd("/home/chong/Dropbox/FSU_research/R21 Grant Application/MWAS/resAna")
library(VennDiagram)

crosRes = readRDS("res_cross_tissue_IGAP1_1000G_SPrediXcan.rds")
single = readRDS("res_signle_tissue_IGAP1_1000G_SPrediXcan.rds")

#crosRes = readRDS("res_cross_tissue_AD_Jansene_1000G_SPrediXcan.rds")
#single = readRDS("res_signle_tissue_AD_Jansene_1000G_SPrediXcan.rds")


dim(single)
crosRes = crosRes[crosRes[, "ntissue"] > 0,]
crosRes[crosRes[, "ntissue"] == 1, "Omni_ACAT_MultiXcan"] =
  crosRes[crosRes[, "ntissue"] == 1, "best_p-value"]

crosRes[crosRes[, "ntissue"] == 1, "MultiXcan"] =
  crosRes[crosRes[, "ntissue"] == 1, "best_p-value"]

cto = crosRes[crosRes[, "Omni_ACAT_MultiXcan"] < 2.1e-6 & !is.na(crosRes[, "Omni_ACAT_MultiXcan"]),]
cto.gene = cto[, "Gene"]
cto.gene = as.character(cto.gene)

list.len = length(cto.gene)

gwax = readRDS("res_cross_tissue_AD_GWAX_1000G_SPrediXcan.rds")

gwax[gwax[, "ntissue"] == 1, "Omni_ACAT_MultiXcan"] =
  gwax[gwax[, "ntissue"] == 1, "best_p-value"]
gwax = gwax[!is.na(gwax[, "Omni_ACAT_MultiXcan"]),]
PopSize = dim(gwax)[1]

list2 = sum(gwax[, "Omni_ACAT_MultiXcan"] < 0.05, na.rm = T)
list = sum(gwax[, "Omni_ACAT_MultiXcan"] < 0.05 / PopSize, na.rm = T)


gwax = gwax[gwax[, "Gene"] %in% cto.gene,]
overlap2 = sum(gwax[, "Omni_ACAT_MultiXcan"] < 0.05, na.rm = T)
overlap = sum(gwax[, "Omni_ACAT_MultiXcan"] < 0.05 / PopSize, na.rm = T)

phyper(overlap - 1, list.len, PopSize - list.len, list, lower.tail = FALSE)
overlap
list

phyper(overlap2 - 1, list.len, PopSize - list.len, list2, lower.tail = FALSE)
overlap2
list2


## replication by GWAS Catalog
res = as.data.frame(matrix(NA, 1, 5))
library(data.table)

gene.inf = fread("gwas_catalog_v1.0-associations_e96_r2019-04-21.tsv")
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


cto.gene[cto.gene %in% gene.list]

phyper(sum(cto.gene %in% gene.list) - 1, list.len, PopSize - list.len, length(gene.list), lower.tail = FALSE)

cto.gene[!cto.gene %in% gene.list]
