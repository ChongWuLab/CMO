setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna")

single = readRDS("res_signle_tissue_IGAP1_1000G_SPrediXcan.rds")


cpgRes = readRDS("res3_cross_CpG_IGAP1.rds")
cpgRes2 = readRDS("res3_cross_CpG_AD_Jansene.rds")

cpgRes = cpgRes[cpgRes[, "cross_CMO"] < 0.05 / dim(cpgRes)[1] & !is.na(cpgRes[, "cross_CMO"]),]
cpgGene = cpgRes[, 2]

cpgRes2 = cpgRes2[cpgRes2[, "gene_ACAT"] < 5e-8 & !is.na(cpgRes2[, "gene_ACAT"]),]
cpgGene2 = cpgRes2[, 2]

length(cpgGene)
length(cpgGene2)
length(intersect(cpgGene, cpgGene2))
intersect(cpgGene, cpgGene2)
colSums(cpgRes < 0.05 / 20000, na.rm = T)



cpgRes = readRDS("res_cross_CpG_IGAP1.rds")
crosRes = readRDS("res_cross_tissue_IGAP1_1000G_SPrediXcan.rds")

colSums(crosRes < 0.05 / 20000, na.rm = T)
cpgRes = cpgRes[cpgRes[, "gene_ACAT"] < 0.05 / 20000 & !is.na(cpgRes[, "gene_ACAT"]),]

crosRes[crosRes[, "ntissue"] == 1, "Omni_ACAT_MultiXcan"] =
  crosRes[crosRes[, "ntissue"] == 1, "best_p-value"]
crosRes = crosRes[crosRes[, "Omni_ACAT_MultiXcan"] < 0.05 / 20000 & !is.na(crosRes[, "Omni_ACAT_MultiXcan"]),]

cpg.gene = cpgRes[, 2]
cros.gene = crosRes[, 3]

length(cpg.gene)
length(cros.gene)
intersect(cpg.gene, cros.gene)



single = readRDS("res_signle_CpG_IGAP1.rds")
single2 = readRDS("res_signle_CpG_AD_Jansene.rds")

dim(single)
single = single[!duplicated(single[, "CpG"]),]
dim(single)
colSums(single < 0.05 / dim(single)[1], na.rm = T)

single = single2
dim(single)
single = single[!duplicated(single[, "CpG"]),]
dim(single)
colSums(single < 0.05 / dim(single)[1], na.rm = T)
