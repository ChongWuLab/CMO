######################################################
# Multiple Imputation for cell type proportion    ####
# Version 1.00                                    ####
# Feb 2, 2015                                     ####
# Author: Chong Wu, Weihua Guan                   ####
######################################################

library(data.table)

dat = fread("/gpfs/research/chongwu/shared/summary_statistics/AD/AD_sumstats_Jansenetal.txt")
dat = as.data.frame(dat)



gene.list = readRDS("/gpfs/research/chongwu/Chong/Multi-tissue/PrediXcan_enet_weight_all_tissues.rds")
gene.list = gene.list[, 3:7]
gene.list = gene.list[!duplicated(gene.list),]

gene.list = gene.list[order(gene.list[, "CHR"], gene.list[, "P0"]),]

gene.list[, 2] = as.character(gene.list[, 2])

final.res = as.data.frame(matrix(NA, dim(gene.list)[1], 7))
final.res[, 1:5] = gene.list

out.i = 1

for (i in 1:22) {

  dat.tmp = dat[dat[, 2] == i,]
  gene.candidate = gene.list[gene.list[, "CHR"] == i,]

  for (j in 1:nrow(gene.candidate)) {
    cis500.tmp = dat.tmp[dat.tmp[, 3] > as.numeric(gene.candidate[j, "P0"]) - 500 * 1000 & dat.tmp[, 3] < as.numeric(gene.candidate[j, "P1"]) + 500 * 1000,]

    cis1000.tmp = dat.tmp[dat.tmp[, 3] > as.numeric(gene.candidate[j, "P0"]) - 1000 * 1000 & dat.tmp[, 3] < as.numeric(gene.candidate[j, "P1"]) + 1000 * 1000,]

    final.res[out.i,] = c(as.character(gene.candidate[j,]), min(cis500.tmp[, 8]), min(cis1000.tmp[, 8]))
    out.i = out.i + 1
  }
}

colnames(final.res) = c("ensembl", "Gene", "CHR", "P0", "P1", "most_sig_500Kb", "most_sig_1000Kb")
saveRDS(final.res, "/gpfs/research/chongwu/Chong/MWAS/resAna/summary_best_cis_MWAS_AD_Jansene.rds")
