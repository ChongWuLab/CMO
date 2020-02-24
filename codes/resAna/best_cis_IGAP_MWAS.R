######################################################
# Multiple Imputation for cell type proportion    ####
# Version 1.00                                    ####
# Feb 2, 2015                                     ####
# Author: Chong Wu, Weihua Guan                   ####
######################################################

library(data.table)

dat = fread("/gpfs/research/chongwu/shared/summary_statistics/AD/IGAP_stage_1.txt")
dat = as.data.frame(dat)



gene.list = readRDS("/gpfs/research/chongwu/Chong/MWAS/processed_gene_CpG_mapped_inf.rds")
gene.list = gene.list[, c(1, 3, 4, 5, 6)]
colnames(gene.list) = c("ensembl", "Gene", "CHR", "P0", "P1")

final.res = as.data.frame(matrix(NA, dim(gene.list)[1], 7))
final.res[, 1:5] = gene.list

out.i = 1
dat[, "Pvalue"] = as.numeric(dat[, "Pvalue"])
dat[, 2] = as.numeric(dat[, 2])
for (i in 1:22) {

  dat.tmp = dat[dat[, "Chromosome"] == i,]
  gene.candidate = gene.list[gene.list[, "CHR"] == i,]

  for (j in 1:nrow(gene.candidate)) {
    cis500.tmp = dat.tmp[dat.tmp[, "Position"] > as.numeric(gene.candidate[j, "P0"]) - 500 * 1000 & dat.tmp[, "Position"] < as.numeric(gene.candidate[j, "P1"]) + 500 * 1000,]

    cis1000.tmp = dat.tmp[dat.tmp[, "Position"] > as.numeric(gene.candidate[j, "P0"]) - 1000 * 1000 & dat.tmp[, "Position"] < as.numeric(gene.candidate[j, "P1"]) + 1000 * 1000,]

    final.res[out.i,] = c(as.character(gene.candidate[j,]), min(cis500.tmp[, "Pvalue"]), min(cis1000.tmp[, "Pvalue"]))
    out.i = out.i + 1
  }
}

colnames(final.res) = c("ensembl", "Gene", "CHR", "P0", "P1", "most_sig_500Kb", "most_sig_1000Kb")
saveRDS(final.res, "/gpfs/research/chongwu/Chong/MWAS/resAna/summary_best_cis_MWAS_IGAP1.rds")
