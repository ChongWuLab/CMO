library(xtable)

read.file2 <- function(file.name) {
  file <- try(readRDS(file.name))
  if (class(file) == "try-error") {
    #cat("Caught an error during fread, trying read.table.\n")
    #cat(file.name)

  }
  file
}

library(data.table)

out.fun <- function(pheno.input) {
  setwd(paste("/gpfs/research/chongwu/Chong/MWAS/res/res_", pheno.input, sep = ""))

  res <- NULL
  single.res <- NULL
  for (index in 1:22) {
    file.name <- paste("out_CHR", index, ".rds", sep = "")

    res.temp <- read.file2(file.name)

    if (!is.null(res.temp) & is.data.frame(res.temp)) {
      res = rbind(res, res.temp)
    }

    file.name <- paste("single_CpG_out_CHR", index, ".rds", sep = "")

    res.temp <- read.file2(file.name)

    if (!is.null(res.temp) & is.data.frame(res.temp)) {
      single.res = rbind(single.res, res.temp)
    }
  }


  out.file = paste("/gpfs/research/chongwu/Chong/MWAS/resAna/res_cross_CpG_", pheno.input, ".rds", sep = "")
  saveRDS(res, out.file)

  out.file = paste("/gpfs/research/chongwu/Chong/MWAS/resAna/res_signle_CpG_", pheno.input, ".rds", sep = "")
  saveRDS(single.res, out.file)

}

out.fun("GWAX")

out.fun("AD_Jansene")
out.fun("IGAP1")



out.fun <- function(pheno.input) {
  setwd(paste("/gpfs/research/chongwu/Chong/Multi-tissue/res/cross_tissue_", pheno.input, "_1000G_SPrediXcan", sep = ""))

  res <- NULL
  single.res <- NULL
  for (index in 1:22) {
    file.name <- paste("out_CHR", index, ".rds", sep = "")

    res.temp <- read.file2(file.name)

    if (!is.null(res.temp) & is.data.frame(res.temp)) {
      res = rbind(res, res.temp)
    }

    file.name <- paste("single_tissue_out_CHR", index, ".rds", sep = "")

    res.temp <- read.file2(file.name)

    if (!is.null(res.temp) & is.data.frame(res.temp)) {
      single.res = rbind(single.res, res.temp)
    }
  }


  out.file = paste("/gpfs/research/chongwu/Chong/MWAS/resAna/res_cross_tissue_", pheno.input, "_1000G_SPrediXcan.rds", sep = "")
  saveRDS(res, out.file)

  out.file = paste("/gpfs/research/chongwu/Chong/MWAS/resAna/res_signle_tissue_", pheno.input, "_1000G_SPrediXcan.rds", sep = "")
  saveRDS(single.res, out.file)

}

out.fun("AD_Jansene")
out.fun("IGAP1")


