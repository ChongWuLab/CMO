
SMWAS.sum <- function(wgtlist0, opt) {
  final.out = NULL
  snp.gene = NULL
  snp.enhancer = NULL

  for (w in 1:nrow(wgtlist0)) {
    wgt.file <- paste(opt$weights_dir, wgtlist0$weights_dir[w], sep = "")
    load(wgt.file)

    info[, 1] = as.character(info[, 1])
    res.save = as.data.frame(matrix(NA, length(NET), 12))
    start.time = proc.time()[3]
    snp.used = NULL
    for (i in 1:length(NET)) {

      ##        select SNPs present in GWAS
      ii = which(!is.na(match(NET[[i]][, 2], gwas_snp)));
      if (length(ii) > 2) {
        #  we only select  models with > 2 SNPs

        NET.used = NET[[i]][ii,]
        LD.used = LD[[i]][ii, ii]

        weight = t(as.matrix(as.numeric(NET.used[, 3]))) # MWAS/TWAS model weights
        weight_diag <- diag(as.vector(weight), nrow = length(weight))

        ii = match(NET.used[, 2], gwas_snp);
        gwasz.used = gwas_z[ii]
        snp.tmp = gwas_snp[ii]
        snp.used = c(snp.used, snp.tmp)

        Zstat.w <- weight_diag %*% gwasz.used
        corSNP.w <- weight_diag %*% LD.used %*% t(weight_diag)
        cur.twas <- (weight %*% gwasz.used) / sqrt(weight %*% LD.used %*% t(weight))

        pSUM <- 2 * (pnorm(abs(cur.twas), lower.tail = F))
        pSSU <- SumSqU(U = Zstat.w, CovS = corSNP.w, method = "Pan")

        gwasp.used = pnorm(abs(gwasz.used),lower.tail=F) * 2
        used.weight = abs(weight)/ sum(abs(weight))
        pACAT = ACAT(gwasp.used,used.weight)
        
        paSumSSU <- ACAT(c(pSUM, pSSU), c(1 / 2, 1 / 2))
        paSumSSUACAT <- ACAT(c(pSUM, pSSU,pACAT), c(1 / 3, 1 / 3, 1/3))
        
        res.save[i, 1:5] <- c(length(gwasz.used), cur.twas, pSUM, pSSU, paSumSSU)
        res.save[i, 7:10] <- info[i, 2:5]
        res.save[i, 6] <- info[i, 1]
        res.save[i,11:12] <- c(pACAT,paSumSSUACAT)
      }
    }
    res.save = res.save[!is.na(res.save[, 1]),]
    res.save$feature = wgtlist0[w, 2]

    if (wgtlist0[w, 2] == "GeneBody") {
      snp.gene = c(snp.gene, snp.used)
    } else {
      snp.enhancer = c(snp.enhancer, snp.used)
    }
    final.out = rbind(final.out, res.save)
  }

  colnames(final.out) <- c("n_infoSNP", "Z_score", "SUM", "SSU", "SUM_plus_SSU", "CpG", "Lambda", "n_SNPs_into_Lasso", "n_SNPs_out_Lasso", "explained_variance", "ACAT","aSumSSUACAT","Feature")
  return(list(final.out = final.out, snp.enhancer = snp.enhancer, snp.gene = snp.gene))
}








