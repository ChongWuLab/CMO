#!/usr/bin/env Rscript
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
chr.id <- as.numeric(slurm_arrayid)

# set working dir
dat.dir <- "/gpfs/research/chongwu/shared/summary_statistics/AD/"
work.dir <- getwd()

setwd(work.dir)

library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
library(mvtnorm)
library(data.table)
#sourceCpp("twas_perm.cpp")
#source("aSPUwscore.R")

source("dist_support.R")
source("mwas_support.R")

suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))

option_list = list(
make_option("--sumstats", action="store", default=NA, type='character',
help="summary statistics (must have SNP and Z column headers; support rds and txt format) [required]"),
make_option("--out", action="store", default=NA, type='character',
help="Path to output files [required]"),
make_option("--chr_id", action="store", default=-1, type='character',
help="The chromosome ID. We recommend parallel the computations by chromosomes")
)


opt = parse_args(OptionParser(option_list=option_list))


gene = readRDS("processed_gene_CpG_mapped_inf.rds") #17296 genes
gene$Feature = "GeneBody"
gene$enhancer_score = NA
gene$geneScore = NA

gene = gene[, c(4, 11, 5, 6, 12, 3, 13, 9, 10, 1)]
colnames(gene) = c("CHR", "Feature", "P0", "P1", "enhancer_score", "geneID", "geneScore", "nCpG", "weights_dir", "ensembl")



enhancer = readRDS("processed_enhancer_CpG_mapped_inf.rds")
enhancer[, 1] = gsub("chr", "", enhancer[, 1])
enhancer = enhancer[enhancer[, "geneID"] %in% gene[, "geneID"],]

enhancer = enhancer[enhancer[, "geneScore"] > 10 & enhancer[, "score"] > 1,]

enhancer = enhancer[order(enhancer[, 1], enhancer[, 4], enhancer[, 5], - enhancer[, 9]),]

enhancer = enhancer[, c(1, 3, 10, 11, 6, 8, 9, 13, 14)]
colnames(enhancer) = c("CHR", "Feature", "P0", "P1", "enhancer_score", "geneID", "geneScore", "nCpG", "weights_dir")

enhancer$ensembl = NA
enhancer = enhancer[enhancer[, "geneID"] %in% gene[, "geneID"],]
enhancer = enhancer[order(enhancer[, 1], enhancer[, 3], enhancer[, 4]),]

wgtlist = rbind(enhancer, gene)
wgtlist = wgtlist[order(wgtlist[, 1], wgtlist[, 3], wgtlist[, 4]),]



outd <- paste(work.dir, "/res/", pfx, "_", pheno, sep = "")
system(paste("mkdir -p ", outd, sep = ""))


if(grepl(".rds",opt$sumstats)) {
    sumstat.orgin = readRDS(opt$sumstats)
} else if (grepl(".txt",opt$sumstats)) {
    sumstat.orgin = fread(opt$sumstats)
    sumstat.orgin = as.data.frame(sumstat.orgin)
} else {
    sumstat.orgin = fread(opt$sumstats)
    sumstat.orgin = as.data.frame(sumstat.orgin)
    warning("We try to read GWAS summary data other than txt and rds format. The data may not read properly")
}

if(opt$chr_id<=22 & opt$chr_id>=1) {
    job = opt$chr_id
    wgtlist <- wgtlist[wgtlist[, 1] == job,]
    sumstat.orgin <- sumstat.orgin[sumstat.orgin$CHR == job,]

} else if( opt$chr_id!= -1) {
    warning("Chromosomes ID must between 1 and 22. We ignore this argument and run the analyses for all the genes.")
}
used.gene <- unique(wgtlist$geneID)


#w = 1
#wgtlist0 = wgtlist[wgtlist[,"geneID"]=="APOE",]
#wgt.file <- paste(opt$weights_dir, wgtlist0$weights_dir[w], sep = "")
#load(wgt.file)
gwas_snp = sumstat.orgin[, "SNP"] # rs ID gwas SNP

a1 = sumstat.orgin[, "A1"]
a2 = sumstat.orgin[, "A2"]

keep = !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
sumstat.orgin = sumstat.orgin[keep,]

load("/gpfs/research/chongwu/Chong/MWAS/WEIGHTS/NET_Methylation_download.RData") # loads: LD, NET, info

NET_all = do.call(rbind, NET) # make matrix from list NET
NET_all = NET_all[which(!duplicated(NET_all[, 2])),] # only keep unique SNPs
m = match(NET_all[, 2], gwas_snp);
NET_all = NET_all[which(!is.na(m)),] # only keep SNPs present in GWAS to speed up matching
m = m[which(!is.na(m))]

sumstat.orgin = sumstat.orgin[m,]


a1 = sumstat.orgin[, "A1"]
a2 = NET_all[, 1] # get gwas allele and allele used in TWAS/MWAS models
a3 = a2;
a3[which(a2 == 'A')] = 'T';
a3[which(a2 == 'G')] = 'C';
a3[which(a2 == 'T')] = 'A';
a3[which(a2 == 'C')] = 'G' ## account for different strand mapping
m1 = (a1 == a2);
m2 = (a1 == a3)
sign = as.numeric(m1 | m2);
sign[which(sign == 0)] = -1

sum(sign == -1)
gwas_z = sumstat.orgin[, "Z"] * sign # transform gwas Z for non matching allele's

gwas_snp = sumstat.orgin[, "SNP"]

# need add the corresponding $p$-values information.

rm(sumstat.orgin)
rm(LD)
rm(NET)
rm(NET_all)


out.res <- as.data.frame(matrix(NA, length(used.gene), 27))
colnames(out.res) <- c("CHR", "geneID", "ensembl", "P0", "P1", "n_enhancer", "best_SUM_Feature", "best_n_SNP", "best_Z_score", "best_p-SUM", "worst_SUM_Feature", "worst_n_SNP", "worst_Z_score", "worst_p-SUM", "best_ACAT_Feature", "best_nSNP_ACAT", "best_p_ACAT", "worst_ACAT_Feature", "worst_nSNP_ACAT", "worst_p_ACAT", "n_enhancer_CpG", "n_gene_CpG", "n_enhancer_SNP"
, "n_gene_SNP", "n_total_SNP", "CMO", "runtime(s)")

res.single <- NULL

for (gene.indx in 1:length(used.gene)) {
  #length(used.gene)
  tryCatch({
    #gene.indx = 10
    # Load in summary stats
    start.time <- proc.time()[3]
    wgtlist1 <- wgtlist[wgtlist[, "geneID"] == used.gene[gene.indx],]
    #wgtlist1 <- wgtlist[wgtlist[, "geneID"] == "AGRN",]

    res <- SMWAS.sum(wgtlist1, opt) #need revise here

    res.save <- res$final.out
    
    n.enhancer.snp = length(unique(res$snp.enhancer))
    n.gene.snp = length(unique(res$snp.gene))
    n.snp = length(unique(c(res$snp.enhancer, res$snp.gene)))

    tmp = wgtlist1[wgtlist1[, "Feature"] == "GeneBody",]
    res.save$geneID <- tmp[1, "geneID"]
    res.save$ensembl <- tmp[1, "ensembl"]

    res.single = rbind(res.single, res.save)

    ###########################
    # combine the results   ###
    ###########################

    # combine Sum Test results
    acat.tmP <- res.save[, "SUM"]
    acat.tmP <- acat.tmP[!is.na(acat.tmP)]
    weight <- rep(1 / length(acat.tmP), length(acat.tmP))
    cross_SUM_p <- ACAT(acat.tmP, weight)

    res.tmp = res.save[res.save[, "Feature"] == "GeneBody",]
    n.gene.cpg = dim(res.tmp)[1]
    acat.tmP <- res.tmp[, "SUM"]
    acat.tmP <- acat.tmP[!is.na(acat.tmP)]
    weight <- rep(1 / length(acat.tmP), length(acat.tmP))
    gene_SUM_p <- ACAT(acat.tmP, weight)

    res.tmp = res.save[res.save[, "Feature"] != "GeneBody",]
    n.enhancer.cpg = dim(res.tmp)[1]
    acat.tmP <- res.tmp[, "SUM"]
    acat.tmP <- acat.tmP[!is.na(acat.tmP)]
    weight <- rep(1 / length(acat.tmP), length(acat.tmP))
    enhancer_SUM_p <- ACAT(acat.tmP, weight)

    # combine adaptive test results
    acat.tmP <- res.save[, "aSumSSUACAT"]
    acat.tmP <- acat.tmP[!is.na(acat.tmP)]
    weight <- rep(1 / length(acat.tmP), length(acat.tmP))
    cross_ACAT_p <- ACAT(acat.tmP, weight)
    
    res.tmp = res.save[res.save[, "Feature"] == "GeneBody",]
    n.gene.cpg = dim(res.tmp)[1]
    acat.tmP <- res.tmp[, "aSumSSUACAT"]
    acat.tmP <- acat.tmP[!is.na(acat.tmP)]
    weight <- rep(1 / length(acat.tmP), length(acat.tmP))
    gene_ACAT_p <- ACAT(acat.tmP, weight)
    
    res.tmp = res.save[res.save[, "Feature"] != "GeneBody",]
    n.enhancer.cpg = dim(res.tmp)[1]
    acat.tmP <- res.tmp[, "aSumSSUACAT"]
    acat.tmP <- acat.tmP[!is.na(acat.tmP)]
    weight <- rep(1 / length(acat.tmP), length(acat.tmP))
    enhancer_ACAT_p <- ACAT(acat.tmP, weight)
    
    ###########################
    ### save information    ###
    ###########################
    tmp = wgtlist1[wgtlist1[, "Feature"] == "GeneBody",]
    n.enhancer = dim(wgtlist1)[1] - 1
    basic.inf <- c(tmp[1, "CHR"], as.character(tmp[1, "geneID"]), as.character(tmp[1, "ensembl"]), tmp[1, "P0"], tmp[1, "P1"], n.enhancer)
    out.res[gene.indx, 1:6] <- basic.inf
    
    tmp <- res.save[res.save[, "SUM"] == min(res.save[, "SUM"]),]
    tmp <- tmp[1,]
    out.res[gene.indx, 7] <- as.character(tmp[, "Feature"])
    out.res[gene.indx, 8:10] <- tmp[c("n_infoSNP", "Z_score", "SUM")]

    tmp <- res.save[res.save[, "SUM"] == max(res.save[, "SUM"]),]
    tmp <- tmp[1,]
    out.res[gene.indx, 11] <- as.character(tmp[, "Feature"])
    out.res[gene.indx, 12:14] <- tmp[c("n_infoSNP", "Z_score", "SUM")]

    tmp <- res.save[res.save[, "aSumSSUACAT"] == min(res.save[, "aSumSSUACAT"]),]
    tmp <- tmp[1,]
    out.res[gene.indx, 15] <- as.character(tmp[, "Feature"])
    out.res[gene.indx, 16:17] <- tmp[c("n_infoSNP", "aSumSSUACAT")]

    tmp <- res.save[res.save[, "aSumSSUACAT"] == max(res.save[, "aSumSSUACAT"]),]
    tmp <- tmp[1,]
    out.res[gene.indx, 18] <- as.character(tmp[, "Feature"])
    out.res[gene.indx, 19:20] <- tmp[c("n_infoSNP", "aSumSSUACAT")]

    out.res[gene.indx, 21:25] = c(n.enhancer.cpg, n.gene.cpg, n.enhancer.snp, n.gene.snp, n.snp)
    out.res[gene.indx, 26] = cross_ACAT_p
    
    cat("Finish Gene",gene.indx,"\t")
    
    end.time <- proc.time()[3]
    run.time <- (end.time - start.time)
    out.res[gene.indx, 27] <- run.time

  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")

    tmp = wgtlist1[wgtlist1[, "Feature"] == "GeneBody",]
    n.enhancer = dim(wgtlist1)[1] - 1
    basic.inf <- c(tmp[1, "CHR"], as.character(wgtlist1[1, "geneID"]), as.character(wgtlist1[1, "ensembl"]), wgtlist1[1, "P0"], wgtlist1[1, "P1"], n.enhancer)
    out.res[gene.indx, 1:6] <- basic.inf
  })
}

out.file <- paste(outd, "/out_CHR", job, ".rds", sep = "")
saveRDS(out.res, out.file)

out.file <- paste(outd, "/single_CpG_out_CHR", job, ".rds", sep = "")

saveRDS(res.single, out.file)



