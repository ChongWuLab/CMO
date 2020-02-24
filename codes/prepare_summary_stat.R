library(data.table)
setwd("/gpfs/research/chongwu/shared/summary_statistics/AD")
dat = fread("AD.gwax.assoc")
dat = as.data.frame(dat)

dat$Z = log(dat$OR) / dat$SE
tmp = (1 - pnorm(abs(dat$Z))) * 2

dat = dat[, c("SNP", "A1", "A2", "Z", "CHR")]

dat = dat[nchar(dat$A1) == 1 & nchar(dat$A2) == 1,]
saveRDS(dat, "AD_GWAX.rds")


dat = fread("IGAP_stage_1.txt")
dat = as.data.frame(dat)

dat[, 8] = as.numeric(dat[, 8])
dat$Z = dat$Beta / dat$SE

tmp = (1 - pnorm(abs(dat$Z))) * 2
dat = dat[, c(3, 4, 5, 9, 1)]

colnames(dat) = c("SNP", "A1", "A2", "Z", "CHR")

saveRDS(dat, "IGAP1.rds")


dat = fread("AD_sumstats_Jansenetal.txt")
dat = as.data.frame(dat)

dat = dat[, c(6, 4, 5, 7, 2)]
colnames(dat) = c("SNP", "A1", "A2", "Z", "CHR")
saveRDS(dat, "AD_Jansene.rds")

