#!/usr/bin/env Rscript
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
chr.id <- as.numeric(slurm_arrayid)

setwd("/gpfs/research/chongwu/Chong/MWAS")
library(data.table)


ann <- readRDS("IlluminaHumanMethylation450kanno.rds")


final.dat = readRDS("processed_enhancer_inf.rds")
final.dat$CpG = NA
final.dat = final.dat[final.dat[,1]==paste("chr",chr.id,sep=""),]

ann = ann[ann[,1]==paste("chr",chr.id,sep=""),]
all.id = unique(final.dat$id)

tmp.dat = final.dat[!duplicated(final.dat$id),]
all.id = tmp.dat$id

for(indx in 1:dim(tmp.dat)[1]) {
    tmp.ann = ann[ann$pos >= tmp.dat[indx,"P0hg19"] & ann$pos <= tmp.dat[indx,"P1hg19"], ]
    if(dim(tmp.ann)[1] > 0) {
        tmp.dat[indx,"CpG"] = paste(tmp.ann[,"Name"],collapse=";")
    }
    if(indx %% 1000 ==0) {
        cat("indx ",indx)
    }
}

rownames(tmp.dat) = tmp.dat$id
tmp.dat2 = tmp.dat[final.dat$id,]
final.dat$CpG = tmp.dat2$CpG

saveRDS(final.dat,paste("processed_enhancer_CpG_inf_CHR",chr.id,".rds",sep=""))


final.dat = NULL
for(chr.id in 1:22) {
    tmp = readRDS(paste("processed_enhancer_CpG_inf_CHR",chr.id,".rds",sep=""))
    final.dat = rbind(final.dat,tmp)
}

saveRDS(final.dat,"processed_enhancer_CpG_inf.rds")