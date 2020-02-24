#!/usr/bin/env Rscript
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
chr.id <- as.numeric(slurm_arrayid)

setwd("/gpfs/research/chongwu/Chong/MWAS")
library(data.table)


ann <- readRDS("IlluminaHumanMethylation450kanno.rds")


final.dat = fread("ENSEMBL_GRch37_gene_list.txt")
final.dat = as.data.frame(final.dat)
final.dat$CpG = NA


final.dat = final.dat[final.dat[,"chromosome_name"]==chr.id,]

ann = ann[ann[,1]==paste("chr",chr.id,sep=""),]
all.id = unique(final.dat$id)


for(indx in 1:dim(final.dat)[1]) {
    tmp.ann = ann[ann$pos >= (final.dat[indx,"start_position"] - 500) & ann$pos <= (final.dat[indx,"end_position"] + 500), ]
    if(dim(tmp.ann)[1] > 0) {
        final.dat[indx,"CpG"] = paste(tmp.ann[,"Name"],collapse=";")
    }
}

saveRDS(final.dat,paste("processed_gene_CpG_inf_CHR",chr.id,".rds",sep=""))


final.dat = NULL
for(chr.id in 1:22) {
    tmp = readRDS(paste("processed_gene_CpG_inf_CHR",chr.id,".rds",sep=""))
    final.dat = rbind(final.dat,tmp)
}

saveRDS(final.dat,"processed_gene_CpG_inf.rds")
