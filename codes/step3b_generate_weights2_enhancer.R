#!/usr/bin/env Rscript
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
chr.id <- as.numeric(slurm_arrayid)

setwd("/gpfs/research/chongwu/Chong/MWAS")
library(data.table)

load("./WEIGHTS/NET_Methylation_download.RData")

info = as.data.frame(info)
info[,2:5] = apply(info[,2:5],2,as.numeric)
info$indx = 1:dim(info)[1]

enhancer = readRDS("processed_enhancer_CpG_inf.rds")

enhancer = enhancer[enhancer[,1]==paste("chr",chr.id,sep=""),]
enhancer = enhancer[!is.na(enhancer$CpG),]

enhancer.out = enhancer

enhancer = enhancer[!duplicated(enhancer$id),]
enhancer$nCpG = -1
enhancer$weights_dir = NA

info2 = info
NET2 = NET
LD2 = LD


rm(LD)
rm(NET)

for(i in 1:dim(enhancer)[1]) {
    tmp = strsplit(enhancer[i,"CpG"],";")
    tmp = unlist(tmp)
    
    info = info2[info2[,1] %in% tmp,]

    used.cpg = info$indx

    NET = NET2[used.cpg]
    LD = LD2[used.cpg]

    weights.dir = paste("enhancer/",enhancer[i,"id"],".RData",sep="")

    enhancer[i,"nCpG"] = length(tmp)
    enhancer[i,"weights_dir"] = weights.dir

    save(info,NET,LD ,file = paste("/gpfs/research/chongwu/Chong/MWAS/WEIGHTS/", weights.dir,sep=""))
}

rownames(enhancer) = enhancer[,"id"]
enhancer = enhancer[enhancer.out[,"id"],]
enhancer.out$nCpG = enhancer[,"nCpG"] 
enhancer.out$weights_dir = enhancer[,"weights_dir"] 


saveRDS(enhancer.out,paste("processed_enhancer_CpG_mapped_CHR",chr.id,".rds",sep=""))

final.dat = NULL
for(chr.id in 1:22) {
    tmp = readRDS(paste("processed_enhancer_CpG_mapped_CHR",chr.id,".rds",sep=""))
    final.dat = rbind(final.dat,tmp)
}

saveRDS(final.dat,"processed_enhancer_CpG_mapped_inf.rds")
