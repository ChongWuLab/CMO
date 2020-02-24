setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/UKBiobankImage_res")

all = as.character(c(paste("0",166:179,sep=""),2649:2655))

# prepare the common set gene list
indx = all[1]
multiXcan = readRDS(paste("res_cross_tissue_IDP_",indx,"_1000G_SPrediXcan.rds",sep=""))
multiXcan = multiXcan[,c(1:3,15)]
multiXcan = multiXcan[!is.na(multiXcan[,4]),]

gene.list = multiXcan[,3]

TWAS = readRDS(paste("res_signle_tissue_IDP_",indx,"_1000G_SPrediXcan.rds",sep=""))

gene.list = gene.list[gene.list %in% TWAS[,"Gene_name"]]

CTO = readRDS(paste("res3_cross_CpG_Z_",indx,".rds",sep=""))
CTO = CTO[,c(1:3,29)]
CTO = CTO[!is.na(CTO[,4]),]

gene.list = gene.list[gene.list %in% CTO[,"geneID"]]
CpG = readRDS(paste("res3_signle_CpG_Z_",indx,".rds",sep=""))

gene.list = gene.list[gene.list %in% CpG[,"geneID"]]
length(gene.list)

final.res = matrix(NA,length(all),8)
out.i = 1

colnames(final.res) = c("MultiXcan","TWAS","CMO","CpG","commonSet_MultiXcan","commonSet_TWAS","commonSet_CMO","commonSet_CpG")
for(indx in all) {
    multiXcan = readRDS(paste("res_cross_tissue_IDP_",indx,"_1000G_SPrediXcan.rds",sep=""))
    multiXcan = multiXcan[,c(1:3,15)]
    multiXcan = multiXcan[!is.na(multiXcan[,4]),]
    
    multiXcan1 = multiXcan[multiXcan[,4]<0.05/dim(multiXcan)[1],]
    final.res[out.i,1] = length(unique(multiXcan1[,"Gene"]))
    
    multiXcan2 = multiXcan[multiXcan[,"Gene"] %in%gene.list,]
    multiXcan2 = multiXcan2[multiXcan2[,4]<0.05/length(gene.list),]
    final.res[out.i,5] = length(unique(multiXcan2[,"Gene"]))

    
    TWAS = readRDS(paste("res_signle_tissue_IDP_",indx,"_1000G_SPrediXcan.rds",sep=""))
    TWAS1 = TWAS[TWAS[,"p-value"]<0.05/dim(TWAS)[1],]
    
    final.res[out.i,2] = length(unique(TWAS1[,"Gene_name"]))
    
    TWAS2 = TWAS[TWAS[,"Gene_name"]%in% gene.list,]
    TWAS2 = TWAS2[TWAS2[,"p-value"]<0.05/dim(TWAS2)[1],]

    final.res[out.i,6] = length(unique(TWAS2[,"Gene_name"]))

    CTO = readRDS(paste("res3_cross_CpG_Z_",indx,".rds",sep=""))
    CTO = CTO[,c(1:3,29)]
    CTO = CTO[!is.na(CTO[,4]),]
    CTO = CTO[CTO[,4]>=0,]
    CTO1 = CTO[CTO[,4]<0.05/dim(CTO)[1],]
    final.res[out.i,3] = length(unique(CTO1[,"geneID"]))
    
    CTO2 = CTO[CTO[,"geneID"]%in% gene.list,]
    CTO2 = CTO2[CTO2[,4]<0.05/length(gene.list),]
    final.res[out.i,7] = length(unique(CTO2[,"geneID"]))
    
    CpG = readRDS(paste("res3_signle_CpG_Z_",indx,".rds",sep=""))
    CpG = CpG[CpG[,"Feature"] %in% c("GeneBody","Promoter"),]
    n.CpG = length(CpG[,"CpG"])
    CpG = CpG[!is.na(CpG[,"SUM"]),]
    CpG = CpG[CpG[,"SUM"]<0.05/n.CpG,]
    final.res[out.i,4] = length(unique(CpG[,"geneID"]))
    
    CpG = CpG[CpG[,"geneID"]%in% gene.list,]
    n.CpG = length(CpG[,"CpG"])
    CpG = CpG[CpG[,"SUM"]<0.05/n.CpG,]
    final.res[out.i,8] = length(unique(CpG[,"geneID"]))
    
    out.i = out.i + 1
}

rownames(final.res) = all
saveRDS(final.res,"summary_res.rds")

final.res = readRDS("summary_res.rds")
colSums(final.res)

colSums(final.res)[3] / colSums(final.res)[1]
colSums(final.res)[3] / colSums(final.res)[2]
colSums(final.res)[3] / colSums(final.res)[4]

out = final.res[,c(2,3)]
out = out[rowSums(out)>0,]
wilcox.test(out[,2] - out[,1], alternative = "greater")

out = final.res[,c(1,3)]
out = out[rowSums(out)>0,]
wilcox.test(out[,2] - out[,1], alternative = "greater")



wilcox.test(out[,1], out[,2], paired = TRUE, alternative = "less")


colnames(final.res) = c("MultiXcan","TWAS","CMO","CpG","commonSet_MultiXcan","commonSet_TWAS","commonSet_CMO","commonSet_CpG")

out = as.data.frame(matrix(NA,2,5))
out1 = final.res[,1:4]
out1 = colSums(out1)
out[1,1:4] = out1
out[,5] = "ALL available genes"
out2 = final.res[,5:8]
out2 = colSums(out2)

out[2,1:4] = out2
out[2,5] = "Common set of 11,708 genes"

colnames(out) = c("MultiXcan","Union of TWAS","CMO","Union of MWAS","setting")


library(reshape)
final.dat3 <- melt(out, id=c("setting"))

library(ggplot2)
library("ggsci")
library("ggplot2")
library("gridExtra")


# generate figure
colnames(final.dat3) <- c("Setting", "Methods", "Power")

final.dat3[, 1] <- as.factor(final.dat3[, 1])
final.dat3[, 2] <- factor(final.dat3[, 2],levels=c("Union of TWAS","MultiXcan","CMO","Union of MWAS") )

# Function to produce summary statistics (mean and +/- sd)
dodge <- position_dodge(width = 1.2)


ggplot(data = final.dat3, aes(x = Setting, y = Power, fill = Methods)) +geom_col(width=0.75,position = position_dodge(width=0.9))+ ylab("Number of significant genes") + xlab("") + facet_grid(. ~ Setting, scales = "free_x") + theme(text = element_text(size=15),strip.background = element_blank(),strip.text.x = element_blank()) +scale_fill_npg() + scale_color_npg()

outd = "/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/"
fig.name = "fig2_UKBiobankImage"
ggsave(paste(outd,fig.name, ".pdf", sep = ""))


