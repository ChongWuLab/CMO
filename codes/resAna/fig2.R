setwd("/home/chong/Dropbox/FSU_research/R21 Grant Application/MWAS/resAna")
setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/")
library(VennDiagram)

crosRes = readRDS("res_cross_tissue_IGAP1_1000G_SPrediXcan.rds")
single = readRDS("res_signle_tissue_IGAP1_1000G_SPrediXcan.rds")

#crosRes = readRDS("res_cross_tissue_AD_Jansene_1000G_SPrediXcan.rds")
#single = readRDS("res_signle_tissue_AD_Jansene_1000G_SPrediXcan.rds")


dim(single)
crosRes = crosRes[crosRes[, "ntissue"] > 0,]
crosRes[crosRes[, "ntissue"] == 1, "Omni_ACAT_MultiXcan"] =
  crosRes[crosRes[, "ntissue"] == 1, "best_p-value"]

crosRes[crosRes[, "ntissue"] == 1, "MultiXcan"] =
  crosRes[crosRes[, "ntissue"] == 1, "best_p-value"]

cto = crosRes[crosRes[, "Omni_ACAT_MultiXcan"] < 2.1e-6 & !is.na(crosRes[, "Omni_ACAT_MultiXcan"]),]
cto.gene = cto[, "Gene"]
cto.gene = as.character(cto.gene)
dim(cto)

multixcan = crosRes[crosRes[, "MultiXcan"] < 2.1e-6 & !is.na(crosRes[, "MultiXcan"]),]
#multixcan = multixcan[multixcan[, "best_p-value"] < 1e-4,]
dim(multixcan)
multixcan.gene = multixcan[, "Gene"]
multixcan.gene = as.character(multixcan.gene)

twas = single[single[, "p-value"] < 2.4e-7,]

twas.gene = unique(twas[, "Gene_name"])
twas.gene = as.character(twas.gene)

length(twas.gene)

cpgRes = readRDS("res2_cross_CpG_IGAP1.rds")
cpgRes = cpgRes[!is.na(cpgRes[,"cross_ACAT"]),]
cpgRes = cpgRes[cpgRes[,"cross_ACAT"]>=0,]
cpgRes = cpgRes[cpgRes[,"cross_ACAT"]<2.9e-6,]
cpgRes.gene = cpgRes[,2]

length(cpgRes.gene)
pdf("fig2_venn.pdf")

venn.plot <- draw.quad.venn(
        area1 = length(cto.gene),
        area2 = length(multixcan.gene),
        area3 = length(twas.gene),
        area4 = length(cpgRes.gene),
        n12 = sum(cto.gene %in% multixcan.gene),
        n13 = sum(cto.gene %in% twas.gene),
        n14 = sum(cto.gene %in% cpgRes.gene),
        n23 = sum(multixcan.gene %in% twas.gene),
        n24 = sum(multixcan.gene %in% cpgRes.gene),
        n34 = sum(twas.gene %in% cpgRes.gene),
        n123 = length(intersect(intersect(cto.gene, multixcan.gene), twas.gene)),
        n124 = length(intersect(intersect(cto.gene, multixcan.gene), cpgRes.gene)),
        n134 = length(intersect(intersect(cto.gene, twas.gene), cpgRes.gene)),
        n234 = length(intersect(intersect(multixcan.gene, twas.gene), cpgRes.gene)),
        n1234 = length(intersect(intersect(intersect(cto.gene, multixcan.gene), twas.gene), cpgRes.gene)),
        category = c("CTO", "MultiXcan", "TWAS(All)", "CMO"),
        lty = "blank",
        cex = 2,
        margin = 0.1,
        fill = c(
            "skyblue", "pink1",
            "mediumorchid", "orange"
        )
    )
dev.off()
