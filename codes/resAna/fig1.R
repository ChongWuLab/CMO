setwd("/home/chong/Dropbox/FSU_research/R21 Grant Application/MWAS/resAna")
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

tmp = single[single[, 1] == "Brain_Cerebellum",]
best.single[i, 3] = dim(tmp)[1]

tmp = tmp[tmp[, "p-value"] < 0.05 / dim(tmp)[1],]
tmp = unique(tmp[, "Gene_name"])
besttissue.gene = as.character(tmp)

pdf("fig1_venn.pdf")

venn.plot <- draw.quad.venn(
        area1 = length(cto.gene),
        area2 = length(multixcan.gene),
        area3 = length(twas.gene),
        area4 = length(besttissue.gene),
        n12 = sum(cto.gene %in% multixcan.gene),
        n13 = sum(cto.gene %in% twas.gene),
        n14 = sum(cto.gene %in% besttissue.gene),
        n23 = sum(multixcan.gene %in% twas.gene),
        n24 = sum(multixcan.gene %in% besttissue.gene),
        n34 = sum(twas.gene %in% besttissue.gene),
        n123 = length(intersect(intersect(cto.gene, multixcan.gene), twas.gene)),
        n124 = length(intersect(intersect(cto.gene, multixcan.gene), besttissue.gene)),
        n134 = length(intersect(intersect(cto.gene, twas.gene), besttissue.gene)),
        n234 = length(intersect(intersect(multixcan.gene, twas.gene), besttissue.gene)),
        n1234 = length(intersect(intersect(intersect(cto.gene, multixcan.gene), twas.gene), besttissue.gene)),
        category = c("CTO", "MultiXcan", "TWAS(All)", "TWAS(Brain Cerebellum)"),
        lty = "blank",
        cex = 2,
        margin = 0.1,
        fill = c(
            "skyblue", "pink1",
            "mediumorchid", "orange"
        )
    )
dev.off()




cis = readRDS("summary_best_cis_IGAP1.rds")
cis[, 3:7] = apply(cis[, 3:7], 2, as.numeric)

cis2 = cis[cis[, "Gene"] %in% cto.gene,]

prop = matrix(NA, 4, 2)
colnames(prop) = c("CTO", "MultiXcan", "TWAS", "TWAS-brain")

cis2 = cis[cis[, "Gene"] %in% cto.gene,]
prop[1, 1] = length(cto.gene)
prop[1, 2] = sum(cis2[, "most_sig_500Kb"] > 5e-8)

cis2 = cis[cis[, "Gene"] %in% multixcan.gene,]
prop[2, 1] = length(multixcan.gene)
prop[2, 2] = sum(cis2[, "most_sig_500Kb"] > 5e-8)

cis2 = cis[cis[, "Gene"] %in% multixcan.gene,]
prop[3, 1] = length(multixcan.gene)
prop[3, 2] = sum(cis2[, "most_sig_500Kb"] > 5e-8)

cis2 = cis[cis[, "Gene"] %in% besttissue.gene,]
prop[4, 1] = length(besttissue.gene)
prop[4, 2] = sum(cis2[, "most_sig_500Kb"] > 5e-8)









panel = unique(single[, 1])
panel = as.character(panel)
best.single = as.data.frame(matrix(NA, length(panel), 3))
#single = readRDS("res_signle_tissue_IGAP1_1000G_SPrediXcan.rds")

for (i in 1:length(panel)) {
  tmp = single[single[, 1] == panel[i],]
  best.single[i, 3] = dim(tmp)[1]

  tmp = tmp[tmp[, "p-value"] < 0.05 / dim(tmp)[1],]
  tmp = unique(tmp[, "Gene_name"])
  best.single[i, 1] = panel[i]
  best.single[i, 2] = length(tmp)
}
