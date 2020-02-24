##ENSEMBL_GRch37_gene_list.txt:start, end : ChR38;  P0hg19    P1hg19: the transformed one, use this 

library(data.table)
setwd("/home/chong/data/MWAS/")
library(biomaRt)


#############################################
# Process ensembl gene data
##############################################
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")
chrY_genes <- getBM(attributes = c('ensembl_gene_id', 'gene_biotype', 'hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'), filters = 'chromosome_name', values = "1", mart = ensembl)

gene.list <- chrY_genes
for (chr in 2:22) {
  chrY_genes <- getBM(attributes = c('ensembl_gene_id', 'gene_biotype', 'hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'), filters = 'chromosome_name', values = chr, mart = ensembl)
  tmp <- chrY_genes #[chrY_genes[,"gene_biotype"] == "protein_coding",]
  gene.list <- rbind(gene.list, tmp)
}
write.table(gene.list, "ENSEMBL_GRch37_gene_list.txt", quote = F, row.names = F)

gene.list2 = gene.list[gene.list[, 2] %in% c("pseudogene", "lincRNA", "protein_coding"),]
gene.list3 = gene.list2[gene.list2[, 3] != "",]
gene.list3$uniqind = paste(gene.list3$chromosome_name, gene.list3$start_position, gene.list3$end_position, sep = ":")
gene.list3 = gene.list3[!duplicated(gene.list3$uniqind),]

write.table(gene.list3, "ENSEMBL_GRch37_gene_list.txt", quote = F, row.names = F)




################################################
## Processed enhancer data
################################################
genehancer <- fread("genehancer.csv", sep = ",")
genehancer <- as.data.frame(genehancer)
#genehancer <- genehancer[,c(1,4,5)]
genehancer$id <- paste("enhancer", 1:dim(genehancer)[1], sep = "")

#write.table(genehancer,"enhnacer_hg38.bed",col.names=F,row.names=F,quote=F)
genehancer19 = fread("enhancer_hg19.bed")
genehancer19 <- as.data.frame(genehancer19)

final.dat <- as.data.frame(matrix(NA, 3e6, 10))


att <- genehancer[, "attributes"]
att <- strsplit(att, ";")
att <- lapply(att, function(x) x[-1])

tmp.len <- lapply(att, function(x) length(x))
tmp.len <- unlist(tmp.len)
genehancer$ngene <- tmp.len / 2

final.dat <- genehancer[rep(row.names(genehancer), genehancer$ngene), c(1:6, 10)]

att <- unlist(att)
att <- gsub("connected_gene=", "", att)
att <- gsub("score=", "", att)

final.dat$geneID = att[1:length(att) %% 2 == 1]
final.dat$geneScore = att[1:length(att) %% 2 == 0]

tmphg19 <- table(genehancer19[, 4])
tmp <- as.data.frame(matrix(NA, dim(tmphg19)[1], 2))
tmp[, 1] = names(tmphg19)
tmp[, 2] = tmphg19
rownames(tmp) = tmp[, 1]
tmp = tmp[tmp[, 2] == 1,]

genehancer19 = genehancer19[genehancer19[, 4] %in% tmp[, 1],]

common.enhancer = intersect(unique(final.dat$id), genehancer19[, 4])

final.dat = final.dat[final.dat$id %in% common.enhancer,]
rownames(genehancer19) = genehancer19[, 4]

genehancer19 = genehancer19[final.dat$id,]

final.dat$P0hg19 = genehancer19[, 2]
final.dat$P1hg19 = genehancer19[, 3]

final.dat[, 9] = as.numeric(final.dat[, 9])
tmp.name = colnames(final.dat)
tmp.name[3] = "feature_name"
colnames(final.dat) = tmp.name

saveRDS(final.dat, "processed_enhancer_inf.rds")
write.table(final.dat, "processed_enhancer_inf.txt", row.names = F, quote = F)

final.dat2 = final.dat[final.dat[, "geneScore"] > 6,] n

final.dat2 = final.dat2[final.dat[, "score"] > 1,]




load("NET_Methylation_download.RData")

info = as.data.frame(info)
info[, 2:5] = apply(info[, 2:5], 2, as.numeric)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(org.Hs.eg.db)
library(limma)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- as.data.frame(ann)

saveRDS(ann, "all_IlluminaHumanMethylation450kanno.rds")
tmp <- ann[ann[, "Name"] %in% info[, 1],]

saveRDS(tmp, "IlluminaHumanMethylation450kanno.rds")
gene.name = tmp[, "UCSC_RefGene_Accession"]
gene.name = strsplit(gene.name, ";")
gene.name = unlist(gene.name)















sum(tmp[, "UCSC_RefGene_Accession"] == "")

tmp[1:100, "UCSC_RefGene_Name"]

methyl.epic <- fread("MethylationEPIC_v-1-0_B4.csv")
methyl.epic <- as.data.frame(methyl.epic)

tmp <- methyl.epic[methyl.epic[, "Name"] %in% info[, 1],]

# data can be downloaded from http://genome.ucsc.edu/cgi-bin/hgTables
