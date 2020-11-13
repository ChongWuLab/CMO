# check if the required packages have been installed
list.of.packages <- c("data.table","optparse","CompQuadForm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(library(data.table))
library(CompQuadForm)

source("dist_support.R")
source("mwas_support.R")

suppressMessages(library("optparse"))


option_list = list(
make_option("--sumstats", action="store", default=NA, type='character',
help="summary statistics (txt format) [required]"),
make_option("--out", action="store", default=NA, type='character',
help="Path to output files [required]"),
make_option("--weights_dir", action="store", default=NA, type='character',
help="Path to output files [required]"),
make_option("--chr_id", action="store", default=-1, type='character',
help="The chromosome ID. We recommend parallel the computations by chromosomes [required]")
)


opt = parse_args(OptionParser(option_list=option_list))

#opt = list(sumstats = "/gpfs/research/chongwu/Chong/MWAS/CMO/ANA_B2_V4_Eur.txt",out = "/gpfs/research/chongwu/Chong/MWAS/CMO/COVID19_B2_ALL", weights_dir = "/gpfs/research/chongwu/Chong/MWAS/WEIGHTS/",chr_id =22)

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



outd <- opt$out
system(paste("mkdir -p ", outd, sep = ""))


if (grepl(".txt",opt$sumstats)) {
    raw = fread(opt$sumstats)
    raw = as.data.frame(raw)
} else {
    raw = fread(opt$sumstats)
    raw = as.data.frame(raw)
    warning("We try to read GWAS summary data other txt format. The data may not read properly")
    
    cat("Following are the first few lines of the dataset:", sep = "\n")
    print(head(raw))
}

header.inner = colnames(raw)
# Initialize the header
header.inner <- tolower(header.inner)

# SNP
try.snp <- c("snp", "markername", "snpid", "rs", "rsid", "rs_number", "snps")
header.inner[header.inner %in% try.snp] <- "SNP"

# A1
try.a1 <- c("a1", "allele1", "allele_1", "effect_allele", "reference_allele", "inc_allele", "ea", "ref", "a1lele1", "al1ele1")
header.inner[header.inner %in% try.a1] <- "A1"

# A2
try.a2 <- c("a2", "allele2", "allele_2", "other_allele", "non_effect_allele", "dec_allele", "nea", "alt", "a0")
header.inner[header.inner %in% try.a2] <- "A2"

# Z-score
try.z <- c("zscore", "z-score", "gc_zscore", "z")
header.inner[header.inner %in% try.z] <- "Z"


try.chromosome <- c("chrom", "ch", "chr", "chromosome","#chr")
header.inner[header.inner %in% try.chromosome] <- "CHR"

# P
try.p <- c("pvalue", "p_value", "pval", "p_val", "gc_pvalue", "p")
header.inner[header.inner %in% try.p] <- "P"


# Beta
try.beta <- c("b", "beta", "effects", "effect","all_inv_var_meta_beta")
header.inner[header.inner %in% try.beta] <- "BETA"

# Odds ratio
try.or <- c("or")
header.inner[header.inner %in% try.or] <- "ODDS_RATIO"

# Log odds
try.logodds <- c("log_odds", "logor", "log_or")
header.inner[header.inner %in% try.logodds] <- "LOG_ODDS"

# Standard error
try.se <- c("se", "sebeta", "beta_se","all_inv_var_meta_sebeta")
header.inner[header.inner %in% try.se] <- "SE"

colnames(raw) <- header.inner


list.coerce <- c("Z", "BETA", "ODDS_RATIO", "LOG_ODDS", "SE")

for (i in 1:length(header.inner)) {
    if (header.inner[i] %in% list.coerce) {
        if (class(raw[, header.inner[i]]) != "numeric") {
            class(raw[, header.inner[i]]) <- "numeric"
            cat(paste0("Column ", header.inner[i], " has wrong class and has been coerced to numeric."), sep = "\n")
            cat("=============================================================================================================", sep = "\n")
        }
    }
}

# Missing z-score?
calculate.z <- FALSE
if (!("Z" %in% header.inner)) {
    warning("No Z score column, we calculate one based on the information we have")
    
    if ("BETA" %in% header.inner & "SE" %in% header.inner) {
        raw['Z'] <- raw$BETA / raw$SE
        calculate.z <- TRUE
    } else if ("ODDS_RATIO" %in% header.inner & "SE" %in% header.inner) {
        raw['Z'] <- log(raw$ODDS_RATIO) / raw$SE
        calculate.z <- TRUE
    } else if ("LOG_ODDS" %in% header.inner & "SE" %in% header.inner) {
        raw['Z'] <- raw$LOG_ODDS / raw$SE
        calculate.z <- TRUE
    } else if ("BETA" %in% header.inner & "P" %in% header.inner) {
        raw['Z'] <- sign(raw$BETA) * abs(qnorm(raw$P / 2))
        calculate.z <- TRUE
    } else if ("ODDS_RATIO" %in% header.inner & "P" %in% header.inner) {
        raw['Z'] <- sign(log(raw$ODDS_RATIO)) * abs(qnorm(raw$P / 2))
        calculate.z <- TRUE
    } else if ("LOG_ODDS" %in% header.inner & "P" %in% header.inner) {
        raw['Z'] <- sign(raw$ODDS_RATIO) * abs(qnorm(raw$P / 2))
        calculate.z <- TRUE
    } else {
        stop("I can't calculate z-score based on the information I have. SAD FACE EMOJI.", sep = "\n")
    }
    
    if (sum(is.na(raw$Z)) != 0) {
        n.start <- nrow(raw)
        
        raw <- raw[!is.na(raw$Z),]
        n.end <- nrow(raw)
        cat(paste0(n.start - n.end, " rows removed for having invalid z-score!"), sep = "\n")
    }
    cat("=============================================================================================================", sep = "\n")
}

if(sum(header.inner %in% c("SNP","A1","A2","CHR")) !=4) {
    stop("We tried our best to match the colnames of the summary data. Please revise the colnames to the following format. You can also report in GitHub to make the lists more comphrensive:\n SNP: snp, markername, snpid, rs, rsid, rs_number, snps\n A1: a1, allele1, allele_1, effect_allele, reference_allele, inc_allele, ea, ref, a1lele1, al1ele1\n A2: a2, allele2, allele_2, other_allele, non_effect_allele, dec_allele, nea, alt, a0\n Z: zscore, z-score, gc_zscore, z\n CHR: chrom, ch, chr, chromosome")
}

sumstat.orgin = raw[,c("CHR","SNP","A1","A2","Z")]
#rm(raw)

if(opt$chr_id<=22 & opt$chr_id>=1) {
    job = opt$chr_id
    wgtlist <- wgtlist[wgtlist[, 1] == job,]
    sumstat.orgin <- sumstat.orgin[sumstat.orgin$CHR == job,]
    
} else if( opt$chr_id!= -1) {
    warning("Chromosomes ID must between 1 and 22. We ignore this argument and run the analyses for all the genes.")
}
used.gene <- unique(wgtlist$geneID)


# Prepare the data
sumstat.orgin$A1 <- toupper(sumstat.orgin$A1)
sumstat.orgin$A2 <- toupper(sumstat.orgin$A2)

a1 = sumstat.orgin[, "A1"]
a2 = sumstat.orgin[, "A2"]

keep = !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))

keep2 = nchar(a1)==1 &nchar(a2)==1

keep = keep & keep2
sumstat.orgin = sumstat.orgin[keep,]
gwas_snp = sumstat.orgin[, "SNP"] # rs ID gwas SNP


load(paste(opt$weights_dir,"/NET_Methylation_download.RData",sep="")) # loads: LD, NET, info

NET_all = do.call(rbind, NET) # make matrix from list NET
NET_all = NET_all[which(!duplicated(NET_all[, 2])),] # only keep unique SNPs
m = match(NET_all[, 2], gwas_snp);
NET_all = NET_all[which(!is.na(m)),] # only keep SNPs present in GWAS to speed up matching
m = m[which(!is.na(m))]

sumstat.orgin = sumstat.orgin[m,]
sum(is.na(sumstat.orgin))

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

rm(sumstat.orgin)
rm(LD)
rm(NET)
rm(NET_all)


out.res <- as.data.frame(matrix(NA, length(used.gene), 13))
colnames(out.res) <- c("CHR", "geneID", "ensembl", "P0", "P1", "n_enhancer", "n_enhancer_CpG", "n_gene_CpG", "n_enhancer_SNP"
, "n_gene_SNP", "n_total_SNP", "CMO", "runtime(s)")

res.single <- NULL

for (gene.indx in 1:length(used.gene)) {
    #length(used.gene)
    tryCatch({
        #gene.indx = 10
        # Load in summary stats
        start.time <- proc.time()[3]
        wgtlist1 <- wgtlist[wgtlist[, "geneID"] == used.gene[gene.indx],]
        
        res <- SMWAS.sum(wgtlist1, opt) #need revise here
        
        res.save <- res$final.out
        
        n.enhancer.snp = length(unique(res$snp.enhancer))
        n.gene.snp = length(unique(res$snp.gene))
        n.snp = length(unique(c(res$snp.enhancer, res$snp.gene)))
        
        tmp = wgtlist1[wgtlist1[, "Feature"] == "GeneBody",]
        res.save$geneID <- tmp[1, "geneID"]
        res.save$ensembl <- tmp[1, "ensembl"]
        
        res.single = rbind(res.single, res.save)
        
        
        res.tmp = res.save[res.save[, "Feature"] == "GeneBody",]
        n.gene.cpg = dim(res.tmp)[1]
        
        res.tmp = res.save[res.save[, "Feature"] != "GeneBody",]
        n.enhancer.cpg = dim(res.tmp)[1]
        
        # combine adaptive test results
        acat.tmP <- res.save[, "aSumSSUACAT"]
        acat.tmP <- acat.tmP[!is.na(acat.tmP)]
        weight <- rep(1 / length(acat.tmP), length(acat.tmP))
        cross_ACAT_p <- ACAT(acat.tmP, weight)
        
        ###########################
        ### save information    ###
        ###########################
        tmp = wgtlist1[wgtlist1[, "Feature"] == "GeneBody",]
        n.enhancer = dim(wgtlist1)[1] - 1
        basic.inf <- c(tmp[1, "CHR"], as.character(tmp[1, "geneID"]), as.character(tmp[1, "ensembl"]), tmp[1, "P0"], tmp[1, "P1"], n.enhancer)
        out.res[gene.indx, 1:6] <- basic.inf
        
        out.res[gene.indx, 7:11] = c(n.enhancer.cpg, n.gene.cpg, n.enhancer.snp, n.gene.snp, n.snp)
        out.res[gene.indx, 12] = cross_ACAT_p
        
        ensembl.id = unique(wgtlist1[,"ensembl"])
        ensembl.id = ensembl.id[!is.na(ensembl.id)]
        
        cat("Finish Gene",gene.indx," ",ensembl.id,"\n")
        
        end.time <- proc.time()[3]
        run.time <- (end.time - start.time)
        out.res[gene.indx, 13] <- run.time
        
    }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
        
        tmp = wgtlist1[wgtlist1[, "Feature"] == "GeneBody",]
        n.enhancer = dim(wgtlist1)[1] - 1
        basic.inf <- c(tmp[1, "CHR"], as.character(wgtlist1[1, "geneID"]), as.character(wgtlist1[1, "ensembl"]), wgtlist1[1, "P0"], wgtlist1[1, "P1"], n.enhancer)
        out.res[gene.indx, 1:6] <- basic.inf
    })
}

out.file <- paste(outd, "/res_CMO_CHR", job, ".txt", sep = "")
write.table(out.res,out.file,quote=FALSE, row.names=FALSE,col.names=TRUE)
#saveRDS(out.res, out.file)

out.file <- paste(outd, "/res_MWAS_CHR", job, ".txt", sep = "")
write.table(res.single,out.file,quote=FALSE, row.names=FALSE,col.names=TRUE)



