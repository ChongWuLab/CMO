library(xtable)
library(cowplot)
library(xtable)
library(qqman)
library(RColorBrewer)
library(ggrepel)
library(lattice)


#########################################################################################
# Code is based on https://github.com/timknut/gg.manhattan/blob/master/gg.manhattan.R
# Author: Chong Wu (wuxx0845@umn.edu)
# Added function: mark some significant genes.
# If you share this codes, please let me know.
##########################################################################################
gg.manhattan <- function(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", annot = "annot", trait = "trait",
cols = c("gray10", "gray60"), plot_title = NULL, chrlabs = NULL,
suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8),
logp = TRUE) {
    
    
    # Require ggplot2, install if not present in labrary.
    if (!require(ggplot2)) install.packages("ggpplot2")
    
    # Not sure why, but package check will warn without this.
    CHR <- BP <- P <- index <- NULL
    
    # Check for sensible dataset
    ## Make sure you have chr, bp and p columns.
    if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
    ## warn if you don't have a snp column
    if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!(trait %in% names(x))) warning(paste("No trait column found. OK unless you're trying to plot multiple traits."))
    ## make sure chr, bp, and p columns are numeric.
    if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
    
    # Create a new data.frame with columns called CHR, BP, and P.
    d <- data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], SNP = x[[snp]], ANNOT = x[[annot]])
    
    # If the input data frame has a SNP column, add it to the new data frame you're creating.
    # if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
    # If the input data frame has a trait column, add it to the new data frame you're creating.
    if (!is.null(x[[trait]])) d <- transform(d, TRAIT = x[[trait]])
    
    # Set positions, ticks, and labels for plotting
    ## Sort and keep only values where is numeric.
    # d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP),]
    # d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
    if (logp) {
        d$logp <- -log10(d$P)
    } else {
        d$logp <- d$P
    }
    d$pos <- NA
    
    
    # Fixes the bug where one chromosome is missing by adding a sequential index column.
    d$index <- NA
    ind <- 0
    for (i in unique(d$CHR)) {
        ind <- ind + 1
        d[d$CHR == i,]$index <- ind
    }
    
    # This section sets up positions and ticks. Ticks should be placed in the
    # middle of a chromosome. The a new pos column is added that keeps a running
    # sum of the positions of each successive chromsome. For example:
    # chr bp pos
    # 1   1  1
    # 1   2  2
    # 2   1  3
    # 2   2  4
    # 3   1  5
    nchr <- length(unique(d$CHR))
    if (nchr == 1) {
        ## For a single chromosome
        ## Uncomment the next line to plot single chr results in Mb
        d$pos <- d$BP / 1e6
        ticks <- floor(length(d$pos)) / 2 + 1
        xlabel <- paste("Chromosome", unique(d$CHR), "position(Mb)")
        labs <- ticks
    } else {
        ## For multiple chromosomes
        lastbase <- 0
        ticks <- NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i,]$pos <- d[d$index == i,]$BP
            } else {
                lastbase <- lastbase + tail(subset(d, index == i - 1)$BP, 1)
                d[d$index == i,]$pos <- d[d$index == i,]$BP + lastbase
            }
            # Old way: assumes SNPs evenly distributed
            # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
            # New way: doesn't make that assumption
            ticks <- c(ticks, (min(d[d$CHR == i,]$pos) + max(d[d$CHR == i,]$pos)) / 2 + 1)
        }
        xlabel <- "Genomic position"
        ## Set x labels to Chromsomes, or provided custom list from chrlabs argument.
        if (!is.null(chrlabs)) {
            labs <- chrlabs
        } else {
            labs <- unique(d$CHR)
        }
    }
    
    # Initialize plot
    xmax <- ceiling(max(d$pos) * 1.01)
    xmin <- floor(max(d$pos) * -0.01)
    ymin <- floor(min(d$logp))
    ymax <- ceiling(max(d$logp + 1))
    # Set ymax to even number.
    if ((ymax %% 2) != 0) {
        ymax <- ymax + 1
    }
    mycols <- rep(cols, nchr / 2 + 1)
    if (nchr == 1) {
        plot <- ggplot(
        data = d, aes(pos, logp), ylab = expression(-log[10](italic(p))),
        xlab = xlabel
        ) + geom_point(size = 0.5)
        if (!is.null(d[[trait]]) && length(unique(d$TRAIT)) >= 2) {
            plot <- plot + geom_point(aes(colour = TRAIT))
        }
        plot <- plot + scale_y_continuous(
        breaks = seq(2, ymax, 2), labels = seq(2, ymax, 2),
        limits = c(ymin - 0.5, ymax), expand = c(0, 0)
        )
        plot <- plot + geom_text_repel(data = subset(d, ANNOT == 1), aes(label = SNP), segment.color = "#3d3d3d", size = 4, arrow = arrow(length = unit(0.05, "npc")) ,    nudge_y      = 0.05,
        direction    = "y")
    } else {
        plot <- ggplot(d, aes(x = pos, y = logp))
        plot <- plot + geom_point(aes(colour = factor(CHR)),size = 0.7)
        plot <- plot + scale_x_continuous(name = xlabel, breaks = ticks, labels = labs)
        plot <- plot + geom_text_repel(data = subset(d, ANNOT == 1), aes(label = SNP), segment.color = "#636363", size = 4, arrow = arrow(length = unit(0.01, "npc")),nudge_y      = 0.01,
        direction    = "x")
        
        if (logp) {
            plot <- plot + scale_y_continuous(
            breaks = seq(2, ymax, 2), labels = seq(2, ymax, 2),
            limits = c(ymin - 0.5, ymax), expand = c(0, 0)
            )
        }
        # expand ensures pretty y-axis
        plot <- plot + scale_colour_manual(values = mycols)
    }
    # if (annotate)   plot=plot + geom_point(data=d.annotate, colour=I("green3"))
    if (!is.null(plot_title)) plot <- plot + ggtitle(plot_title)
    plot <- plot + theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(),
    axis.line.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(vjust = 0.5, angle = -90) # size=14,
    ) + ylab(expression(paste(-log[10], "(P value)")))
    if (suggestiveline) plot <- plot + geom_hline(yintercept = suggestiveline, colour = "black", linetype = "dashed")
    # if (genomewideline) plot <- plot + geom_hline(yintercept = genomewideline, colour = "red")
    plot
}


setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/AD_res/")
library(VennDiagram)

crosRes = readRDS("res3_cross_CpG_AD_Jansene.rds")
single = readRDS("res_signle_tissue_AD_Jansene_1000G_SPrediXcan.rds")

#crosRes = readRDS("res_cross_tissue_AD_Jansene_1000G_SPrediXcan.rds")
#single = readRDS("res_signle_tissue_AD_Jansene_1000G_SPrediXcan.rds")

dim(single)

crosRes = crosRes[!is.na(crosRes[, "cross_CMO"]),]
crosRes = crosRes[crosRes[, "cross_CMO"]>=0,]
crosRes = crosRes[!is.na(crosRes[,1]),]
crosRes[, "P0"] = as.numeric(crosRes[, "P0"])
crosRes[, "P1"] = as.numeric(crosRes[, "P1"])
crosRes[, "P0"] = floor((crosRes[, "P0"] + crosRes[, "P1"]) / 2)
crosRes = crosRes[, c("CHR", "geneID", "P0", "cross_CMO")]

#tmp = crosRes[crosRes[, 4] < 0.05/dim(crosRes)[1],]
#tmp.gene = tmp[, 2]
colnames(crosRes) = c("CHR", "SNP", "BP", "p-value")

crosRes[, 1] = as.numeric(crosRes[, 1])
crosRes[, 3] = as.numeric(crosRes[, 3])
crosRes[, 4] = as.numeric(crosRes[, 4])

sig_level = 0.05/dim(crosRes)[1]

sigres = crosRes[crosRes[, 4] < sig_level,]

chr.id = unique(sigres[, 1])
annot = NULL
for (i in chr.id) {
    tmp = sigres[sigres[, 1] == i,]
    tmp = tmp[order(tmp[,3]),]
    if (dim(tmp)[1] == 1) {
        annot = rbind(annot, tmp)
    } else if (tmp[dim(tmp)[1], 3] - tmp[1, 3] < 500 * 1000) {
        tmp = tmp[tmp[, 4] == min(tmp[, 4]),]
        annot = rbind(annot, tmp)
    } else {
        tmp.save = tmp[1,]
        most.sig = tmp.save[1, 4]
        bp = tmp[1, 3]
        for (j in 2:dim(tmp)[1]) {
            if (tmp[j, 3] < bp + 500 * 1000) {
                if (tmp[j, 4] < most.sig) {
                    most.sig = tmp[j, 4]
                    tmp.save = tmp[j,]
                }
                if (j == dim(tmp)[1]) {
                    annot = rbind(annot, tmp.save)
                }
                bp = tmp[j,3]
                cat(j,"\n")
            } else {
                annot = rbind(annot, tmp.save)
                bp = tmp[j, 3]
                most.sig = tmp[j, 4]
                tmp.save = tmp[j,]
                #if (j == dim(tmp)[1]) {
                #    annot = rbind(annot, tmp.save)
                #}
            }
        }
    }
}

crosRes[crosRes[, 4] < 1e-30, 4] = 1e-30
annot.gene = annot[, 2]

annot.gene = annot.gene[!annot.gene %in% c("PVR","ZNF404")]
annot.gene = c(annot.gene, "APOE")
crosRes$annot = as.numeric(crosRes[, 2] %in% annot.gene)

gg.manhattan(crosRes, cols = c("#3B77AE", "#CD342B"), suggestiveline = -log10(sig_level), p = "p-value", plot_title = "")
ggsave("fig2_manhattan.pdf")


#else if (dim(tmp)[1] == 2) {
#annot = rbind(annot, tmp)
#}
