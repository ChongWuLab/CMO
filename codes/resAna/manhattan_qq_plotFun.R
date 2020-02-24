
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
gg.manhattan <- function(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", trait = "trait",
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
    d <- data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], SNP = x[[snp]])

    # If the input data frame has a SNP column, add it to the new data frame you're creating.
    # if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
    # If the input data frame has a trait column, add it to the new data frame you're creating.
    if (!is.null(x[[trait]])) d <- transform(d, TRAIT = x[[trait]])

    # Set positions, ticks, and labels for plotting
    ## Sort and keep only values where is numeric.
    # d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
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
        d[d$CHR == i, ]$index <- ind
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
    if (nchr == 1) { ## For a single chromosome
        ## Uncomment the next line to plot single chr results in Mb
        d$pos <- d$BP / 1e6
        ticks <- floor(length(d$pos)) / 2 + 1
        xlabel <- paste("Chromosome", unique(d$CHR), "position(Mb)")
        labs <- ticks
    } else { ## For multiple chromosomes
        lastbase <- 0
        ticks <- NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos <- d[d$index == i, ]$BP
            } else {
                lastbase <- lastbase + tail(subset(d, index == i - 1)$BP, 1)
                d[d$index == i, ]$pos <- d[d$index == i, ]$BP + lastbase
            }
            # Old way: assumes SNPs evenly distributed
            # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
            # New way: doesn't make that assumption
            ticks <- c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR == i, ]$pos)) / 2 + 1)
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
    xmax <- ceiling(max(d$pos) * 1.02)
    xmin <- floor(max(d$pos) * -0.02)
    ymin <- floor(min(d$logp))
    ymax <- ceiling(max(d$logp))
    # Set ymax to even number.
    if ((ymax %% 2) != 0) {
        ymax <- ymax + 1
    }
    mycols <- rep(cols, nchr / 2 + 1)
    if (nchr == 1) {
        plot <- ggplot(
            data = d, aes(pos, logp), ylab = expression(-log[10](italic(p))),
            xlab = xlabel
        ) + geom_point(size = 2)
        if (!is.null(d[[trait]]) && length(unique(d$TRAIT)) >= 2) {
            plot <- plot + geom_point(aes(colour = TRAIT))
        }
        plot <- plot + scale_y_continuous(
            breaks = seq(2, ymax, 2), labels = seq(2, ymax, 2),
            limits = c(ymin - 0.5, ymax), expand = c(0, 0)
        )
        plot <- plot + geom_text_repel(data = subset(d, logp > suggestiveline), aes(label = SNP), segment.color = "#3d3d3d", size = 4, arrow = arrow(length = unit(0.01, "npc")))
    } else {
        plot <- ggplot(d, aes(x = pos, y = logp))
        plot <- plot + geom_point(aes(colour = factor(CHR)))
        plot <- plot + scale_x_continuous(name = xlabel, breaks = ticks, labels = labs)
        plot <- plot + geom_text_repel(data = subset(d, logp > suggestiveline), aes(label = SNP), segment.color = "#636363", size = 4, arrow = arrow(length = unit(0.01, "npc")))

        if (logp) {
            plot <- plot + scale_y_continuous(
                breaks = seq(2, ymax, 2), labels = seq(2, ymax, 2),
                limits = c(ymin - 0.5, ymax), expand = c(0, 0)
            )
        } # expand ensures pretty y-axis
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



qqunif.plot <- function(pvalues,
                        should.thin = T, thin.obs.places = 4, thin.exp.places = 4,
                        xlab = expression(paste("Expected (", -log[10], " p-value)")),
                        ylab = expression(paste("Observed (", -log[10], " p-value)")),
                        draw.conf = TRUE, conf.points = 1000, conf.col = "lightgray", conf.alpha = .05,
                        already.transformed = FALSE, pch = 20, ...) {
    pvalues <- my.pvalue.list
    should.thin <- T
    thin.obs.places <- 2
    thin.exp.places <- 2
    xlab <- expression(paste("Expected (", -log[10], " p-value)"))
    ylab <- expression(paste("Observed (", -log[10], " p-value)"))
    already.transformed <- FALSE
    # error checking
    if (length(pvalues) == 0) stop("pvalue vector is empty, can't draw plot")
    if (!(class(pvalues) == "numeric" ||
        (class(pvalues) == "list" && all(sapply(pvalues, class) == "numeric")))) {
        stop("pvalue vector is not numeric, can't draw plot")
    }
    if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
    if (already.transformed == FALSE) {
        if (any(unlist(pvalues) == 0)) stop("pvalue vector contains zeros, can't draw plot")
    } else {
        if (any(unlist(pvalues) < 0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
    }

    grp <- NULL
    n <- 1
    exp.x <- c()
    if (is.list(pvalues)) {
        nn <- sapply(pvalues, length)
        rs <- cumsum(nn)
        re <- rs - nn + 1
        n <- min(nn)
        if (!is.null(names(pvalues))) {
            grp <- factor(rep(names(pvalues), nn), levels = names(pvalues))
            names(pvalues) <- NULL
        } else {
            grp <- factor(rep(1:length(pvalues), nn))
        }
        pvo <- pvalues
        pvalues <- numeric(sum(nn))
        exp.x <- numeric(sum(nn))
        for (i in 1:length(pvo)) {
            if (!already.transformed) {
                pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
                exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method = "first") - .5) / nn[i])
            } else {
                pvalues[rs[i]:re[i]] <- pvo[[i]]
                exp.x[rs[i]:re[i]] <- -log10((nn[i] + 1 - rank(pvo[[i]], ties.method = "first") - .5) / (nn[i] + 1))
            }
        }
    } else {
        n <- length(pvalues) + 1
        if (!already.transformed) {
            exp.x <- -log10((rank(pvalues, ties.method = "first") - .5) / n)
            pvalues <- -log10(pvalues)
        } else {
            exp.x <- -log10((n - rank(pvalues, ties.method = "first") - .5) / n)
        }
    }

    # reduce number of points to plot
    if (should.thin == T) {
        if (!is.null(grp)) {
            thin <- unique(data.frame(
                pvalues = round(pvalues, thin.obs.places),
                exp.x = round(exp.x, thin.exp.places),
                grp = grp
            ))
            grp <- thin$grp
        } else {
            thin <- unique(data.frame(
                pvalues = round(pvalues, thin.obs.places),
                exp.x = round(exp.x, thin.exp.places)
            ))
        }
        pvalues <- thin$pvalues
        exp.x <- thin$exp.x
    }

    base_size <- 12
    base_family <- ""

    input.dat <- as.data.frame(pvalues)
    input.dat$expX <- exp.x
    input.dat$Methods <- grp
    p <- ggplot(input.dat, aes(expX, pvalues))
    p <- p + geom_point(aes(colour = Methods, shape = Methods)) + geom_abline(intercept = 0, slope = 1, alpha = 0.5) + xlab(xlab) + ylab(ylab) +  theme_grey(base_size = base_size, base_family = base_family) %+replace%
        theme(
            legend.position = c(0.2, .618), legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 12), legend.text.align = 0, axis.text = element_text(size = rel(1)), axis.ticks = element_line(colour = "black"),
            legend.key = element_rect(colour = "grey80"), panel.background = element_rect(
                fill = "white",
                colour = NA
            ), panel.border = element_rect(
                fill = NA,
                colour = "grey50"
            ), panel.grid.major = element_line(
                colour = "grey90",
                size = 0.2
            ), panel.grid.minor = element_line(
                colour = "grey98",
                size = 0.5
            ), strip.background = element_rect(
                fill = "grey80",
                colour = "grey50", size = 0.2
            ), plot.title = element_text(hjust = 0.5, size = 16)
        )
    return(p)

    # prepare the data and then use ggplot to draw
}
