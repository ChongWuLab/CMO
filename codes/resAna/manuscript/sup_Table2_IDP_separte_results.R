setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/resAna/UKBiobankImage_res")

all = as.character(c(paste("0",166:179,sep=""),2649:2655))

IDP.name =
c("left thalamus",
"right thalamus",
"left caudate",
"right caudate",
"left putamen",
"right putamen",
"left pallidum",
"right pallidum",
"left hippocampus",
"right hippocampus",
"left amygdala",
"right amygdala",
"left accumbens",
"right accumbens",
"left thalamus plus right thalamus",
"left caudate plus right caudate",
"left putamen plus right putamen",
"left pallidum plus right pallidum",
"left hippocampus plus right hippocampus",
"left amygdala plus right amygdala",
"left accumbens plus right accumbens","Overall")

simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
    sep="", collapse=" ")
}


IDP.name = sapply(IDP.name, simpleCap)

final.res = readRDS("summary_res.rds")
final.res = rbind(final.res,colSums(final.res))
rownames(final.res) = IDP.name

library(xtable)

# Stable2
xtable(final.res[,c(2,1,3,4)],digits = 0)

#Stable 3
xtable(final.res[,c(2+4,1+4,3+4,4+4)],digits = 0)
