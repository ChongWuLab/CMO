setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MWAS-CMO/codes/MWAS/")
library(readr)


tmp <- read_file("job16.sh")

start.id <- 17

for (job in paste("Z_0",166:179,sep="")) {
    tmp2 <- gsub("Z_0166", job, tmp)
    write_file(tmp2, paste("job", start.id, ".sh", sep = ""))
    start.id <- start.id + 1
}


tmp <- read_file("job31.sh")

start.id <- 32

for (job in paste("Z_",2650:2655,sep="")) {
    tmp2 <- gsub("Z_2649", job, tmp)
    write_file(tmp2, paste("job", start.id, ".sh", sep = ""))
    start.id <- start.id + 1
}
