

args <- commandArgs(trailingOnly=T)

inf <- args[1]
infsplit <- gsub("\\.sfs$", "_split.sfs", inf)
outpre <- gsub("\\.sfs$", "", inf)

allsfs <- scan(inf, skip=1, what=.3)
splitsfs <- readLines(infsplit)
# hs are actually indices (of the header lines)
hs <- seq(1, length(splitsfs), 2)
# sfss are the actual sfs values not indices
# get leave-one out sfs
sfss <- t(sapply(splitsfs[hs + 1], function(x) as.numeric(unlist(strsplit(x, split=" ")))))
dimnames(sfss) <- NULL
sfss_levaoneout <- t(sapply(1:nrow(sfss), function(x) allsfs - sfss[x,]))

# hs2 are no longer indices but the text we want to write. not very consistent
hs2 <- sapply(strsplit(splitsfs[hs], split="="), function(x) paste(gsub("<|>", "", gsub("/", " ", x[2])), "unfolded"))

for(i in 1:nrow(sfss)){
    outfile <- file(paste0(outpre, "_", i, "out", ".sfs"))
    sfs <-  paste(sfss_levaoneout[i,], collapse=" ")
    writeLines(c(hs2[i], sfs), outfile)
    close(outfile)
}
