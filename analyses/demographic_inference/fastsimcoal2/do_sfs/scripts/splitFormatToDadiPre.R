

args <- commandArgs(trailingOnly=T)

inf <- args[1]
outpre <- gsub("\\.sfs$", "", inf)

allsfs <- readLines(inf)

# hs are actually indices (of the header lines)
hs <- seq(1, length(allsfs), 2)
# sfss are the actual sfs values not indices
sfss <- allsfs[hs + 1]

# hs2 are no longer indices but the text we want to write. not very consistent
hs2 <- sapply(strsplit(allsfs[hs], split="="), function(x) paste(gsub("<|>", "", gsub("/", " ", x[2])), "unfolded"))

for(i in 1:length(sfss)){
    outfile <- file(paste0(outpre, "_", i, ".sfs"))
    writeLines(c(hs2[i], sfss[i]), outfile)
    close(outfile)
}
