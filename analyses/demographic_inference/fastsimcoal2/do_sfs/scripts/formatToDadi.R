args <- commandArgs(trailingOnly=TRUE)


inf <- args[1]
outf <- gsub("\\.sfs$", "_dadiformat.sfs", inf)

sfsfile <- readLines(inf)

h <- sfsfile[1]
sfs <- sfsfile[2]

h2 <- sapply(strsplit(h, split="="), function(x) paste(gsub("<|>", "", gsub("/", " ", x[2])), "unfolded"))

writeLines(c(h2, sfs), outf)
