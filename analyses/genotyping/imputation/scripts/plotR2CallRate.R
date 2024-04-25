library(data.table)

args <- commandArgs(trailingOnly=T)

indir <- args[1]
outpng <- args[2]


indir <- "/emc/genis/impala/vcf/vcf_goat/impute/"

f <- list.files(indir, ".r2$", full.names=T)

r2 <- numeric(0)
for(i in 1:length(f)) r2 <- c(r2,fread(f[i], data.table=F)$V2)
r2[is.nan(r2)] <- 0

callrate <- function(x) sum(r2 >= x) / length(r2)

probs <- c(1, 0.99, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, 0.825, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5)

cr <- sapply(probs, callrate)

bitmap(outpng, w=6, h=6, res=300)
plot(probs, cr, xlab="r2", ylab="Call rate", type="l", main="Genotype calling rate", col="red")
points(probs, cr, pch="x", col="blue", cex=2)
dev.off()
