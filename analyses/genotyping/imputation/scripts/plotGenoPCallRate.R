library(data.table)

args <- commandArgs(trailingOnly=T)

indir <- args[1]
outpng <- args[2]


f <- list.files(indir, ".gprobs.gz$", full.names=T)

gprobmax <- numeric(0)
for(filename in f){
cat("Using file", filename, "\n")

df <- fread(filename, data.table=F)

n <- (ncol(df) - 3) / 3

cat("Read beagle genotype probability file, there are", n, "samples and", nrow(df), "sites.\nWill now extract maximum probability for each site\n")

idx <- split(1:(n*3), ceiling((1:(n*3))/3))

gprobmax <- c(gprobmax, apply(df[,-c(1,2,3)], 1, function(x) sapply(idx, function(y) max(x[y]))))
}

cat("Done, will now plot\n\n")

callrate <- function(x) sum(gprobmax >= x) / length(gprobmax)

probs <- c(1, 0.9999, 0.999, 0.99, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, 0.825, 0.8)

cr <- sapply(probs, callrate)

bitmap(outpng, w=6, h=6, res=300)
plot(probs, cr, xlab="Posterior probability cutoff", ylab="Call rate", type="l", main="Genotype calling rate", col="red")
points(probs, cr, pch="x", col="blue", cex=2)
dev.off()

cat("Finished, saved plot to", outpng, "\n")
