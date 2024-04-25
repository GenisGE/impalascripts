library(data.table)

args <- commandArgs(trailingOnly=T)

indir <- args[1]
outpng <- args[2]

# load r2 values
f <- list.files(indir, ".r2$", full.names=T)
r2 <- numeric(0)
for(i in 1:length(f)) r2 <- c(r2,fread(f[i], data.table=F)$V2)
r2[is.nan(r2)] <- 0


# load dosages and use it to estimate allele frequency
f <- list.files(indir, ".dose.gz$", full.names=T)
freqs <- numeric(0)
for(i in 1:length(f)){
    dose <- fread(f[i], data.table=F, h=T)
    freqs <- c(freqs, rowMeans(dose[,4:ncol(dose)])/2)
}

# join, order by allele freq and make allele freq bins
dat <- cbind(freqs,r2)
dat <- dat[order(dat[,1]),]

bins <- cut(dat[,1], breaks=20)


means <- tapply(X=dat[,2], INDEX=bins, FUN=mean)
sds <- tapply(X=dat[,2], INDEX=bins, FUN=sd) # add sd to reflect some bins have only a few sites



# THIS PLOT LOOKS NICER BUT IS LESS INFORMATIVE
if(FALSE){
bitmap(outpng, w=6, h=6, res=300)
par(oma=c(4,0,0,0))
plot(means, type="l", xaxt="n", xlab="", ylab="R2", ylim=c(0,1), main="Mean R2 by allele frequency bins")
points(means, pch=16)
axis(at=1:20, labels=levels(bins), side=1, las=2, xpd=NA)
title(xlab="Allele frequency", line=7, xpd=NA)
arrows(x0=1:20, y0=means-sds, x1=1:20, y1=means+sds, length=0.1, code=3, angle=90, lwd=1)
dev.off()
}



bitmap(outpng, w=6, h=6, res=300)
par(oma=c(4,0,0,0))
boxplot(dat[,2] ~ bins, las=2, main="Mean R2 by allele frequency bins", xlab="", ylab="R2", outpch=16)
title(xlab="Allele frequency", line=7, xpd=NA)
dev.off()

