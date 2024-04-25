source("/home/genis/impala/info_files/loadPopInfo2.R")
#library(RcppCNPy)
#library(ape)


indir <- "/home/genis/impala/analyses/impalaMap/pcangsd/results"
outdir <- "/home/genis/impala/analyses/impalaMap/pcangsd/plots"
make_path <- function(file) paste(indir, file, sep="/")

# PCA PLOTS
cov_mat <- as.matrix(read.table(make_path("impalaImpalaMappedRepeatMapInbDepFiltersNo1BadNo5dupse7.cov")))
ei <- eigen(cov_mat)
e <- ei$vectors
v <- ei$values
vars <- v/sum(v)*100

outpng <- "/home/genis/impala/paperplots/figure1/pca.png"
bitmap(outpng, width=4, height=4, res=300)
par(mar=par("mar") + c(0,1.5,-2,-0.5))
plot(-e[,1], e[,2], pch=21, cex=4, bg=impala_colors[info$Locality], xlab=paste0("PC 1 (",round(vars[1], 2),"%)"), ylab=paste0("PC 2 (",round(vars[2], 2),"%)"),cex.lab=2.5, cex.axis=2)
#legend("bottomleft",
#       legend=names(impala_colors), pt.bg=impala_colors,
#       pch=21,cex=1.5, pt.cex=2,bty='n')
dev.off()
