source("/home/genis/impala/info_files/loadPopInfo2.R")

args <- commandArgs(trailingOnly=T)

inpre <- args[1]
outpng <- args[2]

egfile <- paste0(inpre, ".eigenvec")
evfile <- paste0(inpre, ".eigenval")

e <- read.table(egfile)
v <- scan(evfile)
vars <- v/sum(v)*100

names(impala_colors) <- gsub(" ","", names(impala_colors))

#bitmap("/home/genis/impala/analyses/goatMapV2/pcangsd/plots/pca_all_impala12.png", width=6, height=6, res=300)
bitmap(outpng, width=6, height=6, res=300)
plot(-e$V3, e$V4, pch=21, cex=2, bg=impala_colors[as.character(e$V2)], xlab=paste0("PC 1 (",round(vars[1], 2),"%)"), ylab=paste0("PC 2 (",round(vars[2], 2),"%)"),cex.lab=1.2)
legend("bottomleft",
       legend=names(impala_colors), pt.bg=impala_colors,
       pch=21,cex=1.5, pt.cex=2,bty='n')
dev.off()
