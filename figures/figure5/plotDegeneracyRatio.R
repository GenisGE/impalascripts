source("/home/genis/impala/info_files/loadPopInfo2.R")

sfsfolder <- "/home/genis/impala/analyses/goatMapV3/dfe/results/sfs"


pops <- unique(gsub( " ", "", info$Locality))


d0 <- c()
d4 <- c()
for(p in pops){
    if(p == "LakeMburo") p <- "LakeMburu"
    sfs1 <- scan(list.files(sfsfolder, paste0(p, "_", 0, "f_degenerate.sfs"), full.names=T), what=.3)
    sfs2 <- scan(list.files(sfsfolder, paste0(p, "_", 4, "f_degenerate.sfs"), full.names=T), what=.3)

    d0 <- c(d0, sum(sfs1 * 0:(length(sfs1)-1)) / sum(sfs1 * (length(sfs1)-1)))
    d4 <- c(d4, sum(sfs2 * 0:(length(sfs2)-1)) / sum(sfs2 * (length(sfs2)-1)))

}


ord <- sapply(gsub(" ", "", popord), function(x) which(pops==x)) 

groups <- factor(c("Black-faced", "Southern common", "Southern common", "Southern common", "Southern common", "Southern common", "Southern common", "Selous common", "Eastern common", "Eastern common", "Eastern common", "Eastern common", "Eastern common"),
            levels = c("Black-faced", "Southern common", "Selous common", "Eastern common"))

k <- ! pops[ord] %in% c("Shangani")

outpng <- "boxplotPopRelLoad.png"
bitmap(outpng, w=4, h=4, res=300)
par(mar=c(4,7,4,3))
boxplot((d0/d4)[ord][k] ~ groups[k], outline=F, 
    ylab="# derived 4-fold degenerate alleles\nto # derived 0-fold degenerate alleles", 
    cex.lab=1.5, cex.axis=1.25,
     names= rep("", 4)) #gsub(" ", "\n", levels(groups)))
axis(1, at=1:4, gsub(" ", "\n", levels(groups)), line=1, cex.axis=1.25, tick=F)
points((d0/d4)[ord][k] ~ jitter(as.integer(groups[k]), 1),pch=21,bg=impala_colors[popord[k]], cex=3)
dev.off()

