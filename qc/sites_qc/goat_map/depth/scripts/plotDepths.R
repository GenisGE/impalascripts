
highdepth <- "/home/genis/impala/localSitesQC/depth/output/highDepth/%s_highDepth.depthGlobal"
lowdepth <- "/home/genis/impala/localSitesQC/depth/output/lowDepth/%s_lowDepth.depthGlobal"
#alldepth <- "/home/genis/impala/localSitesQC/depth/output/allDepth/%s_allDepth.depthGlobal"

chrs <- scan("/home/genis/impala/info_files/goatChrsAutoandX.txt", what="dl")

ldall <- numeric(5001)
hdall <- numeric(5001)
#adall <- numeric(5001)

for(chr in chrs){
    cat("reading chr ", chr, "\n")
    ld <- scan(sprintf(lowdepth,chr),what=2)
    hd <- scan(sprintf(highdepth,chr),what=2)
#    ad <- scan(sprintf(alldepth,chr),what=2)
     
    ldall <- ld + ldall
    hdall <- hd + hdall
#    adall <- ad + adall

    ld <- c(ld[1:500], sum(ld[501:5001]))
    hd <- c(hd[1:500], sum(hd[501:5001]))
#    ad <- c(ad[1:500], sum(ad[501:5001]))

    yl <- c(0, max(c(max(hd), max(ld)))

    png(paste0("/home/genis/impala/localSitesQC/depth/plots/depthDist_",chr,".png"), width=600, height=600)
    par(mfrow=c(2,1))
#    barplot(ad,names.arg=0:500, ylim=yl,space=0, border=NA, main=paste0("Depth distribution ", chr, " all samples"))
    barplot(ld,names.arg=0:500, ylim=yl, space=0, border=NA,  main=paste0("Depth distribution ", chr, " low depth samples"))
    barplot(hd,names.arg=0:500, ylim=yl,space=0, border=NA, main=paste0("Depth distribution ", chr, " high depth samples"))
    dev.off()
      
}


ld <- c(ldall[1:500], sum(ldall[501:5001]))
hd <- c(hdall[1:500], sum(hdall[501:5001]))
#ad <- c(adall[1:500], sum(adall[501:5001]))

yl <- c(0, max(c(max(hd), max(ld)))

png("/home/genis/impala/localSitesQC/depth/plots/depthDist_allChromosomes.png", width=600, height=600)
par(mfrow=c(2,1))
#barplot(ad,names.arg=0:500, ylim=yl,space=0, main=paste0("Depth distribution all samples"))
barplot(ld,names.arg=0:500, ylim=yl,space=0, main=paste0("Depth distribution low depth samples"))
barplot(hd,names.arg=0:500, ylim=yl,space=0, main=paste0("Depth distribution high depth samples"))
dev.off()
