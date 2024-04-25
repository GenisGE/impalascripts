
library(data.table)

highdepth <- "/home/genis/impala/localSitesQC/depth/output/highDepth"
lowdepth <- "/home/genis/impala/localSitesQC/depth/output/lowDepth"
 
alldepth <- "/home/genis/impala/localSitesQC/depth/output/allDepth"

tl <- scan(paste0(lowdepth,"/allChrs_lowDepth.depthGlobal"), what=4)
th <- scan(paste0(highdepth,"/allChrs_highDepth.depthGlobal"), what=43)
ta <- scan(paste0(alldepth,"/allChrs_allDepth.depthGlobal"), what=4)

mtl <- median(rep.int(0:5000,tl))
mth <- median(rep.int(0:5000,th))
mta <- median(rep.int(0:5000,ta))

tl_plot <- c(tl[1:500], sum(tl[501:5001]))
th_plot <- c(th[1:500], sum(th[501:5001]))
ta_plot <- c(ta[1:500], sum(ta[501:5001]))


# test using median

getQth <- function(d,q){
    m <- median(d)
    thres <- m*q
}

plotallthres <- function(fmin, fmax){
    
    yl <- c(0, max(ta_plot))
    par(mfrow=c(3,1))
    barplot(ta_plot, names.arg=0:500, ylim = yl,
            space=0, col="grey", ,border=NA,
            main=paste("Depth distribution all samples thresholds median times", fmin, "and", fmax))
    abline(v=mta*fmin, col=2, lty=2)
    abline(v=mta, col=1, lty=1)
    abline(v=mta*fmax, col=2, lty=2)
    
    barplot(tl_plot, names.arg=0:500, ylim = yl,
            space=0, col="grey", ,border=NA,
            main=paste("Depth distribution low depth samples thresholds median times", fmin, "and", fmax))
    abline(v=mtl*fmin, col=2, lty=2)
    abline(v=mtl, col=1, lty=1)
    abline(v=mtl*fmax, col=2, lty=2)
    
    barplot(th_plot, names.arg=0:500, ylim = yl,
            space=0, col="grey", ,border=NA,
            main=paste("Depth distribution high depth samples thresholds median times", fmin, "and", fmax))
    abline(v=mth*fmin, col=2, lty=2)
    abline(v=mth, col=1, lty=1)
    abline(v=mth*fmax, col=2, lty=2)
    
}

fmin<-0.5
fmax <- 1.5

png(sprintf("/home/genis/impala/localSitesQC/depth/plots/depthDist_allChromosomes_median%sand%s.png",fmin,fmax), width=600, height=800)
plotallthres(fmin, fmax)
dev.off()



getQth <- function(d, q){

    d <- as.numeric(d)
    cs <- cumsum(d)
    qx <- sum(d) * q
    thres <- min(which(cs > qx))
    thres
}

plotallthres <- function(qmin,qmax){
    yl <- c(0, max(ta_plot))
    par(mfrow=c(3,1))
    barplot(ta_plot, names.arg=0:500, ylim = yl,
            space=0, col="grey", ,border=NA,
            main=paste("Depth distribution all samples quantile thresholds", qmin, "and", qmax))
    abline(v=getQth(ta, qmin), col=2, lty=2)
    abline(v=getQth(ta, qmax), col=2, lty=2)

        barplot(tl_plot, names.arg=0:500, ylim = yl,
            space=0, col="grey", ,border=NA,
            main=paste("Depth distribution low depth samples quantile thresholds", qmin, "and", qmax))
    abline(v=getQth(tl, qmin), col=2, lty=2)
    abline(v=getQth(tl, qmax), col=2, lty=2)
    
    barplot(th_plot, names.arg=0:500, ylim = yl,
            space=0, col="grey", ,border=NA,
            main=paste("Depth distribution high depth samples quantile thresholds", qmin, "and", qmax))
    abline(v=getQth(ta, qmin), col=2, lty=2)
    abline(v=getQth(ta, qmax), col=2, lty=2)

}

qmin <- 0.01
qmax <- 0.99

png(sprintf("/home/genis/impala/localSitesQC/depth/plots/depthDist_allChromosomes_q%sand%s.png",qmin,qmax), width=600, height=800)
plotallthres(qmin, qmax)
dev.off()


qmin <- 0.02
qmax <- 0.98

png(sprintf("/home/genis/impala/localSitesQC/depth/plots/depthDist_allChromosomes_q%sand%s.png",qmin,qmax), width=600, height=800)
plotallthres(qmin, qmax)
dev.off()

qmin <- 0.03
qmax <- 0.99

png(sprintf("/home/genis/impala/localSitesQC/depth/plots/depthDist_allChromosomes_q%sand%s.png",qmin,qmax), width=600, height=800)
plotallthres(qmin, qmax)
dev.off()


qmin <- 0.05
qmax <- 0.99

png(sprintf("/home/genis/impala/localSitesQC/depth/plots/depthDist_allChromosomes_q%sand%s.png",qmin,qmax), width=600, height=800)
plotallthres(qmin, qmax)
dev.off()

q005th <- getQth(ta, 0.05)

q099th <- getQth(ta, 0.99)


qmin <- 0.07
qmax <- 0.99

png(sprintf("/home/genis/impala/localSitesQC/depth/plots/depthDist_allChromosomes_q%sand%s.png",qmin,qmax), width=600, height=800)
plotallthres(qmin, qmax)
dev.off()
