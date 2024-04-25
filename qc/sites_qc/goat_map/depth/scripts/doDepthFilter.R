library(data.table)

lowDepthDir <- "/home/genis/impala/localSitesQC/depth/output/lowDepth/"
highDepthDir <- "/home/genis/impala/localSitesQC/depth/output/highDepth/"

chrs <- scan("/home/genis/impala/info_files/goatChrsAutoandX.txt", what="da")

badlow <- data.frame(chr=character(0), start=integer(0), end=integer(0))
badhigh <- data.frame(chr=character(0), start=integer(0), end=integer(0))

# thresholds are lower 0.5 * median and higher 1.5 * median, separately for low depth and high depth samples
# see graphical distribution in two bottom panels (top panel is wrong) of /home/genis/impala/localSitesQC/depth/plots/depthDist_allChromosomes_median0.5and1.5.png

# low depth median 256
minl <- 128
maxl <- 384

# high depth median 136
minh <- 68
maxh <- 204

for(chr in chrs){
    # make low depth filts
    x <- fread(paste0(lowDepthDir, chr, "_lowDepth.pos.gz"), h=T, data.table=F)
    
    i <- 1
    
    check <- TRUE
    maxpos <- nrow(x)
    while(check){
        
        if(x$pos[i] < minl | x$pos[i] > maxl) start <- pos[i]
       
        while(x$pos[i] < minl | x$pos[i] > maxl){
            end <- pos[i]
            i <- i + 1
            cat(i, "/n")
            if(i > maxpos){
                badlow$chr <- c(badlow$chr, chr)
                badlow$start <- c(badlow$start, start)
                badlow$chr <- c(badlow$end, end)
                break
            }
        }

        badlow$chr <- c(badlow$chr, chr)
        badlow$start <- c(badlow$start, start)
        badlow$chr <- c(badlow$end, end)
        
        i <- i + 1
        cat(i,"/n")
        if(i > maxpos) break
            
    }
}

/n8 126 754
