
l <- lapply(list.files("/home/genis/impala/localSitesQC/depth/output/lowDepth/",".depthSample", full.names=T), function(x) as.matrix(read.table(x, colClasses="numeric")))
h <- lapply(list.files("/home/genis/impala/localSitesQC/depth/output/highDepth/",".depthSample", full.names=T), function(x) as.matrix(read.table(x, colClasses="numeric")))

ltot <- matrix(data=0, nrow=nrow(l[[1]]), ncol=ncol(l[[1]]))
htot <- matrix(data=0, nrow=nrow(h[[1]]), ncol=ncol(h[[1]]))

for(c in l){
ltot <- ltot + c
}

for(c in h){
htot <- htot + c
}

meanDepth <- function(x) sum(x* 0:(length(x)-1))/sum(x)
medianDepth <- function(x) median(rep.int(0:(length(x)-1), x))

lmean <- apply(ltot,1,meanDepth)
lmedian <- apply(ltot,1,medianDepth)

hmean <- apply(htot,1,meanDepth)
hmedian <- apply(htot,1,medianDepth)
