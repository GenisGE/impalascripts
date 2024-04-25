

makeCol <- function( nHalf=10){
     color_palette=c("#001260", "#EAEDE9", "#601200")
   
    ## Make vector of colors for values below threshold
    rc1 <- colorRampPalette(colors = color_palette[1:2], space="Lab")(nHalf)    
    ## Make vector of colors for values above threshold
    rc2 <- colorRampPalette(colors = color_palette[2:3], space="Lab")(nHalf)
    rampcols <- c(rc1, rc2)
    rampcols[c(nHalf, nHalf+1)] <- rgb(t(col2rgb(color_palette[2])), maxColorValue=256)
    return(rampcols)
}

getMean <- function(mat,id){
    uid <- unique(id)
    m <- matrix(NA,length(uid),length(uid),dimnames=list(uid,uid))
    for(i in uid)
        for(j in uid)
            m[i,j] <- mean(mat[id==i,id==j],na.rm=T)
    m
}



getColor <- function(val,corColors,maxCor=0.25){
   
    Nc <- length(corColors)
    val <- ifelse(val < -maxCor,-maxCor,val)
    val <- ifelse(val > maxCor,maxCor,val)    
    corColors[cut(val,seq(-maxCor,maxCor,length.out=Nc),include.lowest=T)]
}
addKey<- function(from,to,N,maxCor=0.1,corColors){
    if(missing(corColors))
        corColors <- c("#00125F","#2A266E","#443C7D","#5C528C","#746A9B","#8B83AB","#A29CBA","#BAB6C9","#D2D1D9","#E9ECE8","#E9ECE8","#DDD3CC","#CFB9B0","#C1A194","#B3887A","#A37160","#935948","#834230","#712B1A","#5F1200")

    y0 <- (to-from)*0.6+from
    y1 <- (to-from)*0.75+from
    y2 <- (to-from)*0.8+from
    y3 <- (to-from)*0.9+from
    hx <- N/6
    rasterImage(as.raster(matrix(corColors, nrow=1)),0,y0,hx,y1,xpd=NA)
    text(0,y2,-maxCor,xpd=NA, cex=1.4)
    text(hx,y2,maxCor,xpd=NA, cex=1.4)
    text(hx,y3,"Correlation of residuals",xpd=NA,adj=1, cex=1.5)
}

addCor <- function(corMat,from=1,to=1.2,popID,withinOnly= FALSE,maxCor=0.25,lines=0,x,corColors,meanCor=FALSE){
    if(missing(corColors))
        corColors <- c("#00125F","#2A266E","#443C7D","#5C528C","#746A9B","#8B83AB","#A29CBA","#BAB6C9","#D2D1D9","#E9ECE8","#E9ECE8","#DDD3CC","#CFB9B0","#C1A194","#B3887A","#A37160","#935948","#834230","#712B1A","#5F1200")
 
    if(missing(x))
        x <- 1:nrow(corMat)-0.5
    if((withinOnly | lines | meanCor) &  missing(popID))
        stop("must sypply popID")
    N <- ncol(corMat)
    y <- from + (1:N/N)*(to - from)
    ySize <- diff(y[1:2])/2
    xSize <- diff(x[1:2])/2

    if(meanCor){
        mCor <- getMean(corMat,popID)
        intID <- match(popID,unique(popID))
        for(i in 1:N-1)
            for(j in i:N)
                corMat[i,j] <- mCor[intID[i],intID[j]]
    }
    ## plot individuals pairs
    for(i in 2:N-1){
        cat("\r",i)
        for(j in (i+1):N){
            if(withinOnly)
                if(popID[i]!=popID[j])
                    next
            xc <- (x[j]-x[i])/2 + x[i]
            yc <- y[j-i]
            polygon(
                x=c(0,-1,0,1,0)*xSize+xc ,
                y=c(-1,0,1,0,-1)*ySize+yc,
                col=getColor(corMat[i,j],maxCor=maxCor,corColors),
                xpd=NA,
                border=getColor(corMat[i,j],maxCor=maxCor,corColors)
                #,lty=0
            )
        }
    }
    
    ## add lines seperating pops
    if(lines>0){
        popSep <- c(0,cumsum(table(popID)[unique(popID)]))
        popSepR <- rev(popSep)
        N <- max(popSep)
        for(i in 2:length(popSep)-1){
            lines( c( x[popSep[i]+1],x[(N-popSep[i])/2 + popSep[i]]+xSize) , y[c( 1, N-popSep[i])] , xpd=NA,lwd=lines )
            lines( c(x[popSepR[i]],(x[popSepR[i]]+xSize)/2) , y[c( 1, popSepR[i])] , xpd=NA,lwd=lines )
            
        }
    }
 }


if(FALSE){ #starts here
    source("/home/albrecht/github/evalAdmix/visFuns.R")
    source("/home/albrecht/github/evalAdmix/aaplot.R")    

# read population labels and estimated admixture proportions
    pop<-read.table("admixTjeck2.fam",stringsAsFactors=F)
    q<-read.table("admixTjeck2.3.Q")
    corMat <-as.matrix(read.table("output.corres.txt"))


# order according to population and plot the ADMIXTURE reults
    popID <- pop[,2]
    table(popID)
    ord<-orderInds(pop = as.vector(popID), q = q)
    q <- q[ord,]
    popID <- popID[ord]
    corMat <- corMat[ord,ord]


## 
Q <- t(q)
N <- length(popID)
popNames <- unique(popID)
Npop <- length(popNames)
popIDint <- rep(1:Npop,table(popID)[popNames])
popSep <- c(0,cumsum(table(popIDint)))
meanPopX <- popSep[-1]-table(popIDint)/2


    
## colors for admixture
colorpal= c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
# colors for eval
corColors <-makeCol(10)

## plot individual cor above
bitmap("tmp.png",w=10,h=10,res=200)
par(mar=c(5.1,4.1,4.1+6,2.1))

x <- barplot(Q, col=colorpal, space=0, border=NA, cex.axis=1.2,cex.lab=1.8,axisnames=FALSE,
             ylab="Admixture proportions", xlab="", main="", cex.main=1.5,xpd=NA)
    #make population labels
    text(meanPopX,rep(-0.05,length(meanPopX)),unique(popID),xpd=T,font=2,cex=1.5)
    #lines to between populations
    abline(v=popSep)
    text(mean(x),1.35,"add correlation above",font=2,cex=2,xpd=T)
    ## from and to is the placement on the y-axis. 0-1 is where the admixture plot it so 1-1.3 will place evalAdmix above
    addKey(from=1,to=1.3,N=ncol(Q)) ## from 1 to 1.3 change to make evalAdmix hgher/lower
    addCor(corMat,from=1,to=1.3,maxCor=0.1,lines=0.2,popID=popID)

dev.off()




## plot individual cor below
bitmap("tmp2.png",w=10,h=10,res=200)
par(mar=c(5.1+6,4.1,4.1,2.1))

x <- barplot(Q, col=colorpal, space=0, border=NA, cex.axis=1.2,cex.lab=1.8,axisnames=FALSE,
             ylab="Admixture proportions", xlab="", main="", cex.main=1.5,xpd=NA)
text(meanPopX,rep(1.05,length(meanPopX)),unique(popID),xpd=T,font=2,cex=1.5)

abline(v=popSep)
 addKey(from=0,to=-.3,N=ncol(Q))
addCor(corMat,from=0,to=-0.3,maxCor=0.1,lines=0.2,popID=popID)

dev.off()

## both but below as mean
bitmap("tmp3.png",w=10,h=10,res=200)
par(mar=c(5.1+4,4.1,4.1+4,2.1))

x <- barplot(Q, col=colorpal, space=0, border=NA, cex.axis=1.2,cex.lab=1.8,axisnames=FALSE,
             ylab="Admixture proportions", xlab="", main="", cex.main=1.5,xpd=NA)
abline(v=popSep)
addKey(from=1,to=1.3,N=ncol(Q))
addCor(corMat,from=1,to=1.3,maxCor=0.1,lines=0.2,popID=popID)
addCor(corMat,from=-0.1,to=-0.4,maxCor=0.1,lines=0.2,popID=popID,mean=T)
text(0,-0.2,"Mean\n correlation",adj=0,xpd=T)
text(meanPopX,rep(-0.05,length(meanPopX)),unique(popID),xpd=T,font=2,cex=1.5)

dev.off()

}
 







