source("/home/genis/impala/analyses/goatMapV3/roh/roh_call/rohrfuns/rohfuns.R")

source("/home/genis/impala/info_files/loadPopInfo2.R")

#args <- commandArgs(trailingOnly=T)

#f <- args[1]
#outpng<- args[2]
#mergedist <- as.numeric(args[3])

mergedist <- 5e5
g <- "/davidData/data/genis/impala/ref/goat/goat_ref_renamed_autosomes.bed"

f <- "/home/genis/impala/analyses/goatMapV3/roh/roh_call_imputed/roh_miss20_het5_den100_kb500_win50.hom"
autosome <- read.table(g)
autosome_len <- sum(as.numeric(autosome$V3))

roh <- read.table(f, h=T)

samples <- as.character(unique(roh$IID))
chrs <- sort(unique(roh$CHR))


getROHproportion <- function(roh, s, autosome_len){

    x <- roh[roh$IID==s,]
    d <- c(
        sum(x$KB[x$KB >= 1e3 &x$KB < 2.5e3] * 1e3),
        sum(x$KB[x$KB >= 2.5e3 & x$KB < 5e3] * 1e3),
        sum(x$KB[x$KB >= 5e3 & x$KB < 10e3] * 1e3),
        sum(x$KB[x$KB >= 10e3] * 1e3)
    )

    return(d/autosome_len)
}



rohcolslegend <- c("#cccccc", "#a6a6a6", "#8f8f8f", "#5d5d5d")


l <- apply(as.matrix(expand.grid(chrs, samples)),1, function(x) combineROH(roh[roh$CHR==x[1] & roh$IID == x[2],], maxdist=mergedist))
roh <- do.call("rbind", l)


famfile <- "/home/genis/impala/analyses/goatMapV3/call_genos/results/impute/call_imputed_all.fam"

fam <- read.table(famfile, h=F, stringsAsFactors=F)
ids <- fam$V1

pop <- fam$V2
#pop <- info$Locality
#names(pop) <- info$Serial_number
#pop <- pop[as.character(ids)]
pop <- gsub("LakeMburu", "LakeMburo", pop)

names(pop) <- as.character(ids)
superpop <- info$Superpop
names(superpop) <- as.character(info$Serial_number)

m <- sapply(ids, getROHproportion, roh=roh, autosome_len=autosome_len)

colnames(m) <- ids
popord2 <- gsub(" ", "", popord)
ord <- unlist(sapply(popord2, function(x) which(pop==x)))
m <- m[,ord]

npop <- length(popord)
m2 <- do.call(rbind, replicate(npop, m, simplify=FALSE)) # where m is your matrix
#m2 <- m2[,ord]
m2 <- matrix(0, ncol=ncol(m), nrow=4 * npop)

for(p in popord2){

    i <- which(popord2==p)
    start <- (i - 1) * 4 + 1
    end <- start + 3
    w <- which(pop[ord] == p)

    m2[ (start:end) , w] <- m[,w]
    
}


impala_colors <- c("Etosha" = "#8A0E83", "Ovita" = "#F72B69", "Chobe" = "#B81E22", "Shangani" = "#FF927C", "Mana Pools" = "#DB9139", "Kafue" = "#E0D152", "Luangwa" = "#C18E8E", "Selous" = "#698065", "Ugalla" = "#86ABCB", "Masai Mara" = "#26C6C9", "Tsavo" = "#94D9DC", "Samburu" = "#2B597D", "Lake Mburo" = "#0559E8")
library(colorspace)

roh_cols <- c(rbind(
    lighten(impala_colors, 0.8),
    lighten(impala_colors, 0.4),
    darken(impala_colors, 0.1),
    darken(impala_colors, 0.6)
))

legend_cols <- c(
    lighten("grey", 0.8),
    lighten("grey", 0.4),
    darken("grey", 0.1),
    darken("grey", 0.6)
)
legend_cols <- c("#E6E6E6", "#C3C3C3", "#969696", "#787878")


outpng <- "/home/genis/impala/paperplots/figure6/impalaRohsnice.png"

bitmap(outpng, h=4, w=12, res=300)
par(oma=c(0,1,1.3,6.5))
barplot(m2,
        #names.arg=ids[ord],
        las=2, col=roh_cols, #ylab="Genome fraction",
        main="",
        cex.lab=1.5, cex.axis=1.2, cex.main=1.4, space=0, border=NA)#, lwd=0.1)
title(ylab="Genome fraction", cex.lab=2, line=2, xpd=NA)

ybot <- seq(max(colSums(m)) * 0.2, max(colSums(m)) * 0.5, length.out=4) 
ytop <- seq(max(colSums(m)) * 0.2, max(colSums(m)) * 0.5, length.out=4) + max(colSums(m)) * 0.1
#ytop <- seq(0.1,0.2,length.out=4) + 0.1/3
ymed <- (ybot + ytop) / 2
rect(xleft=ncol(m) *1.035, xright=ncol(m) *1.055, ybottom=ybot,
     ytop=ytop,
     col=legend_cols,xpd=NA)


text(x=ncol(m) *1.09, y=ymed, labels=c("1-2.5 Mbp", "2.5-5Mbp", "5-10 Mbp", ">10 Mbp"), xpd=NA,cex=1.2)
text(font=2,x=ncol(m) *1.055,y=ytop[4] * 1.1, labels="ROH length", xpd=NA,cex=1.4)

spop_labs <- c("Black-faced\nimpala", "Southern common\nimpala", "Selous\ncommmon\nimpala", "Eastern common\nimpala")
text(x=sort(tapply(1:length(superpop),superpop[ord],mean)),y=max(colSums(m)) * 1.1,labels=spop_labs,xpd=NA, cex=2, srt=0, font=2)
abline(v=rev(cumsum(sapply(unique(superpop[ord]),function(x){sum(superpop[ord]==x)})))[-1],col=1,lwd=1, lty=2)

text(sort(tapply(1:length(pop),pop[ord],mean)),-0.018,unique(pop[ord]),xpd=NA, cex=1.7, srt=30)
#abline(v=cumsum(sapply(unique(pop[ord]),function(x){sum(pop[ord]==x)})),col=1,lwd=2)

dev.off()
