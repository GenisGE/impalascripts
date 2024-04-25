source("/home/genis/impala/info_files/loadPopInfo2.R")
source("/home/genis/software/evalAdmix/visFuns.R")
source("/home/genis/impala/paperplots/figure2/aaplot_gge2.R")


qfile <- "/home/genis/impala/analyses/impalaMap/ngsadmix/results/12/admixResultImpalaMap.12.44.qopt_conv"
rfile <- "/home/genis/impala/analyses/impalaMap/ngsadmix/evaladmix/evalAdmixResultsImpalaMap_K12.corres"


refinds <- c("Etosha" = 2077, "Ovita" = 973,  "Chobe" = 7830, "Shangani" = 1511, "Mana Pools" = 1619, "Luangwa" = 4555, "Selous" = 8517, "Ugalla" = 5459, "Masai Mara" = 5, "Tsavo" = 55, "Lake Mburu" = 4660, "Samburu" = 37)

kpopord <- c("Ovita", "Samburu", "Etosha", "Lake Mburo", "Selous", "Tsavo", "Kafue", "Ugalla", "Chobe", "Shangani", "Masai Mara", "Mana Pools")


# try to do proper plot with all k plots
pop <- info$Locality


ord <- orderInds(pop=info$Locality, popord=popord)

q <- as.matrix(read.table(qfile))
r <- as.matrix(read.table(rfile))
k <- 12
k_ord <- orderK(q=q, refpops=kpopord[1:k], pop=info$Locality)


outpdf <- "/home/genis/impala/paperplots/figure1/admixeval.pdf"
pdf(outpdf, h=8,w=16)
par(mar=c(6,4,16.5,2))
barplot(t(q[ord,k_ord]), col=impala_colors[kpopord[1:k]], space=0, border=NA, cex.axis=1.2,cex.lab=1.5,
        ylab="Admixture proportions", xlab="",xpd=NA)
abline(v=1:nrow(q), col="white", lwd=0.2)

text(x=117, y=0.5, labels="K = 12", xpd=NA,cex=1.8)

abline(v=cumsum(sapply(unique(pop[ord])[-length(unique(pop))],function(x){sum(pop[ord]==x)})),col=1,lwd=1.2)

addKey(from=1., to=1.9, N=nrow(q), maxCor=0.15)
addCor(r[ord,ord], popID = pop[ord],from = 1, to = 1.9, maxCor = 0.15, lines=1)

text(sort(tapply(1:length(pop),pop[ord],mean)),-0.16,gsub(" ", "\n", unique(pop[ord])),xpd=NA, srt=90, cex=1.5)

dev.off()

if(FALSE){
outpng <- "/home/genis/impala/paperplots/figure1/admixeval.png"
bitmap(outpng, h=4,w=8, res=300)
par(mar=c(6,4,16.5,2))
barplot(t(q[ord,k_ord]), col=impala_colors[kpopord[1:k]], space=0, border=NA, cex.axis=1.2,cex.lab=1.5,
        ylab="Admixture proportions", xlab="",xpd=NA)
abline(v=1:nrow(q), col="white", lwd=0.2)

text(x=117, y=0.5, labels="K = 12", xpd=NA,cex=1.8)

abline(v=cumsum(sapply(unique(pop[ord])[-length(unique(pop))],function(x){sum(pop[ord]==x)})),col=1,lwd=1.2)

addKey(from=1., to=1.9, N=nrow(q), maxCor=0.15)
addCor(r[ord,ord], popID = pop[ord],from = 1, to = 1.9, maxCor = 0.15, lines=1)

text(sort(tapply(1:length(pop),pop[ord],mean)),-0.16,gsub(" ", "\n", unique(pop[ord])),xpd=NA, srt=90, cex=1.5)

dev.off()
}
