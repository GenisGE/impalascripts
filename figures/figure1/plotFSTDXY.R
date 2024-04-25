source("/home/genis/impala/info_files/loadPopInfo2.R")

infile <- "/home/genis/impala/analyses/impalaMap/sfss_fst/results2/fst/all_fst.txt"
infile <- "/kellyData/home/genis/impala/analyses/speciescompare/fstvsdxy/results/impala/impala_all_fst_dxy.tsv"

fstdf <- read.table(infile, header=T, stringsAsFactors=F)

fstdf$pop1 <- sapply(strsplit(fstdf$pop_pair, split="_"), function(x) x[1])
fstdf$pop2 <- sapply(strsplit(fstdf$pop_pair, split="_"), function(x) x[2])

fstdf$pop1 <- gsub(pattern="MasaiMara", replacement="Masai Mara", x=fstdf$pop1)
fstdf$pop2 <- gsub(pattern="MasaiMara", replacement="Masai Mara", x=fstdf$pop2)

fstdf$pop1 <- gsub(pattern="LakeMburu", replacement="Lake Mburo", x=fstdf$pop1)
fstdf$pop2 <- gsub(pattern="LakeMburu", replacement="Lake Mburo", x=fstdf$pop2)

fstdf$pop1 <- gsub(pattern="ManaPools", replacement="Mana Pools", x=fstdf$pop1)
fstdf$pop2 <- gsub(pattern="ManaPools", replacement="Mana Pools", x=fstdf$pop2)


colpal_fst <- colorRampPalette(colors = c("#F7FBFF", "#08306B"), space="Lab")(50)
colpal_dxy <- colorRampPalette(colors = c("#FFF9F7", "#6B1508"), space="Lab")(50)


fstmat <- matrix(nrow=12, ncol=12, dimnames=list(popord[-length(popord)], popord[-1]))
fstmat <- matrix(nrow=13, ncol=13, dimnames=list(popord, popord))

for(i in 1:(length(popord)-1)) for(j in (i+1):(length(popord))) fstmat[popord[i], popord[j]] <- fstdf$fst[(fstdf$pop1==popord[i] & fstdf$pop2==popord[j]) | (fstdf$pop1==popord[j] & fstdf$pop2==popord[i])]

scols <- sapply(colnames(fstmat), function(x) info$Superpop[which(info$Locality==x)[1]])
srows <- sapply(rownames(fstmat), function(x) info$Superpop[which(info$Locality==x)[1]])

dxymat <- matrix(nrow=12, ncol=12, dimnames=list(popord[-length(popord)], popord[-1]))
dxymat <- matrix(nrow=13, ncol=13, dimnames=list(popord, popord))

for(i in 1:(length(popord)-1)) for(j in (i+1):(length(popord))) dxymat[popord[i], popord[j]] <- fstdf$dxy[(fstdf$pop1==popord[i] & fstdf$pop2==popord[j]) | (fstdf$pop1==popord[j] & fstdf$pop2==popord[i])]




regcols <- c(impala_colors["Etosha"], "darkred", impala_colors["Selous"],  "darkblue")

outpng <- "/home/genis/impala/paperplots/figure1/fstsdxy.png"
bitmap(outpng, width=6, height=5.5, res=300)
# https://stackoverflow.com/questions/3789549/display-a-matrix-including-the-values-as-a-heatmap (second answer)
par(oma=c(2,2.7,0,6.5), mar=par("mar") + c(0, 0, -2.5, 0))
image(1:ncol(fstmat), 1:nrow(fstmat), t(fstmat), col = colpal_fst, axes = FALSE, ylab="", xlab="")
#axis(1, 1:ncol(fstmat), colnames(fstmat), las=2, tick=FALSE, xpd=NA, cex=1.5)
#axis(2, 1:nrow(fstmat), rownames(fstmat), las=2, tick=FALSE, xpd=NA, cex=1.5)


text(y=-0.89, x=1:length(popord), popord,
     col=impala_colors, font=2, srt=90, xpd=NA, cex=1.4)
text(x= -1.1, y=1:length(popord), popord,
     col=impala_colors, font=2, xpd=NA, cex=1.4)



image(1:ncol(dxymat), 1:nrow(dxymat), dxymat, col = colpal_dxy, axes = FALSE, ylab="", xlab="", add=T)


abline(v=rev(cumsum(sapply(unique(scols),function(x){sum(scols==x)})))[-1] + 0.5, lwd=2.5, col="black")
abline(h=rev(cumsum(sapply(unique(srows),function(x){sum(srows==x)})))[-1] + 0.5, lwd=2.5, col="black")

segments(x0=0.5,
          y0=c(0, cumsum(sapply(unique(scols),function(x){sum(scols==x)}))[-4]) + 0.5,
          y1=cumsum(sapply(unique(scols),function(x){sum(scols==x)})) + 0.5,
          lwd=5, col=regcols, xpd=NA)

segments(y0=0.5,
          x0=c(0, cumsum(sapply(unique(srows),function(x){sum(srows==x)}))[-4]) + 0.5,
          x1=cumsum(sapply(unique(srows),function(x){sum(srows==x)})) + 0.5,
          lwd=5, col=regcols, xpd=NA)

#text(x= -0.5 - c(-0.4,0.2,-0.4, 0.2),y=sort(tapply(1:length(srows),srows,mean)),
#     gsub("impala", "", unique(srows)),xpd=NA, srt=90, cex=1.2,font=2)

#text(y=12 + (2.3+c(-0.5,0.1,-0.5, 0.1)), x=sort(tapply(1:length(scols),scols,mean)),
#     gsub("impala", "", unique(scols)), xpd=NA, srt=0, cex=1.2, font=2)

# plot legends
fst_legend <- rev(as.raster(matrix(colpal_fst, ncol=1)))

rasterImage(fst_legend, xleft=14, ybottom=1, xright = 15, ytop = 5, xpd=NA)
rect(xleft=14, ybottom=1, xright = 15, ytop = 5, xpd=NA)
text(x=16, y = seq(1 + 0.25,5 - 0.25,length.out=5),
     labels = signif(seq(min(fstmat, na.rm=T), max(fstmat, na.rm=T), length.out=5), 2),
     cex=1.2,xpd=NA)
text(x=14.5, y=5.5, cex=2, font=2, labels=expression("F"[ST]), xpd=NA)

dxy_legend <- rev(as.raster(matrix(colpal_dxy, ncol=1)))

rasterImage(dxy_legend, xleft=14, ybottom=8, xright = 15, ytop = 12, xpd=NA)
rect(xleft=14, ybottom=8, xright = 15, ytop = 12, xpd=NA)
text(x=16, y = seq(8 + 0.25, 12 - 0.25,length.out=5),
     labels = signif(seq(min(dxymat, na.rm=T), max(dxymat, na.rm=T), length.out=5), 2),
     cex=1.2,xpd=NA)
text(x=14.5, y=12.5, cex=2, font=2, labels=expression("D"[xy]), xpd=NA)

dev.off()


outpdf <- "/home/genis/impala/paperplots/figure1/fstsdxy.pdf"
pdf(outpdf, width=6 * 1.5, height=5.5 * 1.5)
# https://stackoverflow.com/questions/3789549/display-a-matrix-including-the-values-as-a-heatmap (second answer)
par(oma=c(2,2.8,0,6.5), mar=par("mar") + c(0, 0, -2.5, 0))
image(1:ncol(fstmat), 1:nrow(fstmat), t(fstmat), col = colpal_fst, axes = FALSE, ylab="", xlab="")
#axis(1, 1:ncol(fstmat), colnames(fstmat), las=2, tick=FALSE, xpd=NA, cex=1.5)
#axis(2, 1:nrow(fstmat), rownames(fstmat), las=2, tick=FALSE, xpd=NA, cex=1.5)


text(y=-0.89, x=1:length(popord), popord,
     col=impala_colors, font=2, srt=90, xpd=NA, cex=1.4)
text(x= -1.1, y=1:length(popord), popord,
     col=impala_colors, font=2, xpd=NA, cex=1.4)



image(1:ncol(dxymat), 1:nrow(dxymat), dxymat, col = colpal_dxy, axes = FALSE, ylab="", xlab="", add=T)


abline(v=rev(cumsum(sapply(unique(scols),function(x){sum(scols==x)})))[-1] + 0.5, lwd=2.5, col="black")
abline(h=rev(cumsum(sapply(unique(srows),function(x){sum(srows==x)})))[-1] + 0.5, lwd=2.5, col="black")

segments(x0=0.5,
          y0=c(0, cumsum(sapply(unique(scols),function(x){sum(scols==x)}))[-4]) + 0.5,
          y1=cumsum(sapply(unique(scols),function(x){sum(scols==x)})) + 0.5,
          lwd=5, col=regcols, xpd=NA)

segments(y0=0.5,
          x0=c(0, cumsum(sapply(unique(srows),function(x){sum(srows==x)}))[-4]) + 0.5,
          x1=cumsum(sapply(unique(srows),function(x){sum(srows==x)})) + 0.5,
          lwd=5, col=regcols, xpd=NA)

#text(x= -0.5 - c(-0.4,0.2,-0.4, 0.2),y=sort(tapply(1:length(srows),srows,mean)),
#     gsub("impala", "", unique(srows)),xpd=NA, srt=90, cex=1.2,font=2)

#text(y=12 + (2.3+c(-0.5,0.1,-0.5, 0.1)), x=sort(tapply(1:length(scols),scols,mean)),
#     gsub("impala", "", unique(scols)), xpd=NA, srt=0, cex=1.2, font=2)

# plot legends
fst_legend <- rev(as.raster(matrix(colpal_fst, ncol=1)))

rasterImage(fst_legend, xleft=14, ybottom=1, xright = 15, ytop = 5, xpd=NA)
rect(xleft=14, ybottom=1, xright = 15, ytop = 5, xpd=NA)
text(x=16, y = seq(1 + 0.25,5 - 0.25,length.out=5),
     labels = signif(seq(min(fstmat, na.rm=T), max(fstmat, na.rm=T), length.out=5), 2),
     cex=1.2,xpd=NA)
text(x=14.5, y=5.5, cex=2, font=2, labels=expression("F"[ST]), xpd=NA)

dxy_legend <- rev(as.raster(matrix(colpal_dxy, ncol=1)))

rasterImage(dxy_legend, xleft=14, ybottom=8, xright = 15, ytop = 12, xpd=NA)
rect(xleft=14, ybottom=8, xright = 15, ytop = 12, xpd=NA)
text(x=16, y = seq(8 + 0.25, 12 - 0.25,length.out=5),
     labels = signif(seq(min(dxymat, na.rm=T), max(dxymat, na.rm=T), length.out=5), 2),
     cex=1.2,xpd=NA)
text(x=14.5, y=12.5, cex=2, font=2, labels=expression("D"[xy]), xpd=NA)

dev.off()


