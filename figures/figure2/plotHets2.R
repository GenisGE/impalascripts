source("/home/genis/impala/info_files/loadPopInfo2.R")
library(data.table)

pop <- df$my_locality2[is.na(df$Exclude_reason)]
names(pop) <- df$`Serial number`[is.na(df$Exclude_reason)]

#dat <- read.table("/home/genis/impala/analyses/goatMapV2/relatedness_heterozygosities/results/sfs/collected.txt",h=T)

dat <- read.table("/home/genis/impala/analyses/impalaMap/relatedness_heterozygosities/results/sfs/collected.txt", h=T)
dat <- dat[!dat$id %in% df$`Serial number`[!is.na(df$Exclude_reason)],]

dat$locality <- factor(pop[as.character(dat$id)], levels=popord)
errfile <- "/home/genis/impala/analyses/goatMapV3/errorRates/results/angsdErrorEst.txt"
errfile <- "/home/genis/impala/analyses/goatMapV3/errorRates/errrates2/results/angsdErrorEst.txt"

err <- fread(cmd=sprintf("grep %% %s | cut -f2 -d' '", errfile),
                 data.table=F, h=F)$V1/100

inds <- gsub(x=read.table("/home/genis/impala/info_files/bamsImpalaGoatMappedCollapsedV0.list", sep="/")$V8, pattern=".bam", replacement="")
errdat <- as.data.frame(cbind(inds,err)[is.na(df$Exclude_reason),], stringsAsFactors=FALSE)
errdat$err <- as.numeric(errdat$err)



dat <- merge(x=dat, y=errdat, by.x="id", by.y="inds")
hd <- as.character(dat$id) %in% as.character(info$Serial_number[info$Depth == "hi"])

k <- (dat$err < 0.0002) | hd

outpdf <- "/home/genis/impala/paperplots/figure2/hets.pdf"

pdf(outpdf, width=12, height=7)
par(oma=c(5,5,0,0))
boxplot(dat$het[k] ~ dat$locality[k], col=impala_colors[popord], outline=F,las=2, xpd=NA, ylim=c(0.0007, 0.0023), cex.lab=2, cex.axis=2, xaxt="n")
axis(1, at=1:length(popord), labels=F)
axis(1, at=1:length(popord), tick=F, labels= popord, cex.axis=1.7, las=2, line=1)
title(ylab="Heterozygosity", line=7, xpd=NA, cex.lab=2)
points(dat$het[k] ~ jitter(as.integer(dat$locality[k]),1),pch=21,bg=impala_colors[dat$locality[k]], cex=1.5)
dev.off()




if(FALSE){

bitmap("/home/genis/impala/paperplots/figure2/heterozygosities_err002orhd.png", width=7.5, height=6,res=300)
par(oma=c(5,5,0,0))
boxplot(dat$het[k] ~ dat$locality[k], col=impala_colors[popord], outline=F,las=2, xpd=NA, ylim=c(0.0007, 0.0023), cex.lab=2, cex.axis=2, xaxt="n")
axis(1, at=1:length(popord), labels=F)
axis(1, at=1:length(popord), tick=F, labels= popord, cex.axis=1.7, las=2, line=1)
title(ylab="Heterozygosity", line=7, xpd=NA, cex.lab=2)
points(dat$het[k] ~ jitter(as.integer(dat$locality[k]),1),pch=21,bg=impala_colors[dat$locality[k]], cex=1.5)
dev.off()



bitmap("/home/genis/impala/paperplots/figure2/heterozygosities_err002orhd_wide.png", width=8, height=4,res=300)
par(oma=c(4,5,0,0))
boxplot(dat$het[k] ~ dat$locality[k], col=impala_colors[popord], outline=F,las=2, xpd=NA, ylim=c(0.0007, 0.0023), cex.lab=2, cex.axis=2, xaxt="n")
axis(1, at=1:length(popord), labels=F)
axis(1, at=1:length(popord), tick=F, labels=gsub(" ", "\n", popord), cex.axis=2, las=2, line=1)
title(ylab="Heterozygosity", line=7, xpd=NA, cex.lab=2)
points(dat$het[k] ~ jitter(as.integer(dat$locality[k]),1),pch=21,bg=impala_colors[dat$locality[k]], cex=2)
dev.off()

}
