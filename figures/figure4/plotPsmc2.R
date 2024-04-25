library(scales)

indir <- "/home/genis/impala/analyses/impalaMap/psmc/results/psmc_output/all_filters/2" #opt$indir
bootdir <- "/home/genis/impala/analyses/impalaMap/psmc/results/psmc_output/all_filters/2/bootstrap" #opt$bootdir
mu <- 1.41e-8 #opt$mu
g <- 5.7 #opt$generation


# read pscm output psmc.result() adapted from function in https://datadryad.org/stash/dataset/doi:10.5061/dryad.0618v

##-------Rescale the ith iteration result of PSMC, and make ready for plotting
# file: result file from PSMC
# i.iteration: the ith iteration
# mu: mutation rate
# s: bin size
# g: years per generation

psmc.result<-function(file,i.iteration=25,mu=1e-8,s=100,g=1){
       #cat("will reescale file", file,"\n")
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	
	START<-grep("^RD",X)
	END<-grep("^//",X)

        if(length(START) < (i.iteration+1)) i.iteration <- length(START) - 1
	X<-X[START[i.iteration+1]:END[i.iteration+1]]
	
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	write(TR,"temp.psmc.result")
	theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
	N0<-theta0/4/mu/s 
	
	write(RS,"temp.psmc.result")
	a<-read.table("temp.psmc.result")
	Generation<-as.numeric(2*N0*a[,3]) # a[,3] is t_k
	Ne<-as.numeric(N0*a[,4]) #a[,4] is lambda_k
	
	file.remove("temp.psmc.result")
	
	n.points<-length(Ne)
	YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])*g
	Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
		Ne[n.points])

       #cat("finished reescale file", file,"\n\n")	
       data.frame(YearsAgo,Ne)

}



files <-list.files(indir,".psmc",full.names=T)

res <- lapply(files, psmc.result,mu=mu,g=g)
names(res) <- gsub(".psmc", "",basename(files))


# THIS IS GENERALIZABLE, BUT WILL FAIL IF THERE ARE MORE THAN 9 INDIVIDUALS TO PLOT.
#inds <- names(res)
#cols <- RColorBrewer::brewer.pal("Set1", n=length(inds))


#### ALL THIS IS SPECIFIC TO IMPALA. CHANGE WHEN DOING WITH OTHER DATASET ####
###### WOULD BE NICE TO DO A GENERALIZABLE FUNCTION SOME DAY #######
source("/home/genis/impala/info_files/loadPopInfo2.R")

inds <- as.character(info$Serial_number[info$Depth=='hi'])
pop <- gsub("Mburu", "Mburo", info$Locality[info$Depth=='hi'])
inds <- c(inds, "rgp")

pop <- c(pop, "Laikipia")
names(pop) <- inds


impala_colors["Laikipia"] <-  "#1d0989" 
cols <- impala_colors
#cols["Laikipia"] <-  "#1d0989" 

#res2 <- lapply(res, function(x) x[-((nrow(x)-8):nrow(x)),])


rm_last <- function(x){ ## FUNCITON TO CLEAN DATA, IS A BIT ARBITRARY MIGHT BE GOOD TO CHECK AGAIN SOME DAY?

    nes <- unique(x$Ne)
    rmv <- nes[(length(nes)-1):length(nes)]
    x[!x$Ne%in%rmv,]
}

#res2 <- lapply(res,rm_last)
res2 <- res

bf <- c("Etosha")
southern <- c("Ovita", "Chobe", "Shangani", "Mana Pools", "Kafue", "Luangwa")
selous <- c("Selous")
eastern <- c("Ugalla", "Masai Mara", "Tsavo", "Samburu", "Lake Mburo", "Laikipia")


# get path to all bootstrapls and format names to only keep sample name
bootfiles <- list.files(bootdir, ".psmc", full.names=T)
names(bootfiles) <- gsub("_[a-zA-Z0-9]*","", gsub(".psmc","",basename(bootfiles)))

# load, reescale and clean bootrap
res_boot <- lapply(bootfiles, psmc.result,mu=mu,g=g)
#res_boot <- lapply(res_boot, function(x) x[-((nrow(x)-8):nrow(x)),])
#res_boot <- lapply(res_boot, rm_last)
#outname <- "/home/genis/impala/paperplots/figure5/impalaPsmcNice2.png" #opt$outfile
#outname <- "/home/genis/impala/paperplots/figure5/impalaPsmcNice2b.png" #opt$outfile



doplot <- function(){

ymax <- 6e4 # ARBITRARY THRESHOLD BECAUSE PLOTS USUALY LOOK NICE WITH THIS
par(mar=c(3,5,0,15)+0.1, mfrow=c(4,1), oma=c(2,0,1.5,0))

###############################################
############# 1 plot blackfaced ###############
###############################################
plot(type='n', x=1, log='x',
     xlab="",sprintf("Years ago (mu=%.2e, g=%.1f)",mu,g),
     ylab="Ne",cex.lab=1.85, xaxt='n', yaxt='n',
     xlim=xlims, ylim=c(0, ymax), cex.axis=1.5, xpd=NA)
xpos <- c(0.5e4, seq(1e4,1e5,by=1e4), seq(1e5, 1e6, by=1e5))
xlabs <- c(0.5e4, 1e4, 2e4, 1e5, 5e5,1e6)

axis(2, at=seq(0, ymax, 1e4), labels=c("0",paste0(as.integer(seq(1e4, ymax, 1e4)/1e3), "k")), cex.axis=1.5)
#axis(2, at=seq(0, ymax, 1e4), labels=c("0",paste0(as.integer(seq(10, ymax, 1e4)/1e3), "k"), cex.axis=1.5)

#axis(2, at=xlabs, tick=F, labels= paste(as.integer(xlabs/1e3), "kya"), cex.axis=1.5)

axis(1, at=xpos, labels=F)
axis(1, at=xlabs, tick=F, labels= paste(as.integer(xlabs/1e3), "kya"), cex.axis=1.5)

dogray <-   which(!pop %in% bf)
docol <-  which(pop %in% bf)
doboot <- which(names(res_boot) %in% inds[docol])

for(i in dogray) lines(x=res2[[inds[i]]]$YearsAgo, y=res2[[inds[i]]]$Ne,  col = "lightgray",lwd=3)
for(i in doboot) lines(x=res_boot[[i]]$YearsAgo, y=res_boot[[i]]$Ne,  col = alpha(cols[pop[names(res_boot)[i]]], 0.1),lwd=0.5,lty=1.5)

for(i in docol) lines(x=res2[[inds[i]]]$YearsAgo, y=res2[[inds[i]]]$Ne,  col = cols[pop[inds[i]]],lwd=3)

legend(x=xlims[2] * 1.4, y=ymax,
       title="Black-faced", legend=bf, col=impala_colors[bf],
       lty=1,lwd=3,border=NA, bty="n",cex=1.5, xpd=NA)

###############################################
############# 2 plot southern ###############
###############################################
plot(type='n', x=1, log='x',
     xlab="",#sprintf("Years ago (mu=%.2e, g=%.1f)",mu,g),
     ylab="Ne",cex.lab=1.85, xaxt='n', yaxt='n',
     xlim=xlims, ylim=c(0, ymax), cex.axis=1.5, xpd=NA)
xpos <- c(0.5e4, seq(1e4,1e5,by=1e4), seq(1e5, 1e6, by=1e5))
xlabs <- c(0.5e4, 1e4, 2e4, 1e5, 5e5,1e6)

axis(2, at=seq(0, ymax, 1e4), labels=c("0",paste0(as.integer(seq(1e4, ymax, 1e4)/1e3), "k")), cex.axis=1.5)

axis(1, at=xpos, labels=F)
axis(1, at=xlabs, tick=F, labels= paste(as.integer(xlabs/1e3), "kya"), cex.axis=1.5)

dogray <-   which(!pop %in% southern)
docol <-  which(pop %in% southern)
doboot <- which(names(res_boot) %in% inds[docol])

for(i in dogray) lines(x=res2[[inds[i]]]$YearsAgo, y=res2[[inds[i]]]$Ne,  col = "lightgray",lwd=3)
for(i in doboot) lines(x=res_boot[[i]]$YearsAgo, y=res_boot[[i]]$Ne,  col = alpha(cols[pop[names(res_boot)[i]]], 0.1),lwd=0.5,lty=1.5)

for(i in docol) lines(x=res2[[inds[i]]]$YearsAgo, y=res2[[inds[i]]]$Ne,  col = cols[pop[inds[i]]],lwd=3)

legend(x=xlims[2] * 1.4, y=ymax *1.25,
       title="Southern common", legend=southern, col=impala_colors[southern],
       lty=1,lwd=3,border=NA, bty="n",cex=1.5, xpd=NA)

###############################################
############# 3 plot selous ###############
###############################################
plot(type='n', x=1, log='x',
     xlab="",#sprintf("Years ago (mu=%.2e, g=%.1f)",mu,g),
     ylab="Ne",cex.lab=1.85, xaxt='n', yaxt='n',
     xlim=xlims, ylim=c(0, ymax), cex.axis=1.5, xpd=NA)
xpos <- c(0.5e4, seq(1e4,1e5,by=1e4), seq(1e5, 1e6, by=1e5))
xlabs <- c(0.5e4, 1e4, 2e4, 1e5, 5e5,1e6)

axis(2, at=seq(0, ymax, 1e4), labels=c("0",paste0(as.integer(seq(1e4, ymax, 1e4)/1e3), "k")), cex.axis=1.5)

axis(1, at=xpos, labels=F)
axis(1, at=xlabs, tick=F, labels= paste(as.integer(xlabs/1e3), "kya"), cex.axis=1.5)

dogray <-   which(!pop %in% selous)
docol <-  which(pop %in% selous)
doboot <- which(names(res_boot) %in% inds[docol])

for(i in dogray) lines(x=res2[[inds[i]]]$YearsAgo, y=res2[[inds[i]]]$Ne,  col = "lightgray",lwd=3)
for(i in doboot) lines(x=res_boot[[i]]$YearsAgo, y=res_boot[[i]]$Ne,  col = alpha(cols[pop[names(res_boot)[i]]], 0.1),lwd=0.5,lty=1.5)

for(i in docol) lines(x=res2[[inds[i]]]$YearsAgo, y=res2[[inds[i]]]$Ne,  col = cols[pop[inds[i]]],lwd=3)

legend(x=xlims[2] * 1.4, y=ymax,
       title="Selous common", legend=selous, col=impala_colors[selous],
       lty=1,lwd=3,border=NA, bty="n",cex=1.5, xpd=NA)


###############################################
############# 4 plot eastern ###############
###############################################
plot(type='n', x=1, log='x',
     xlab=sprintf("Years ago (mu=%.2e, g=%.1f)",mu,g),
     ylab="Ne",cex.lab=1.85, xaxt='n', yaxt='n',
     xlim=xlims, ylim=c(0, ymax), cex.axis=1.5, xpd=NA)
xpos <- c(0.5e4, seq(1e4,1e5,by=1e4), seq(1e5, 1e6, by=1e5))
xlabs <- c(0.5e4, 1e4, 2e4, 1e5, 5e5,1e6)

axis(2, at=seq(0, ymax, 1e4), labels=c("0",paste0(as.integer(seq(1e4, ymax, 1e4)/1e3), "k")), cex.axis=1.5)

axis(1, at=xpos, labels=F)
axis(1, at=xlabs, tick=F, labels= paste(as.integer(xlabs/1e3), "kya"), cex.axis=1.5)

dogray <-   which(!pop %in% eastern)
docol <-  which(pop %in% eastern)
doboot <- which(names(res_boot) %in% inds[docol])

for(i in dogray) lines(x=res2[[inds[i]]]$YearsAgo, y=res2[[inds[i]]]$Ne,  col = "lightgray",lwd=3)
for(i in doboot) lines(x=res_boot[[i]]$YearsAgo, y=res_boot[[i]]$Ne,  col = alpha(cols[pop[names(res_boot)[i]]], 0.1),lwd=0.5,lty=1.5)

for(i in docol) lines(x=res2[[inds[i]]]$YearsAgo, y=res2[[inds[i]]]$Ne,  col = cols[pop[inds[i]]],lwd=3)

legend(x=xlims[2] * 1.4, y=ymax * 1.25,
       title="Eastern common", legend=eastern, col=impala_colors[eastern],
       lty=1,lwd=3,border=NA, bty="n",cex=1.5, xpd=NA)


}



xlims <- c(0.4 * 10^4, 1.4 * 10^6)
outname <- "impalaPsmcsNice_gen57.png"
png(outname,width=10,height=6,res=300, units="in")
#ymax <- max(sapply(res2,function(x)max(x$Ne)))
doplot()
dev.off()


xlims <- c(0.4 * 10^4, 1.4 * 10^6)
outname <- "impalaPsmcsNice_gen57.pdf"
pdf(outname,width=10,height=6)
#ymax <- max(sapply(res2,function(x)max(x$Ne)))
doplot()
dev.off()



###### RE-DO PLOTS WITH GENERATION TIME OF 4.03
g <- 4.03 #5.7 #4.03 #5.7 #opt$generation

files <-list.files(indir,".psmc",full.names=T)

res <- lapply(files, psmc.result,mu=mu,g=g)
names(res) <- gsub(".psmc", "",basename(files))


inds <- as.character(info$Serial_number[info$Depth=='hi'])
pop <- gsub("Mburu", "Mburo", info$Locality[info$Depth=='hi'])
inds <- c(inds, "rgp")

pop <- c(pop, "Laikipia")
names(pop) <- inds


impala_colors["Laikipia"] <-  "#1d0989" 
cols <- impala_colors
#cols["Laikipia"] <-  "#1d0989" 

#res2 <- lapply(res, function(x) x[-((nrow(x)-8):nrow(x)),])


rm_last <- function(x){ ## FUNCITON TO CLEAN DATA, IS A BIT ARBITRARY MIGHT BE GOOD TO CHECK AGAIN SOME DAY?

    nes <- unique(x$Ne)
    rmv <- nes[(length(nes)-1):length(nes)]
    x[!x$Ne%in%rmv,]
}

#res2 <- lapply(res,rm_last)
res2 <- res

bf <- c("Etosha")
southern <- c("Ovita", "Chobe", "Shangani", "Mana Pools", "Kafue", "Luangwa")
selous <- c("Selous")
eastern <- c("Ugalla", "Masai Mara", "Tsavo", "Samburu", "Lake Mburo", "Laikipia")


# get path to all bootstrapls and format names to only keep sample name
bootfiles <- list.files(bootdir, ".psmc", full.names=T)
names(bootfiles) <- gsub("_[a-zA-Z0-9]*","", gsub(".psmc","",basename(bootfiles)))

# load, reescale and clean bootrap
res_boot <- lapply(bootfiles, psmc.result,mu=mu,g=g)

xlims <- c(0.4 * 10^4, 1. * 10^6)
outname <- "impalaPsmcsNice_gen4.png"
png(outname,width=10,height=6,res=300, units="in")
#ymax <- max(sapply(res2,function(x)max(x$Ne)))
doplot()
dev.off()


xlims <- c(0.4 * 10^4, 1. * 10^6)
outname <- "impalaPsmcsNice_gen4.pdf"
pdf(outname,width=10,height=6)
#ymax <- max(sapply(res2,function(x)max(x$Ne)))
doplot()
dev.off()
