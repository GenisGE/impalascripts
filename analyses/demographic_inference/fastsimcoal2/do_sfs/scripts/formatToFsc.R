args <- commandArgs(trailingOnly=TRUE)

insfs <- args[1]
outsfs <- args[2]
p1 <- args[3]
p2 <- args[4]

popord <- c("Etosha", "Chobe", "Shangani", "MasaiMara")

idx1 <- which(popord==p1) - 1
idx2 <- which(popord==p2) - 1

# here, p1 will be rows and p2 columns
sfs <- as.matrix(read.table(insfs))

rownames(sfs) <- paste0("d",idx1,"_",0:(nrow(sfs)-1))
colnames(sfs) <- paste0("d",idx2,"_",0:(ncol(sfs)-1))


if(idx2>idx1) sfs <- t(sfs)


cat("1 observation", file=outsfs)
write.table(sfs, outsfs, quote=F, sep="\t", col.names=NA, append=T)
