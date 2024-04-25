
pairs <- c("Chobe-Etosha", 
            "Etosha-Shangani",
            "Etosha-MasaiMara", 
            "Chobe-Shangani",
            "Chobe-MasaiMara", 
            "MasaiMara-Shangani")

sfsdir <- "../do_sfs/results/2dsfs/"

readSumGlobal <- function(p)
    sum(scan(paste0(sfsdir, p, ".sfs"), skip=1))

readSumSplits <- function(p){

    a <- readLines((paste0(sfsdir, p, "_split.sfs")))[seq(2, 100,2)]
    res <- sapply(a, function(x) sum(as.numeric(unlist(strsplit(x, " ")))))
    names(res) <- NULL
    return(res)
}


get_jk_est <- function(est_global, est_splits, n, ms, g=50)
    g * est_global - sum((1 - ms / n) * est_splits)

get_jk_sd <- function(est_global, est_splits, n, ms, g=50){

    hs <- n / ms
    jk_est <- get_jk_est(est_global, est_splits, n, ms, g)
    jk_var <- (1/g) * sum((1/(hs - 1)) * ((est_splits - jk_est) ^ 2))
    jk_se <- sqrt(jk_var)
    return(jk_se)
}
 
# number of sites acorss all 6 2dsfs
n <- sum(sapply(pairs, readSumGlobal))

# number of sites removes for each split across all 6 2dsfs
m_js <- colSums(do.call("rbind", lapply(pairs, readSumSplits)))

##################################
########## DO MODEL 2 ############
##################################
maxlike <- read.table("../runfsc/model2/max_like_run/5/5.bestlhoods", h=T)
jk_ests <- read.table("model2_jacknife50.params", h=T)

sds <- sapply(1:(ncol(jk_ests) - 2),
              function(x) get_jk_sd(est_global=maxlike[,x], est_splits=jk_ests[,x], n=n, ms=m_js, g=50)) 
sds2 <- sapply(1:(ncol(jk_ests) - 2),
              function(x) get_jk_sd2(est_global=maxlike[,x], est_splits=jk_ests[,x], n=n, ms=m_js, g=50)) 
 
res <- t(rbind(maxlike[,1:(ncol(jk_ests) - 2)], sds))

outtsv <- "model2_jackknife_se.tsv"
write.table(res, outtsv, col.names=c("maxres", "se"), quote=F, row.names=T, sep="\t")


##################################
########## DO MODEL 4 ############
##################################
maxlike <- read.table("../runfsc/model4/max_like_run/5/5.bestlhoods", h=T)
jk_ests <- read.table("model4_jacknife50.params", h=T)

sds <- sapply(1:(ncol(jk_ests) - 2),
              function(x) get_jk_sd(est_global=maxlike[,x], est_splits=jk_ests[,x], n=n, m=m_js, g=50)) 

res <- t(rbind(maxlike[,1:(ncol(jk_ests) - 2)], sds))

outtsv <- "model4_jackknife_se.tsv"
write.table(res, outtsv, col.names=c("maxres", "se"), quote=F, row.names=T, sep="\t")










if(FALSE){
    
    get_jk_sd2 <- function(est_global, est_splits, n, ms, g=50){

        hs <- n / ms
        jk_est <- get_jk_est(est_global, est_splits, n, ms, g)
        jk_var <- (1 / g) * sum((1 / (1-hs)) * (hs * est_global - (hs - 1) * est_splits - g * est_global + sum((1-(ms/n)) * est_splits)) ^ 2)
        #jk_var <- (1/g) * sum((1/(hs - 1)) * ((est_splits - jk_est) ^ 2))
        jk_se <- sqrt(jk_var)
        return(jk_se)
    }

    est_global <- maxlike[,1]
    est_splits <- jk_ests[,1]
    m=m_js
    g=50

    sds2 <- sapply(1:(ncol(jk_ests) - 2),
                   function(x) get_jk_sd2(est_global=maxlike[,x], est_splits=jk_ests[,x], n=n, ms=m_js, g=50)) 

}
