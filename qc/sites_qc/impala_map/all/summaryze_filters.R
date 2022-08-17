library(data.table)


clean <- function(x) x <- x[x$V1 %in% allregions$V1,]
nox <- function(x) x <-  x[x$V1 %in% autosome$V1,]
len <- function(x) sum(as.numeric(x$V3 - x$V2))


allregions <- fread("/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/impalaRef100kbScaffolds.bed", h=F, data.table=F)
autosome <- fread("/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/autosomeScaffoldsImpala.bed", h=F, data.table=F)


allfilts <- fread("autosome_map_inb_rep_dep.bed",h=F,data.table=F)


inbreed <- fread("/home/genis/impala/localSitesQC/impalaMap/inbreedSites/results/keeplistInbreedSites_impalaRef.BED",h=F,data.table=F)
repmask <- fread("/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/impalaNonRepeat.bed",h=F,data.table=F)
repmask <- clean(repmask)
map <- fread("/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/mappability_ImpalaRef_K100_E2_good.bed",h=F, data.table=F)
map <- clean(map)
depth <- fread("/home/genis/impala/localSitesQC/impalaMap/depth/bedFiles/keep_depth.bed", h=F, data.table=F)


### check here general if a feel like. then keep only autosomes





#### keep only autosomes
inbreed <- nox(inbreed)
repmask <- nox(repmask)
map <- nox(map)
depth <- nox(depth)

kept <- c(len(inbreed), len(repmask), len(map), len(depth), len(allfilts))/len(autosome)
names(kept) <- c("Inbreeding", "RepeatMasker", "Mappability", "Depth", "All")


# this is last one when removing F < -0.95
#  Inbreeding RepeatMasker  Mappability        Depth          All 
#   0.9596420    0.6147814    0.9060962    0.9169945    0.5421016 



# previous when inbreeding removed F < -0.5
#  Inbreeding RepeatMasker  Mappability        Depth          All 
#   0.8604423    0.6147814    0.9060962    0.9169945    0.4905401 



# Same with goat for comparison from:

#  Inbreeding RepeatMasker  Mappability        Depth          All 
#   0.6588008    0.6121615    0.9239331    0.7684793    0.3557530 
# len(allfilts)
# 877355094


### if I feel like look properly at intersection between different filters
