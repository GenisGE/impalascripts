library(data.table)


clean <- function(x) x <- x[x$V1 %in% allregions$V1,]
nox <- function(x) x <-  x[x$V1 %in% autosome$V1,]
len <- function(x) sum(as.numeric(x$V3 - x$V2))


allregions <- fread("allGoatAutoXchrs.bed", h=F, data.table=F)
autosome <- fread("allGoatAutosomes.bed", h=F, data.table=F)


allfilts <- fread("autosome_map_inb_rep_dep.bed",h=F,data.table=F)


inbreed <- fread("keeplistInbreedSites_remapped.BED",h=F,data.table=F)
repmask <- fread("impalaNonRepeat.bed",h=F,data.table=F)
repmask <- clean(repmask)
map <- fread("mappability_V1Goat_K100_E2_mappability1.bed",h=F, data.table=F)
map <- clean(map)
depth <- fread("keep_depth.bed", h=F, data.table=F)


### check here general if a feel like. then keep only autosomes





#### keep only autosomes
inbreed <- nox(inbreed)
repmask <- nox(repmask)
map <- nox(map)
depth <- nox(depth)

kept <- c(len(inbreed), len(repmask), len(map), len(depth), len(allfilts))/len(autosome)
names(kept) <- c("Inbreeding", "RepeatMasker", "Mappability", "Depth", "All")

#> kept
#  Inbreeding RepeatMasker  Mappability        Depth          All 
#   0.7983774    0.6121615    0.9239331    0.7684793    0.4371494 


### if I feel like look properly at intersection between different filters
