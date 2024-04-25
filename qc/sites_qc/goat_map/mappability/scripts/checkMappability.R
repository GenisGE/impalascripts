library(data.table)

mp <- fread("/home/genis/impala/localSitesQC/mappability/mappability_V1Goat_K100_E2.bedgraph", data.table=F,h=F)

> colnames(mp) <- c("chr", "start", "end", "mappability")
> head(mp)
     chr start end mappability
1 CHR_01     0  66    1.000000
2 CHR_01    66  84    0.500000
3 CHR_01    84 176    1.000000
4 CHR_01   176 261    0.500000
5 CHR_01   261 267    0.333333
6 CHR_01   267 268    0.500000
> len
Error: object 'len' not found
> len <- function(x) sum(as.numeric(x$end-x$start))
> k <- mp$mappability == 1
> sum(k)
[1] 4218979
> sum(k)/ncol(mp)
[1] 1054745
> sum(k)/nrow(mp)
[1] 0.01783125
> nrow(mp[k,])
[1] 4218979
> len(mp[k,])/len(mp)
[1] 0.854158
> len(mp)
[1] 2919751498
> len(mp[k,])
[1] 2493929246
