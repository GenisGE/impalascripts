# combine all filters

Remove short scaffolds (<1e5 bp):

```
awk '$2 > 1e5 {print $1"\t"0"\t"$2}' /davidData/genis/impala/ref/impala_draft/Impala.scaf.sorted.gz.fai > impalaRef100kbScaffolds.bed
```

1. make bed and angsd files with regions to keep `/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/scripts/intersectFilters.sh`

2. check how much is removed by different filters, and how much is kept overall: `summaryze_filters.R`
