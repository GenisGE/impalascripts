# Combine the bed files of differnt filters


## make bed with all regions in autosome and X chromosomes

```
awk '{print $1"\t"0"\t"$2}' /davidData/genis/impala/ref/goat/goat_ref_renamed.fa.gz.fai | head -31 > /home/genis/impala/localSitesQC/allRegionsFiles/allGoatAutoXchrs.bed

head -29 /home/genis/impala/localSitesQC/allRegionsFiles/allGoatAutoXchrs.bed > allGoatAutosomes.bed
```

## intersect the different filters

```
bash intersectFilters.sh
```

## convert from bed format to format used by angsd sites

```
awk '{print $1"\t"$2+1"\t"$3}' autosome_map_inb_rep_dep.bed > autosome_map_inb_rep_dep.region
angsd sites index autosome_map_inb_rep_dep.region


awk '{print $1"\t"$2+1"\t"$3}' all_map_inb_rep_dep.bed > all_map_inb_rep_dep.region
angsd sites index all_map_inb_rep_dep.region
```

## sumarize filters

Run interatively `summaryze_filters.R`. 
