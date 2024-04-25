BEDTOOLS=/home/genis/software/bedtools2/bin/bedtools

bedall=/home/genis/impala/localSitesQC/goatMap/allRegionsFiles/allGoatAutoXchrs.bed

exclude=/home/genis/impala/localSitesQC/goatMap/inbreedSites2/results/excludelistInbreedSites_minF-0.9_rm20000
keep=/home/genis/impala/localSitesQC/goatMap/inbreedSites2/results/keeplistInbreedSites_minF-0.9_rm20000k

$BEDTOOLS sort -i ${exclude}_unmerged_unsorted.BED > ${exclude}_unmerged.BED
$BEDTOOLS merge -i ${exclude}_unmerged.BED > ${exclude}.BED

$BEDTOOLS subtract -a $bedall -b $exclude.BED > $keep.BED

