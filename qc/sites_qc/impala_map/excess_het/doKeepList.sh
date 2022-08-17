BEDTOOLS=/home/genis/software/bedtools2/bin/bedtools

bedall=/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/impalaRef100kbScaffolds.bed

exclude=/home/genis/impala/localSitesQC/impalaMap/inbreedSites/results/excludelistInbreedSites_minF-0.5_rm20000k_remapped
keep=/home/genis/impala/localSitesQC/impalaMap/inbreedSites/results/keeplistInbreedSites_impalaRef_oldminf05

$BEDTOOLS sort -i ${exclude}_unmerged_unsorted.BED > ${exclude}_unmerged.BED
$BEDTOOLS merge -i ${exclude}_unmerged.BED > ${exclude}.BED

$BEDTOOLS subtract -a $bedall -b $exclude.BED > $keep.BED

