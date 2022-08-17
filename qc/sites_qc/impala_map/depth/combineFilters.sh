BEDTOOLS=/home/genis/software/bedtools2/bin/bedtools

bedall=/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/impalaRef100kbScaffolds.bed

badlow=/home/genis/impala/localSitesQC/impalaMap/depth/bedFiles/remove_depth_lowDepth
badhigh=/home/genis/impala/localSitesQC/impalaMap/depth/bedFiles/remove_depth_highDepth

badall=/home/genis/impala/localSitesQC/impalaMap/depth/bedFiles/remove_depth_all

keep=/home/genis/impala/localSitesQC/impalaMap/depth/bedFiles/keep_depth.bed

# make sure everything is fine
$BEDTOOLS sort -i ${badlow}.bed > ${badlow}_sorted.bed
$BEDTOOLS merge -i ${badlow}_sorted.bed > ${badlow}_sorted_merged.bed

$BEDTOOLS sort -i ${badhigh}.bed > ${badhigh}_sorted.bed
$BEDTOOLS merge -i ${badhigh}_sorted.bed > ${badhigh}_sorted_merged.bed


# make list of bad sites combining both
cat $badlow.bed $badhigh.bed > ${badall}_unmerged_unsorted.bed

$BEDTOOLS sort -i ${badall}_unmerged_unsorted.bed > ${badall}_unmerged.bed
$BEDTOOLS merge -i ${badall}_unmerged.bed > ${badall}.bed

# make list of good sites to keep
$BEDTOOLS subtract -a $bedall -b $badall.bed > $keep

