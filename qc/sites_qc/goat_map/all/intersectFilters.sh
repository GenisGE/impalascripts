BEDTOOLS=/home/genis/software/bedtools2/bin/bedtools

all=allGoatAutoXchrs.bed
auto=allGoatAutosomes.bed
map=mappability_V1Goat_K100_E2_mappability1.bed
rep=impalaNonRepeat.bed
inb=keeplistInbreedSites_remapped.BED
dep=keep_depth.bed



$BEDTOOLS intersect -a $all -b $map > all_map.bed
$BEDTOOLS intersect -a all_map.bed -b $inb > all_map_inb.bed
$BEDTOOLS intersect -a all_map_inb.bed -b $rep > all_map_inb_rep.bed
$BEDTOOLS intersect -a all_map_inb_rep.bed -b $dep > all_map_inb_rep_dep.bed


$BEDTOOLS intersect -a $auto -b $rep > autosome_rep.bed
$BEDTOOLS intersect -a $auto -b $inb > autosome_inb.bed
$BEDTOOLS intersect -a $auto -b $map > autosome_map.bed
$BEDTOOLS intersect -a $auto -b $dep > autosome_dep.bed
$BEDTOOLS intersect -a autosome_inb.bed -b $dep > autosome_inb_dep.bed
$BEDTOOLS intersect -a autosome_inb_dep.bed -b $rep > autosome_inb_dep_rep.bed
$BEDTOOLS intersect -a autosome_map.bed -b $inb > autosome_map_inb.bed
$BEDTOOLS intersect -a autosome_map.bed -b $dep > autosome_map_dep.bed
$BEDTOOLS intersect -a autosome_map_inb.bed -b $dep > autosome_map_inb_dep.bed
$BEDTOOLS intersect -a autosome_map_dep.bed -b $rep > autosome_map_dep_rep.bed
$BEDTOOLS intersect -a autosome_map_inb.bed -b $rep > autosome_map_inb_rep.bed
$BEDTOOLS intersect -a autosome_map_inb_rep.bed -b $dep > autosome_map_inb_rep_dep.bed


$BEDTOOLS subtract -a $all -b $map > remove/map_remove.bed
$BEDTOOLS subtract -a $all -b $rep > remove/rep_remove.bed
$BEDTOOLS subtract -a $all -b $inb > remove/inb_remove.bed

$BEDTOOLS intersect -a  remove/map_remove.bed -b remove/rep_remove.bed > remove/intersect_map_rep_remove.bed
$BEDTOOLS intersect -a  remove/map_remove.bed -b remove/inb_remove.bed > remove/intersect_map_inb_remove.bed
$BEDTOOLS intersect -a  remove/rep_remove.bed -b remove/inb_remove.bed > remove/intersect_rep_inb_remove.bed
$BEDTOOLS intersect -a  remove/intersect_map_rep_remove.bed -b remove/inb_remove.bed > remove/intersect_map_rep_inb_remove.bed

