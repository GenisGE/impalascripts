# adapted from script made by Kristian for leoaprd in /home/leopard/users/krishang/bed_region/get_nonrepeatmask.sh

ref=/davidData/genis/impala/ref/goat/goat_ref_renamed.fa
bed=/home/genis/impala/localSitesQC/repeatMasker/impalaNonRepeat.bed

/home/krishang/software/vir_python36/bin/python /home/leopard/users/krishang/bed_region/scripts/get_upper.py $ref  |  /home/krishang/software/bedops/bin/bedops -m - > $bed

#awk '{print $1,0,$2}' /davidData/genis/impala/ref/goat/goat_ref_renamed.fa.gz.fai > ENTIRE_GENOME.bed

#~/software/bedops/bin/bedops -d <(head -n1 ENTIRE_GENOME.bed ) nonrepeat.bed > repeat_NW_017619845.1.bed
