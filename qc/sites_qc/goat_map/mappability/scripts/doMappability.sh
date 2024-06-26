GENMAP=/home/genis/software/genmap-build/bin/genmap

refdir=/davidData/genis/impala/ref/goat
ref=$refdir/goat_ref_renamed.fa
index=$refdir/index
out=/home/genis/impala/localSitesQC/mappability/mappability_V1Goat_K100_E2

#gzip -d $ref.gz
#export TMPDIR=/davidData/genis/tmp
# Do genmap reference index
#$GENMAP index -F $ref -I $index -A skew

# DO mappability
#$GENMAP map -K 100 -E 2 -I $index -O $out -t -w -bg -T 40

awk '$4 == 1 {print $1"\t"$2"\t"$3}' $out.bedgraph > ${out}_good.bed
