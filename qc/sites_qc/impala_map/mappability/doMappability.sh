GENMAP=/home/genis/software/genmap-build/bin/genmap

refdir=/davidData/genis/impala/ref/impala_draft
ref=$refdir/Impala.scaf.sorted
index=$refdir/index
out=/home/genis/impala/localSitesQC/impalaMap/mappability/mappability_ImpalaRef_K100_E2

#gzip -d $ref.gz
ln -s /davidData/genis/impala/ref/impala_draft/Impala.scaf.sorted /davidData/genis/impala/ref/impala_draft/Impala.scaf.sorted.fasta
export TMPDIR=/davidData/genis/tmp
# Do genmap reference index
$GENMAP index -F $ref.fasta -I $index -A skew

# DO mappability
$GENMAP map -K 100 -E 2 -I $index -O $out -t -w -bg -T 40

#awk '$4 == 1 {print $1"\t"$2"\t"$3}' $out.bedgraph > ${out}_good.bed
