ANGSD=/home/genis/software/angsd/angsd

dir=/home/genis/impala/analyses/impalaMap/beagle
chrs=/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/impalaRefAutosomeScaffoldsUsed.txt
bgl=$dir/impalaV3No1badNo5dupsImpalaMap
bams=/home/genis/impala/info_files/impalaMap/bamsImpalaMapV2No1badNo5dups.list
sites=/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/autosome_map_inb_rep_dep.regions

nscaff=`wc -l $chrs | cut -f1 -d" "`


for n in `seq 1 $nscaff`
do
    chr=`head -$n $chrs | tail -1`
    bgl=tmpbgl_$n
    echo "$ANGSD -GL 2 -out $bgl -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -sites $sites -bam $bams -minmapQ 30 -minQ 30 -r $chr -minmaf 0.05"
done | xargs -L1 -I CMD -P20 nice bash -c CMD


echo "Finished estiamting genotype likelihoods going to concatenate all scaffolds in $bgl.beagle"

zcat tmpbgl_1.beagle.gz | head -1 > $bgl.beagle
for i in `seq 1 $nscaff`
do
    zcat tmpbgl_$i.beagle.gz | sed '1d' >> $bgl.beagle
    cat tmpbgl_${i}.arg >> $bgl.arg
    echo "done scaffold $i"
done

gzip $bgl.beagle

rm tmpbgl_*
