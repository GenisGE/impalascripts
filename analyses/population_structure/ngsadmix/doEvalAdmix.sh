EVALADMIX=/home/genis/software/evalAdmix/evalAdmix

bgl=/home/genis/impala/analyses/impalaMap/beagle/impalaV3No1badNo5dupsImpalaMap.beagle.gz

resdir=/home/genis/impala/analyses/impalaMap/ngsadmix/results
outdir=/home/genis/impala/analyses/impalaMap/ngsadmix/evaladmix
outname=evalAdmixResultsImpalaMap

for k in `seq 2 12`
do
    qfile=`find $resdir/$k | grep qopt_conv`
    ffile=`find $resdir/$k | grep fopt_conv`
    out=$outdir/${outname}_K${k}
    $EVALADMIX -beagle $bgl -qname $qfile -fname $ffile -P 20 -o $out.corres &> $out.log
done
