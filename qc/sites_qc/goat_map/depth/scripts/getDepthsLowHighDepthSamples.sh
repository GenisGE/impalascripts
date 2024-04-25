# modified from frederik scripts /home/leopard/users/ffs/depth/scripts/02_lowdepthsamples_dept_pr_scaffold.sh and /home/leopard/users/ffs/depth/scripts/03_highdepthsamples_dept_pr_scaffold.sh
ANGSD=/home/genis/software/angsd/angsd

BAMSLOW=/home/genis/impala/info_files/bamsLowDepthImpalaGoatMappedV2Nobad.txt
BAMSHIGH=/home/genis/impala/info_files/bamsHighDepthImpalaGoatMappedV2Nobad.txt
OUT=/home/genis/impala/localSitesQC/depth/output
CHRS=/home/genis/impala/info_files/goatChrsAutoandX.txt


cat $CHRS | while read chr
do
    echo $ANGSD -minQ 30 -minMapQ 30 -doCounts 1 -doDepth 1 -dumpCounts 1 -maxdepth 5000 -b $BAMSLOW -out ${OUT}/lowDepth/${chr}_lowDepth -r $chr
    echo $ANGSD -minQ 30 -minMapQ 30 -doCounts 1 -doDepth 1 -dumpCounts 1 -maxdepth 5000 -b $BAMSHIGH -out ${OUT}/highDepth/${chr}_highDepth -r $chr
done | xargs -L 1 -I CMD -P 10 bash -c CMD

