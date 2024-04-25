PCANGSD=/home/genis/software/pcangsd2/pcangsd/pcangsd.py

bgl=/home/genis/impala/localSitesQC/goatMap/inbreedSites2/results/ImpalaGoatMappedCollapsedNoQCFiltersNo1BadNo5dupsOnlyCommonImpala.beagle.gz # there are 26472465 snps
#out=/home/genis/impala/localSitesQC/goatMap/inbreedSites2/results/impalaGoatMappedNo1BadNo5DupsOnlyCommonImpalaE6


#python3 $PCANGSD -beagle $bgl -sites_save -inbreedSites -o $out -threads 20 -e 6 > $out.log

out=/home/genis/impala/localSitesQC/goatMap/inbreedSites2/results/impalaGoatMappedNo1BadNo5DupsOnlyCommonImpalaE11

python3 $PCANGSD -beagle $bgl -sites_save -inbreedSites -o $out -threads 20 -e 11 > $out.log
