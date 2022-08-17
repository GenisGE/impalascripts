PCANGSD=/home/genis/software/pcangsd2/pcangsd/pcangsd.py

bgl=/home/genis/impala/localSitesQC/impalaMap/inbreedSites/results/ImpalaImpalaMappedCollapsedNoQCFiltersNo1BadNo5dupsOnlyCommonImpala.beagle.gz # there are 15392266 snps
out=/home/genis/impala/localSitesQC/impalaMap/inbreedSites/results/impalaImpalaMappedCollapsedNoQCFiltersNo1BadNo5dupsOnlyCommonImpalae6

python3 $PCANGSD -beagle $bgl -sites_save -inbreedSites -o $out -threads 20 -e 6 > $out.log
