PCANGSD=/home/genis/software/pcangsd2/pcangsd/pcangsd.py

bgl=/home/genis/impala/analyses/impalaMap/beagle/impalaV2No1badNo5dupsImpalaMap.beagle.gz # there are 6221853 snps
out=/home/genis/impala/analyses/impalaMap/pcangsd/results/impalaImpalaMappedRepeatMapInbDepFiltersNo1BadNo5dupse7


#python3 $PCANGSD -beagle $bgl -sites_save -inbreedSites -kinship -inbreed 1 -tree -admix -admix_save -admix_fst -admix_tree -o $out -threads 20 -e 7 > $out.log


#for e in `seq 2 6`
#for e in 1
#do   
#    out=/home/genis/impala/analyses/goatMapV2/pcangsd/results/impalaV2GoatMappedRepeatMapInbDepFiltersNo1Bad_testKin_e$e
#    python3 $PCANGSD -beagle $bgl -kinship -o $out -threads 20 -e $e > $out.log
#done

bgl=/home/genis/impala/analyses/impalaMap/beagle/impalaV2No1badNo5dupsImpalaMap.beagle.gz
out=/home/genis/impala/analyses/impalaMap/pcangsd/results/impalaImpalaMappedRepeatMapInbDepFiltersNo1BadNo5dupse10


#python3 $PCANGSD -beagle $bgl -sites_save -inbreedSites -kinship -inbreed 1 -tree -admix -admix_save -admix_fst -admix_tree -o $out -threads 20 -e 10 > $out.log


bgl=/home/genis/impala/analyses/impalaMap/beagle/impalaV2No1badNo5dupsImpalaMap.beagle.gz
out=/home/genis/impala/analyses/impalaMap/pcangsd/results/impalaImpalaMappedRepeatMapInbDepFiltersNo1BadNo5dupse12


python3 $PCANGSD -beagle $bgl -sites_save -inbreedSites -kinship -inbreed 1 -tree -admix -admix_save -admix_fst -admix_tree -o $out -threads 20 -e 12 > $out.log
