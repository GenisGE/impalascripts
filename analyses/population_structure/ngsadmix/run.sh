NGSA=NGSadmixConv.sh

#echo 13 12 `seq 2 11` | xargs -n1 -P2 $NGSA
echo `seq 2 11` | xargs -n1 -P2 $NGSA

# rm /home/genis/impala/analyses/impalaMap/ngsadmix/results/*/admixResultImpalaMap.*fopt.gz
