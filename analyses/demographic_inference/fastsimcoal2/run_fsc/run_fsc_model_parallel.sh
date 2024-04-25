
model=$1  ####model name
nrunsparallel=$2 ##### number of runs to do in parallel, each run uses 8 threads so total cpus used will be nrunsparallel * 8

input=../do_sfs/results/2dsfs_fscformat
nruns=50

# prepare folders to run each split
cp $input/Chobe-Etosha_folded_fscformat.sfs $model/_jointMAFpop1_0.obs
cp $input/Etosha-Shangani_folded_fscformat.sfs $model/_jointMAFpop2_0.obs
cp $input/Etosha-MasaiMara_folded_fscformat.sfs $model/_jointMAFpop3_0.obs
cp $input/Chobe-Shangani_folded_fscformat.sfs $model/_jointMAFpop2_1.obs
cp $input/Chobe-MasaiMara_folded_fscformat.sfs $model/_jointMAFpop3_1.obs
cp $input/MasaiMara-Shangani_folded_fscformat.sfs $model/_jointMAFpop3_2.obs

#cp ${model}/${model}_bestlikes.pv $model/split_${split}/${model}_bestlikes.pv


# run splits in parallel
seq 1 $nruns | xargs -n 1 -P $nrunsparallel bash run_fsc.sh $model
