
model=$1  ####model name

input=../do_sfs/results/2dsfs_splits_fscformat
nruns=50
nrunsparallel=10

# prepare folders to run each split
for split in `seq 1 $nruns`
do
  mkdir $model/split_${split}
  cp $input/split_${split}/Chobe-Etosha_split_${split}_folded_fscformat.sfs $model/split_${split}/_jointMAFpop1_0.obs
  cp $input/split_${split}/Etosha-Shangani_split_${split}_folded_fscformat.sfs $model/split_${split}/_jointMAFpop2_0.obs
  cp $input/split_${split}/Etosha-MasaiMara_split_${split}_folded_fscformat.sfs $model/split_${split}/_jointMAFpop3_0.obs
  cp $input/split_${split}/Chobe-Shangani_split_${split}_folded_fscformat.sfs $model/split_${split}/_jointMAFpop2_1.obs
  cp $input/split_${split}/Chobe-MasaiMara_split_${split}_folded_fscformat.sfs $model/split_${split}/_jointMAFpop3_1.obs
  cp $input/split_${split}/MasaiMara-Shangani_split_${split}_folded_fscformat.sfs $model/split_${split}/_jointMAFpop3_2.obs
  cp ${model}/${model}_bestlikes.pv $model/split_${split}/${model}_bestlikes.pv
done

# run splits in parallel
seq 1 $nruns | xargs -n 1 -P $nrunsparallel bash run_fsc_split.sh $model
