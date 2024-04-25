
FSC27=/kellyData/home/genis/software/fsc27_linux64/fsc2702

model=$1
split=$2
nits=5
ncores=8

for i in `seq 1 $nits`
  do
  mkdir $model/split_${split}/$i
  cp $model/$model.est $model/split_${split}/$i/$i.est
  cp $model/$model.tpl $model/split_${split}/$i/$i.tpl
  cp $model/split_${split}/_jointMAFpop1_0.obs $model/split_${split}/$i/${i}_jointMAFpop1_0.obs
  cp $model/split_${split}/_jointMAFpop2_0.obs $model/split_${split}/$i/${i}_jointMAFpop2_0.obs
  cp $model/split_${split}/_jointMAFpop3_0.obs $model/split_${split}/$i/${i}_jointMAFpop3_0.obs
  cp $model/split_${split}/_jointMAFpop2_1.obs $model/split_${split}/$i/${i}_jointMAFpop2_1.obs
  cp $model/split_${split}/_jointMAFpop3_1.obs $model/split_${split}/$i/${i}_jointMAFpop3_1.obs
  cp $model/split_${split}/_jointMAFpop3_2.obs $model/split_${split}/$i/${i}_jointMAFpop3_2.obs
  cp $model/split_${split}/${model}_bestlikes.pv $model/split_${split}/$i/${model}_bestlikes.pv

  export DIR=$model/split_${split}/$i
  cd $DIR
  echo "start running split $split for model $model run $i"
  $FSC27 -t $i.tpl -n500000 --initValute ${model}_bestlikes.pv -e $i.est -M -L100 -m --cores $ncores --numBatches $ncores -C100 --seed ${i}
  rm $i.est
  rm $i.tpl
  rm ${i}_jointMAFpop1_0.obs
  rm ${i}_jointMAFpop2_0.obs
  rm ${i}_jointMAFpop3_0.obs
  rm ${i}_jointMAFpop2_1.obs
  rm ${i}_jointMAFpop3_1.obs
  rm ${i}_jointMAFpop3_2.obs
  cd ../../..
done
