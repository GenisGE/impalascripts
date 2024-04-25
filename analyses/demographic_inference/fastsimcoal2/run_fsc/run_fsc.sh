
FSC27=/kellyData/home/genis/software/fsc27_linux64/fsc2702

model=$1
i=$2
ncores=8

mkdir $model/$i
cp $model/$model.est $model/$i/$i.est
cp $model/$model.tpl $model/$i/$i.tpl
cp $model/_jointMAFpop1_0.obs $model/$i/${i}_jointMAFpop1_0.obs
cp $model/_jointMAFpop2_0.obs $model/$i/${i}_jointMAFpop2_0.obs
cp $model/_jointMAFpop3_0.obs $model/$i/${i}_jointMAFpop3_0.obs
cp $model/_jointMAFpop2_1.obs $model/$i/${i}_jointMAFpop2_1.obs
cp $model/_jointMAFpop3_1.obs $model/$i/${i}_jointMAFpop3_1.obs
cp $model/_jointMAFpop3_2.obs $model/$i/${i}_jointMAFpop3_2.obs

cd $model/$i
echo "start running model $model run $i"
$FSC27 -t $i.tpl -n500000 -e $i.est -M -L100 -m --cores $ncores --numBatches $ncores -C100 --seed ${i}
cd ../../
