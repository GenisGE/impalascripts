# Demographic inference of joint demographic history based on pariwise 2DSFS between populations


## SFS estiamtion

We inferred the joint demographic history of 4 populations using the pairwise 2DSFS between all of them (so 6 2DSFS in total). In the folder `do_sfs` you can find the pipeline used to infer the 2DSFS from bam file using ANGSD + winsfs. We inferred both genome-wide 2DSFS and 50 splits for each to do jackknife.

## Demographic inference

We tested 4 different models (1 to 4), see manuscript for details. The folder `do_sfs` contains the input 2DSFSs (the six files `_jointMAFpop{i}_{j}.obs` for 2dSFS between pop i and j) and the model specificiation and paramter search range files (`model{k}.tpl` and `model{k}.est` for model k from 1 to 4).

## Uncertainty estimation

We used a block-jacknife approach to estimate uncertainty on the parameter estiamtes fitted to models 2 and 4. In `do_sfs` pipeline we use winsfs split to estiamte sfs in 50 different block splits, and generate corresponding 50 leave one block out 2DSFS estimates. In the folder `jackknife` we refit the models 50 times using the genbome-wide maximum likelihood estimates as initial guesses (`model{k}_bestlikes.pv`) and use the jackknife parameter estiamtes to estimate standard errors.
