# Pipelines for genotyping individauls based on different approahces

## Genotype likelihoods

Done jointly on all inidivuals, estimate genotype likelihoods with ANGSD

## Genotype calling

Done only on medium and high depth samples, called gentoypes with bcftools. We used the `a1kg gentoype calling` pipeline.

## Imputation

Done on all samples, used bcftools to call variants and estimate genotype likelihoods and then BEAGLE for imputing common variants.