# Piepline to estiamte D-statistics from called genotypes

Takes a vcf file as input, cleans it and does format conversions. It also converts the reference sequence to gentoype format (both binary plink and vcf; will be homozygous at all positions since it is from the fasta) to be used as outgroup from the Dstats.

It produces a vcf file that is used to estimate D-statistics, summaryize as fbranch and plot it wiht Dsuit, and also to bed format which is used as input for qpgraph ADMIXTOOLS2 (see qpgraph folder).