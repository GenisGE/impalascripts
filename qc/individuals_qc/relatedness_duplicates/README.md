# Detect duplciates using ibsRelate

Implements allele frequency free ratios of relatedness from ibsRelate. Based on estimating 2dSFS between individiual samples (with AGNSD + realSFS), then estiamte three ratios (R0, R1 and KING) which are all functions of the sample 2dSFS. They were used to idenfiy duplciate smaples.

The pipeline also estimates per sample heterzogysoity since the steps are largely shared. The heterozygosities were not sued here (see `analyses/genetic_diversity`).
