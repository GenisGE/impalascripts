# Pipeline for estimating SFS of 4-fold and 0-fold degenerate sites in coding regions from low-depth sequencing data

Estimates the SFS in different popualtions in the 4-fold and 0-fold degenerate sites, which is used then to estiamte the ration of the total number of segregating variation that is used as a measure of the relative efficiency of selection across populations.

First uses the fasta reference genome and a gff annotation to divide coding sequence in 0-fold, 2-fold and 4-fold degenerate. From there just use normal low-depth SFS estimation pipeline with ANGSD saf and realSFS to estiamte the SFS for the different sites. 

The ratio reported in the manuscript is calculated from the SFSs in the plotting stage (see figures/figure6).