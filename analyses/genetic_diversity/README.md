# Estimate heterozygosity (and relatedness ratios)

This pipeline estimates per sample genome-wide heterozygosity by estiamting the sample 1dSFS with ANGSD + realSFS. Heterozygosity is then the counts in heterozygous category divided by the total count.

It also estiamtes relatedness ratios from ibsRelate. It is done jointly because the intermediate files and steps are common. Realtedness results from ibsRelate were used to detect duplciated samples and exlude one duplicate from analyses.