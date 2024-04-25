# Pipeline used to mapping the data to reference genome

The snakemake file contains all steps from qc of the fastq to mapping, fitlering and QC the mapped bam files.

The data was mapped to two genomes, a scaffold-level impala genome and a chromosome-level goat genome. 

- Impala genome (IMP): https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_006408695.1/
- Goat Genome (ARS1): https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001704415.2/

The pipelines used for mapping are identical for both genomes, just need to set the `REF` variable in the snakemake file to the coresponding genome.