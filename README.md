# Scripts used for the paper "Extensive population structure highlights an apparent paradox of stasis in the impala (Aepyceros melampus)"

The data used in the paper is whole genome sequencing of 119 impala samples (105 samples were low-depth sequenced at 2-4X and 14 samples were sequenced to medium-high depth at 7-18X). The samples are from 13 different populations across its range in sub Saharian Africa and is publicly available in SRA with project ID PRJNA862915. We used several population genomic analyses to characterize impala's population structure, genetic diversity, and the phylogeographic and evoltuionary factors that have shaped it.

This repository contains the code used to perform analysis. It is divided in 5 different folders:

- `mapping` contains the pipeline used for mapping the short read data to the reference genome.
- `qc` contains different pipeliens used for filtering the reference genome to mask regions more likely to contain mapping and gentoyping errors.
- `resources` contains several scripts and tables with information such as the masks resulting from the qc pipeline and sample metadata such as location, coordinates...
- `analyses` contains code and pipelines for the different analyses.
- `figures` contains code used to generate the figures.

See README within each subfolder for further information.

## Reference

Extensive population structure highlights an apparent paradox of stasis in the impala (Aepyceros melampus)
Genís Garcia-Erill, Xi Wang, Malthe S. Rasmussen, Liam Quinn, Anubhab Khan, Laura D. Bertola, Cindy G. Santander, Renzo F. Balboa, Joseph O. Ogutu, Patricia Pecnerova, Kristian Hanghoej, Josiah Kuja, Casia Nursyifa, Charles Masembe, Vincent Muwanika, Faysal Bibi, Ida Moltke, Hans R. Siegismund, Anders Albrechtsen, Rasmus Heller
bioRxiv 2024.04.19.590257; doi: https://doi.org/10.1101/2024.04.19.590257
