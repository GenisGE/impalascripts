# remove outlier regions regarding depth across all samples (too low or too high).

A more reproducible version of this pipeline can be found here:

https://github.com/GenisGE/warthogscripts/tree/main/analyses/qc/sites_qc/depth



Steps:

1. Get depth for each site, separately for low depth and high depth samples, in parallel across chromosomes: `getDepthsLowHighDepthSamples.sh`

2. Combine global depth distribution across chromosomes: `combineDepths.R`

3. Plot global depth distribution and fitler thresholds (based on medina globla depth): `plotDepths.R`

4. Make bed files with regions to remove based on min and max depth thresholds: `doDepthFilter.py`

5. Combine differnet groups and produce final bed file with sites to keep after depth fitlering: `combineFilters.sh`
 
