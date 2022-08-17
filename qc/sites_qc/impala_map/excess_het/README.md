

# filter to remove regions with extreme excess of heterozygosity, usually caused by mapping errors due to paralogs regions. Steps:

1. make beagle file from all the genome: `makeBeagleRaw.sh`

2. Do site specific HWE test accounting for popualtion structure: `doHWEpcangsd.sh`

3. Use pcangsd hwe test output to remove windows showing evidence for extreme excess of heterozygosity: `doHWEfilter.R`

4. Make bed file with regions to keep using output of 3: `doKeepList.sh`