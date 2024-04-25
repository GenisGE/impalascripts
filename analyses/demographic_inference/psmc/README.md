# Population size trajectoreis inferred with PMSC

We use PSMC to estiamte the trajecotires of population size throught time for all medium to high depth samples plus the high depth sample from the Ruminant Genome Project.

For each sample we call genotypes indivdiually and use them to run PSMC after doing some cleaning up. All steps (genotype calling, cleaning up, formatting, running psmc...) are included in the snakemake file (`run_psmc.snakefile`).
