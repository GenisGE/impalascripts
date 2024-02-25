# stolen and modified from kristian snakemake in /home/leopard/users/krishang/old/relatedness/run_2dsfs_cleanref.snakefile

import itertools as it
ANGSD="/home/genis/software/angsd"

REF="/davidData/genis/impala/ref/impala_draft/Impala.scaf.sorted.gz"
BASEQ=30

BAM_IN_DIR="/davidData/genis/impala/qc_mapping_collapse/impalaRef/bam"

SAMPLES = glob_wildcards(os.path.join(BAM_IN_DIR, "{s}.bam")).s

MY_SITES = "/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/autosome_map_inb_rep_dep.regions"
CHROMS="/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/impalaRefAutosomeScaffoldsUsed.txt"

configfile: "config_sfshetrel.yaml"

OUTMAIN="results"
#print(list(it.combinations(SAMPLES, 2)))

wildcard_constraints:
    pop = "|".join(config.keys())

rule all:
    input:
        expand(os.path.join(OUTMAIN, "safs", "{s1}.saf.idx"), s1=SAMPLES),
        #os.path.join(OUTMAIN, "sfs_2d", "collected.txt"),
        #os.path.join(OUTMAIN, "sfs", "collected.txt"),
    #    os.path.join(OUTMAIN,"finished2DSts.txt")



rule all_saf_hets:
     input:
        os.path.join(OUTMAIN, "sfs", "collected.txt")

rule per_sample_saf:
    input:
        bam=os.path.join(BAM_IN_DIR, "{s}.bam"),
    output:
        saf_idx = os.path.join(OUTMAIN, "safs", "{s}.saf.idx"),
        saf = os.path.join(OUTMAIN, "safs", "{s}.saf.gz"),
        saf_pos = os.path.join(OUTMAIN, "safs", "{s}.saf.pos.gz"),
        arg = os.path.join(OUTMAIN, "safs", "{s}.arg"),
    params:
        outbase = lambda wildcards, output: output.saf_idx.replace(".saf.idx", "")
    threads: 3
    shell:
        "{ANGSD}/angsd -i {input.bam} -out {params.outbase} -minQ {BASEQ} -minMapQ 30 -dosaf 1 -sites {MY_SITES} -rf {CHROMS} -anc {REF} -GL 2"

rule sfs:
    input:
        saf_idx1 = os.path.join(OUTMAIN, "safs", "{s1}.saf.idx"),
    output:
        os.path.join(OUTMAIN, "sfs", "{s1}.sfs")
    threads: 15
    log:
        os.path.join(OUTMAIN, "sfs", "{s1}.log")
    shell:
        "{ANGSD}/misc/realSFS -P {threads} {input.saf_idx1} > {output} 2> {log}"



rule collect_sfs:
    input:
        expand(os.path.join(OUTMAIN, "sfs", "{s1}.sfs"), s1=SAMPLES)
    output:
        f=os.path.join(OUTMAIN, "sfs", "collected.txt")
    run:
        import os
        import pandas as pd
        data = []
        names = []
        for x in input:
            name = os.path.basename(x).replace(".sfs", "")
            names.append(name)
            with open(x, 'r') as fh:
                t = fh.read()
                data.append([float(x) for x in t.rstrip().split()])
        a = pd.DataFrame(data, index=names, columns=["aa","ad","dd"])
        a["het"] = a["ad"] / a.sum(1)
        a.to_csv(output.f, index=True, header=True, index_label="id", sep=" ")



rule sfs_2d_bypop:
    input:
        saf_idx1 = os.path.join(OUTMAIN, "safs", "{s1}.saf.idx"),
        saf_idx2 = os.path.join(OUTMAIN, "safs", "{s2}.saf.idx"),
    output:
        os.path.join(OUTMAIN, "sfs_2d", "{s1}_{s2}.sfs")
    threads: 15
    log:
        os.path.join(OUTMAIN, "sfs_2d", "{s1}_{s2}.pre_log")
    shell:
        "{ANGSD}/misc/realSFS -P {threads} {input.saf_idx1} {input.saf_idx2} > {output} 2> {log}"

rule clean_sfs_log:
    input:
        log = os.path.join(OUTMAIN, "sfs_2d", "{s1}_{s2}.pre_log")
    output:
        log = os.path.join(OUTMAIN, "sfs_2d", "{s1}_{s2}.log")
    shell: "tail -100 {input} > {output}"
        

rule collect_2dsfs_bypop:
    input:
       # lambda wildcards: [OUTMAIN+"sfs_2d/{s[0]}_{s[1]}.sfs"]
        lambda wildcards: expand(os.path.join(OUTMAIN, "sfs_2d", "{s[0]}_{s[1]}.sfs"), s=it.combinations(config[wildcards.pop], 2)),
        lambda wildcards: expand(os.path.join(OUTMAIN, "sfs_2d", "{s[0]}_{s[1]}.log"), s=it.combinations(config[wildcards.pop], 2))
    output:
        f=os.path.join(OUTMAIN, "sfs_2d", "{pop}_collected.txt")
    run:
        import os
        import pandas as pd
        data = []
        names = []
        for x in input:
            name = os.path.basename(x).replace(".sfs", "")
            names.append(name)
            with open(x, 'r') as fh:
                t = fh.read()
                data.append([float(x) for x in t.rstrip().split()])
        a = pd.DataFrame(data, index=names, columns=["aaAA","aaAD","aaDD","adAA","adAD","adDD","ddAA","ddAD","ddDD"])

        a["r0"] = (a["aaDD"]+a["ddAA"])/a["adAD"]
        a["r1"] = a["adAD"] / (a.iloc[:,[1,2,3,5,6,7]].sum(1))
        a["king"] = (a["adAD"] - 2*a[["aaDD", "ddAA"]].sum(1)) / (a.iloc[:,[1,3,5,7]].sum(1) + 2*a["adAD"])
        a.to_csv(output.f, index=True, header=True, index_label="id", sep=" ")




rule check_2dsfs:
     input:
        lambda wildcards: expand(os.path.join(OUTMAIN, "sfs_2d", "{p}_collected.txt"), p=config.keys())
     output:
        os.path.join(OUTMAIN,"finished2DSts.txt")
     shell: "echo done > {output}"
