import pandas as pd
import itertools as it


BEDTOOLS="/home/genis/software/bedtools2/bin/bedtools"
ANGSDIR="/home/genis/software/angsd"
ANGSD=os.path.join(ANGSDIR, "angsd")
REALSFS=os.path.join(ANGSDIR, "misc", "realSFS")
WINSFS="winsfs"

SCRIPTSDIR=config["scriptsdir"]
DEGENERATE=os.path.join(SCRIPTSDIR, "make_degeneracy_bed.py")
PLOTSFS=os.path.join(SCRIPTSDIR, "plotDegenerateSFSs.R")

PYTHON="python3"
R="Rscript"

OUTMAIN = config["outmain"]

# laod inds paths to get individual heterozygositeis
INDS={}
for p in config['pops'].keys():
    with open(config['pops'][p]) as f:
        INDS[p] = {}
        for b in [x.rstrip() for x in f.readlines()]:
            ind = os.path.basename(b).replace('.bam','')
            INDS[p][ind] = b 
        


rename_chrs = []
with open(config["rename_chrs"], "r") as fh:
    for line in fh.readlines():
        rename_chrs.append(line.strip().split())


        
rule all:
    input:
        #os.path.join(OUTMAIN, "sfs_indivs", "collected_heterozygosities_all.txt"),
        expand(os.path.join(OUTMAIN, "2dsfs_winsfs", "{p[0]}_{p[1]}_{d}f_degenerate_sites.sfs"), p = it.combinations(list(config["pops"].keys()), 2), d=[0,2,4]),
        #os.path.join(OUTMAIN, "fst", "all_fst.txt"),
        #expand(os.path.join(OUTMAIN, "sfs_indivs", "collected_heterozygosities_{p}.txt"), p=INDS.keys()),
        expand(os.path.join(OUTMAIN, "sfs_winsfs", "{p}_{d}f_degenerate.sfs"), p = config["pops"].keys(), d=[0,2,4]),
        #expand(os.path.join(OUTMAIN, "beds", "good_sites_{d}f_degenerate.sites.idx"), d=[0,2,4]),
        expand(os.path.join(OUTMAIN, "plot", "{p}_degenerate_sfs{w}.png"), p=config["pops"].keys(),
               w=["", "_nofixedderived"]),
        expand(os.path.join(OUTMAIN, "plot_winsfs", "{p}_degenerate_sfs{w}.png"), p=config["pops"].keys(),
               w=["", "_nofixedderived"])


rule select_gff_cds:
    input:
        gff = config["gff"]
    output:
        gff = os.path.join(OUTMAIN, "inputs", "cds.gff")
    params:
        db = config["db"]
    shell: """
    grep -v "^#" {input.gff} | awk '$3 == "CDS" && $2 == "{params.db}" ' > {output.gff}
"""



rule select_fasta_cds:
    input:
        ref = config["ref"],
        gff = os.path.join(OUTMAIN, "inputs", "cds.gff")
    output:
        fas = os.path.join(OUTMAIN, "inputs", "cds.fas")
    shell: """
    {BEDTOOLS} getfasta -fi {input.ref} -bed {input.gff} -s -fo {output.fas}
"""



rule degenerate_bed:
    input:
        fas = os.path.join(OUTMAIN, "inputs", "cds.fas"),
        gff = os.path.join(OUTMAIN, "inputs", "cds.gff")
    output:
        bed = os.path.join(OUTMAIN, "beds", "degeneracy.bed")
    shell:"""
    {PYTHON} {DEGENERATE} -g {input.gff} -f {input.fas} -o {output.bed}
"""


rule rename_chrs:
    input:
        bed = os.path.join(OUTMAIN, "beds", "degeneracy.bed")
    output:
        bed = os.path.join(OUTMAIN, "beds", "degeneracy_chrrenamed.bed")
    params:
        renamer = ";".join(expand("s/{r[0]}/{r[1]}/g", r = rename_chrs))
    shell: """
    sed '{params.renamer}' {input.bed} > {output.bed}
"""



rule split_by_degeneracy:
    input:
        bed = os.path.join(OUTMAIN, "beds", "degeneracy_chrrenamed.bed")
    output:
        bed4f = os.path.join(OUTMAIN, "beds", "sites_4f_degenerate.bed"),
        bed2f = os.path.join(OUTMAIN, "beds", "sites_2f_degenerate.bed"),
        bed0f = os.path.join(OUTMAIN, "beds", "sites_0f_degenerate.bed")
    shell:"""
    awk '$4 == 0 {{print $1"\t"$2"\t"$3}}' {input.bed} | {BEDTOOLS} sort | uniq > {output.bed0f}
    awk '$4 == 2 {{print $1"\t"$2"\t"$3}}' {input.bed} | {BEDTOOLS} sort | uniq > {output.bed2f}
    awk '$4 == 4 {{print $1"\t"$2"\t"$3}}' {input.bed} | {BEDTOOLS} sort | uniq > {output.bed4f}
"""



rule intersect_good_regions:
    input:
        bed = os.path.join(OUTMAIN, "beds", "sites_{d}f_degenerate.bed")
    output:
        bed = os.path.join(OUTMAIN, "beds", "good_sites_{d}f_degenerate.bed")
    params:
        mask = config["mask"]
    shell: """
    {BEDTOOLS} intersect -a {params.mask} -b {input.bed} | {BEDTOOLS} sort > {output.bed}
"""



rule bed_to_angsd_sites:
    input:
        bed = os.path.join(OUTMAIN, "beds", "good_sites_{d}f_degenerate.bed")
    output:
        sites = os.path.join(OUTMAIN, "beds", "good_sites_{d}f_degenerate.sites")
    shell:"""
    cut -f1,3 {input.bed} > {output.sites}
"""


rule index_angsd_sites:
    input:
        sites = os.path.join(OUTMAIN, "beds", "good_sites_{d}f_degenerate.sites")
    output:
        multiext(os.path.join(OUTMAIN, "beds", "good_sites_{d}f_degenerate.sites"), ".idx", ".bin")
    shell: """
    {ANGSD} sites index {input.sites}
"""


rule do_saf_degenerate_individual:
    input:
        bam = lambda wildcards: INDS[wildcards.p][wildcards.ind],
        sites = os.path.join(OUTMAIN, "beds", "good_sites_{d}f_degenerate.sites"),
        sites_idx = os.path.join(OUTMAIN, "beds", "good_sites_{d}f_degenerate.sites.idx")
    output:
        saf = os.path.join(OUTMAIN, "safs_indivs", "{p}_{ind}_{d}f_degenerate_sites.saf.gz"),
        saf_idx = os.path.join(OUTMAIN, "safs_indivs", "{p}_{ind}_{d}f_degenerate_sites.saf.idx")
    params:
        rf = config["chromlist"],
        outprefix = lambda wildcards, output: output.saf_idx.replace(".saf.idx", ""),
        minQ = 30,
        minMapQ = 30,
        ref = config["ref_renamed"]
    log: os.path.join(OUTMAIN, "safs", "{p}_{ind}_{d}f_degenerate_sites.arg")
    shell: """
    {ANGSD} -i {input.bam} -gl 2 -doSaf 1 -anc {params.ref} -out {params.outprefix} -rf {params.rf} -sites {input.sites} -minQ {params.minQ} -minMapQ {params.minMapQ}
"""


rule do_sfs_degenerate_individual:
    input:
        saf = os.path.join(OUTMAIN, "safs_indivs", "{p}_{ind}_{d}f_degenerate_sites.saf.idx")
    output:
        sfs  = os.path.join(OUTMAIN, "sfs_indivs", "{p}_{ind}_{d}f_degenerate.sfs")
    log: os.path.join(OUTMAIN, "sfs_indivs", "{p}_{ind}_{d}f_degenerate.log")
    threads: 15
    shell:"""
    {REALSFS} {input.saf} > {output.sfs} 2> {log}
"""

    

rule collect_pop_hets:
    input:
        sfs = lambda wildcards: expand(os.path.join(OUTMAIN, "sfs_indivs", "{{p}}_{ind}_{d}f_degenerate.sfs"), ind = INDS[wildcards.p].keys(), d = [0,2,4])
    output:
        hets = os.path.join(OUTMAIN, "sfs_indivs", "collected_heterozygosities_{p}.txt"),
    run:
        import os
        import pandas as pd
        data = []
        names = []
        for x in input:
            #name = os.path.basename(x).replace("degenerate.sfs", "")
            name = os.path.basename(x).split("_")[1]
            names.append(name)
            pop = os.path.basename(x).split("_")[0]
            degeneracy = os.path.basename(x).split("_")[2]
            with open(x, 'r') as fh:
                t = fh.read()
                data.append([pop, degeneracy] + [float(x) for x in t.rstrip().split()])
                
        a = pd.DataFrame(data, index=names, columns=["population", "degeneracy", "aa","ad","dd"])
        a["het"] = a["ad"] / a[["aa","ad","dd"]].sum(1)
        a.to_csv(output.hets, index=True, header=True, index_label="id", sep=" ")



rule collect_hets:
    input:
        expand(os.path.join(OUTMAIN, "sfs_indivs", "collected_heterozygosities_{p}.txt"), p=config["pops"].keys())
    output:
        hets = os.path.join(OUTMAIN, "sfs_indivs", "collected_heterozygosities_all.txt")
    run:
        dfs = []
        for f in input:
            dfs.append(pd.read_table(f))
        df = pd.concat(dfs)
        df.to_csv(output.hets, index=False, header=True)
        


rule do_safs_degeneate:
    input:
        bamlist = lambda wildcards: config["pops"][wildcards.p],
        sites = os.path.join(OUTMAIN, "beds", "good_sites_{d}f_degenerate.sites"),
        sites_idx = os.path.join(OUTMAIN, "beds", "good_sites_{d}f_degenerate.sites.idx")
    output:
        saf = os.path.join(OUTMAIN, "safs", "{p}_{d}f_degenerate_sites.saf.gz"),
        saf_idx = os.path.join(OUTMAIN, "safs", "{p}_{d}f_degenerate_sites.saf.idx")
    params:
        rf = config["chromlist"],
        outprefix = os.path.join(OUTMAIN, "safs", "{p}_{d}f_degenerate_sites"),
        minQ = 30,
        minMapQ = 30,
        ref = config["ref_renamed"]
    log: os.path.join(OUTMAIN, "safs", "{p}_{d}f_degenerate_sites.arg")
    shell: """
    {ANGSD} -b {input.bamlist} -gl 2 -doSaf 1 -anc {params.ref} -out {params.outprefix} -rf {params.rf} -sites {input.sites} -minQ {params.minQ} -minMapQ {params.minMapQ}
"""



rule do_sfs_degenerate:
    input:
        saf = os.path.join(OUTMAIN, "safs", "{p}_{d}f_degenerate_sites.saf.idx")
    output:
        sfs = os.path.join(OUTMAIN, "sfs", "{p}_{d}f_degenerate.sfs")
    log: os.path.join(OUTMAIN, "sfs", "{p}_{d}f_degenerate.log")
    threads: 15
    shell:"""
    {REALSFS} {input.saf} > {output.sfs} 2> {log}
"""




rule plot_sfs_degenerate:
    input:
        expand(os.path.join(OUTMAIN, "sfs", "{{p}}_{d}f_degenerate.sfs"), d=[0,2,4])
    output:
        png1 = os.path.join(OUTMAIN, "plot", "{p}_degenerate_sfs.png"),
        png2 = os.path.join(OUTMAIN, "plot", "{p}_degenerate_sfs_nofixedderived.png"),
    params:
        indir = os.path.join(OUTMAIN, "sfs"),
        p = "{p}"
    shell:"""
    {R} {PLOTSFS} {params.indir} {params.p} {output.png1}
"""
