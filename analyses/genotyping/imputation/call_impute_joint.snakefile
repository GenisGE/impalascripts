# snakemake to do joint genotype calling with bcftools plus beagle imputation.
# config needs:
#     ref: path to reference fasta file
#     outmain: path to main results folder
#     bamlist: path to file with list of bam files
#     chroms: path to file with chromosomes to call genotypes on
#     samples: path to file with sample names, in same order as in bamlist
#     grouplist: path to file with population assignment for each sample (tab delimited two columns, first column sample name and second a
#ssigned popualtion). joint calling is done within populations
#     filters: ditionary filter for genotype calling
#                 bed: bedf ile with regions to include
#                 minq: minimum base quality to include base
#                 min_mapq: minimum mapping quality to include read
#     imputed_thres: thershold to keep imputed sites/genotype
#                 minr2: minimum r2 value to include site (sites below are exclude)
#                 mingp: minimum genotype probability to call a genotype (genotypes below are set to missing)


BCFTOOLS="/home/genis/software/bcftools/bcftools"
JAVA="java"
R="Rscript"
BEAGLE3="/home/genis/software/beagle.jar"
ANGSD="/home/genis/software/angsd/angsd"
PLINK="/home/genis/software/plink"
PLOTCALLRATE="scripts/plotR2CallRate.R"
PLOTGENOCALLRATE="scripts/plotGenoPCallRate.R"
PLOTR2VSMAF="scripts/plotR2vsMAF.R"
PLOTPCA="scripts/plotPlinkPCA.R"



REF=config["ref"]
OUTMAIN=config["outmain"]
OUTBIG=config["outbig"]
BAMLIST=config["bamlist"]
CHROMFILE=config["chroms"]
SAMPLES=config["samplelist"]
GROUPS=config["grouplist"]

with open(CHROMFILE, "r") as fh:
    CHROMS=[x.rstrip() for x in fh.readlines()]

BED=config["filters"]["bed"]
MIN_BQ=config["filters"]["minq"]
MIN_MQ=config["filters"]["min_mapq"]

MINR2=config["imputed_thres"]["minr2"]
MINGP=config["imputed_thres"]["minr2"]

wildcard_constraints:
    chrom="|".join(CHROMS),


    
rule all:
    input:
        os.path.join(OUTMAIN, "vcf_joint", "joint_call_filtered.bcf.gz"),
        os.path.join(OUTMAIN, "impute", "r2vsMAF.png"),
        multiext(os.path.join(OUTBIG, "call_imputed", "call_imputed_all"),".bed", ".bim", ".fam"),
        expand(os.path.join(OUTMAIN, "impute", "{s}callrate.png"), s=["r2", "gp"]),
        multiext(os.path.join(OUTMAIN, "impute", "call_imputed_all_check"), ".eigenvec", ".eigenval", ".lmiss.gz", ".imiss.gz"),
        os.path.join(OUTMAIN, "impute", "call_imputed_pca.png"),
        os.path.join(OUTBIG, "call_imputed", "call_imputed_all.vcf")
        

        
        
rule gen_bcftools_genome_wide_joint:
    input:
        BAMLIST
    output:
        temp(os.path.join(OUTBIG, "vcf_joint", "{chrom}.bcf.gz"))
    params:
        samples = SAMPLES, # file for bcftools -samples-file argumenthttp://samtools.github.io/bcftools/bcftools.html#common_options
        groups = GROUPS # file for bcftools call --group-samples http://samtools.github.io/bcftools/bcftools.html#call
    threads: 1
    shell: """
    {BCFTOOLS} mpileup -r {wildcards.chrom}  -B -Q {MIN_BQ}  -q {MIN_MQ} --threads {threads} -O u --fasta-ref {REF} -b {input} --per-sample-mF -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  | {BCFTOOLS} call -Ob -o {output} --threads {threads} --multiallelic-caller --variants-only --samples-file {params.samples} --group-samples {params.groups}
    #{BCFTOOLS} index {output}
"""


    
rule concat_chroms_joint:
    input:
        expand(os.path.join(OUTBIG, "vcf_joint", "{chrom}.bcf.gz"), chrom = CHROMS)
    output:
        vcf=os.path.join(OUTBIG, "vcf_joint", "joint_call.bcf.gz"),
        idx=os.path.join(OUTBIG, "vcf_joint", "joint_call.bcf.gz.csi")
    shell: """
    {BCFTOOLS} concat --naive -Ob -o {output.vcf} {input}
    {BCFTOOLS} index {output.vcf}
    """



rule filter_vcf:
    input:
        vcf = rules.concat_chroms_joint.output.vcf,
        idx = rules.concat_chroms_joint.output.idx,
    output:
        vcf = os.path.join(OUTBIG, "vcf_joint", "joint_call_filtered.bcf.gz"),
        idx = os.path.join(OUTBIG, "vcf_joint", "joint_call_filtered.bcf.gz.csi")
    params:
        bed = BED,
    threads: 3
    shell: """
    {BCFTOOLS} view --threads {threads} -T {params.bed} {input.vcf} | {BCFTOOLS} filter --threads {threads} -g 10:indel,other | {BCFTOOLS} view --threads {threads} --types snps -m2 -M2 | {BCFTOOLS} annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' -Ob -o {output.vcf}
    {BCFTOOLS} index {output.vcf}
"""



rule link_vcf:
    input:
        vcf = rules.filter_vcf.output.vcf,
        idx = rules.filter_vcf.output.idx
    output:
        vcf=os.path.join(OUTMAIN, "vcf_joint", "joint_call_filtered.bcf.gz"),
        idx=os.path.join(OUTMAIN, "vcf_joint", "joint_call_filtered.bcf.gz.csi")
    shell: """
    ln -s {input.vcf} {output.vcf}
    ln -s {input.idx} {output.idx}
"""



rule do_bgl_gl:
    input:
        vcf = rules.filter_vcf.output.vcf,
    output:
        bgl_gl = temp(os.path.join(OUTBIG, "beagle", "{chrom}.beagle.gz")),
    params:
        chrom = "{chrom}",
        minmaf = config["minmaf"],
        outprefix = os.path.join(OUTBIG, "beagle", "{chrom}"),
    log: os.path.join(OUTBIG, "beagle", "{chrom}.arg")
    shell: """
    {ANGSD} -vcf-pl {input.vcf} -doGlf 2 -r {params.chrom} -doMaf 1 -doMajorMinor 1 -minMaf {params.minmaf} -out {params.outprefix}
"""



rule impute:
    input:
        bgl_gl = rules.do_bgl_gl.output.bgl_gl,
    output:
        gprob = os.path.join(OUTBIG, "impute", "impute_{chrom}.{chrom}.beagle.gz.gprobs.gz"),
        r2 = os.path.join(OUTBIG, "impute", "impute_{chrom}.{chrom}.beagle.gz.r2"),
    log: os.path.join(OUTBIG, "impute", "impute_{chrom}.log")
    params:
        outprefix = os.path.join(OUTBIG, "impute", "impute_{chrom}")
    shell: """
    {JAVA} -jar {BEAGLE3} like={input.bgl_gl} out={params.outprefix}
"""



rule plot_r2_callrate:
    input:
        r2 = expand(os.path.join(OUTBIG, "impute", "impute_{chrom}.{chrom}.beagle.gz.r2"),chrom=CHROMS),
    output:
        os.path.join(OUTMAIN, "impute", "r2callrate.png")
    params:
        indir=os.path.join(OUTBIG, "impute")
    shell: """
    {R} {PLOTCALLRATE} {params.indir} {output}
"""


rule plot_r2_vs_maf:
    input:
        r2 = expand(os.path.join(OUTBIG, "impute", "impute_{chrom}.{chrom}.beagle.gz.r2"),chrom=CHROMS),
        dosages = expand(os.path.join(OUTBIG, "impute", "impute_{chrom}.{chrom}.beagle.gz.dose.gz"),chrom=CHROMS)
    output:
        os.path.join(OUTMAIN, "impute", "r2vsMAF.png")
    params:
        indir = os.path.join(OUTBIG, "impute")
    shell:"""
    {R} {PLOTR2VSMAF} {params.indir} {output}
"""


rule plot_gp_callrate:
    input:
        gp = expand(os.path.join(OUTBIG, "impute", "impute_{chrom}.{chrom}.beagle.gz.gprobs.gz"),chrom=CHROMS),
    output:
        os.path.join(OUTMAIN, "impute", "gpcallrate.png")
    params:
        indir=os.path.join(OUTBIG, "impute")
    shell: """
    {R} {PLOTGPCALLRATE} {params.indir} {output}
"""





rule select_r2_sites:
    input:
        r2 = os.path.join(OUTBIG, "impute", "impute_{chrom}.{chrom}.beagle.gz.r2"),
        r2_plot = os.path.join(OUTMAIN, "impute", "r2callrate.png")
    output:
        sites = temp(os.path.join(OUTBIG, "call_imputed", "{chrom}_use.sites")),
    params:
        minr2 = config["imputed_thres"]["minr2"],
        chrom = "{chrom}"
    shell: """
    awk '$2 > {params.minr2} {{print $0}}' {input.r2} | awk -F'[_ \t]' '{{print "chr"$2"\t"$3}}' | sed 's/CHR_/chr/g' > {output.sites}
"""



rule index_sites:
    input:
        sites =  os.path.join(OUTBIG, "call_imputed", "{chrom}_use.sites"),
    output:
        idx = temp(os.path.join(OUTBIG, "call_imputed", "{chrom}_use.sites.idx")),
        binn = temp(os.path.join(OUTBIG, "call_imputed", "{chrom}_use.sites.bin"))
    shell:"""
    {ANGSD} sites index {input.sites}
"""



rule rename_chr_fai:
    """ Stupid rule needed because chromosome names have underscores which messes with angsd, so need to change that
    in fai file and in .gprobs.gz file. This will probably create mess downstream but will woryy about that later"""
    input:
        fai= REF+".fai"
    output:
        fai = os.path.join(OUTMAIN, "fai", "ref_renamed_chrs.fai")
    shell: """
    sed 's/CHR_/chr/g' {input.fai} > {output.fai}
"""



rule rename_chr_gprobs:
    """ Stupid rule needed because chromosome names have underscores which messes with angsd, so need to change that
    in fai file and in .gprobs.gz file"""
    input:
        gprobs = os.path.join(OUTBIG, "impute", "impute_{chrom}.{chrom}.beagle.gz.gprobs.gz"),
    output:
        gprobs = temp(os.path.join(OUTBIG, "impute", "imputed_renamed_{chrom}.gprobs.gz"))
    shell: """
    zcat {input.gprobs} | sed 's/CHR_/chr/g' | gzip -c > {output.gprobs}
"""



rule call_imputed:
    input:
        gprobs = os.path.join(OUTBIG, "impute", "imputed_renamed_{chrom}.gprobs.gz"),
        fai = os.path.join(OUTMAIN, "fai", "ref_renamed_chrs.fai"),
        sites = os.path.join(OUTBIG, "call_imputed", "{chrom}_use.sites"),
        sites_idx = os.path.join(OUTBIG, "call_imputed", "{chrom}_use.sites.idx")
    output:
        tped = temp(os.path.join(OUTBIG, "call_imputed", "{chrom}.tped")),
        tfam = temp(os.path.join(OUTBIG, "call_imputed", "{chrom}.tfam"))
    params:
        outprefix = os.path.join(OUTBIG, "call_imputed", "{chrom}"),
        mingp = MINGP
    log: os.path.join(OUTBIG, "imputed_called", "{chrom}.args")
    shell: """
    {ANGSD} -beagle {input.gprobs} -fai {input.fai} -doPlink 2 -doGeno -1 -postCutoff {params.mingp} -sites {input.sites} -out {params.outprefix}
"""



rule concat_chr_call_imputed_tped:
    input:
        tpeds = expand(os.path.join(OUTBIG, "call_imputed", "{chrom}.tped"), chrom=CHROMS),
    output:
        tped = os.path.join(OUTBIG, "call_imputed", "call_imputed_all.tped"),
    shell: """
    cat {input.tpeds} > {output.tped}
"""



rule make_good_tfam:
    input:
        tfams = expand(os.path.join(OUTBIG, "call_imputed", "{chrom}.tfam"), chrom=CHROMS),
        popinf = GROUPS
    output:
        tfam = os.path.join(OUTBIG, "call_imputed", "call_imputed_all.tfam"),
    shell:"""
    awk '{{print $1" "$2" 0 0 0 -9"}}' {input.popinf} > {output.tfam}
"""



rule tped_to_bed:
    input:
        tped = os.path.join(OUTBIG, "call_imputed", "call_imputed_all.tped"),
        tfam = os.path.join(OUTBIG, "call_imputed", "call_imputed_all.tfam")
    output:
        bed = os.path.join(OUTBIG, "call_imputed", "call_imputed_all.bed"),
        bim = os.path.join(OUTBIG, "call_imputed", "call_imputed_all.bim"),
        fam = os.path.join(OUTBIG, "call_imputed", "call_imputed_all.fam")
    params:
        prefix = os.path.join(OUTBIG, "call_imputed", "call_imputed_all")
    log: os.path.join(OUTBIG, "call_imputed", "call_imputed_all.log")
    shell:"""
    {PLINK} --tfile {params.prefix} --make-bed --chr-set 29 --out {params.prefix}
"""


rule bed_to_vcf:
    input:
        multiext(os.path.join(OUTBIG, "call_imputed", "call_imputed_all"), ".bed", ".bim", ".fam")
    output:
        vcf = os.path.join(OUTBIG, "call_imputed", "call_imputed_all.vcf")
    params:
        prefix = os.path.join(OUTBIG, "call_imputed", "call_imputed_all")
    shell:"""
    {PLINK} --bfile {params.prefix} --chr-set 29 --recode vcf --out {params.prefix}
"""


rule do_pca_check:
    """rule to do a pca to make sure all is more or less ok"""
    input:
        multiext(os.path.join(OUTBIG, "call_imputed", "call_imputed_all"), ".bed", ".bim", ".fam")
    output:
        multiext(os.path.join(OUTMAIN, "impute", "call_imputed_all_check"), ".eigenvec", ".eigenval", ".lmiss.gz", ".imiss.gz")
    log: os.path.join(OUTMAIN, "impute", "call_imputed_all_check.log")
    params:
        inprefix = os.path.join(OUTBIG, "call_imputed", "call_imputed_all"),
        outprefix = os.path.join(OUTMAIN, "impute", "call_imputed_all_check")
    shell:"""
    {PLINK} --bfile {params.inprefix} --pca --missing gz --out {params.outprefix} --chr-set 29
"""



rule plot_pca_check:
    """rule to plot pca to make sure all i smore or less ok """
    input:
        multiext(os.path.join(OUTMAIN, "impute", "call_imputed_all_check"), ".eigenvec", ".eigenval")
    output:
        png = os.path.join(OUTMAIN, "impute", "call_imputed_pca.png")
    params:
        inprefix = os.path.join(OUTMAIN, "impute", "call_imputed_all_check")
    shell:"""
    {R} {PLOTPCA} {params.inprefix} {output.png}
"""
