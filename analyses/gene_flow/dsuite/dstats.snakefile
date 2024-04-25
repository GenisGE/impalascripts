
import pandas as pd

BCFTOOLS="/home/genis/software/bcftools/bcftools"
VCFTOOLS="/home/genis/software/vcftools/src/cpp/vcftools"
VCF2EIGENSTRAT="python /home/genis/software/gdc/vcf2eigenstrat.py"
PLINK="/home/genis/software/plink"
REFFINDER="/home/genis/software/refFinder/refFinder"
DSUITE="/home/genis/software/Dsuite/Build/Dsuite"
DTOOLS="/home/genis/software/Dsuite/utils/dtools.py"
BGZIP="bgzip"


OUTBIG=config["outbig"]
OUTMAIN=config["outmain"]
REF=config["ref"]
CHROMLIST=config["chroms"]


with open(CHROMLIST, "r") as fh:
    CHROMS=[x.rstrip() for x in fh.readlines()]


wildcard_constraints:
    sample = "|".join(config["samples"].keys()),
    chrom="|".join(CHROMS)


rule all:
    input:
        os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt_annotated_variants_conservation.vcf.gz"),
        os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt_annotated_variants_conservation.vcf.gz.csi"),
        multiext(os.path.join(OUTMAIN, "dsuite", "fbranch"), ".png", ".svg"),


rule clean_vcf:
    """
    removes sites where all samples are homozyogus derived from bcf
"""
    input:
        vcf = config["vcf"]
    output:
        vcf = os.path.join(OUTBIG, "vcf", "merged_snps_filtered.bcf.gz")
    shell:"""
    {BCFTOOLS} view -i 'GT=="0/0" || GT=="0/1"' -Ob -o {output.vcf} {input.vcf}
"""


## convert to plink for admixtools2 and to add ref as sample
rule bcf_to_plink:
    input:
        bcf = os.path.join(OUTBIG, "vcf", "merged_snps_filtered.bcf.gz")
    output:
        plink = multiext(os.path.join(OUTBIG, "plink", "merged_snps_filtered"), ".bed", ".bim", ".fam")
    params:
        prefix = os.path.join(OUTBIG, "plink", "merged_snps_filtered"),
        nchr = len(config["chroms"])
    shell:
        """
        {BCFTOOLS} annotate -Ov -x ID -I +'%CHROM:%POS:%REF:%ALT' {input.bcf} | {PLINK} --vcf /dev/stdin --make-bed --allow-extra-chr --chr-set {params.nchr} --out {params.prefix}
"""


rule bed_to_tped:
    input:
        multiext(os.path.join(OUTBIG, "plink", "merged_snps_filtered"), ".bed", ".bim", ".fam")
    output:
        temp(multiext(os.path.join(OUTBIG, "plink", "merged_snps_filtered"), ".tped", ".tfam"))
    params:
        prefix = os.path.join(OUTBIG, "plink", "merged_snps_filtered")
    shell: """
    {PLINK} --bfile {params.prefix} --recode transpose --out {params.prefix} --allow-extra-chr
"""


rule get_ref_tped:
    input:
        bim = os.path.join(OUTBIG, "plink", "merged_snps_filtered.bim"),
        ref = REF
    output:
        tped = temp(os.path.join(OUTBIG, "plink", "ref.tped")),
        tfam = temp(os.path.join(OUTBIG, "plink", "ref.tfam"))
    params:
        os.path.join(OUTBIG, "plink", "ref")
    shell: """
    cat {input.bim} | {REFFINDER} {input.ref} plink > {output.tped}
    echo "ref Goat 0 0 0 -9" > {output.tfam}
"""


rule add_ref:
    input:
        my_tped = os.path.join(OUTBIG, "plink", "merged_snps_filtered.tped"),
        my_tfam = os.path.join(OUTBIG, "plink", "merged_snps_filtered.tfam"),
        ref_tped = os.path.join(OUTBIG, "plink", "ref.tped"),
        ref_tfam = os.path.join(OUTBIG, "plink", "ref.tfam"),
    output:
        tped = temp(os.path.join(OUTBIG, "plink", "merged_snps_filtered_withref.tped")),
        tfam = temp(os.path.join(OUTBIG, "plink", "merged_snps_filtered_withref.tfam")),
    shell: """
    paste -d" " {input.my_tped} {input.ref_tped} > {output.tped}
    cat {input.my_tfam} {input.ref_tfam} > {output.tfam}
"""


rule tped_to_bed_withref:
    input:
        multiext(os.path.join(OUTBIG, "plink", "merged_snps_filtered_withref"),".tped", ".tfam")  
    output:
        multiext(os.path.join(OUTBIG, "plink", "merged_snps_filtered_withref"),".bed", ".fam", ".bim")
    params:
        prefix = os.path.join(OUTBIG, "plink", "merged_snps_filtered_withref")
    shell:"""
    {PLINK} -tfile {params.prefix} --make-bed --out {params.prefix} --allow-extra-chr
"""



rule bed_to_vcf_withref:
    input:
        multiext(os.path.join(OUTBIG, "plink", "merged_snps_filtered_withref"),".bed", ".fam", ".bim"),
    output:
        temp(os.path.join(OUTBIG, "vcf", "merged_snps_filtered_withref.vcf"))
    params:
        inprefix = os.path.join(OUTBIG, "plink", "merged_snps_filtered_withref"),
        outprefix = os.path.join(OUTBIG, "vcf", "merged_snps_filtered_withref")
    shell:"""
    {PLINK} -bfile {params.inprefix} --recode vcf --out {params.outprefix} --allow-extra-chr
"""


rule fix_ref_alt_vcf:
    input:
        vcf = os.path.join(OUTBIG, "vcf", "merged_snps_filtered_withref.vcf")
    output:
        vcf = os.path.join(OUTBIG, "vcf", "merged_snps_filtered_withref_fixrefalt.vcf.gz")
    params:
        ref=config["ref"]
    shell:"""
    {BCFTOOLS} norm -c s -f {params.ref} -N -Oz -o {output.vcf} {input.vcf}
"""




### DSUITE RULES
rule do_pop_file_to_dsuite:
    input:
        fam = os.path.join(OUTBIG, "plink", "merged_snps_filtered_withref.fam")
    output:
        popfile = os.path.join(OUTMAIN, "dsuite", "pops_withref.list")
    params:
        outgroup = "Goat"
    shell: """
    awk '{{print $1"_"$2"\t"$2}}' {input.fam} | sed 's/{params.outgroup}$/Outgroup/g' > {output.popfile}
"""


rule do_dsuite_dstats:
    input:
        vcf = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt_annotated_variants_conservation.vcf.gz"),
        popfile = os.path.join(OUTMAIN, "dsuite", "pops_withref.list")
    output:
        multiext(os.path.join(OUTMAIN, "dsuite", "dstatsall"), "_BBAA.txt", "_Dmin.txt")
    params:
        outprefix = os.path.join(OUTMAIN, "dsuite", "dstatsall"),
        n_blocks = 450 # for mammalian normal sized genome this makes approx blocks of 5 mb
    shell: """
    {DSUITE} Dtrios -k {params.n_blocks} -o {params.outprefix} {input.vcf} {input.popfile}
"""

    
rule do_dsuite_dstats_withtree:
    input:
        vcf = os.path.join(OUTMAIN, "vcf", "merged_snps_filtered_withref_fixrefalt_annotated_variants_conservation.vcf.gz"),
        popfile = os.path.join(OUTMAIN, "dsuite", "pops_withref.list"),
        tree = config["tree"]
    output:
        multiext(os.path.join(OUTMAIN, "dsuite", "dstatsall_withtree"), "_BBAA.txt", "_Dmin.txt", "_tree.txt"),
    params:
        outprefix = os.path.join(OUTMAIN, "dsuite", "dstatsall_withtree"),
        n_blocks = 450 # for mammalian normal sized genome this makes approx blocks of 5 mb
    shell: """
    {DSUITE} Dtrios -k {params.n_blocks} -t {input.tree} -o {params.outprefix} {input.vcf} {input.popfile}
"""
        

rule do_fbranch:
    input:
        tree=config["tree"],
        fvals = os.path.join(OUTMAIN, "dsuite", "dstatsall_withtree_tree.txt")
    output:
        fbranch = os.path.join(OUTMAIN, "dsuite", "fbranch.txt")
    shell:"""
    {DSUITE} Fbranch -p 0.001 {input.tree} {input.fvals} > {output.fbranch}
"""


rule plot_fbranch:
    input:
        fbranch=os.path.join(OUTMAIN, "dsuite", "fbranch.txt"),
        tree = config["tree"]
    output:
        multiext(os.path.join(OUTMAIN, "dsuite", "fbranch"), ".png", ".svg"),
    params:
        outprefix = os.path.join(OUTMAIN, "dsuite", "fbranch")
    shell:"""
    {DTOOLS} -n {params.outprefix} {input.fbranch} {input.tree}
"""
