# snakemake to estimate 2dsfs, sfs and fst between popualtions pairs.
# version 2 intention is to do a more efficient version which will use downsampled saf



import pandas as pd
import itertools as it

ANGSD="/home/genis/software/angsd/angsd"
REALSFS="/home/genis/software/angsd/misc/realSFS"


OUTMAIN="results2"


# ancestral used is impala reference genome. not really ancestral but it does not matter since we cannot polarize anyways
#configfile:"config_sfs_fst.yaml"
configfile:"config_sfs_fst_v2.yaml"

infofile=config['info']

info = pd.read_table(infofile)
pops = list(info.my_locality2.unique()) # here the pop column in info file is called my_locality2. it is very likely in future uses we want to change that
#pops = ["Chobe",  "Ugalla"]


rule all:
    input:
        expand(os.path.join(OUTMAIN, "sfs", "{pop}.sfs"), pop=pops),
        os.path.join(OUTMAIN, "fst", "all_fst.txt")
#    shell: "ln -s -f {OUTMAIN} ."


        
rule do_bam_lists:
    input:
        config['info']
    output:
        expand(os.path.join(OUTMAIN, "bamlists", "bams_{pop}.list"), pop = pops),
        expand(os.path.join(OUTMAIN, "bamlists", "samples_{pop}.list"), pop = pops)
    run:
        for pop in pops:
            bams = info.loc[info.my_locality2==pop].bams_impala_map
            outfile1 = os.path.join(OUTMAIN, "bamlists", "bams_{}.list".format(pop))
            ids = info.loc[info.my_locality2==pop].Serial_number
            outfile2 = os.path.join(OUTMAIN, "bamlists", "samples_{}.list".format(pop))
            bams.to_csv(outfile1, header=False, index = False)
            ids.to_csv(outfile2, header=False, index = False)



rule do_saf_full:
    input:
        bamlist = ancient(os.path.join(OUTMAIN, "bamlists", "bams_{pop}.list"))
    output:
        saf = os.path.join(OUTMAIN, "safs", "{pop}.saf.gz"),
        saf_idx = os.path.join(OUTMAIN, "safs", "{pop}.saf.idx"),
        saf_pos = os.path.join(OUTMAIN, "safs", "{pop}.saf.pos.gz")
    params:
        anc = config['anc'],
        outprefix = os.path.join(OUTMAIN, "safs", "{pop}"),
        sites = config['sites'],
        rf = config['chroms'],
        minQ = 30,
        minMapQ = 30
    log: os.path.join(OUTMAIN, "safs", "{pop}.arg")
    threads: 5
    shell: "{ANGSD} -b {input.bamlist} -gl 2 -dosaf 1 -anc {params.anc} -out {params.outprefix} -P {threads} -rf {params.rf} -sites {params.sites} -minQ {params.minQ} -minMapQ {params.minMapQ} 2> /dev/null"

           
rule do_saf_small:
    input:
        bamlist = os.path.join(OUTMAIN, "bamlists", "bams_{pop}.list")
    output:
        saf = os.path.join(OUTMAIN, "small_safs", "{pop}.saf.gz"),
        saf_idx = os.path.join(OUTMAIN, "small_safs", "{pop}.saf.idx"),
        saf_pos = os.path.join(OUTMAIN, "small_safs", "{pop}.saf.pos.gz")
    log: os.path.join(OUTMAIN, "small_safs", "{pop}.arg")         
    params:
        sites = config['small_sites'],
        rf = config['small_chroms'], #small chroms are actually big chroms (>1mb)
        anc = config['anc'],
        outprefix = os.path.join(OUTMAIN, "small_safs", "{pop}"),
        minQ = 30,
        minMapQ = 30
    threads: 5
    shell: "{ANGSD} -b {input.bamlist} -gl 2 -dosaf 1 -anc {params.anc} -out {params.outprefix} -P {threads} -rf {params.rf} -sites {params.sites} -minQ {params.minQ} -minMapQ {params.minMapQ} 2> /dev/null"
        
    
rule do_sfs:
    input:
        idx = os.path.join(OUTMAIN, "small_safs", "{pop}.saf.idx")
    output:
        sfs = os.path.join(OUTMAIN, "sfs", "{pop}.sfs")
    threads: 10
    log: os.path.join(OUTMAIN, "sfs", "{pop}.log")
    shell: "{REALSFS} {input.idx} -P {threads} > {output.sfs} 2> {log}"



rule do_2dsfs:
    input:
        idx1 = os.path.join(OUTMAIN, "small_safs", "{pop1}.saf.idx"),
        idx2 = os.path.join(OUTMAIN, "small_safs", "{pop2}.saf.idx")
    output:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{pop1}_{pop2}.sfs")
    log: os.path.join(OUTMAIN, "2dsfs", "{pop1}_{pop2}.log")
    threads: 10
    shell: "{REALSFS} {input.idx1} {input.idx2} -P {threads} > {output.sfs} 2> {log}"


           
rule do_fst_idx:
    input:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{pop1}_{pop2}.sfs"),
        idx1 = os.path.join(OUTMAIN, "safs", "{pop1}.saf.idx"),
        idx2 = os.path.join(OUTMAIN, "safs", "{pop2}.saf.idx")
    output:
        fst_idx = temp(os.path.join(OUTMAIN, "fst", "{pop1}_{pop2}.fst.idx")),
        fst_gz = temp(os.path.join(OUTMAIN, "fst", "{pop1}_{pop2}.fst.gz"))
    params:
        outprefix = os.path.join(OUTMAIN, "fst","{pop1}_{pop2}"),
        whichfst = 1
    shell: "{REALSFS} fst index {input.idx1} {input.idx2} -sfs {input.sfs} -fstout {params.outprefix} -whichFst {params.whichfst} 2> /dev/null"



rule do_fst:
    input:
        fst_idx = os.path.join(OUTMAIN, "fst", "{pop1}_{pop2}.fst.idx"),
        fst_gz = os.path.join(OUTMAIN, "fst", "{pop1}_{pop2}.fst.gz")
    output:
        fst = os.path.join(OUTMAIN, "fst", "{pop1}_{pop2}.fst")
    shell: "{REALSFS} fst stats {input.fst_idx} > {output.fst}"

           

rule collect_fst:
    input:
        expand(os.path.join(OUTMAIN, "fst", "{p[0]}_{p[1]}.fst"), p = it.combinations(pops, 2))
    output:
        f = os.path.join(OUTMAIN, "fst", "all_fst.txt")
    run:
        data = []
        for x in input:
            pops = os.path.basename(x).replace(".sfs", "").split("_")
#            popss.append(pops)
            with open(x, 'r') as fh:
                ff = fh.read()
                fst = [float(x) for x in ff.rstrip().split()]
            data.append(pops+fst)
        out = pd.DataFrame(data, columns = ['pop1', 'pop2', 'fst_unweight', 'fst_weight'])
        out.to_csv(output.f, index=False, header=True, sep = "\t")
