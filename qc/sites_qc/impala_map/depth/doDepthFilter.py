import gzip
# adpated from frederik ffs script in /home/leopard/users/ffs/depth/scripts/06_find_ranges.py

def is_within_thresholds(x, min_thres, max_thres):
    
    if x == 'NA\n':
        return False

    x = int(x)
    if x < min_thres or x > max_thres:
        return False
    
    return True


def get_bad_ranges(infile, out, min_thres, max_thres):

    #out = open(outfile, 'w+')
    is_out = False
    with gzip.open(infile, 'rt') as f:
        first_line = True
        for line in f:
            if first_line:
                first_line = False
                continue
            chrr, pos, d = line.split('\t')
            if not is_out:
                if is_within_thresholds(d, min_thres, max_thres):
                    continue
                elif not is_within_thresholds(d, min_thres, max_thres):
                    start = int(pos) - 1 # -1 bc bed format is 0 indexed. APPARENTLY THIS IS WRONG
                    is_out = True
            elif is_out:
                if is_within_thresholds(d, min_thres, max_thres):
                    end = int(pos) - 1 # -1 bc this is good but previous was bad, bed format does not include end position in range. APPARENTLY THIS IS WRONG
                    out.write(f'{chrr}\t{start}\t{end}\n')
                    is_out = False

                elif not is_within_thresholds(d, min_thres, max_thres):
                    continue

        if is_out:
             end = int(pos)
             out.write(f'{chrr}\t{start}\t{end}\n')
             
    #out.close()


    
# depths by pos across scaffolds done by casia in /home/casia16/projects/impala/results/depth
lowDepthDir = "/home/casia16/projects/impala/results/depth/lowDepth/"
highDepthDir = "/home/casia16/projects/impala/results/depth/highDepth/"

outDir = "/home/genis/impala/localSitesQC/impalaMap/depth/bedFiles/"

outfile_low = outDir + 'remove_depth_lowDepth.bed'
outfile_high = outDir + 'remove_depth_highDepth.bed'

chr_file = "/home/genis/impala/localSitesQC/impalaMap/allRegionsFiles/impalaRef100kbScaffolds.txt"


# thresholds are lower 0.5 * median and higher 1.5 * median, separately for low depth and high depth samples
# see graphical distribution in two bottom panels /home/casia16/projects/impala/results/depth/plots/depthDist_allChromosomes.png
# low depth thresholds
#print(c(0.5*mld,mld,1.5*mld))
#133.5 267.0 400.5
 
minl = 133.5
maxl = 400.5

# high depth thresholds based on casia
#print(c(0.5*mhd,mhd,1.5*mhd))
#71.5 143.0 214.5

minh = 71.5
maxh = 214.5

outlow = open(outfile_low, 'w+')
outhigh = open(outfile_high, 'w+')

with open(chr_file, 'rt') as f:
    for line in f:
        
        chrr = line.strip('\n')
        infile_low = lowDepthDir + chrr + '_lowDepth.pos.gz'
        infile_high = highDepthDir + chrr + '_highDepth.pos.gz'
        
        get_bad_ranges(infile_low, outlow, minl, maxl)
        get_bad_ranges(infile_high, outhigh, minh, maxh)


outlow.close()
outhigh.close()
        
    
