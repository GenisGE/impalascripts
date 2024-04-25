import gzip
# adpated from frederik ffs script in /home/leopard/users/ffs/depth/scripts/06_find_ranges.py

def is_within_thresholds(x, min_thres, max_thres):
    
    if x == 'NA\n':
        return False

    x = int(x)
    if x < min_thres or x > max_thres:
        return False
    
    return True


def get_bad_ranges(infile, outfile, min_thres, max_thres):

    out = open(outfile, 'w+')
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
                    start = int(pos) - 1 # -1 bc bed format is 0 indexed
                    is_out = True
            elif is_out:
                if is_within_thresholds(d, min_thres, max_thres):
                    end = int(pos) - 1 # -1 bc this is good but previous was bad, bed format does not include end position in range
                    out.write(f'{chrr}\t{start}\t{end}\n')
                    is_out = False

                elif not is_within_thresholds(d, min_thres, max_thres):
                    continue

        if is_out:
             end = int(pos)
             out.write(f'{chrr}\t{start}\t{end}\n')
             
    out.close()


    

lowDepthDir = "/home/genis/impala/localSitesQC/depth/output/lowDepth/"
highDepthDir = "/home/genis/impala/localSitesQC/depth/output/highDepth/"

outDir = '/home/genis/impala/localSitesQC/depth/bedFiles/'

outfile_low = outDir + 'remove_depth_lowDepth.bed'
outfile_high = outDir + 'remove_depth_highDepth.bed'

chr_file = "/home/genis/impala/info_files/goatChrsAutoandX.txt"


# thresholds are lower 0.5 * median and higher 1.5 * median, separately for low depth and high depth samples
# see graphical distribution in two bottom panels (top panel is wrong) of /home/genis/impala/localSitesQC/depth/plots/depthDist_allChromosomes_median0.5and1.5.png

# low depth median 256
minl = 128
maxl = 384

# high depth median 136
minh = 68
maxh = 204


chrr = 'CHR_29'
infile_low = lowDepthDir + chrr + '_lowDepth.pos.gz'
infile_high = highDepthDir + chrr + '_highDepth.pos.gz'

get_bad_ranges(infile_low, outfile_low, minl, maxl)
get_bad_ranges(infile_high, outfile_high, minh, maxh)
