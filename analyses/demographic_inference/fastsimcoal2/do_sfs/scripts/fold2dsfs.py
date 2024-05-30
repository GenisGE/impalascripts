import dadi
import sys
import numpy as np

insfs=sys.argv[1]
outsfs=sys.argv[2]
#n1=sys.argv[3]
#n2=sys.argv[4]



fs = dadi.Spectrum.from_file(insfs)
fs = fs.fold()

with open(outsfs, "w+") as f:
    for line in fs.data:
        np.savetxt(f, line, fmt="%.2f", delimiter=" ",  newline=" ")
        f.write("\n")

