import pdb
import sys
import os
import numpy as np
import cooler

cooler_file   = sys.argv[1]
contact_file  = sys.argv[2]
bins_file     = sys.argv[3]

res  = 1000000
cool = cooler.Cooler(str(cooler_file)+"::/resolutions/"+str(res))
bins_full = np.array(cool.bins()[:])
bin_idx   = list(range(0, bins_full.shape[0]))

def getcoord(pos1, pos2):
    chr1 = bins_full[pos1][0]
    bin1 = bins_full[pos1][1]
    chr2 = bins_full[pos2][0]
    bin2 = bins_full[pos2][1]
    return chr1, bin1, chr2, bin2

vfunc             = np.vectorize(getcoord)
idx1s             = cool.pixels()[:]['bin1_id']
idx2s             = cool.pixels()[:]['bin2_id']
chr1s, bin1s, chr2s, bin2s = vfunc(idx1s, idx2s)
vals              = cool.pixels()[:]['count']
contact_map_array = np.vstack([chr1s, bin1s, chr2s, bin2s, vals]).transpose()
bins_array        = np.vstack([bins_full[:,0],
                                bins_full[:,1],
                                bins_full[:,2],
                                bins_full[:,1]]).transpose()
np.savetxt(contact_file,
        contact_map_array,
        fmt="%s",
        delimiter='\t')

np.savetxt(bins_file,
        bins_array,
        fmt="%s",
        delimiter='\t')

os.system("gzip "+str(bins_file))
os.system("gzip "+str(contact_file))
