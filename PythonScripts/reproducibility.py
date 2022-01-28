'''
This file implements
HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient.
'''
import numpy as np
from hicrep.utils import readMcool
from hicrep import hicrepSCC
import itertools

from PythonScripts.utlis import plot_rep

samples = ['OF', 'OM', 'OMF', 'YF', 'YM', 'YMF']
pair_order_list = list(itertools.combinations(samples, 2))
resolution_bin = [(1000000, 500000), (1000000, 1000000)]
chros = ['chr{}'.format(i) for i in list(range(1, 20))] + ['chrX', 'chrY']
h = 1
bDownSample = False
result = {}
datatype = 'norm'

for res in resolution_bin:
    binSize = res[1]
    dBPMax = res[0]
    for chrom in chros:
        for i in pair_order_list:
            if datatype == 'raw':
                fmcool1 = "../MuSC_HiC_files/HiC_CPB_normalized/{}_{}.cool".format(i[0], binSize)
                fmcool2 = "../MuSC_HiC_files/HiC_CPB_normalized/{}_{}.cool".format(i[1], binSize)
            elif datatype == 'norm':
                fmcool1 = "../MuSC_HiC_files/HiC_CPB_normalized/normalized_{}_{}.cool".format(i[0].lower(), binSize)
                fmcool2 = "../MuSC_HiC_files/HiC_CPB_normalized/normalized_{}_{}.cool".format(i[1].lower(), binSize)
            elif datatype == 'not norm':
                fmcool1 = "../MuSC_HiC_files/HiC_CPB_normalized/not_normalized_{}_{}.cool".format(i[0].lower(), binSize)
                fmcool2 = "../MuSC_HiC_files/HiC_CPB_normalized/not_normalized_{}_{}.cool".format(i[1].lower(), binSize)


            cool1, binSize1 = readMcool(fmcool1, -1)
            cool2, binSize2 = readMcool(fmcool2, -1)

            result['_'.join(i)] = sccSub = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample, np.array([chrom], dtype=str))[0]
        plot_rep(result, chrom, res, False)



