'''
This file implements
HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient.
'''
import pickle

import numpy as np
import pandas as pd
from hicrep.utils import readMcool
from hicrep import hicrepSCC
import itertools

# url = 'https://raw.githubusercontent.com/holtzy/The-Python-Graph-Gallery/master/static/data/mtcars.csv'
# df = pd.read_csv(url)
# df = df.set_index('model')
# print(df)
# exit()
from PythonScripts.utlis import plot_rep, line_rep

samples = ['OF', 'OM', 'YF', 'YM']
pair_order_list = list(itertools.combinations_with_replacement(samples, 2))

resolution_bin = [(1000000, 500000), (1000000, 1000000)]
chros = ['chr{}'.format(i) for i in list(range(1, 20))] #+ ['chrX']#, 'chrY']
# chros = ['chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrM', 'chrX', 'chrY']


h = 1
bDownSample = False
result = {}
datatype = 'norm'


def pickle_save(data, filename):
    with open('{}.pickle'.format(filename), 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)


def pickle_load(filename):
    with open('{}.pickle'.format(filename), 'rb') as handle:
        return pickle.load(handle)


def compute_rep():
    for res in resolution_bin:
        binSize = res[1]
        dBPMax = res[0]
        for i in pair_order_list:
            pair_name = '_'.join(i)
            print("SCC between {} at resolution {}".format(pair_name, res))
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

            sc = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample, np.array(chros, dtype=str))

            # result['_'.join(i)] = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample, np.array(chros, dtype=str))[0]
            result[pair_name] = sc

        line_rep(result, chros, res[1], True)
compute_rep()
#
# data = pickle_load("data")
# line_rep(data, chros, 500000, True)