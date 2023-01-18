import os

import cooler
import pandas as pd
from matplotlib import pyplot as plt

samples = ['of', 'om', 'yf', 'ym', 'omf', 'ymf']

jar_3dmax = "/home/3DMaxInput/3DMax.jar"
params_file = "/home/3Dmax_params/parameters_{}.txt"


# Generate Interaction Frequencies
def eigen_pvalue(datatype, resolution, save_fig=False):
    # Using 1mb for AB compartment
    chros = ['chr{}'.format(i) for i in list(range(1, 20))]
    for chro in chros:
        print("Chromosome is {}".format(chro))
        for pos, i in enumerate(samples[4:]):
            if datatype == 'raw':
                clr = cooler.Cooler('../MuSC_HiC_files/raw/{}_{}.cool'.format(i.upper(), resolution))
            elif datatype == 'norm':
                clr = cooler.Cooler('../MuSC_HiC_files/HiC_CPB_normalized/normalized_{}_{}.cool'.format(i, resolution))
            elif datatype == 'not norm':
                clr = cooler.Cooler(
                    '../MuSC_HiC_files/HiC_not_normalized/not_normalized_{}_{}.cool'.format(i, resolution))

            # print(clr.pixels()[:].columns)
            # print(clr.bins()[:].columns)
            # print(set(clr.bins()[:]['chrom'])

            # c.matrix(balance=False, as_pixels=True, join=True)[1000:1005, 1000:1005]
            pi = clr.matrix(as_pixels=True, balance=True, join=True).fetch(('{}'.format(chro)))
            # klp = pi[['start1', 'end1', 'start2', 'end2', 'count', 'balanced']]


            start_regions = pi.set_index(['start1', 'end1'])
            end_regions = pi.set_index(['start2', 'end2'])
            regions = sorted(set(start_regions.index.to_list()).union(set(end_regions.index.to_list())))

            regions = {k: v + 1 for v, k in enumerate(regions)}


            pi['start_region'] = pi.set_index(['start1', 'end1']).index.map(regions)
            pi['end_region'] = pi.set_index(['start2', 'end2']).index.map(regions)

            pi.to_csv('../Results/IFs_all/{}.csv'.format(chro), index=False)

            pi[['start_region', 'end_region', 'count']] \
                .to_csv('../Results/IFs/{}.csv'
                        .format(chro), sep='\t', encoding='utf-8', index=False, header=False)

            pi[['start_region', 'end_region', 'balanced']].dropna() \
                .to_csv('../Results/IFs_balanced/{}.csv'
                        .format(chro), sep='\t', encoding='utf-8', index=False, header=False)

# run with windows intepreter
def generate_3d_structure(chroms=None):
    # No argument all chromosomes
    # list for specific chromosomes
    if chroms is None:
        chros = ['chr{}'.format(i) for i in list(range(1, 20))]
    else:
        chros = ['chr{}'.format(i) for i in chroms]

    # Generate 3D structure
    for chro in chros:
        command = " ".join(["java -jar",
                            jar_3dmax,
                            params_file.format(chro)])
        os.system(command)


eigen_pvalue(datatype='norm', resolution=1000000, save_fig=False)

# generate_3d_structure()