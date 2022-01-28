# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import cooltools
import cooler
from cooltools import insulation
# from cooltools.insulation import calculate_insulation_score, find_boundaries
from matplotlib.lines import Line2D

from utlis import coverage, insulation_plot, boundaries, balance_matrix, difference_matrix, get_windows, ims, ab_plot, \
    ab_plot_overlay, convert_bytes
from pathlib import Path
import itertools
import bioframe
import os, subprocess


DATA_PATH = '../MuSC_HiC_files/HiC_not_normalized/not_normalized_{}_{}.cool'
DATA_PATH = '../MuSC_HiC_files/HiC_CPB_normalized/normalized_{}_{}.cool'
# DATA_PATH = '../MuSC_HiC_files/raw/{}_{}.cool'
INSULATION_PATH = '../Results/Insulation/insulation_score_{}_{}.csv'

paths = {
    'of': 'OF_new_CPBnorm_inter_30.hic',
    'om': 'OM_new_CPBnorm_inter_30.hic',
    'omf': 'OMF_new_CPBnorm_inter_30.hic',
    'yf': 'YF_new_CPBnorm_inter_30.hic',
    'ym': 'YM_new_CPBnorm_inter_30.hic',
    'ymf': 'OM_new_CPBnorm_inter_30.hic'
}

resolutions = [('50kb', 50000), ('100kb', 100000), ('500kb', 500000), ('1mb', 1000000)]
samples = ['of', 'om', 'omf', 'yf', 'ym', 'ymf']
# samples = ['of', 'om', 'yf', 'ym']


pairs = [('of', 'of'),  ('of', 'om'),  ('of', 'omf'),  ('of', 'yf'),  ('of', 'ym'),  ('of', 'ymf'),
         ('om', 'of'),  ('om', 'om'),  ('om', 'omf'),  ('om', 'yf'),  ('om', 'ym'),  ('om', 'ymf'),
         ('omf', 'of'), ('omf', 'om'), ('omf', 'omf'), ('omf', 'yf'), ('omf', 'ym'), ('omf', 'ymf'),
         ('yf', 'of'),  ('yf', 'om'),  ('yf', 'omf'),  ('yf', 'yf'),  ('yf', 'ym'),  ('yf', 'ymf'),
         ('ym', 'of'),  ('ym', 'om'),  ('ym', 'omf'),  ('ym', 'yf'),  ('ym', 'ym'),  ('ym', 'ymf'),
         ('ymf', 'of'), ('ymf', 'om'), ('ymf', 'omf'), ('ymf', 'yf'), ('ymf', 'ym'), ('ymf', 'ymf'),]

# pairs = [('of', 'of'),  ('of', 'om'),  ('of', 'yf'),  ('of', 'ym'),
#          ('om', 'of'),  ('om', 'om'),  ('om', 'yf'),  ('om', 'ym'),
#          ('yf', 'of'),  ('yf', 'om'),  ('yf', 'yf'),  ('yf', 'ym'),
#          ('ym', 'of'),  ('ym', 'om'),  ('ym', 'yf'),  ('ym', 'ym'),]

regions = [('Section F', 'chr8', 101000000, 103000000), ('Section G', 'chr8', 15000000, 19000000),
           ('Section 15-A', 'chr15', 15000000, 20000000), ('Section 15-AB', 'chr15', 15000000, 26000000),
           ('Section 15-B', 'chr15', 21000000, 26000000),
           ('Section 1-A', 'chr1', 25000000, 35000000), ('Section 1-B', 'chr1', 90000000, 112500000),
           ('Section 1-C', 'chr1', 140000000, 150000000), ('Section 1-C', 'chr1', 140000000, 150000000),
           ('Section 1-C', 'chr1', 140000000, 150000000),
           ('Section 2-A', 'chr2', 25000000, 30000000),
           ('Section 14-A', 'chr14', 11000000, 13000000),  ('Section 14-B', 'chr14', 38000000, 40000000),
           ('Section 14-C', 'chr14', 96000000, 98000000),  ('Section 14-D', 'chr14', 108000000, 110800000),
           ('Section 3-A', 'chr3', 10000000, 15000000), ('Section 3-B', 'chr3', 41000000, 45000000),
           ('Section 6-A', 'chr6', 11000000, 12500000), ('Section 6-B', 'chr6', 19000000, 21000000),
           ('Section 6-C', 'chr6', 79000000, 81000000), ('Section 6-D', 'chr6', 101000000, 107000000),
           ('Section 6-E', 'chr6', 132000000, 136000000),
           ('Section 12-A', 'chr12', 43000000, 45000000),
           ('Section 7-A', 'chr7', 53000000, 56000000), ('Section 7-B', 'chr7', 85000000, 86000000),
           ('Section 17-A', 'chr17', 18000000, 20000000), ('Section 17-B', 'chr17', 38000000, 41000000),
           ('Section 17-C', 'chr17', 53000000, 56000000), ('Section 17-D', 'chr17', 58000000, 62000000),
           ('Section 19-A', 'chr19', 59000000, 60000000),
           ('Section 18-A', 'chr18', 16000000, 20000000),
           ('Section 11-A', 'chr11', 10000000, 11000000), ('Section 11-B', 'chr11', 36000000, 42000000),
           ('Section 10-A', 'chr10', 122000000, 124000000),]

tad_regions = [('Section G', 'chr8', 5000000, 40000000), ('Section F', 'chr8', 10000000, 30000000)]


def compartmentalization(datatype, resolution, save_fig=True):
    # fasta sequence is required for calculating binned profile of GC conent
    if not os.path.isfile('../MuSC_HiC_files/fasta/mm10.fa'):
        # note downloading a ~1Gb file can take a minute
        subprocess.call(
            'wget -O ../MuSC_HiC_files/fasta/mm10.fa.gz https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz',
            shell=True)
        subprocess.call('gunzip ../MuSC_HiC_files/fasta/mm10.fa.gz', shell=True)

    # Using 1mb for AB compartment
    chros = ['chr{}'.format(i) for i in list(range(1, 20))]
    for chro in chros:
        fig, ax = plt.subplots(6, 1, figsize=(20, 16))
        for pos, i in enumerate(samples):
            if datatype == 'raw':
                clr = cooler.Cooler('../MuSC_HiC_files/raw/{}_{}.cool'.format(i.upper(), resolution))
            elif datatype == 'norm':
                clr = cooler.Cooler('../MuSC_HiC_files/HiC_CPB_normalized/normalized_{}_{}.cool'.format(i, resolution))
            elif datatype == 'not norm':
                clr = cooler.Cooler('../MuSC_HiC_files/HiC_not_normalized/not_normalized_{}_{}.cool'.format(i, resolution))

            # Check if fraction of GC file exist
            GC_PATH = '../Results/GC/{}_mm10_gc_cov_{}_{}.tsv'.format(datatype, i, resolution)
            if Path(GC_PATH).is_file():
                gc_cov = pd.read_csv(GC_PATH, sep='\t')
            else:
                bins = clr.bins()[:]
                mm10_genome = bioframe.load_fasta('../MuSC_HiC_files/fasta/mm10.fa');
                #  note the next command may require installing pysam
                gc_cov = bioframe.frac_gc(bins, mm10_genome)
                gc_cov.to_csv(GC_PATH, index=False, sep='\t')

            view_df = pd.DataFrame({'chrom': clr.chromnames, 'start': 0, 'end': clr.chromsizes.values, 'name': clr.chromnames})
            cis_eigs = cooltools.eigs_cis(clr, gc_cov, view_df=view_df, n_eigs=3, phasing_track_col='GC',)

            # cis_eigs[0] returns eigenvalues, here we focus on eigenvectors
            eigenvector_track = cis_eigs[1][['chrom', 'start', 'end', 'E1']]
            eigenvector_track = eigenvector_track.loc[eigenvector_track['chrom'] == chro]
            eigenvector_track.reset_index(drop=True, inplace=True)

            ab_plot(eigenvector_track, ax[pos], i)

            ax[-1].set_xlabel("genomic loci")
            ax[0].set_title("A-B Compartment Chromosome --- " + str(chro), fontsize=18)
            cmap = plt.cm.coolwarm
            custom_lines = [Line2D([0], [0], color=cmap(0.0), lw=4),
                            Line2D([0], [0], color=cmap(1.0), lw=4)]
            ax[0].legend(custom_lines, ['B', 'A'], loc="upper right", fontsize=10)

        if save_fig:
            plt.savefig("../plots/AB/{}_{}_{}.png".format(chro, datatype, resolution))
        else:
            plt.show()



def compartment_overlay(datatype, resolution, save_fig=False):
    # fasta sequence is required for calculating binned profile of GC content
    if not os.path.isfile('../MuSC_HiC_files/fasta/mm10.fa'):
        # note downloading a ~1Gb file can take a minute
        subprocess.call(
            'wget -O ../MuSC_HiC_files/fasta/mm10.fa.gz https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz',
            shell=True)
        subprocess.call('gunzip ../MuSC_HiC_files/fasta/mm10.fa.gz', shell=True)

    # Using 1mb for AB compartment
    chros = ['chr{}'.format(i) for i in list(range(1, 20))]
    for chro in chros:
        # fig, ax = plt.subplots(figsize=(20, 16))
        fig, ax = plt.subplots(figsize=(12, 6), dpi=120)
        for pos, i in enumerate(samples):
            if datatype == 'raw':
                clr = cooler.Cooler('../MuSC_HiC_files/raw/{}_{}.cool'.format(i.upper(), resolution))
            elif datatype == 'norm':
                clr = cooler.Cooler('../MuSC_HiC_files/HiC_CPB_normalized/normalized_{}_{}.cool'.format(i, resolution))
            elif datatype == 'not norm':
                clr = cooler.Cooler('../MuSC_HiC_files/HiC_not_normalized/not_normalized_{}_{}.cool'.format(i, resolution))

            # Check if fraction of GC file exist
            GC_PATH = '../Results/GC/{}_mm10_gc_cov_{}_{}.tsv'.format(datatype, i, resolution)
            if Path(GC_PATH).is_file():
                gc_cov = pd.read_csv(GC_PATH, sep='\t')
            else:
                bins = clr.bins()[:]
                mm10_genome = bioframe.load_fasta('../MuSC_HiC_files/fasta/mm10.fa');
                #  note the next command may require installing pysam
                gc_cov = bioframe.frac_gc(bins, mm10_genome)
                gc_cov.to_csv(GC_PATH, index=False, sep='\t')

            view_df = pd.DataFrame({'chrom': clr.chromnames, 'start': 0, 'end': clr.chromsizes.values, 'name': clr.chromnames})
            cis_eigs = cooltools.eigs_cis(clr, gc_cov, view_df=view_df, n_eigs=3, phasing_track_col='GC',)

            # cis_eigs[0] returns eigenvalues, here we focus on eigenvectors
            eigenvector_track = cis_eigs[1][['chrom', 'start', 'end', 'E1']]
            eigenvector_track = eigenvector_track.loc[eigenvector_track['chrom'] == chro]
            eigenvector_track.reset_index(drop=True, inplace=True)



            ab_plot_overlay(eigenvector_track, ax, (pos, i))
        ax.set_ylabel('Sample')
        plt.legend(loc="upper right")
        plt.title('AB Compartment overlay --- {} --- {}'.format(convert_bytes(resolution), chro))

        if save_fig:
            plt.savefig("../plots/AB_Overlay/{}_{}_{}.png".format(chro, datatype, resolution))
        else:
            plt.show()



def generate_balance(save_fig=True):
    for region in regions[2:]:
        for resolution in resolutions:
            loaded_coolers = {}
            for sample in samples:
                print("generating for {} {}".format(sample, resolution))
                print('{}'.format(DATA_PATH.format(sample, resolution[1])))
                tmp = cooler.Cooler('{}'.format(DATA_PATH.format(sample, resolution[1]))).matrix().fetch(region[1:])
                loaded_coolers[sample] = tmp
            balance_matrix(loaded_coolers, region=region, res=resolution[0], save_fig=save_fig, data_type='norm')


def generate_difference(save_fig=True):
    for region in regions[28:29]:
        print(region)
        for resolution in resolutions[1:]:
            loaded_coolers = {}
            for pair in pairs:
                tmp1 = cooler.Cooler('{}'.format(DATA_PATH.format(pair[0], resolution[1]))).matrix().fetch(region[1:])
                tmp1 = tmp1 - np.eye(tmp1.shape[0]) * tmp1

                tmp2 = cooler.Cooler('{}'.format(DATA_PATH.format(pair[1], resolution[1]))).matrix().fetch(region[1:])
                tmp2 = tmp2 - np.eye(tmp2.shape[0]) * tmp2

                loaded_coolers[pair[0]+'_'+pair[1]] = tmp1 - tmp2
            print("ploting for {}".format(resolution))
            difference_matrix(loaded_coolers, region=region, res=resolution[0], save_fig=save_fig, data_type='norm')


def generate_insulation_scores():
    for sample in samples:
        for resolution in resolutions[1:]:
            file_name = INSULATION_PATH.format(sample, resolution[1])
            if Path(file_name).is_file():
                 pass
            else:
                 clr = cooler.Cooler('{}'.format(DATA_PATH.format(sample, resolution[1])))
                 windows = get_windows(resolution[1])
                 insulation_df = insulation(clr, windows, verbose=True)
                 insulation_df.index.name = 'ID'
                 insulation_df.to_csv(file_name)



def generate_insulation(save_fig=True):
    for tad_region in tad_regions[1:]:
        for resolution in resolutions:
            loaded_coolers = {}
            for sample in samples:
                clr = cooler.Cooler('{}'.format(DATA_PATH.format(sample, resolution[1])))
                loaded_coolers[sample + '_' + str(resolution[1])] = [clr, tad_region]

            insulation_plot(loaded_coolers, save_fig=save_fig)
            boundaries(loaded_coolers, save_fig=save_fig)

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
# compartmentalization(datatype='norm', resolution=100000)
# compartment_overlay(datatype='raw', resolution=100000, save_fig=False)
# generate_balance(save_fig=False)
generate_difference(save_fig=True)
# generate_insulation_scores()
# generate_insulation(save_fig=False)
