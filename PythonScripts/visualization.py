# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cooltools
import cooler
from cooltools import insulation
from matplotlib.lines import Line2D

from utlis import insulation_plot, boundaries, balance_matrix, difference_matrix, get_windows, ims, ab_plot, \
    ab_plot_overlay, convert_bytes
from pathlib import Path
import bioframe
import os, subprocess
from scipy.stats import zscore


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
samples = ['omf', 'ymf']


pairs = [('of', 'of'),  ('of', 'om'),  ('of', 'omf'),  ('of', 'yf'),  ('of', 'ym'),  ('of', 'ymf'),
         ('om', 'of'),  ('om', 'om'),  ('om', 'omf'),  ('om', 'yf'),  ('om', 'ym'),  ('om', 'ymf'),
         ('omf', 'of'), ('omf', 'om'), ('omf', 'omf'), ('omf', 'yf'), ('omf', 'ym'), ('omf', 'ymf'),
         ('yf', 'of'),  ('yf', 'om'),  ('yf', 'omf'),  ('yf', 'yf'),  ('yf', 'ym'),  ('yf', 'ymf'),
         ('ym', 'of'),  ('ym', 'om'),  ('ym', 'omf'),  ('ym', 'yf'),  ('ym', 'ym'),  ('ym', 'ymf'),
         ('ymf', 'of'), ('ymf', 'om'), ('ymf', 'omf'), ('ymf', 'yf'), ('ymf', 'ym'), ('ymf', 'ymf'),]

pairs = [('omf', 'ymf'), ('ymf', 'omf'),]


regions = [('Section-Zscore--3.05', 'chr1', 30500000, 35500000),
           ('Section-Zscore--2.59', 'chr1', 149500000, 154500000),
           ('Section-Zscore--2.20', 'chr1', 113000000, 118000000),
           ('Section-Zscore-- -1.85', 'chr1', 85000000, 90000000),]

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
            plt.savefig("../Plots/AB/{}_{}_{}.png".format(chro, datatype, resolution))
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
            plt.savefig("../Plots/AB_Overlay/{}_{}_{}.png".format(chro, datatype, resolution))
        else:
            plt.show()

def difference(x):
    return abs(x['omf'] - x['ymf'])


def filter(x):
    covered_regions = []

    print(x['z_score'])

    exit()

    if x[(x['z_score'] < 1.5) | (x['z_score'] > -1.5)]:
        return False
    # else:
    #     for i in covered_regions:
    #         covered = False
    #         if x[(x['start'] > i[1]) | (x['end1'] < i[0])]:
    #             covered_regions.append((x['start'], x['end1']))
    #             return True




def eigen_pvalue(datatype, resolution, save_fig=False, window=50):
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
        df = pd.DataFrame()
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

            # eigenvector_track.reset_index(drop=True, inplace=True

            if pos == 0:
                df = eigenvector_track
                df.rename({'E1': i}, axis=1, inplace=True)
            else:
                df = pd.merge(df, eigenvector_track, how='inner', left_on=['chrom', 'start', 'end'],
                              right_on=['chrom', 'start', 'end'])
                df.rename({'E1': i}, axis=1, inplace=True)

        df['diff_omf_ymf'] = difference(df)
        df = df[['chrom', 'start', 'end', 'diff_omf_ymf']]

        df['diff_omf_ymf'] = df['diff_omf_ymf'].rolling(window).mean()
        diff = df.loc[0].at['end']
        df['end1'] = df['end'] + (window-1)*diff

        # df = df.dropna()
        df = df.fillna(0)

        df['z_score'] = zscore(df['diff_omf_ymf'])

        df.sort_values(by=['z_score'], inplace=True, key=abs, ascending=False)


        ## Not so good coding practice here.
        start = df['start'].to_list()
        end1 = df['end1'].to_list()
        z_scores = df['z_score'].to_list()

        assert len(start) == len(end1) == len(z_scores)

        covered_regions = []
        actual_regions = []

        for i, j, k in zip(start, end1, z_scores):
            if k < -1.5 or k > 1.5:
                identified = False
                for x in actual_regions:
                    if i < x[0] and j > x[0]:
                        identified = True
                    elif  i < x[1] and j > x[1]:
                        identified = True
                if identified == False:
                    covered_regions.append(True)
                    actual_regions.append((i, j))
                else:
                    covered_regions.append(False)
            else:
                covered_regions.append(False)

        assert len(covered_regions) == len(z_scores)
        df['filtered'] = covered_regions

        df.sort_values(by=['filtered', 'z_score'], inplace=True, key=abs, ascending=False)

        df.to_csv('../Results/Pvalues/{}.csv'.format(chro), sep='\t', encoding='utf-8', index=False)



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
    for region in regions:
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
# compartment_overlay(datatype='norm', resolution=1000000, save_fig=True)
# eigen_pvalue(datatype='norm', resolution=500000, save_fig=True, window=10)
# generate_balance(save_fig=True)
generate_difference(save_fig=True)
# generate_insulation_scores()
# generate_insulation(save_fig=False)
