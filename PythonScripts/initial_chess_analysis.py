import os, subprocess

import fanc
import numpy as np
import pandas as pd

from PythonScripts.utlis import convert_bytes, plot_signal_noise_ratio

pairs = [('of', 'om'),  ('of', 'omf'), ('of', 'yf'),   ('of', 'ym'),  ('of', 'ymf'),
         ('om', 'of'),  ('om', 'omf'), ('om', 'yf'),   ('om', 'ym'),  ('om', 'ymf'),
         ('omf', 'of'), ('omf', 'om'), ('omf', 'yf'),  ('omf', 'ym'), ('omf', 'ymf'),
         ('yf', 'of'),  ('yf', 'om'),  ('yf', 'omf'),  ('yf', 'ym'),  ('yf', 'ymf'),
         ('ym', 'of'),  ('ym', 'om'),  ('ym', 'omf'),  ('ym', 'yf'),  ('ym', 'ymf'),
         ('ymf', 'of'), ('ymf', 'om'), ('ymf', 'omf'), ('ymf', 'yf'), ('ymf', 'ym'), ]

def create_bed(chrom='all', window_size=50000000, step_size=10000):
    '''
        Generating a pairs input file
        The window size should not be smaller than 20x the bin size
        of the data that will be compared with chess sim.
        :param chrom: 'all', 1-19
        '''
    bed_path = '../Results/Bed/mm10_{}_{}_win_{}_step.bed'
    if chrom == 'all':
        bed_path = bed_path.format('full', convert_bytes(window_size).replace(" ", ""), convert_bytes(step_size).replace(" ", ""))
    else:
        bed_path = bed_path.format(chrom, convert_bytes(window_size).replace(" ", ""), convert_bytes(step_size).replace(" ", ""))

    if not os.path.isfile(bed_path):
        print('Creating Bed file {} '.format(bed_path))

        if chrom == 'all':
            subprocess.call('chess pairs mm10 {} {} {}'.format(window_size, step_size, bed_path), shell=True)
        else:
            subprocess.call('chess pairs mm10  {} {} {} --chromosome {}'.\
                            format(window_size, step_size, bed_path, chrom), shell=True)

    else:
        print("Bed file {} exist".format(bed_path))


def search_similarity(window_size, step_size, pair_1, pair_2):
    print(window_size)
    print(step_size)
    print(pair_1)
    print(pair_2)
    exit()
    '''
        run the search with the chess sim subcommand:
    :return:
    '''
    chros = ['chr{}'.format(i) for i in list(range(1, 20))] # + ['all', 'chrX', 'chrY']
    for i in chros:
        create_bed(i, window_size, step_size)
        chess_sim_results = '../Results/SIM/{}_{}_{}_vs_{}_chess_results.tsv'.format(i, convert_bytes(window_size).replace(" ", ""), pair_1, pair_2)
        if not os.path.isfile(chess_sim_results):
            print('Generating similarity for {} {} and {}'.format(i, pair_1, pair_2))
            subprocess.call('chess sim \
                                        ../MuSC_HiC_files/HiC_CPB_normalized/{}_new_CPBnorm_inter_30.hic \
                                        ../MuSC_HiC_files/HiC_CPB_normalized/{}_new_CPBnorm_inter_30.hic \
                                        ../Results/Bed/mm10_{}_{}_win_{}_step.bed \
                                        {} \
                                        -p 4'.format(pair_1, pair_2, i, convert_bytes(window_size).replace(" ", ""),
                                                     convert_bytes(step_size).replace(" ", ""), chess_sim_results),
                            shell=True)
        else:
            print("Similarity file {} exist".format(chess_sim_results))
        chess_extract_results = '../Results/Features/{}_{}_{}'.format(i, pair_1, pair_2)
        if not os.path.isfile(chess_extract_results):
            os.mkdir(chess_extract_results)
            print('Extracting features for {} {} and {}'.format(i, pair_1, pair_2))
            print(chess_extract_results)
            subprocess.call('chess extract \
                                        ../Results/Bed/mm10_{}_{}_win_{}_step.bed \
                                        ../MuSC_HiC_files/HiC_CPB_normalized/{}_new_CPBnorm_inter_30.hic \
                                        ../MuSC_HiC_files/HiC_CPB_normalized/{}_new_CPBnorm_inter_30.hic \
                                        {} \
                                        '.format(i, convert_bytes(window_size).replace(" ", ""),
                                                     convert_bytes(step_size).replace(" ", ""),
                                                     pair_1, pair_2, chess_extract_results),
                            shell=True)
        else:
            print("Features already extracted".format(chess_extract_results))



# window_size=50000000, step_size=10000
def examine_similarity(window_size=5000000, step_size=10000):
    for i in pairs:
        search_similarity(window_size, step_size, i[0].upper(), i[1].upper())


def vis_sim_results(window_size=50000000, step_size=10000, resolution=1000000, features=False, save_fig=False):
    sn_thr = 0.5
    sim_field = "z_ssim"
    sn_thr = sn_thr
    zsim_thr = -1

    chros = ['chr{}'.format(i) for i in list(range(1, 20))]  # + ['all', 'chrX', 'chrY']
    for chrom in chros:
        for pair in pairs:
            regions = pd.read_csv('../Results/Bed/mm10_{}_{}_win_{}_step.bed'. \
                                  format(chrom, convert_bytes(window_size).replace(" ", ""), convert_bytes(step_size).replace(" ", "")), sep='\t', header=None)
            similarities = pd.read_csv('../Results/SIM/{}_{}_{}_vs_{}_chess_results.tsv'. \
                                  format(chrom, convert_bytes(window_size).replace(" ", ""), pair[0].upper(), pair[1].upper()), sep='\t', index_col=0)

            sub_sim = similarities[(similarities["SN"] >= sn_thr) & (similarities[sim_field] <= zsim_thr)]

            pair_A_hic = fanc.load("../MuSC_HiC_files/HiC_CPB_normalized/{}_new_CPBnorm_inter_30.hic@{}".format(pair[0].upper(), resolution))
            pair_B_hic = fanc.load("../MuSC_HiC_files/HiC_CPB_normalized/{}_new_CPBnorm_inter_30.hic@{}".format(pair[1].upper(), resolution))

            # Get the region with the highest dissimilarity
            region_id = similarities[(similarities["SN"] >= sn_thr)].sort_values("z_ssim").index.values[0]
            window_start, window_end = regions.loc[region_id][1:3]

            region_string = "{}:{}-{}".format(chrom, window_start, window_end)
            pair_A_sub = pair_A_hic[region_string, region_string].data
            pair_B_sub = pair_B_hic[region_string, region_string].data

            min_v = min(
                [
                    np.min(np.extract(pair_A_sub > 0, pair_A_sub)),
                    np.min(np.extract(pair_B_sub > 0, pair_B_sub))
                ]
            )

            pair_A_sub += min_v
            pair_B_sub += min_v
            l2fcm = np.log2(pair_A_sub / pair_B_sub)


            # Feature Extraction

            ## obtaining regions of interest
            # regions2compare = regions.loc[sub_sim.index]
            # regions2compare.to_csv('filtered_regions_chr2_{}_100kb.tsv'.format(winsize), '\t', index=False, header=False)

            ## load gained and lost features
            if features:
                gained = pd.read_csv('../Results/Features/{}_{}_{}/gained_features.tsv'.format(chrom, pair[0].upper(), pair[1].upper()), delimiter=',', usecols=[0, 1, 2, 3, 4, 5],
                                     header=None, index_col=[0])
                lost = pd.read_csv('../Results/Features/{}_{}_{}/lost_features.tsv'.format(chrom, pair[0].upper(), pair[1].upper()), delimiter=',', usecols=[0, 1, 2, 3, 4, 5],
                                   header=None, index_col=[0])

                plot_signal_noise_ratio(sub_sim=sub_sim,
                                        sim_field=sim_field,
                                        zsim_thr=zsim_thr,
                                        region_id=region_id,
                                        chrom=chrom,
                                        data1=(regions, similarities),
                                        data2=(pair_A_sub, pair_B_sub, l2fcm),
                                        pairs=pair,
                                        save_fig=save_fig,
                                        gained=gained,
                                        lost=lost,
                                        add_features=features)
            else:
                plot_signal_noise_ratio(sub_sim=sub_sim,
                                        sim_field=sim_field,
                                        zsim_thr=zsim_thr,
                                        region_id=region_id,
                                        chrom=chrom,
                                        data1=(regions, similarities),
                                        data2=(pair_A_sub, pair_B_sub, l2fcm),
                                        pairs=pair,
                                        save_fig=save_fig,
                                        add_features=features)


# search_similarity()
examine_similarity(window_size=3000000, step_size=10000)
# vis_sim_results(features=False, save_fig=True)
# import matplotlib.pyplot as plt
# fig, axes = plt.subplots(1, 3, figsize=(16, 8))
# print(axes)

