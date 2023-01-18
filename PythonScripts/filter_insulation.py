import bioframe
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

INSULATION_PATH = '../Results/Insulation/insulation_score_{}_{}.csv'
resolutions = [('50kb', 50000), ('100kb', 100000), ('500kb', 500000), ('1mb', 1000000)]
samples = ['of', 'om', 'omf', 'yf', 'ym', 'ymf']


for resolution in resolutions:
    for sample in samples:
        file_name = INSULATION_PATH.format(sample, resolution[1])
        df = pd.read_csv(file_name)
        # df = df[(df.is_bad_bin == False)]

        insul_region = bioframe.select(df, ('chr10', 18628620, 32129231))

        insul_region = insul_region[(insul_region.chrom == 'chr10')]

        insul_region.rename(columns={'log2_insulation_score_150000': 'log2',
                                     'boundary_strength_150000': 'bds',
                                     'is_boundary_150000': 'isbd'}, inplace=True)

        print(insul_region.head(100)[['chrom', 'start', 'end', 'region',
                                     'log2',
                                     'bds',
                                     'isbd',
                                     'n_valid_pixels_150000']])
        exit()



# for resolution in resolutions[0:1]:
#
#
# clr = cooler.Cooler('{}'.format(DATA_PATH.format(sample, resolution[1])))
# windows = get_windows(resolution[1])
# insulation_df = insulation(clr, windows, verbose=True)
# insulation_df.index.name = 'ID'
# insulation_df.to_csv(file_name)