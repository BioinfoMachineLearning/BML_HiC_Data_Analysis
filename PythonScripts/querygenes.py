import pandas as pd
from pathlib import Path
import cooler
# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)


# ranges are defined as (start, end)

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


def getGene(df, index=1):
    x = df[8].split()[index].replace(";", "").replace('"', "")
    return x

def get_relevant_info():
    GTF_PATH = '../Results/GTF/gtf.csv'
    if Path(GTF_PATH).is_file():
        df = pd.read_csv(GTF_PATH, sep='\t')
    else:
        RAW_GTF_PATH = '../MuSC_HiC_files/fasta/mm10.refGene.gtf'
        df = pd.read_csv(RAW_GTF_PATH, sep='\t', header=None)
        df = df[(df[2] == 'transcript')]
        df.drop([1, 2, 5, 6, 7], axis=1, inplace=True)
        df['gene_id'] = df.apply(getGene, axis=1)
        df['transcript_id'] = df.apply(getGene, index=3, axis=1)
        df['gene_name'] = df.apply(getGene, index=5, axis=1)
        df = df.rename({0: 'Chromosome', 3: 'Start', 4: 'End'}, axis=1)
        df.drop([8], axis=1, inplace=True)
        df.to_csv(GTF_PATH, index=False, sep='\t')
    return df


def search():
    final_df = pd.DataFrame()
    GTF_PATH = '../Results/GTF/gtf.csv'
    if Path(GTF_PATH).is_file():
        df = pd.read_csv(GTF_PATH, sep='\t')
    for i in regions:
        print(i[0:2])
        x = df.loc[(df['Chromosome'] == i[1]) & (df['Start'] >= i[2]) & (df['End'] <= i[3])]
        print(x[['Start', 'End', 'transcript_id', 'gene_name']])
        final_df = final_df.append(x, sort=False)
    final_df.to_csv("../Results/GeneMapping/genes.csv", index=False)

# df = get_relevant_info()
# print(df.dtypes)
search()