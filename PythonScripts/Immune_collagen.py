import pandas as pd

collagen_data = pd.read_excel("../MuSC_HiC_files/immune_collagen_genelist.xlsx", sheet_name="collagen genes")
GTF_PATH = '../MuSC_HiC_files/fasta/mm10.refGene.gtf'

def getGene(df, annot='gene_id'):
    if annot == 'gene_id':
        x = df[8].split()[1].replace(";", "").replace('"', "")
    return x

immune_data = pd.read_excel("../MuSC_HiC_files/immune_collagen_genelist.xlsx", sheet_name="immune gene",
                            header=None, names=['Gene'])
immune_data = set(immune_data['Gene'].tolist())
print(immune_data)

chros = ['chr{}'.format(i) for i in list(range(1, 20))] + ["chrX"]
colnames = {0: 'Gene_Chrom', 2: 'Annotation', 3: 'Gene_Start', 4: 'Gene_End'}
gtf_df = pd.read_csv(GTF_PATH, sep='\t', header=None)
gtf_df = gtf_df[(gtf_df[2] == 'transcript') & (gtf_df[0].isin(chros))]
gtf_df['Gene'] = gtf_df.apply(getGene, annot='gene_id', axis=1)
gtf_df = gtf_df[[0, 2, 3, 4, 'Gene']]
gtf_df.rename(index=str, columns=colnames, inplace=True)
gtf_df = gtf_df[(gtf_df['Gene'].isin(immune_data))]

gtf_df.drop_duplicates(subset=['Gene_Chrom', 'Annotation',  'Gene_Start',
                           'Gene_End', 'Gene'], inplace=True)

# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)
# print(gtf_df)
# print(len(immune_data))
#
# exit()

resolutions = [1000000, 500000, 100000, 50000]
print(len(immune_data))
for i in resolutions:
    identified_genes = pd.read_csv("../Results/HiCCompare_OM_OF/IdentifiedGenes/{}.csv".format(i))
    identified_genes = set(identified_genes['Gene'].tolist())
    print(len(identified_genes.intersection(immune_data)), identified_genes.intersection(immune_data))

resolutions = [1000000, 500000, 100000, 50000]
print(len(immune_data))
for i in resolutions:
    identified_genes = pd.read_csv("../Results/HiCCompare_YM_YF/IdentifiedGenes/{}.csv".format(i))
    identified_genes = set(identified_genes['Gene'].tolist())
    print(len(identified_genes.intersection(immune_data)), identified_genes.intersection(immune_data))