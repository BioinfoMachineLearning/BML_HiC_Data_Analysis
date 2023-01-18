import bioframe as bf
import numpy as np
import pandas as pd
import os
from pathlib import Path

resolutions = ["50000", "100000", "500000", "1000000"]
chros = ['chr{}'.format(i) for i in list(range(1, 19))] + ['chrX']
groups = ["Old_Young", "Male_Female"]
identified_gene_path = "../{}/"
GTF_PATH = '../MuSC_HiC_files/fasta/mm10.refGene.gtf'
bed_pth = "../Results/{}/Results/HiCCompare_{}/{}/{}/{}_top_pair_bed.csv"
save_pth = "../Results/{}/Results/HiCCompare_{}/IdentifiedGenes/{}"


def make_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)


def getGene(df, annot='gene_id'):
    if annot == 'gene_id':
        x = df[8].split()[1].replace(";", "").replace('"', "")
    return x


def is_unique(s):
    a = s.to_numpy()
    return (a[0] == a).all()


chros = ['chr{}'.format(i) for i in list(range(1, 20))] + ["chrX"]
colnames = {0: 'Gene_Chrom', 2: 'Annotation', 3: 'Gene_Start', 4: 'Gene_End'}
gtf_df = pd.read_csv(GTF_PATH, sep='\t', header=None)
gtf_df = gtf_df[(gtf_df[2] == 'transcript') & (gtf_df[0].isin(chros))]
gtf_df['Gene'] = gtf_df.apply(getGene, annot='gene_id', axis=1)
gtf_df = gtf_df[[0, 2, 3, 4, 'Gene']]
gtf_df.rename(index=str, columns=colnames, inplace=True)
sheets = ["H3K27ac_unique_enhancer_list", "H3K27ac_increased_in_old", "H3K27ac_decreased_in_old",
          "H3K4me1_decreased_in_old", "H3K4me1_increase_in_old"]

for sheet in sheets:

    enhancer_file = "../MuSC_HiC_files/H3K4me1_H3K27ac_differential_peaks.xlsx"
    enhancer_df = pd.read_excel(enhancer_file, sheet_name=sheet)
    colsnames = {"seqnames": "CHROM", "start": "START", "end": "STOP", 4: 'Gene_End'}
    enhancer_df.rename(index=str, columns=colsnames, inplace=True)

    writer = pd.ExcelWriter('../MuSC_HiC_files/{}.xlsx'.format(sheet), engine='xlsxwriter')
    for group in groups:
        # make_directory(identified_gene_path.format(group))
        for resolution in resolutions:
            print("Generating for sheet {} resolution {}".format(sheet, resolution))
            # Path(save_pth.format(group, "")).mkdir(parents=True, exist_ok=True)
            dfs = []
            for chr in chros:
                # print("Generating for Chromosome {}".format(chr))
                ch_df = pd.DataFrame()
                current_bed = bed_pth.format(group, group, resolution, chr, chr)
                if os.path.isfile(current_bed):
                    ch_df = pd.read_csv(current_bed)
                    if chr == 'chrX':
                        if not ch_df.empty:
                            ch_df['chr1'] = 'chrX'
                            ch_df['chr2'] = 'chrX'
                    else:
                        if not ch_df.empty:
                            assert is_unique(ch_df['chr1'])
                            assert is_unique(ch_df['chr2'])

                if not ch_df.empty:
                    colnames = {"CHROM": "Enhancer Chromosome", "START": "Enhancer Start", "STOP": "Enhancer End",
                                "Gene_Chrom": "Loop Chromosome", "Gene_Start": "Loop Start", "Gene_End": "Loop End",
                                "Gene": "Loop Gene"}

                    filtered_gtf = gtf_df[(gtf_df['Gene_Chrom'] == chr)]
                    filtered_gtf = filtered_gtf.astype({"Gene_Start": np.int64, "Gene_End": np.int64})

                    # Left side
                    # Join Unique enhancer to differential regions
                    left_temp = bf.overlap(ch_df, enhancer_df, how="inner",
                                           cols1=('chr1', 'start1', 'end1'),
                                           cols2=('CHROM', 'START', 'STOP'),
                                           suffixes=("", ""))
                    # Join gtf
                    left_temp = bf.overlap(left_temp, filtered_gtf, how="inner",
                                           cols1=('chr2', 'start2', 'end2'),
                                           cols2=('Gene_Chrom', 'Gene_Start', 'Gene_End'),
                                           suffixes=("", ""))

                    left_temp.rename(index=str, columns=colnames, inplace=True)
                    left_temp = left_temp.drop_duplicates()
                    dfs.append(left_temp)

                    # Right Side
                    right_temp = bf.overlap(ch_df, enhancer_df, how="inner",
                                            cols1=('chr2', 'start2', 'end2'),
                                            cols2=('CHROM', 'START', 'STOP'),
                                            suffixes=("", ""))

                    right_temp = bf.overlap(right_temp, filtered_gtf, how="inner",
                                            cols1=('chr1', 'start1', 'end1'),
                                            cols2=('Gene_Chrom', 'Gene_Start', 'Gene_End'),
                                            suffixes=("", ""))
                    right_temp.rename(index=str, columns=colnames, inplace=True)
                    right_temp = right_temp.drop_duplicates()
                    dfs.append(right_temp)

            if len(dfs) > 0:
                df = pd.concat(dfs, ignore_index=True)
                df = df.drop_duplicates(["chr1", "start1", "end1", "Loop Gene"])
                df = df.drop_duplicates(["chr2", "start2", "end2", "Loop Gene"])
                df = df.drop_duplicates(["Loop Chromosome", "Loop Start", "Loop End", "Loop Gene"])
                cols = list(df.columns)
                cols = cols[11: 24] + cols[0: 11]
                df = df[cols]
            df.to_excel(writer, sheet_name=group + "_" + resolution)
    writer.save()
