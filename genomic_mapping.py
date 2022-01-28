import pandas as pd
import cooler
# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)


def getGene(df):
    x = df[8].split()[1].replace(";", "").replace('"', "")
    return x


_df = pd.read_csv('/storage/htc/bdm/Frimpong/ref_geneomes/mm39/refGene.gtf', sep='\t', header=None)
_df[9] = _df.apply(getGene,  axis=1)
_df = _df[(_df[2] == 'transcript')]

cooler_files = ["27950", "27951", "27989", "28033"]
for cooler_file in cooler_files:
    file_name = cooler_file+"/output_1mb.cool"
    c = cooler.Cooler(file_name)
    chrom_regions = c.chromsizes.to_dict()
    binsize = c.binsize
    _mat_df = c.matrix(as_pixels=True, join=True)[:]

    writer = pd.ExcelWriter(cooler_file+"/"+cooler_file+"_1mb_genomic.xlsx", engine='xlsxwriter')
    regions_of_interest = ['Mir155', 'Mrpl39', 'Zbtb21']
    for i in regions_of_interest:

        df = _df[(_df[9] == i)].drop_duplicates(subset=[0, 3, 4, 9], keep="first")
        start = int(df[3].values[0]/binsize)*binsize
        end = int(df[4].values[0]/binsize)*binsize + binsize
        print(start, end)

        mat_df = _mat_df[((_mat_df['start1'] >= start) & (_mat_df['end1'] <= end) & (_mat_df['chrom1'] == 'chr16')) |
                         ((_mat_df['start2'] >= start) & (_mat_df['end2'] <= end) & (_mat_df['chrom2'] == 'chr16')) ]

        mat_df.to_excel(writer, sheet_name=i)
    writer.save()
