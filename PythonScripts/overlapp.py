import pandas as pd

resolutions = ["50000", "100000", "500000", "1000000"]

identified_gene_path = "../Results/HiCCompare_{}/IdentifiedGenes/{}.csv"
enhancer_file = "../MuSC_HiC_files/H3K27ac_unique_enhancer_list.xlsx"
groups = ["OM_OF", "YM_YF"]

identified_genes = set()
for group in groups:
    for i in resolutions:
        tmp = set(pd.read_csv(identified_gene_path.format(group, i), sep=',')['Gene'].\
                  to_list())
        enhancer_df = set(pd.read_excel(enhancer_file)['SYMBOL'].to_list())
        print(enhancer_df.intersection(tmp))
        # print(tmp)
        # exit()
        # identified_genes = identified_genes.union(tmp)

    # sheets = ['O_H3K122ac_enhancer_unique', 'O_H3K27ac_enhancer_unique']
    # for sheet in sheets:
    #     enhancer_df = set(pd.read_excel(enhancer_file, sheet_name=sheet)['SYMBOL'].to_list())
    #     print(enhancer_df.intersection(identified_genes))


    # enhancer_df = set(pd.read_excel(enhancer_file, sheet_name=sheet)['SYMBOL'].to_list())
    # print(enhancer_df.intersection(identified_genes))