import cooler
import gzip


samples = ['MuSC_HiC_files/HiC_CPB_normalized', 'MuSC_HiC_files/HiC_not_normalized']
groups = ['omf', 'ymf']
resolutions = [1000000, 500000,  100000, 50000,  10000, 5000]
for i in samples:
    for j in groups:
        for k in resolutions:
            c = cooler.Cooler("{}/normalized_{}_{}.cool".format(i, j, k))
            #c = cooler.Cooler("{}/not_normalized_{}_{}.cool".format(i, j, k))


            print(c.bins())
            for i in c.bins():
                print(i)
            exit()


            print("Generating matrix for {}".format(c))
            mat_df = c.matrix(as_pixels=True, join=True)[:, :]

            mat_df = mat_df.drop(['end1', 'end2', 'balanced'], axis=1)
            out_file = "repqc/{}_{}_mat_{}.csv".format(i, j, k)
            mat_df.to_csv(out_file, index=False, sep='\t', header=False)

            # with open(out_file, 'rb') as src, gzip.open(out_file+".gz", 'wb') as dst:
            #     dst.writelines(src)


            print("Generating Bins for {} {}".format(c))
            bins_df = c.bins()[:]
            bins_df = bins_df.drop('weight', axis=1)
            bins_df['name'] = bins_df['start']
            out_file = "repqc/{}_{}_bin_{}.csv".format(i, j, k)
            bins_df.to_csv(out_file, index=False, sep='\t', header=False)

            # with open(out_file, 'rb') as src, gzip.open(out_file+".gz", 'wb') as dst:
            #     dst.writelines(src)
