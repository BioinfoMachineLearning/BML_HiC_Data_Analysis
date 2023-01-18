import subprocess


samples = ['om', 'omf', 'yf', 'ym', 'ymf']
resolutions = [('50kb', 50000), ('100kb', 100000), ('500kb', 500000), ('1mb', 1000000)]
resolutions = [('50kb', 50000)]

# Extract data for multihic compare
for resolution in resolutions:
    for sample in samples:
        print("Resolution:{} --- Sample: {}".format(resolution[0], sample))
        print("normalized")
        exit()
        subprocess.call('cooler dump --join ../MuSC_HiC_files/HiC_CPB_normalized/normalized_{}_{}.cool > '
                        '../MuSC_HiC_files/HiC_CPB_normalized/normalized_{}_{}.txt'.\
                        format(sample, resolution[1], sample, resolution[1]), shell=True)

        print("unnormalized")
        subprocess.call('cooler dump --join ../MuSC_HiC_files/HiC_not_normalized/not_normalized_{}_{}.cool > '
                        '../MuSC_HiC_files/HiC_not_normalized/not_normalized_{}_{}.txt'. \
                        format(sample, resolution[1], sample, resolution[1]), shell=True)