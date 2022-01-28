# Convert HiC to cool files and analyze.

from hic2cool import hic2cool_convert

resolutions = [0, 1000000, 500000,  100000, 50000]# ,  10000, 5000]

for i in resolutions[0:1]:
    # Normalised

    hic2cool_convert('MuSC_HiC_files/HiC_CPB_normalized/OMF_new_CPBnorm_inter_30.hic',
                     'MuSC_HiC_files/HiC_CPB_normalized/normalized_omf_{}'.format(i), resolution=i)
    hic2cool_convert('MuSC_HiC_files/HiC_CPB_normalized/YMF_new_CPBnorm_inter_30.hic',
                     'MuSC_HiC_files/HiC_CPB_normalized/normalized_ymf_{}'.format(i), resolution=i)
    # hic2cool_convert('MuSC_HiC_files/HiC_CPB_normalized/OM_new_CPBnorm_inter_30.hic',
                      #'MuSC_HiC_files/HiC_CPB_normalized/normalized_om_{}'.format(i), resolution=i)
    hic2cool_convert('MuSC_HiC_files/HiC_CPB_normalized/OF_new_CPBnorm_inter_30.hic',
                     'MuSC_HiC_files/HiC_CPB_normalized/normalized_of_{}'.format(i), resolution=i)
    hic2cool_convert('MuSC_HiC_files/HiC_CPB_normalized/YM_new_CPBnorm_inter_30.hic',
                     'MuSC_HiC_files/HiC_CPB_normalized/normalized_ym_{}'.format(i), resolution=i)
    hic2cool_convert('MuSC_HiC_files/HiC_CPB_normalized/YF_new_CPBnorm_inter_30.hic',
                     'MuSC_HiC_files/HiC_CPB_normalized/normalized_yf_{}'.format(i), resolution=i)


    # Unormalized
    hic2cool_convert('MuSC_HiC_files/HiC_not_normalized/OMF_new_inter_30.hic',
                     'MuSC_HiC_files/HiC_not_normalized/normalized_omf_{}'.format(i), resolution=i)
    hic2cool_convert('MuSC_HiC_files/HiC_not_normalized/YMF_new_inter_30.hic',
                     'MuSC_HiC_files/HiC_not_normalized/normalized_ymf_{}'.format(i), resolution=i)
    hic2cool_convert('MuSC_HiC_files/HiC_not_normalized/OM_new_inter_30.hic',
                     'MuSC_HiC_files/HiC_not_normalized/not_normalized_om_{}'.format(i), resolution=i)
    hic2cool_convert('MuSC_HiC_files/HiC_not_normalized/OF_new_inter_30.hic',
                     'MuSC_HiC_files/HiC_not_normalized/not_normalized_of_{}'.format(i), resolution=i)
    hic2cool_convert('MuSC_HiC_files/HiC_not_normalized/YM_new_inter_30.hic',
                     'MuSC_HiC_files/HiC_not_normalized/not_normalized_ym_{}'.format(i), resolution=i)
    hic2cool_convert('MuSC_HiC_files/HiC_not_normalized/YF_new_inter_30.hic',
                     'MuSC_HiC_files/HiC_not_normalized/not_normalized_yf_{}'.format(i), resolution=i)

