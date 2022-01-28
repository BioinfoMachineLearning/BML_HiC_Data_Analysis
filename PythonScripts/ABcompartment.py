import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

eigens = {}

chros = ['chr{}'.format(i) for i in list(range(1, 20))]
# print(chros)
res = 1000000
samples = [('OF_new_CPBnorm_inter_30.hic', 'OF'), ('OM_new_CPBnorm_inter_30.hic', 'OM'),
           ('OMF_new_CPBnorm_inter_30.hic', 'OMF'), ('YF_new_CPBnorm_inter_30.hic', 'YF'),
           ('YM_new_CPBnorm_inter_30.hic', 'YM'), ('YMF_new_CPBnorm_inter_30.hic', 'YMF')]
juicer_file = "../juicer_tools_1.22.01.jar"


#### AB Compartment #######
def plotAB(main_ax,
           ax_num,
           data,
          sample):
    if data[-5] >0:
        data = data*-1
    dataA = data.copy()
    dataB = data.copy()
    dataA[dataA<0]=0
    dataB[dataB>0]=0
    xcoord = list(range(0, len(dataA)))
    main_ax[ax_num].bar(xcoord, dataA, color="red")
    main_ax[ax_num].bar(xcoord,dataB, color="blue")
    main_ax[ax_num].spines['top'].set_visible(False)
    main_ax[ax_num].spines['right'].set_visible(False)
    main_ax[ax_num].set_ylabel(sample)

# for chro in chros:
#     for sample in samples:
#         command = " ".join(["java -jar",
#                 juicer_file,
#                 "eigenvector KR",
#                 "../MuSC_HiC_files/HiC_CPB_normalized/{}".format(sample[0]),
#                 str(chro),
#                 "BP",
#                 str(res),
#                 "Eigens/eigen_"+sample[1]+"_"+str(chro)+"_"+str(res)+".txt"])
#         os.system(command)


for chro in chros:
    fig, ax = plt.subplots(6, 1, figsize=(16, 12))
    for s, sample in enumerate(samples):
        eigen_name = "../Results/Eigens/norm/eigen_" + sample[1] + "_" + str(chro) + "_" + str(res) + ".txt"
        ab_vec = np.loadtxt(eigen_name)
        eigens[sample[1], chro] = ab_vec
        print(len(ab_vec))
        plotAB(ax, s, ab_vec, sample[1])
    ax[-1].set_xlabel("genomic loci")
    ax[0].set_title("A-B Compartment Chromosome --- " + str(chro), fontsize=18)
    cmap = plt.cm.coolwarm
    custom_lines = [Line2D([0], [0], color=cmap(0.0), lw=4),
                    Line2D([0], [0], color=cmap(1.0), lw=4)]
    ax[-1].legend(custom_lines, ['B', 'A'], loc="upper right", fontsize=10)
    plt.savefig('../plots/ABCompartment/{}.jpg'.format(chro))


        # save plot

exit()
###################### PCA ##############################
chro =1
for chro in chros:
    data = np.stack([
                 eigens[samples[0][1],chro],
                 eigens[samples[1][1],chro],
                 eigens[samples[2][1],chro],
                 eigens[samples[3][1],chro],
                 eigens[samples[4][1],chro],
                 eigens[samples[5][1],chro]])

    data = np.nan_to_num(data, nan=0)
    pca = PCA(n_components=2)
    transform = pca.fit_transform(data)

    fig, ax = plt.subplots(figsize=(16, 12))
    ax.scatter(transform[:,0], transform[:,1])
    for i, stri in enumerate(samples):
        ax.text(transform[i,0], transform[i,1], stri[1])
    ax.set_xlabel("PC1 ("+"{:.2f}".format(100*pca.explained_variance_ratio_[0])+"%)")
    ax.set_ylabel("PC2 ("+"{:.2f}".format(100*pca.explained_variance_ratio_[1])+"%)")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title("Chro "+str(chro))
    plt.savefig('../plots/PCA/{}.jpg'.format(chro))



fig, ax = plt.subplots(6,6, figsize=(16, 12))
for chro in chros:
    for i, di in enumerate(samples):
        for j, dj in enumerate(samples):
            vec1 = eigens[di[1], chro]
            vec1 = np.nan_to_num(vec1, nan=0)
            vec2 = eigens[dj[1], chro]
            vec2 = np.nan_to_num(vec2, nan=0)
            p = spearmanr(vec1, vec2)[0]
            if p<0:
                vec2 = vec2*-1
            ax[i,j].scatter(vec1, vec2, s=1)
plt.show()