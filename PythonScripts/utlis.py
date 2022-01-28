import statistics

import pandas as pd
import cooltools
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.ticker import EngFormatter
import cooltools.lib.plotting
import fanc
import fanc.plotting as fancplot
from fanc.architecture.comparisons import hic_pca
import itertools
from matplotlib.gridspec import GridSpec
import bioframe
import matplotlib.patches as patches
from scipy.ndimage import zoom
from scipy import ndimage as ndi


bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x', rotation=45)


def get_positions(nrows, ncols):
    rows = range(nrows)
    cols = range(ncols)
    return (list(itertools.product(rows, cols)))


def balance_matrix(clrs, region, res, save_fig=True, data_type='raw'):
    """
           Method to plot the difference between the balanced ICE normalised matrix
           Parameters:
           arg1 (clr): cooler object
           arg2 (regions): regions to filter on ICE matrix

           Returns:
           void:
           """
    nrows = 2
    ncols = 3
    fig, axs = plt.subplots(
        nrows, ncols, figsize=(20, 16),
    )
    norm = LogNorm(vmax=0.1)
    positions = get_positions(nrows, ncols)
    for i, j in enumerate(clrs):
        pos = positions[i]
        im = axs[pos[0], pos[1]].matshow(
            clrs[j],
            norm=norm,
            cmap='fall',
            extent=(region[2], region[3], region[3], region[2])
        );
        axs[pos[0], pos[1]].set_xlim(region[2], region[3]);
        divider = make_axes_locatable(axs[pos[0], pos[1]])
        cax = divider.append_axes("right", size="5%", pad=0.1)
        axs[pos[0], pos[1]].set_title(j, fontsize=18, fontweight="bold")
        plt.colorbar(im, cax=cax, label='corrected frequencies');

    plt.suptitle(f'{region[0]}:{region[1]}-{region[2]:,}--{region[3]:,} at {res}\n',
                 fontsize=24, ha='center', fontweight="bold")
    plt.tight_layout()
    if save_fig:
        plt.savefig("../plots/Balanced/{}/{}_{}.png".format(data_type, region[0], res))
    else:
        plt.show()


def difference_matrix(clrs, region, res, save_fig=True, data_type='raw'):
    """
           Method to plot the difference between the balanced ICE normalised matrix
           Parameters:
           arg1 (clr): cooler object
           arg2 (regions): regions to filter on ICE matrix

           Returns:
           void:
           """
    nrows = 6
    ncols = 6
    fig, axs = plt.subplots(
        nrows, ncols, figsize=(20, 16),
    )
    positions = get_positions(nrows, ncols)
    for i, j in enumerate(clrs):
        pos = positions[i]
        im = axs[pos[0], pos[1]].matshow(
            clrs[j],
            vmin=-0.01,
            vmax=.01,
            cmap='coolwarm',
            extent=(region[2], region[3], region[3], region[2])
        );
        axs[pos[0], pos[1]].set_xlim(region[2], region[3]);
        divider = make_axes_locatable(axs[pos[0], pos[1]])
        cax = divider.append_axes("right", size="5%", pad=0.1)
        axs[pos[0], pos[1]].set_title(j, fontsize=18, fontweight="bold")
        plt.colorbar(im, cax=cax, label='corrected frequencies');

    plt.suptitle(f'{region[0]}:{region[1]}-{region[2]:,}--{region[3]:,} at {res}\n',
                 fontsize=24, ha='center', fontweight="bold")
    plt.tight_layout()
    if save_fig:
        plt.savefig("../plots/Difference/{}/{}_{}.png".format(data_type, region[0], res))
    else:
        plt.show()


def ims(clrs, res, gene):
    import matplotlib.pyplot as plt

    normalized_contacts = list(clrs.values())

    fig, ax = plt.subplots(4, 4, figsize=(20, 18))
    for s1, sample1 in enumerate(normalized_contacts):
        for s2, sample2 in enumerate(normalized_contacts):
            sampleA = sample1 - np.eye(sample1.shape[0])*sample1
            sampleB = sample2 - np.eye(sample2.shape[0]) * sample2
            im = ax[s1, s2].imshow(sampleA - sampleB, vmin=-0.01, vmax=.01, cmap='RdBu')

    fig.colorbar(im)
    plt.title("Gene: {} -- ")
    plt.savefig("{}_{}.png".format(res, gene))


def get_windows(m):
    scale = [3, 5, 7, 10, 15, 25]
    return [i*m for i in scale]


def coverage(clr):
    cis_coverage, tot_coverage = cooltools.coverage(clr)

    f, ax = plt.subplots(
        figsize=(15, 10),
    )

    norm = LogNorm(vmax=0.1)

    im = ax.matshow(
        clr.matrix()[:],
        norm=norm,
        cmap='fall'
    );
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax, label='corrected frequencies');
    ax.set_title('full matrix')
    ax.xaxis.set_visible(True)

    ax1 = divider.append_axes("bottom", size="25%", pad=0.1, sharex=ax)
    weights = clr.bins()[:]['weight'].values
    ax1.plot(cis_coverage, label='cis')
    ax1.plot(tot_coverage, label='total')
    ax1.set_xlim([0, len(clr.bins()[:])])
    ax1.set_ylabel('Coverage')
    ax1.legend()
    # ax1.set_xticks([])

    ax2 = divider.append_axes("bottom", size="25%", pad=0.1, sharex=ax)
    ax2.plot(cis_coverage / tot_coverage)
    ax2.set_xlim([0, len(clr.bins()[:])])
    ax2.set_ylabel('coverage ratio')
    plt.savefig("coverage.png")


def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start + resolution * i for i in range(len(matrix_c) + 1)]
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im


def convert_bytes(num):
    """
    this function will convert bytes to MB.... GB... etc
    """
    step_unit = 1000.0  # 1024 bad the size

    for x in ['b', 'kb', 'mb', 'gb', 'tb']:
        if num < step_unit:
            return "%3.0f %s" % (num, x)
        num /= step_unit


def plot_rep(data, chrom, res, save_fig):
    fig, ax = plt.subplots(figsize=(20, 16))

    keys, values = zip(*data.items())

    s = pd.Series(
        values,
        index=keys
    )

    # Set descriptions:
    plt.title("Stratum-adjusted correlation coefficient plot between samples -- {} -- {}".format(chrom, convert_bytes(res[1])))
    plt.ylabel('SCC score')
    plt.xlabel('Samples')

    minimun = min(data.values()) - statistics.stdev(data.values())
    maximum = max(data.values()) + statistics.stdev(data.values())

    ax.set_ylim(minimun, maximum)

    # Plot the data:
    my_colors = list('bbrrrbrrrrrrkkk')

    s.plot(
        kind='bar',
        color=my_colors,
    )

    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='b', lw=4),
                    Line2D([0], [0], color='k', lw=4)]
    ax.legend(custom_lines, ['OLD vs YOUNG', 'OLD vs OLD', 'YOUNG vs YOUNG'], loc="upper right", fontsize=10)

    if save_fig:
        plt.savefig("../plots/REP/{}_{}.png".format(chrom, convert_bytes(res[1])))
    else:
        plt.show()

def ab_plot(eigenvector_track, ax, sample):

    data = eigenvector_track['E1'].values


    dataA = data.copy()
    dataB = data.copy()

    dataA[dataA < 0] = 0
    dataB[dataB > 0] = 0

    xcoord = list(range(0, len(dataA)))

    ax.bar(xcoord, dataA, color="red")
    ax.bar(xcoord, dataB, color="blue")
    ax.plot(eigenvector_track['E1'].values, label='E1')
    ax.plot([0, len(dataA)], [0, 0], 'k', lw=0.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel(sample)


def ab_plot_overlay(eigenvector_track, ax, sample):

    colors = {
        0: 'red',
        1: 'orange',
        2: 'yellow',
        3: 'blue',
        4: 'indigo',
        5: 'violet'
    }
    data = eigenvector_track['E1'].values


    dataA = data.copy()
    dataB = data.copy()

    dataA[dataA < 0] = 0
    dataB[dataB > 0] = 0

    xcoord = list(range(0, len(dataA)))

    #ax.bar(xcoord, dataA, color="red")
    #ax.bar(xcoord, dataB, color="blue")
    ax.plot(eigenvector_track['E1'].values, label=sample[1], color=colors[sample[0]])
    ax.plot([0, len(dataA)], [0, 0], 'k', lw=0.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)




def insulation_plot(data, save_fig=True):
    nrows = 6
    ncols = 1
    fig, axs = plt.subplots(
        nrows, ncols, figsize=(20, 16),
    )
    norm = LogNorm(vmax=0.1, vmin=0.001)
    positions = get_positions(nrows, ncols)

    for i, j in enumerate(data):
        resolution = int(j.split('_')[1])

        window = get_windows(resolution)
        region = data[j][1]

        route = (region[1], region[2], region[2]+ 80*window[0])
        # clr = data[j][0].matrix(balance=True).fetch(route)
        clr = data[j][0].matrix(balance=True).fetch(region[1:])

        sample = j.split('_')[0]
        pos = positions[i]
        im = pcolormesh_45deg(axs[pos[0]], clr, start=region[2], resolution=resolution, norm=norm, cmap='fall')
        # axs[pos[0]].set_aspect(0.5)
        axs[pos[0]].set_ylim(0, window[0] * 10)
        format_ticks(axs[pos[0]], rotate=False)
        axs[pos[0]].xaxis.set_visible(False)
        axs[pos[0]].set_title(sample, fontsize=14)

        divider = make_axes_locatable(axs[pos[0]])
        cax = divider.append_axes("right", size="1%", pad=0.1)
        plt.colorbar(im, cax=cax)

        INSULATION_PATH = '../Results/Insulation/insulation_score_{}_{}.csv'
        file_name = INSULATION_PATH.format(sample, resolution)
        df = pd.read_csv(file_name)

        insul_region = bioframe.select(df, region[1:])
        ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=axs[pos[0]])
        ins_ax.set_prop_cycle(plt.cycler("color", plt.cm.plasma(np.linspace(0, 1, 5))))
        ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
                    insul_region[f'log2_insulation_score_{window[0]}'],
                    label='{}'.format(convert_bytes(window[0])))
        ins_ax.legend(bbox_to_anchor=(1.100, 1.05), loc='upper right')

        format_ticks(ins_ax, y=False, rotate=False)
        axs[pos[0]].set_xlim(region[2], region[3])

        for res in window[1:]:
            ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
                        insul_region[f'log2_insulation_score_{res}'],
                        label=f'{convert_bytes(res)}')
        ins_ax.legend(bbox_to_anchor=(1.06, 0.998), loc='upper right')

        plt.suptitle(f'TADs -- Gene: {region[0]} -- chromosome: {region[1]} -- Region {region[2]:,}-{region[3]:,} -- Resolution: {convert_bytes(resolution)}\n',
                     fontsize=18, ha='center', fontweight="bold")
        plt.tight_layout()

    if save_fig:
        plt.savefig("../plots/Insulation/{}_{}.png".format(region[0], resolution))
    else:
        plt.show()


def boundaries(data, save_fig=True):
    nrows = 6
    ncols = 1
    fig, axs = plt.subplots(
        nrows, ncols, figsize=(20, 16),
    )
    norm = LogNorm(vmax=0.1, vmin=0.001)
    positions = get_positions(nrows, ncols)

    for i, j in enumerate(data):
        resolution = int(j.split('_')[1])
        window = get_windows(resolution)

        region = data[j][1]
        route = (region[1], region[2], region[2] + 80 * window[0])
        # clr = data[j][0].matrix(balance=True).fetch(route)
        clr = data[j][0].matrix(balance=True).fetch(region[1:])

        sample = j.split('_')[0]
        pos = positions[i]
        im = pcolormesh_45deg(axs[pos[0]], clr, start=region[2], resolution=resolution, norm=norm, cmap='fall')
        # axs[pos[0]].set_aspect(0.5)
        axs[pos[0]].set_ylim(0, window[0] * 10)
        format_ticks(axs[pos[0]], rotate=False)
        axs[pos[0]].xaxis.set_visible(False)
        axs[pos[0]].set_title(sample, fontsize=14)

        divider = make_axes_locatable(axs[pos[0]])
        cax = divider.append_axes("right", size="1%", pad=0.1)
        plt.colorbar(im, cax=cax)

        INSULATION_PATH = '../Results/Insulation/insulation_score_{}_{}.csv'
        file_name = INSULATION_PATH.format(sample, resolution)
        df = pd.read_csv(file_name)
        # df = find_boundaries(df)

        insul_region = bioframe.select(df, region[1:])
        ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=axs[pos[0]])
        # ins_ax.set_prop_cycle(plt.cycler("color", plt.cm.plasma(np.linspace(0, 1, 5))))
        ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
                    insul_region[f'log2_insulation_score_{window[0]}'],
                    label='{}'.format(convert_bytes(window[0])))
        ins_ax.legend(bbox_to_anchor=(1.100, 1.05), loc='upper right')

        boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{window[0]}'])]
        weak_boundaries = boundaries[~boundaries[f'is_boundary_{window[0]}']]
        strong_boundaries = boundaries[boundaries[f'is_boundary_{window[0]}']]

        ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
                       weak_boundaries[f'log2_insulation_score_{window[0]}'], label='Weak boundaries')
        ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
                       strong_boundaries[f'log2_insulation_score_{window[0]}'], label='Strong boundaries')

        ins_ax.legend(bbox_to_anchor=(1.125, 1.05), loc='upper right')

        format_ticks(ins_ax, y=False, rotate=False)
        axs[pos[0]].set_xlim(region[2], region[3])

        plt.suptitle(
            f'Annotated valleys: {region[0]} -- chromosome: {region[1]} -- Region {region[2]:,}-{region[3]:,} -- Resolution: {convert_bytes(resolution)}\n',
            fontsize=18, ha='center', fontweight="bold")
        plt.tight_layout()

    if save_fig:
        plt.savefig("../plots/Boundary/{}_{}.png".format(region[0], resolution))
    else:
        plt.show()


def clipped_zoom(img, zoom_factor, **kwargs):
    h, w = img.shape[:2]
    zoom_tuple = (zoom_factor,) * 2 + (1,) * (img.ndim - 2)
    if zoom_factor < 1:
        zh = int(np.round(h * zoom_factor))
        zw = int(np.round(w * zoom_factor))
        top = (h - zh) // 2
        left = (w - zw) // 2
        out = np.zeros_like(img)
        out[top:top+zh, left:left+zw] = zoom(img, zoom_tuple, **kwargs)
    elif zoom_factor > 1:
        zh = int(np.round(h / zoom_factor))
        zw = int(np.round(w / zoom_factor))
        top = (h - zh) // 2
        left = (w - zw) // 2
        out = zoom(img[top:top+zh, left:left+zw], zoom_tuple, **kwargs)
        trim_top = ((out.shape[0] - h) // 2)
        trim_left = ((out.shape[1] - w) // 2)
        out = out[trim_top:trim_top+h, trim_left:trim_left+w]
    else:
        out = img
    return out


def highlight_features(dataframe, region, color, a, axes):
    try:
        features = dataframe.loc[region].values.tolist()
        if type(features[0]) == int:
            _, x_min, x_max, y_min, y_max = features
            rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, linewidth=1.2,
                                     edgecolor=color, facecolor='none')
            axes[a].add_patch(rect)
        else:
            for f in features:
                _, x_min, x_max, y_min, y_max = f
                rect = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, linewidth=1.2,
                                         edgecolor=color, facecolor='none')
                axes[a].add_patch(rect)

    except KeyError:
        next


def plot_signal_noise_ratio(sub_sim, sim_field, zsim_thr, region_id, data1, data2, chrom,
                            pairs, gained=None, lost=None, save_fig=False, add_features=False):
    regions = data1[0]
    similarities = data1[1]

    pair_A_sub = data2[0]
    pair_B_sub = data2[1]
    l2fcm = data2[2]

    # plt.rcParams['text.usetex'] = True
    # plt.style.use(['dark_background'])
    # plt.style.use('classic')
    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    SMALL_SIZE = 13
    MEDIUM_SIZE = 15
    BIGGER_SIZE = 20

    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fig = plt.figure(figsize=(20, 16), constrained_layout=True, dpi=120)
    if add_features:
        gs = GridSpec(3, 3, figure=fig)
        ax1 = fig.add_subplot(gs[0, :])
        ax2 = fig.add_subplot(gs[1, 0:1])
        ax3 = fig.add_subplot(gs[1, 1:2])
        ax4 = fig.add_subplot(gs[1, 2:3])
        ax5 = fig.add_subplot(gs[2, 0:1])
        ax6 = fig.add_subplot(gs[2, 1:2])
        ax7 = fig.add_subplot(gs[2, 2:3])
    else:
        gs = GridSpec(2, 3, figure=fig)
        ax1 = fig.add_subplot(gs[0, :])
        ax2 = fig.add_subplot(gs[1, 0:1])
        ax3 = fig.add_subplot(gs[1, 1:2])
        ax4 = fig.add_subplot(gs[1, 2:3])


    all_X = regions.loc[similarities.index, 1:2].mean(axis=1).values / 10 ** 6
    X = regions.loc[sub_sim.index, 1:2].mean(axis=1).values / 10 ** 6
    S = sub_sim[sim_field]
    SN = sub_sim["SN"]


    ax1.plot(all_X, similarities[sim_field], ":", alpha=0.4)
    ax1.hlines(zsim_thr, 0, max(all_X), linestyle=":", color="red")
    ax1.scatter(all_X, similarities[sim_field], facecolors='none', edgecolors='grey', alpha=0.1, s=20)
    ax1.scatter(X, S, c=SN, marker='.')
    ax1.set_ylabel(sim_field.replace("_", "-"))
    ax1.set_xlabel("window midpoint [Mb]")

    sc = fig.axes[0].collections[0]

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='1%', pad=0.3)
    fig.colorbar(sc, cax=cax, orientation='vertical', label='Signal/Noise')

    ax2.set_title('{}'.format(pairs[0].upper()))
    ax3.set_title('{}'.format(pairs[1].upper()))
    ax4.set_title('$log_2$({} / {})'.format(pairs[0].upper(), pairs[1].upper()))

    m1 = ax2.imshow(pair_A_sub, norm=LogNorm(), cmap='germany')
    m2 = ax3.imshow(pair_B_sub, norm=LogNorm(), cmap='germany')
    m3 = ax4.imshow(l2fcm, cmap='seismic', vmax=5, vmin=-5)

    for m, ax in zip([m1, m2, m3], [ax2, ax3, ax4]):
        ax.axis('off')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('bottom', size='5%', pad=0.05)
        fig.colorbar(m, cax=cax, orientation='horizontal')

    if add_features:
        ax5.set_title('{}'.format(pairs[0].upper()))
        ax6.set_title('{}'.format(pairs[1].upper()))
        ax7.set_title('$log_2$({} / {})'.format(pairs[0].upper(), pairs[1].upper()))
        zml2 = clipped_zoom(l2fcm, 0.7)
        rot_l2 = ndi.rotate(zml2, 45, reshape=False)

        # clipped zoom and rotate patient and control and keep only half-matrix
        zm1 = clipped_zoom(pair_A_sub, 0.7)
        rot_pair_A= ndi.rotate(zm1, 45, reshape=False)

        zm2 = clipped_zoom(pair_B_sub, 0.7)
        rot_pair_B = ndi.rotate(zm2, 45, reshape=False)

        middle = int(np.shape(rot_pair_B)[1] / 2.)

        m1 = ax5.imshow(rot_pair_A[:middle, :], vmin=0, vmax=0.03, cmap='germany')
        m2 = ax6.imshow(rot_pair_B[:middle, :], vmin=0, vmax=0.03, cmap='germany')

        axes = [ax5, ax6, ax7]

        # # per region check if identified features, to highlight
        highlight_features(gained, region_id, 'crimson', 0, axes)
        highlight_features(lost, region_id, 'royalblue', 1, axes)

        m3 = ax7.imshow(rot_l2[:middle, :], cmap='seismic', vmax=5, vmin=-5)

        for m, ax in zip([m1, m2, m3], [ax5, ax6, ax7]):
            ax.axis('off')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('bottom', size='5%', pad=0.05)
            fig.colorbar(m, cax=cax, orientation='horizontal')

    fig.suptitle("Similarity Index & Interesting Regions {} {}--{}".format(chrom, pairs[0], pairs[1]), fontsize=16)
    plt.tight_layout()

    if save_fig:
        plt.savefig("../plots/Chess/{}_{}_{}.png".format(chrom, pairs[0], pairs[1]))
    else:
        plt.show()

    #



