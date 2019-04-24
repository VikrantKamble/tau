import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Circle
from tau import barcode


def main():
    # barcode computation
    file = "../Data/real-original-cl5/real-original2_bars.json"
    shape = (24, 24, 151)

    mybar = barcode.Barcode(file)

    threshold = 1.2
    pos, sizes = mybar.get_voids_at_threshold(threshold)

    # load in the mapped data file
    map_data = np.fromfile("../Data/real-original-cl5/map-real-original2.bin")
    map_data = map_data.reshape(shape)

    # sort and plot the first three
    _co_ord = [(x, y) for y, x in sorted(zip(sizes, pos), reverse=True)]

    # plotting
    fig, ax = plt.subplots(ncols=3, figsize=(12, 4))

    for i, ele in enumerate(_co_ord[:3]):
        # radius of the patch
        radius = ele[1]

        # location of the patch
        ixs = np.unravel_index(ele[0], shape)
        loc = np.mean(ixs, 1).astype(int)

        cbar = ax[i].imshow(map_data[:, :, loc[2]].T, origin='lower',
                            cmap=plt.cm.magma, vmin=-2, vmax=2)

        circ_patch = Circle(loc[:2], radius,
                            facecolor='None', edgecolor='k', linewidth=2,
                            linestyle='dashed')
        ax[i].add_patch(circ_patch)

        ax[i].set_xlabel(r'$x_\perp$')
        ax[i].set_ylabel(r'$x_{||}$')

    fig.colorbar(cbar, ax=ax.ravel().tolist(), orientation='horizontal')
    plt.savefig('voids_bar.pdf')

    # plt.show()


if __name__ == '__main__':
    main()
