import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
from mpl_toolkits.mplot3d import Axes3D


def tau_scatter(x, y, z, v, limits=None):
    if len(x) > 10 ** 4:
        a = input("Do you want to proceed?")
        print(a)

        if a == 'n' or a == 'N':
            sys.exit(-1)

    # 3d visualization
    if limits is None:
        limits = [v.min(), v.max()]

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    cbar1 = ax.scatter(x, y, z, c=v, vmin=limits[0], vmax=limits[1])

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    plt.colorbar(cbar1, ax=ax)

    # 1d visualization
    fig, ax = plt.subplots(ncols=3)
    ax[0].hist(x, bins=50)
    ax[1].hist(y, bins=50)
    ax[2].hist(z, bins=50)

    plt.show()


def visualizeMap(map_file, cfg_file, cuts, dir=2, data_file=None, width=2,
                 prefix='test', plot_skewer=True, plot_bin=False, limits=None):
    """
    Plot the mapped data

    Parameters:
    ----------
    map_file : binary file containing mapped data
    cfg_file : config file use in Dachshund
    cuts : [0, 1, 2], indices along the cut direction
    pixel : file containing the pixel data
    sk_idx : skewer index of data-points
    width : width over which to avg the points in h^-1 Mpc

    Returns: None
    """

    # read data file
    if data_file is not None:
        np_data = np.loadtxt(data_file)
        xx, yy, zz, _, vals, sk_idx = np_data.T

    # read map data
    map_data = np.fromfile(map_file)

    # read metadata from the config file
    with open(cfg_file) as f:
        fields = f.readlines()
        cfg = {}
        for field in fields:
            data = field.split('=')
            key = data[0].strip()
            val = float(data[1].strip())
            cfg[key] = val

    # shape of the mapped data grid
    shape = (int(cfg['map_nx']), int(cfg['map_ny']), int(cfg['map_nz']))
    map_data = map_data.reshape(shape)

    # size of the mapped data cube
    lx, ly, lz = float(cfg['lx']), float(cfg['ly']), float(cfg['lz'])

    # bin edges of the mapped grid
    if dir == 0:
        edges = np.linspace(0, lx, shape[0] + 1)
        l0, l1 = ly, lz
        data_ax0, data_ax1, data_fix = yy, zz, xx
    elif dir == 1:
        edges = np.linspace(0, ly, shape[1] + 1)
        l0, l1 = lx, lz
        data_ax0, data_ax1, data_fix = xx, zz, yy
    else:
        edges = np.linspace(0, lz, shape[2] + 1)
        l0, l1 = lx, ly
        data_ax0, data_ax1, data_fix = xx, yy, zz

    centers = (edges[:-1] + edges[1:]) / 2.

    # set limits for colorbar
    if limits is None:
        limits = [map_data.min(), map_data.max()]

    # number of plots to plot
    n_plots = len(cuts)

    fig, ax = plt.subplots(nrows=int(np.ceil(n_plots / 3)), ncols=3)
    ax = np.atleast_2d(ax)

    for ct, ele in enumerate(cuts):
        ii = (ct // 3, ct % 3)

        # plot the mapped data
        plot_data = np.take(map_data, ele, axis=dir)
        cbar = ax[ii].imshow(plot_data.T, origin="lower", vmin=limits[0],
                             vmax=limits[1], extent=(0, l0, 0, l1),
                             cmap=plt.cm.jet)

        # plot the raw data points
        if data_file is not None:
            ixs = ((data_fix > centers[ele] - width) &
                   (data_fix <= centers[ele] + width))

            if plot_skewer:
                agg = np.array([data_ax0[ixs], data_ax1[ixs], vals[ixs]]).T
                df = pd.DataFrame(agg, sk_idx[ixs], ['x', 'y', 'v'])

                # average over the pixels along the fix axis of a given width
                agg_df = df.groupby(df.index).agg(np.mean)

                # overplot the scatter on the mapped plot
                ax[ii].scatter(agg_df['x'], agg_df['y'], c=agg_df['v'],
                               vmin=limits[0], vmax=limits[1], edgecolor='k',
                               cmap=plt.cm.jet)

            if plot_bin:
                # average over data points in a grid same as Weiner map
                # df = binned_statistic_2d(xx[ixs], yy[ixs], vals[ixs],
                #                          bins=[xedges, yedges]).statistic

                # cbar = ax[ii].imshow(df.T, origin="lower", vmin=vmin, vmax=vmax,
                #                      extent=(0, lx, 0, ly), cmap=plt.cm.jet)

                # ax[ii].scatter(np.ravel(map_data[:, :, ele].T), np.ravel(df.T),
                #                alpha=0.1, c='b')
                # ax[ii].plot([-2, 2], [-2, 2])
                pass

        ax[ii].set_title(r"$z = %.1f [h^{-1}\ Mpc]$" % centers[ele])

        # plot labels
        if ii[1] == 0:
            ax[ii].set_ylabel(r'$y\ [h^{-1} Mpc]$')

        if ii[0] == int(np.ceil(n_plots / 3)) - 1:
            ax[ii].set_xlabel(r'$x\ [h^{-1} Mpc]$')

    # Common colorbar for all the subplots
    fig.colorbar(cbar, ax=ax.ravel().tolist(), orientation='horizontal')

    plt.show()
    # plt.savefig(prefix + 'plot.pdf')
