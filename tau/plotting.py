import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


from tau import barcode


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


def visualizeMap(map_file, cfg_file, cuts, dir=2, json_file=None,
                 data_file=None, all_qso=None, threshold=0.5,
                 width=2, prefix='test', plot_skewer=True, limits=None):
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

    # read file that contains the pixel data
    if data_file is not None:
        if data_file[-3:] == 'npy':
            np_data = np.load(data_file)
        else:
            np_data = np.loadtxt(data_file)
        xx, yy, zz, _, vals, sk_idx = np_data.T

        pos = np.array([xx, yy, zz])

    # read metadata from the config file
    cfg_data = np.loadtxt(cfg_file, delimiter="=", usecols=1)
    lp = cfg_data[:3]
    ls = cfg_data[4:7].astype(int)

    # read map data
    map_data = np.fromfile(map_file)
    map_data = map_data.reshape(ls)

    # read json data and make void catalog
    if json_file is not None:
        mybar = barcode.Barcode(json_file)
        vs, vc = mybar.get_voids_at_threshold(threshold)

    # the center value of the slices for the mapped data
    dx = lp[dir] / ls[dir]
    centers = (np.arange(ls[dir]) + 0.5) * ls[dir]

    lx, ly, lz = lp
    if dir == 0:
        l0, l1 = ly, lz
        if data_file is not None:
            data_ax0, data_ax1, data_fix = yy, zz, xx
    elif dir == 1:
        l0, l1 = lx, lz
        if data_file is not None:
            data_ax0, data_ax1, data_fix = xx, zz, yy
    else:
        l0, l1 = lx, ly
        if data_file is not None:
            data_ax0, data_ax1, data_fix = xx, yy, zz

    n_plots = len(cuts)
    ll = [0, 1, 2]
    ll.remove(dir)

    fig, ax = plt.subplots(ncols=int(np.ceil(n_plots / 3)), nrows=3)
    ax = np.atleast_2d(ax)

    # set limits for colorbar
    if limits is None:
        limits = [map_data.min(), map_data.max()]

    for ct, ele in enumerate(cuts):
        ii = (ct // 3, ct % 3)

        # plot the mapped data
        plot_data = np.take(map_data, ele, axis=dir)
        cbar = ax[ii].imshow(plot_data, origin="lower", vmin=limits[0],
                             vmax=limits[1], extent=(0, l1, 0, l0),
                             cmap=plt.cm.seismic, interpolation='gaussian')

        # plot the voids if json_file
        if json_file is not None:
            # indices where the voids reside
            vcs = vc[:, dir].astype(int)
            v_ixs = np.where(vcs == ele)[0]

            for ixs in v_ixs:
                # radius = (3 * vs[ixs] / (4 * np.pi))
                radius = 10
                cs = ((vc[ixs][ll[1]] + 0.5) * dx, (vc[ixs][ll[0]] + 0.5) * dx)
                circ = Circle(cs, radius,
                              facecolor='None', edgecolor='k', linewidth=2,
                              linestyle='dashed')
                ax[ii].add_patch(circ)

        # plot the raw data points
        if data_file is not None:
            ixs = ((data_fix > centers[ele] - width) &
                   (data_fix <= centers[ele] + width))

            agg = np.array([data_ax0[ixs], data_ax1[ixs], vals[ixs]]).T
            df = pd.DataFrame(agg, sk_idx[ixs], ['x', 'y', 'v'])

            # average over the pixels along the fix axis of a given width
            agg_df = df.groupby(df.index).agg(np.mean)

            # overplot the scatter on the mapped plot
            if dir == 2:
                ax[ii].scatter(agg_df['x'], agg_df['y'], c=agg_df['v'],
                               vmin=limits[0], vmax=limits[1], edgecolor='k',
                               cmap=plt.cm.coolwarm)
            else:
                ax[ii].scatter(data_ax0[ixs][::10], data_ax1[ixs][::10], c=vals[ixs][::10],
                               vmin=limits[0], vmax=limits[1], edgecolor='k',
                               cmap=plt.cm.coolwarm)

        if all_qso is not None:
            rpar = np.take(all_qso, dir, axis=1)
            ixs = np.where((rpar > centers[ele] - width) &
                           (rpar <= centers[ele] + width))[0]

            # overplot the scatter on the mapped plot
            rperp = [0, 1, 2]
            rperp.remove(dir)

            if len(ixs) > 0:
                ax[ii].scatter(*np.take(all_qso, rperp, axis=1)[ixs], c='k')

        ax[ii].set_title(r"$z = %.1f [h^{-1}\ {\rm Mpc}]$" % centers[ele])

        # plot labels
        if ii[1] == 1:
            ax[ii].set_ylabel(r'$y\ [h^{-1} {\rm Mpc}]$')

        if ii[0] == int(np.ceil(n_plots / 3)) - 1:
            ax[ii].set_xlabel(r'$x\ [h^{-1} {\rm Mpc}]$')

    # Common colorbar for all the subplots
    fig.colorbar(cbar, ax=ax.ravel().tolist(), orientation='vertical')

    plt.show()
    # plt.savefig(prefix + 'plot.pdf')
