from tau import tau_class
from mayavi import mlab
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import matplotlib.pyplot as plt


def skyplot(tau_obj):
    plt.figure()
    plt.scatter(tau_obj.ra, tau_obj.dec, c=tau_obj.zq)

    # customizations
    plt.xlabel(r'$RA$')
    plt.ylabel(r'$Dec$')

    plt.colorbar()
    plt.title(r'$\# = %d$' % len(tau_obj.zq))

    plt.tight_layout()
    plt.show()

    # add other functionalities here


def plot_3d(tau_obj):
    x, y, z = tau_obj.pixel_data[:, :3].T
    value = tau_obj.pixel_data[:, 4]
    # uncer = tau_obj.pixel_data[:, 3]

    # plt.errorbar(x, value, uncer, fmt='o', capsize=2)
    # plt.show()
    # # plt.hist(value, bins=np.linspace(-2, 2, 100))
    # # plt.show()

    # 3d visualization
    s_min = value.min()
    s_max = value.max()

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    # ax2 = fig.add_subplot(122, projection='3d')

    sp = ax1.scatter(x, y, z, c=value, vmin=-1, vmax=2)
    plt.colorbar(sp, ax=ax1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')

    plt.show()


if __name__ == '__main__':
    file = '../../data/delta-2.fits.gz'
    my_obj = tau_class.TauClass(file)
    # skyplot(my_obj)

    my_obj.get_data(ylabel='DELTA', skewers_perc=1.,
                    los_sampling=50, noise=0.1)

    # my_obj.process_box(dir_path, [-50, 50], [-50, 50], [])
    plot_3d(my_obj)
