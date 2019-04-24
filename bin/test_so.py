import numpy as np
import matplotlib.pyplot as plt

from tau.so import get_voids
from matplotlib.patches import Circle


def two_voids():
    x, y = np.mgrid[0:150, 0:150]

    # two gaussian voids
    z = 10 * np.exp(- 0.01 * ((x - 0) ** 2 + (y - 0) ** 2)) + \
        10 * np.exp(- 0.001 * ((x - 42) ** 2 + (y - 80) ** 2))

    return z


def grf():
    map_data = np.loadtxt("../data/test_voids2.dat")

    return map_data


def main(type=0):
    if type == 0:
        z = two_voids()
        threshold, cutoff = 9, 8
    elif type == 1:
        z = grf()
        threshold, cutoff = 1.2, 1

    # run the SO finder
    pos = get_voids(z.ravel(), z.shape, threshold, cutoff)

    # plotting
    fig, ax = plt.subplots()
    sp = ax.imshow(z.T, origin='lower', cmap=plt.cm.magma)

    for (pt, radius) in pos:
        circ_patch = Circle(pt, radius,
                            facecolor='None', edgecolor='k', linewidth=2,
                            linestyle='dashed')
        ax.add_patch(circ_patch)

    plt.xlabel(r'$x_\perp$')
    plt.ylabel(r'$x_{||}$')

    plt.colorbar(sp, ax=ax)

    plt.savefig('voids_so.pdf')
    plt.show()


if __name__ == "__main__":
    main(1)
