# implementing a simple spherical overdensity model
# to identify voids in a 3d map

import numpy as np
from scipy import spatial


def get_voids(data_vec, shape, threshold, cutoff, init_rad=4):
    # build a regular grid
    xx, yy = np.mgrid[0:shape[0], 0:shape[1]]

    # create a spatial tree for fast NN lookup
    # THERE MUST BE A EFFICIENT WAY TO DO THIS CAUSE THE DATA IS
    # ON A REGULAR GRID!
    points = np.c_[xx.ravel(), yy.ravel()]
    tree = spatial.KDTree(points)

    mask = data_vec > threshold
    pts_above_threshold = points[mask]

    print("Found {} points above the threshold value. Growing "
          "spheres around these now!".format(len(pts_above_threshold)))

    # Grow spheres around each of these points and store the radius
    # when the mean density falls below cutoff
    radii = []

    for j, pt in enumerate(pts_above_threshold):
        if not j % 100:
            print("On point:", j)

        done = False
        rad = init_rad

        while not done:
            # select all points within a ball
            nn_pts = tree.query_ball_point(pt, rad)

            # mean density within the given ball
            ball_den = data_vec[nn_pts].mean()

            # check for second threshold
            if ball_den < cutoff:
                done = True
                radii.append(rad)
            else:
                rad += 2

    # create mapping from points to radius
    my_dict = dict(zip(range(len(pts_above_threshold)), radii))

    # join balls that overlap
    n_pts = len(pts_above_threshold)

    for i in range(n_pts):
        for j in range(i + 1, n_pts):
            sep = pts_above_threshold[i] - pts_above_threshold[j]

            mag = sep.dot(sep)

            if mag <= (radii[i] + radii[j]) ** 2:  # they overlap
                if radii[i] > radii[j]:
                    my_dict.pop(j, None)
                else:
                    my_dict.pop(i, None)

    positions = []
    for (key, value) in my_dict.items():
        positions.append([pts_above_threshold[key], value])

    return positions



