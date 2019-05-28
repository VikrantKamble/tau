import numpy as np
from astropy.coordinates import cartesian_to_spherical

# make a realization of field in 3D
data = np.loadtxt('../voids_real/pixel_data.dat')
loc = data[:, :3]
ixs = data[:, -1]

skewers = np.unique(ixs)

# put back to the correct locations
shift = np.loadtxt('../voids_real/shift_xyz.dat')
rot = np.loadtxt('../voids_real/rot_matrix.dat')

loc += shift
loc = np.dot(np.linalg.inv(rot), loc.T).T

# convert to spherical coordinates
sphere = cartesian_to_spherical(*loc.T)

r = sphere[0].value
ra = sphere[1].value


# sample a few points out of it as signal

# add noise to the points

