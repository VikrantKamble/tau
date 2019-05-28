import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Circle
from tau import barcode


yulong_path = '/uufs/chpc.utah.edu/common/home/u1143816/'
sim_json = yulong_path + 'voidSimulation/v6.0.0_original.json'
map_file = yulong_path + 'voidSimulation/v6.0.0_original_map.bin'
cfg_file = yulong_path + 'voidSimulation/v6.0.0_original_void.cfg'
cfg_data = np.loadtxt(cfg_file, delimiter="=", usecols=1)
cfg_data = cfg_data[4:7].astype(int)
map_data = np.fromfile(map_file)
map_data = map_data.reshape(cfg_data)

mybar = barcode.Barcode(sim_json)
threshold = 0.5
loc = mybar.get_voids_at_threshold(threshold)

# get void catalog
sizes, centers = [], []
for lst in loc:
    sizes.append(len(lst))
    centers.append(np.mean(lst, 0).astype(int))

sizes = np.array(sizes)
centers = np.array(centers)

arg_srt = np.argsort(sizes)
sizes = sizes[arg_srt[::-1]]
centers = centers[arg_srt[::-1]]

# plotting
fig, ax = plt.subplots(ncols=3, figsize=(12, 4))

for i in range(3):
    radius = np.sqrt(sizes[i])
    zs = centers[i][2]

    print(radius, centers[i][:2])

    cbar = ax[i].imshow(map_data[:, :, zs].T, origin='lower',
                        cmap=plt.cm.seismic, vmin=-1.5, vmax=1.5,
                        interpolation='gaussian')

    circ_patch = Circle(centers[i][:2], radius,
                        facecolor='None', edgecolor='k', linewidth=2,
                        linestyle='dashed')
    ax[i].add_patch(circ_patch)

    ax[i].set_xlabel(r'$x_\perp$')
    ax[i].set_ylabel(r'$x_{||}$')

fig.colorbar(cbar, ax=ax.ravel().tolist(), orientation='horizontal')
# plt.savefig('voids_bar.pdf')
plt.show()
