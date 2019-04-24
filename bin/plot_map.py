from mayavi import mlab
from tvtk.util import ctf
from matplotlib.pyplot import cm

import numpy as np

# load data
file = "/Users/vikrant/Work/tau/data/unreal.bin"

map_data = np.fromfile(file)
map_data = map_data.reshape((150, 150, 150))


# mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(map_data),
#                                  plane_orientation='z_axes',
#                                  slice_index=10,
#                                  colormap='plasma'
#                                  )

volume = mlab.pipeline.volume(mlab.pipeline.scalar_field(map_data))

mlab.outline()

mlab.show()
