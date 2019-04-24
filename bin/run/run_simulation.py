#! /usr/env/python

"""
Script to run the analysis on the data
"""

from tau import tau_class, plotting
# import matplotlib.pyplot as plt


filepath = '../../data/delta_simulation.fits'
tau_obj = tau_class.TauClass(filepath, in_rad=False)

# kwgs = {'limits': [-1, 1]}
# tau_obj.get_data(xlabel="RCOMOV", ylabel="DELTA_T", los_sampling=4,
#                  plotit=False, xtype='rcomov', **kwgs)


# # # # tests
# # # # print(tau_obj.pixel_data.shape)

# # # # data = tau_obj.pixel_data

# # # # plt.figure()
# # # # plt.plot(data[:, 0], data[:, 4], '-k')
# # # # plt.show()

# dir_path = "./run_simulation"
# tau_obj.process_box(dir_path, [-30, 30], [-30, 30], [3600, 3630])


map_file = "./run_simulation/map.bin"
cfg_file = "./run_simulation/void.cfg"
data_file = "./run_simulation/pixel_data.dat"

plotting.visualizeMap(map_file, cfg_file, [1, 2, 3], data_file=data_file)
