#! /usr/env/python

"""
Script to run the analysis on the data
"""

from tau import tau_class, plotting
# import matplotlib.pyplot as plt


filepath = '../../data/delta_data.fits.gz'
tau_obj = tau_class.TauClass(filepath, spacing=10)

kwgs = {'limits': [-1, 1]}
tau_obj.get_data(xlabel="LOGLAM", ylabel="DELTA", nlabel="WEIGHT", los_sampling=3,
                 plotit=False, **kwgs)


# # # # tests
# # # # print(tau_obj.pixel_data.shape)

# # # # data = tau_obj.pixel_data

# # # # plt.figure()
# # # # plt.plot(data[:, 0], data[:, 4], '-k')
# # # # plt.show()

dir_path = "./run_data"
tau_obj.process_box(dir_path)


map_file = "./run_data/map.bin"
cfg_file = "./run_data/void.cfg"
data_file = "./run_data/pixel_data.dat"

plotting.visualizeMap(map_file, cfg_file, [0, 1, 2], data_file=data_file)
