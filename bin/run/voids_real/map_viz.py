from tau import plotting

map_file = 'map.bin'
cfg_file = 'void.cfg'
data_file = 'pixel_data.dat'

# real_file = helion_path + 'picca_run_delta/Backup/Delta/data_delta_transmission_RMplate.fits'
# rd = tau_class.TauClass(real_file)
# all_qso = all_quasars(rd.rot_matrix)

plotting.visualizeMap(map_file, cfg_file, cuts=[10, 11, 12], dir=0, data_file=data_file,
                      width=5, limits=[-1.2, 1.2])
