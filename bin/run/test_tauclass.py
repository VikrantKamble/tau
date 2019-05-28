from tau import tau_class

# real data
helion_path = ('/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/'
               'kdawson/hdumasde/Run_programs/igmhub/Void_analysis/')

# Real data -------------------------------------------------------------------
infile = helion_path + ('picca_run_delta/Backup/'
                        'Delta/data_delta_transmission_RMplate.fits')
tauobj = tau_class.TauClass(infile)
tauobj.get_data(los_sampling=1)

tauobj.process_box("./voids_real", spacing=10)

# London mocks ----------------------------------------------------------------
infile = helion_path + ('lya_forest/london/v6.0/'
                        'v6.0.0_delta_transmission_RM.fits.gz')
tauobj = tau_class.TauClass(infile)
tauobj.get_data(los_sampling=1)

tauobj.process_box("./voids_london", spacing=10)

# Saclay mocks ----------------------------------------------------------------
infile = helion_path + ('lya_forest/saclay/v4.4/'
                        'v4.4.0_delta_transmission_RM.fits.gz')

tauobj = tau_class.TauClass(infile)
tauobj.get_data(los_sampling=1)

tauobj.process_box("./voids_saclay", spacing=10)