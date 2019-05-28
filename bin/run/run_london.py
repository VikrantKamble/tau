from tau import tau_class

# real data
helion_path = ('/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/'
               'kdawson/hdumasde/Run_programs/igmhub/Void_analysis/')

# Saclay mocks ----------------------------------------------------------------
infile = helion_path + ('lya_forest/london/v6.0/'
                        'v6.0.0_delta_transmission_RM.fits.gz')

tauobj = tau_class.TauClass(infile)
tauobj.get_data(los_sampling=1, skewers_perc=0.1)

tauobj.process_box("./voids_london", spacing=5)