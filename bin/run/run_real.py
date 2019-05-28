from tau import tau_class

# real data
helion_path = ('/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/'
               'kdawson/hdumasde/Run_programs/igmhub/Void_analysis/')

# Real data -------------------------------------------------------------------
infile = helion_path + 'lya_forest/data/data_delta_transmission_RMplate.fits'
tauobj = tau_class.TauClass(infile)
tauobj.get_data(los_sampling=1)

tauobj.process_box("./voids_real", spacing=5)