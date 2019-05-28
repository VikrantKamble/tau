from tau import tau_class

# real data
helion_path = ('/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/'
               'kdawson/hdumasde/Run_programs/igmhub/Void_analysis/')

# Saclay mocks ----------------------------------------------------------------

infile = helion_path + ('lya_forest/saclay/v4.4/'
                        'v4.4.0_delta_transmission_RM.fits.gz')

tauobj = tau_class.TauClass(infile)
tauobj.get_data(los_sampling=3)

tauobj.process_box("./voids_saclay", spacing=2, xcuts=[0, 50], ycuts=[0, 50], zcuts=[3800, 3850])