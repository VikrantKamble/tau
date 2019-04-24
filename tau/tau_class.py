#!/usr/bin/env python
import sys
import os
import subprocess
import numpy as np

from tau.data_wrangle import celestial_rot_matrix
from .plotting import tau_scatter
from astropy.io import fits
from astropy.cosmology import Planck15


class TauClass:
    def __init__(self, filepath=None, spacing=4., preprocess=True, in_rad=True):
        if filepath is None:
            filepath = sys.argv[1]
        self.data_file = fits.open(filepath)

        self.q_loc = np.zeros((len(self), 3))

        self.rot_matrix = None
        self.pixel_data = None

        self.spacing = spacing
        self.in_rad = in_rad

        # get locations on the celestial sphere
        if preprocess:
            self.get_locations()
            self.get_matrix()

    def __repr__(self):
        repr_str = "***** Interpolate Object ***** \n" + \
            "Number of skewers are {}".format(len(self.data_file) - 1)
        return repr_str

    def __len__(self):
        return len(self.data_file) - 1

    def get_locations(self):
        """ Get the angular positions and the redshift
        of the quasars -> q_loc = [ra, dec, z_q]
        """
        for i in range(len(self)):
            hdr = self.data_file[i+1].header

            self.q_loc[i] = [hdr['RA'], hdr['DEC'], hdr['Z']]

    def get_matrix(self, center=None):
        """ Rotate the celestial sphere such that the polar axis
        points towards the central ra and dec of the data patch
        """
        print("****** CALCULATING ROTATION MATRIX ******** \n")

        if center is None:
            print("+++ Using calculated RA and DEC mean +++")
            center = [self.q_loc[:, 0].mean(), self.q_loc[:, 1].mean()]

        self.rot_matrix = celestial_rot_matrix(*center, is_rad=self.in_rad)

    def get_data(self, xlabel='LOGLAM', ylabel='DELTA_T', nlabel=None,
                 skewers_num=None, skewers_perc=1., los_sampling=3,
                 xtype='loglam', plotit=True, **kwargs):
        """ Load skewer data and convert to cartesian cordinates

        Parameters:
        -----------
        xlabel : str
            the label of the x-data in the Table HDUs
        ylabel : str
            the label of the field data in the Table HDUs
        skewers_num : integer
            Number of skewers to select
        skewers_perc : float [0, 1]
            Percentage of skewers to select
        log_sampling : integer
            Selecting every n'th pixel along the Line-Of-Sight
        noise : measurement noise

        Returns:
        --------
            None
        """

        # subsampling the skewers
        if skewers_num is None:
            skewers_num = int(len(self) * skewers_perc)

        ixs = np.random.choice(np.arange(1, len(self) + 1), size=skewers_num,
                               replace=False)

        # variables that hold information for Dachshund
        x_sig, y_sig, z_sig = [], [], []
        n_sig, v_sig, sx_sig = [], [], []

        # check if given xlabel and ylabel exist in HDUs
        foo = self.data_file[1].data
        col_names = foo.dtype.names

        if xlabel not in col_names:
            raise AttributeError('X-data not found')
        if ylabel not in col_names:
            raise AttributeError('Y-data not found')

        print('****** READING ****** \n')
        for i in ixs:
            hdu = self.data_file[i].data
            xdata = hdu[xlabel]

            # thin the data along the Line-of-Sight
            xdata = xdata[::los_sampling]

            if xtype == 'loglam':
                # convert from wavelength to redshift
                lam = 10 ** xdata
                z_abs = lam / 1215.67 - 1

                # convert from redshift to comoving distance
                rcomov = Planck15.comoving_distance(z_abs).value * Planck15.h
            else:
                rcomov = xdata

            # convert spherical to cartesian cordinates
            ra, dec = self.q_loc[i-1, :2]

            if not self.in_rad:
                ra, dec = np.deg2rad([ra, dec])

            x = np.sin(dec) * np.cos(ra) * rcomov
            y = np.sin(dec) * np.sin(ra) * rcomov
            z = np.cos(dec) * rcomov

            # get value of the field
            value = hdu[ylabel][::los_sampling]

            # get meaasurement uncertainity
            if nlabel is None:
                noise = 0.1 * np.ones_like(x)
            else:
                noise = 1.0 / np.sqrt(hdu[nlabel][::los_sampling])

            # skewer index
            sx = i * np.ones_like(x)

            # concatenating together
            x_sig = np.hstack((x_sig, x))
            y_sig = np.hstack((y_sig, y))
            z_sig = np.hstack((z_sig, z))
            n_sig = np.hstack((n_sig, noise))
            v_sig = np.hstack((v_sig, value))
            sx_sig = np.hstack((sx_sig, sx))

        # rotating the points
        x_sig, y_sig, z_sig = np.dot(self.rot_matrix, np.array([x_sig, y_sig, z_sig]))

        # final pixel data -  xloc, yloc, zloc, noise, field, skewer_number
        self.pixel_data = np.vstack([x_sig, y_sig, z_sig, n_sig, v_sig, sx_sig]).T
        print('Read {} points from the file'.format(len(x_sig)))

        # diagnostic plot
        if plotit:
            tau_scatter(x_sig, y_sig, z_sig, v_sig, **kwargs)

    def process_box(self, dir_path, xcuts=None, ycuts=None, zcuts=None):
        """ Run Dachsund on the data with given cuts along the different
        directions

        Parameters:
        ----------
            dir_path : the directory to dump the results in
            xcuts, ycuts, zcuts : cuts along x, y and z axes respectively
        """
        xx, yy, zz = self.pixel_data[:, :3].T

        if xcuts is None:
            xcuts = [xx.min(), xx.max()]
        if ycuts is None:
            ycuts = [yy.min(), yy.max()]
        if zcuts is None:
            zcuts = [zz.min(), zz.max()]

        # select data in a cube
        ixs = (xx >= xcuts[0]) & (xx < xcuts[1]) &\
              (yy >= ycuts[0]) & (yy < ycuts[1]) &\
              (zz >= zcuts[0]) & (zz < zcuts[1])
        print("Number of points", len(ixs))

        # data to send to Dachsund - remove the skewer index column
        chunk = self.pixel_data[ixs]

        # shift the cordinates to the edge of the cube
        chunk[:, 0] -= chunk[:, 0].min()
        chunk[:, 1] -= chunk[:, 1].min()
        chunk[:, 2] -= chunk[:, 2].min()

        # make a directory and carry all operations there
        if os.path.isdir(dir_path):
            a = input("directory already exists!, do you want to proceed")
            if a == 'n' or 'N':
                sys.exit(-1)
            else:
                print('Overwriting.....')
        else:
            os.makedirs(dir_path)

        cwd = os.getcwd()
        try:
            os.chdir(dir_path)

            # save file needed for reconstruction
            chunk[:, :-1].tofile('pixel_data.bin')

            # save pixel_data file
            np.savetxt('pixel_data.dat', chunk)

            # run the Weiner reconstruction on the chunk
            self.dachshund(chunk)
        except Exception:
            raise
        finally:
            os.chdir(cwd)

    def dachshund(self, chunk):
        xx, yy, zz = chunk[:, :3].T

        # length along each direction - starts at 0
        self.lx = xx.max()
        self.ly = yy.max()
        self.lz = zz.max()

        # Number of parts in which to divide the length along each direction
        self.npx_x = int(self.lx // self.spacing)
        self.npx_y = int(self.ly // self.spacing)
        self.npx_z = int(self.lz // self.spacing)
        print("Grid dimensions are {} {} {}".format(self.npx_x,
                                                    self.npx_y, self.npx_z))

        # write config file
        with open('void.cfg', 'w') as cf:
            cf.write("lx = %f\n" % self.lx)
            cf.write("ly = %f\n" % self.ly)
            cf.write("lz = %f\n" % self.lz)
            cf.write("num_pixels = %i\n" % len(xx))
            cf.write("map_nx = %i\n" % self.npx_x)
            cf.write("map_ny = %i\n" % self.npx_y)
            cf.write("map_nz = %i\n" % self.npx_z)
            cf.write("corr_var_s = 0.23\n")
            cf.write("corr_l_perp = 5\n")
            cf.write("corr_l_para = 5\n")
            cf.write("pcg_tol = 1.0e-5\n")
            cf.write("pcg_max_iter = 1000\n")

        # run the program on the data
        executable = "../../../dachshund/dachshund.exe"
        message1 = subprocess.run([executable, 'void.cfg'],
                                  stdout=sys.stdout)
        print(message1)
