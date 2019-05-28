#!/usr/bin/env python
import sys
import os
import subprocess
import numpy as np

from tau.data_wrangle import celestial_rot_matrix
from tau.plotting import tau_scatter
from astropy.io import fits


class TauClass:
    def __init__(self, filepath=None, preprocess=True):
        if filepath is None:
            filepath = sys.argv[1]
        self.data_file = fits.open(filepath)

        self.q_loc = np.zeros((len(self), 3))

        self.rot_matrix = None
        self.pixel_data = None

        # get locations on the celestial sphere
        if preprocess:
            self.get_locations()
            self.get_matrix()
            self.set_config()

        print("The file contains a total of {} skewers".format(len(self)))

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

    def get_matrix(self, center=[213.704, 53.083]):
        """ Rotate the celestial sphere such that the polar axis
        points towards the central ra and dec of the data patch
        """
        print("****** CALCULATING ROTATION MATRIX ******** \n")

        if center is None:
            print("+++ Using calculated RA and DEC mean +++")
            center = [self.q_loc[:, 0].mean(), self.q_loc[:, 1].mean()]

        self.rot_matrix = celestial_rot_matrix(*center, is_rad=False)

    def set_config(self, var_s=0.23, l_para=5, l_perp=5):
        """
        Set the parameters of the config file -
        Currently only supports ExpSquared kernel
        """
        self.var_s = var_s
        self.l_para = l_para
        self.l_perp = l_perp

    def get_data(self, skewers_num=None, skewers_perc=1., los_sampling=3,
                 plotit=False, **kwargs):
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
        """

        # subsampling the skewers
        if skewers_num is None:
            skewers_num = int(len(self) * skewers_perc)

        ixs = np.random.choice(np.arange(1, len(self) + 1), size=skewers_num,
                               replace=False)

        # variables that hold information for Dachshund
        x_sig, y_sig, z_sig = [], [], []
        n_sig, v_sig, sx_sig = [], [], []

        print('****** READING {} SKEWERS ****** \n'.format(len(ixs)))
        for tag, i in enumerate(ixs):
            tb = self.data_file[i].data
            rcomov = tb['RCOMOV'][::los_sampling]

            # convert spherical to cartesian cordinates
            ra, dec = np.deg2rad(self.q_loc[i-1, :2])

            x = np.cos(dec) * np.cos(ra) * rcomov
            y = np.cos(dec) * np.sin(ra) * rcomov
            z = np.sin(dec) * rcomov

            # get value of the field
            value = tb['DELTA_T'][::los_sampling]

            # get meaasurement uncertainity
            if 'WEIGHTS' in tb.dtype.names:
                noise = 1.0 / np.sqrt(tb['WEIGHTS'][::los_sampling])
            else:
                noise = 0.1 * np.ones_like(x)

            # skewer index
            sx = tag * np.ones_like(z)

            # concatenating together
            x_sig = np.hstack((x_sig, x))
            y_sig = np.hstack((y_sig, y))
            z_sig = np.hstack((z_sig, z))
            n_sig = np.hstack((n_sig, noise))
            v_sig = np.hstack((v_sig, value))
            sx_sig = np.hstack((sx_sig, sx))

        print("Number of skewers that actually contain"
              "data are {}".format(len(np.unique(sx_sig))))

        # rotating the points
        x_sig, y_sig, z_sig = np.dot(self.rot_matrix,
                                     np.array([x_sig, y_sig, z_sig]))

        # final pixel data -  xloc, yloc, zloc, noise, field, skewer_number
        self.pixel_data = np.vstack([x_sig, y_sig, z_sig, n_sig, v_sig, sx_sig]).T
        print('The file contains a total of {} data-points'.format(len(x_sig)))

        # diagnostic plot
        if plotit:
            tau_scatter(x_sig, y_sig, z_sig, v_sig, **kwargs)

    def process_box(self, dir_path, xcuts=None, ycuts=None, zcuts=None,
                    spacing=10, reconstruct=False):
        """ Run Dachshund on the data with given cuts along the different
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
        print("Points used for Weiner reconstruction: ", ixs.sum())

        # data to send to Dachshund
        chunk = self.pixel_data[ixs]

        # shift the cordinates to the edge of the cube
        self.xshift = chunk[:, 0].min()
        self.yshift = chunk[:, 1].min()
        self.zshift = chunk[:, 2].min()

        chunk[:, 0] -= self.xshift
        chunk[:, 1] -= self.yshift
        chunk[:, 2] -= self.zshift

        xx, yy, zz = chunk[:, :3].T

        # length along each direction - starts at 0
        self.lx = xx.max()
        self.ly = yy.max()
        self.lz = zz.max()

        # Number of parts in which to divide the length along each direction
        self.npx_x = int(self.lx // spacing)
        self.npx_y = int(self.ly // spacing)
        self.npx_z = int(self.lz // spacing)
        print("Grid dimensions are {} {} {}".format(self.npx_x,
                                                    self.npx_y, self.npx_z))

        # make a directory and carry all operations there
        if os.path.isdir(dir_path):
            a = input("directory already exists!, do you want to proceed")
            if (a == 'n') or (a == 'N'):
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

            # save pixel_data file and other variables needed later
            np.savetxt('pixel_data.dat', chunk)
            np.savetxt('rot_matrix.dat', self.rot_matrix)
            np.savetxt('shift_xyz.dat', [self.xshift, self.yshift, self.zshift])

            # write config file
            with open('void.cfg', 'w') as cf:
                cf.write("lx = %f\n" % self.lx)
                cf.write("ly = %f\n" % self.ly)
                cf.write("lz = %f\n" % self.lz)
                cf.write("num_pixels = %i\n" % len(xx))
                cf.write("map_nx = %i\n" % self.npx_x)
                cf.write("map_ny = %i\n" % self.npx_y)
                cf.write("map_nz = %i\n" % self.npx_z)
                cf.write("corr_var_s = %f\n" % self.var_s)
                cf.write("corr_l_perp = %f\n" % self.l_perp)
                cf.write("corr_l_para = %f\n" % self.l_para)
                cf.write("pcg_tol = 1.0e-5\n")
                cf.write("pcg_max_iter = 5000\n")

            # run the Weiner reconstruction on the chunk
            if reconstruct:
                print(os.getcwd())
                executable = "../../../dachshund/dachshund.exe"
                message1 = subprocess.run([executable, 'void.cfg'],
                                          stdout=sys.stdout)
                print(message1)
        except Exception:
            raise
        finally:
            os.chdir(cwd)
