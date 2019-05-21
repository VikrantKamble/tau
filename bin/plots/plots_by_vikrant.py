import scipy as sp
import matplotlib.pyplot as plt
import fitsio
import plot_style


from matplotlib.patches import Circle
from tau import tau_class

helion_path = ('/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/'
               'kdawson/hdumasde/Run_programs/igmhub/Void_analysis/')


def plot_angular_distribution():
    """
    Plot the distribution of quasars in ra and dec for real data and
    for simulated data
    """
    ra_cen = 213.704
    dec_cen = 53.083

    # real and simulated data
    real_file = helion_path + 'picca_run_delta/data_delta_transmission_RMplate.fits'

    simul_file = helion_path + 'lya_forest/london/v6.0/v6.0.0_delta_transmission_RM.fits.gz'
    # simul_file = helion_path + 'lya_forest/saclay/v4.4/v4.4.0_delta_transmission_RM.fits.gz'

    real = tau_class.TauClass(real_file)
    simul = tau_class.TauClass(simul_file)

    tb = fitsio.read(helion_path + 'picca_run_delta/Catalogs/RM-qso.fits', 1)

    fig, ax = plt.subplots(1, figsize=(7, 6.7))
    stretch = sp.cos(sp.deg2rad(dec_cen))

    ax.scatter((tb['RA'] - ra_cen) * stretch, tb['DEC'] - dec_cen,
               marker='+', c='k', alpha=0.6, label=r'$\mathrm{All\ RM\ QSOs}$')

    ax.scatter((real.q_loc[:, 0] - ra_cen) * stretch,
               real.q_loc[:, 1] - dec_cen,
               marker='o', c='r', alpha=0.6, label=r'$\mathrm{QSOs\ used}$')

    ax.scatter((simul.q_loc[:, 0] - ra_cen) * stretch,
               simul.q_loc[:, 1] - dec_cen,
               marker='o', c='g', alpha=0.6, label=r'$\mathrm{Simulation\ 0}$')

    circ = Circle([0, 0], 1.5, facecolor='none', edgecolor='k')
    ax.add_patch(circ)

    title_str = r'$\mathrm{RA}_{cen} = %.2f\ \mathrm{deg},' + \
                r'\mathrm{DEC}_{cen} = %.2f\ \mathrm{deg}$'
    plt.title(title_str % (ra_cen, dec_cen))
    plt.xlabel(r'$\mathrm{RA\ offset\ [deg]}$')
    plt.ylabel(r'$\mathrm{DEC\ offset\ [deg]}$')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('ra_dec_dist.pdf')
    plt.show()


def get_avg_sep():
    real = tau_class.TauClass('../../data/delta-2.fits.gz')
    real.get_data(skewers_perc=1., ylabel='DELTA')

    avg_r = real.pixel_data[:, 2].mean()
    print("Average distance is {}".format(avg_r))

    ra, dec = real.q_loc[:, 0], real.q_loc[:, 1]

    delta_ra = sp.abs(ra - ra[:, None])
    delta_dec = sp.abs(dec - dec[:, None])

    angle = 2 * sp.arcsin(sp.sqrt(sp.sin(delta_dec / 2.) ** 2 +
                                  sp.cos(dec) * sp.cos(dec[:, None]) *
                                  sp.sin(delta_ra / 2.) ** 2))

    # remove self-distances
    sp.fill_diagonal(angle, sp.inf)

    min_angle = angle.min(0) * 180 / sp.pi
    avg_angle = min_angle.mean()

    p_cmd = "Average distance is {} h^-1 Mpc \n".format(avg_r)
    p_cmd += "Average angle is {} degree \n".format(avg_angle)

    avg_sep = avg_angle * avg_r * sp.pi / 180
    p_cmd += "Average transverse separation is {} h^-1 Mpc \n".format(avg_sep)

    print(p_cmd)
    # *************************************************************************

    simul = tau_class.TauClass('../../data/simulations/'
                               'v6.0.1_delta_transmission_RMplate.fits')
    simul.get_data(skewers_perc=1.)

    avg_r = simul.pixel_data[:, 2].mean()
    print("Average distance is {}".format(avg_r))

    ra, dec = sp.deg2rad(simul.q_loc[:, 0]), sp.deg2rad(simul.q_loc[:, 1])

    delta_ra = sp.abs(ra - ra[:, None])
    delta_dec = sp.abs(dec - dec[:, None])

    angle = 2 * sp.arcsin(sp.sqrt(sp.sin(delta_dec / 2.) ** 2 +
                                  sp.cos(dec) * sp.cos(dec[:, None]) *
                                  sp.sin(delta_ra / 2.) ** 2))

    # remove self-distances
    sp.fill_diagonal(angle, sp.inf)

    min_angle = angle.min(0) * 180 / sp.pi
    avg_angle = min_angle.mean()

    p_cmd = "Average distance is {} h^-1 Mpc \n".format(avg_r)
    p_cmd += "Average angle is {} degree \n".format(avg_angle)

    avg_sep = avg_angle * avg_r * sp.pi / 180
    p_cmd += "Average transverse separation is {} h^-1 Mpc \n".format(avg_sep)

    print(p_cmd)


def plot_zq_hist():
    # real and simulated data
    real_file = helion_path + 'picca_run_delta/Backup/Delta/data_delta_transmission_RMplate.fits'

    simul_file = helion_path + 'lya_forest/london/v6.0/v6.0.0_delta_transmission_RM.fits.gz'
    # simul_file = helion_path + 'lya_forest/saclay/v4.4/v4.4.0_delta_transmission_RM.fits.gz'

    real = tau_class.TauClass(real_file)
    simul = tau_class.TauClass(simul_file)

    tb = fitsio.read(helion_path + 'picca_run_delta/Catalogs/RM-qso.fits', 1)

    fig, ax = plt.subplots(1, figsize=(7, 6.7))

    bins = sp.linspace(0, 5, 50)

    ax.hist(tb['Z'], bins=bins, histtype='step', color='k', alpha=0.4, label=r'$\mathrm{All\ RM\ QSOs}$')
    ax.hist(real.q_loc[:, 2], bins=bins, histtype='step', color='r', label=r'$\mathrm{QSOs\ used}$')
    ax.hist(simul.q_loc[:, 2], bins=bins, histtype='step', color='g', label=r'$\mathrm{Simulation\ 0}$')

    plt.xlabel(r'$z_q$')
    plt.ylabel(r'$\#$')

    plt.legend()
    plt.tight_layout()
    plt.savefig('red_dist.pdf')

    plt.show()


def plot_delta_pdf():
    real_file = helion_path + 'picca_run_delta/Backup/Delta/data_delta_transmission_RMplate.fits'

    # simul_file = helion_path + 'lya_forest/london/v6.0/v6.0.0_delta_transmission_RM.fits.gz'
    simul_file = helion_path + 'lya_forest/saclay/v4.4/v4.4.0_delta_transmission_RM.fits.gz'

    real = fitsio.FITS(real_file)
    simul = fitsio.FITS(simul_file)

    real_red, real_delta = [], []
    for i in range(1, len(real)):
        tb = real[i].read()
        real_red.append(tb['LAMBDA'])
        real_delta.append(tb['DELTA_T'])

    simul_red, simul_delta = [], []
    for i in range(1, len(simul)):
        tb = simul[i].read()
        simul_red.append(tb['LAMBDA'])
        simul_delta.append(tb['DELTA_T'])

    real_delta = sp.hstack(real_delta)
    simul_delta = sp.hstack(simul_delta)
    # print(sp.hstack(redshift).shape)

    plt.figure(figsize=(7, 5.6))
    plt.hist(real_delta, bins=sp.linspace(-2.5, 2, 150), density=True,
             color='r', alpha=0.4, label=r'$\mathrm{Real\ data}$')
    plt.hist(simul_delta, bins=sp.linspace(-2.5, 2, 150), density=True,
             color='g', alpha=0.4, label=r'$\mathrm{Simulation\ 0}$')
    plt.xlabel(r'$\delta_F$')
    plt.ylabel(r'$\#$')
    plt.xlim(-1.5, 1.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig('delta_pdf.pdf')
    plt.show()


if __name__ == "__main__":
    plot_delta_pdf()
    # plot_angular_distribution()
    # plot_zq_hist()