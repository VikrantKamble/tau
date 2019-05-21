import fitsio
import scipy as sp
import matplotlib.pyplot as plt
import subprocess
import json
import glob


def GetHisto(data,nbBin,weights=None):
    '''
        Get an histo
    '''

    if weights is None:
        hist, axisX = sp.histogram(data,bins=nbBin)
    else:
        hist, axisX = sp.histogram(data,bins=nbBin, weights=weights)
    xxx = sp.array([ axisX[i]+(axisX[i+1]-axisX[i])/2. for i in range(axisX.size-1) ])

    return sp.asarray(list(zip(xxx,hist)))
def plot_one_forest():

    p_delta = '$WORKDIR/Run_programs/igmhub/Void_analysis/picca_run_delta/Delta_LYA/Delta/delta-2.fits.gz'

    h = fitsio.FITS(p_delta)
    for hh in h[1:]:
        if len(hh['LOGLAM'][:])<622: continue
        print(len(hh['LOGLAM'][:]))
        head = hh.read_header()
        print(head)

        el = [head['PLATE'],head['MJD'],head['FIBERID'],head['THING_ID']]

        PLATE = el[0]
        MJD = el[1]
        FIBERID = el[2]
        THINGID = el[3]
        print(PLATE, MJD, FIBERID)
        cmd = 'python /uufs/astro.utah.edu/common/home/u6011908/Articles/void_paper/Work/Codes/do_plotSpec.py '
        cmd += ' --plate '+str(PLATE)
        cmd += ' --mjd '+str(MJD)
        cmd += ' --fiberid '+str(FIBERID)
        cmd += ' --thingid '+str(THINGID)
        cmd += ' --nside 2'
        cmd += ' --drq $WORKDIR/Run_programs/igmhub/Void_analysis/picca_run_delta/Catalogs/RM-qso.fits'
        cmd += ' --in-dir $BOSS_SPECTRO_REDUX/v5_13_0/'
        cmd += ' --mask-file $WORKDIR/Run_programs/igmhub/picca_DR16/picca_DR16_paper_analysis/dr16-line-sky-mask.txt'
        cmd += ' --flux-calib $WORKDIR/Run_programs/igmhub/picca_DR16/picca_DR16_paper_analysis/Delta_calibration/Log/delta_attributes.fits.gz'
        cmd += ' --ivar-calib $WORKDIR/Run_programs/igmhub/picca_DR16/picca_DR16_paper_analysis/Delta_calibration2/Log/delta_attributes.fits.gz'
        cmd += ' --rebin 1'
        cmd += ' --mode spplate'
        cmd += ' --lambda-min 3600.0'
        cmd += ' --lambda-max 7235.0'
        cmd += ' --lambda-rest-min 10.0'
        cmd += ' --lambda-rest-max 10000.0'
        print(cmd)
        subprocess.call(cmd, shell=True)

    h.close()

    return
def plot_survey():

    ra_mid = 0.#213.704
    dec_mid = 0.#53.083
    p_delta = '$WORKDIR/Run_programs/igmhub/Void_analysis/picca_run_delta/Delta_LYA/Delta/delta-2.fits.gz'
    p_all = '$WORKDIR/Run_programs/igmhub/Void_analysis/picca_run_delta/Catalogs/RM-qso.fits'
    p_allQSO = '$EBOSS_ROOT/qso/DR14Q/DR14Q_v3_1.fits'

    h = fitsio.FITS(p_allQSO)
    ra = h[1]['RA'][:]
    dec = h[1]['DEC'][:]
    z = h[1]['Z'][:]
    w = z>2.
    plt.errorbar(ra[w]-ra_mid,dec[w]-dec_mid,fmt='o')

    h = fitsio.FITS(p_all)
    ra = h[1]['RA'][:]
    dec = h[1]['DEC'][:]
    plt.errorbar(ra-ra_mid,dec-dec_mid,fmt='o')

    h = fitsio.FITS(p_delta)
    ra = sp.array([hh.read_header()['RA'] for hh in h[1:]])
    dec = sp.array([hh.read_header()['DEC'] for hh in h[1:]])
    h.close()
    ra *= 180./sp.pi
    dec *= 180./sp.pi
    plt.errorbar(ra-ra_mid,dec-dec_mid,fmt='o')

    radius = 1.5
    theta = sp.arange(0.,2.*sp.pi,0.01)
    x = radius*sp.cos(theta)/sp.absolute(sp.sin(ra_mid*sp.pi/180.))
    y = radius*sp.sin(theta)
    plt.plot(x,y,color='black')

    #plt.xlim([-3,3.])
    #plt.ylim([-3,3.])
    plt.grid()
    plt.show()

    return


def with_cut(cut=1.2):

    ### Data
    pd = '/uufs/chpc.utah.edu/common/home/u1143816/voidReal/*.json'
    ps = '/uufs/chpc.utah.edu/common/home/u1143816/voidRealShuffle/*.json'
    ### Mocks
    #pd = '/uufs/chpc.utah.edu/common/home/u1143816/voidSimulation/*.json'
    #ps = '/uufs/chpc.utah.edu/common/home/u1143816/voidSimulationShuffle/*.json'

    dic = {}
    fs = sp.sort(sp.array(glob.glob(pd)))
    for i,f in enumerate(fs):
        dic['Data'+str(i)] = {'PATH':f, 'NB':None}
    fs = sp.sort(sp.array(glob.glob(ps)))
    for i,f in enumerate(fs):
        dic['Shuffle'+str(i)] = {'PATH':f, 'NB':None}

    for name,valueDic in dic.items():
        nb = []
        with open(valueDic['PATH']) as json_file:
            data = json.load(json_file)
            for k,v in list(data.items()):
                step = sp.array(sorted(v['history'].keys()),dtype=sp.float64)
                stepKey = sp.array(sorted(v['history'].keys()))
                w = step>cut
                if sp.any(w):
                    keys = stepKey[w]
                    nb += [ sp.vstack([ sp.array(v['history'][key]) for key in keys ]).shape[0] ]
        valueDic['NB'] = sp.array(nb)
        print(set(valueDic['NB']))

    ###
    for name,v in dic.items():
        v['HIST'] = GetHisto(v['NB'],sp.arange(0,1000,1))

    ###
    newDic = {}
    for el in ['Data','Shuffle']:
        newDic[el] = {'XXX':None,'HIST':None, 'VAR':None}
        newDic[el]['XXX'] = sp.vstack([ v['HIST'][:,0] for name,v in dic.items() if el in name ]).mean(axis=0)
        newDic[el]['HIST'] = sp.vstack([ v['HIST'][:,1] for name,v in dic.items() if el in name ]).mean(axis=0)
        newDic[el]['VAR'] = sp.vstack([ v['HIST'][:,1] for name,v in dic.items() if el in name ]).std(axis=0)

    for name,v in newDic.items():
        #plt.step(v['XXX'],v['HIST'],alpha=0.8,linewidth=2,label=name)
        p = plt.step(v['XXX'],v['HIST'], where='pre',label=r'$\mathrm{'+name+'}$')
        plt.fill_between(v['XXX'],v['HIST'],v['HIST']+v['VAR'], facecolor=p[0].get_color(), step='pre',alpha=0.5)
        plt.fill_between(v['XXX'],v['HIST'],v['HIST']-v['VAR'], facecolor=p[0].get_color(), step='pre',alpha=0.5)
    plt.yscale('log')
    plt.xlabel(r'$\mathrm{Size \, of \, void \, at \, dT \, = \, '+str(cut)+'}$')
    plt.ylabel(r'$\#$')
    plt.xlim([0.,50.])
    plt.ylim([0.,200.])
    plt.legend()
    plt.grid()
    plt.show()

    return
#plot_one_forest()
#plot_survey()
#with_cut(cut=1.2)
with_cut(cut=0.5)
