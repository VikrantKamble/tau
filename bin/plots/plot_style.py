import matplotlib as mpl

cfile = "./config.ini"

mpl.rcParams['text.usetex'] = True
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['figure.figsize'] = (6.4, 4.6)
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18

# mpl.rcParams['xtick.direction'] = 'in'
# mpl.rcParams['ytick.direction'] = 'in'

mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['axes.titlesize'] = 20

mpl.rcParams['lines.linewidth'] = 0.8
mpl.rcParams['font.size'] = 14

mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['legend.fontsize'] = 18