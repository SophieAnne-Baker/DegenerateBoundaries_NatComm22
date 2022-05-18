import pandas as pd
import matplotlib
from matplotlib import rc
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
import matplotlib.pyplot as plt
import numpy as np

# Reading into dataframes
line = pd.read_csv('dataT5_line.txt', delimiter=',')
par1 = pd.read_csv('dataT2_wcs.txt', delimiter=',')
par2 = pd.read_csv('dataT3_t.txt', delimiter=',')
par3 = pd.read_csv('dataT4_a.txt', delimiter=',')

# Colour maps
c_par2 = plt.get_cmap('plasma')
c_par3 = plt.get_cmap('viridis')
c_line = plt.get_cmap('gist_gray')
c_par1 = matplotlib.colors.ListedColormap(['dimgrey', 'darkgrey', 'lightgrey'])


# Plotting
fig = plt.figure()
ax1 = fig.add_axes([0.15,0.15,0.55,0.75])

plt.ylim(None, 17)
plt.yticks(np.arange(0, 17, step=5))
plt.xlim(None, 0.5)
ax1.tick_params(labelsize=18)
plt.xlabel('mean error', fontsize=26)
plt.ylabel('mean decision time', fontsize=26)
plt.title(r'$\theta$', fontsize=28)

line = line.sort_values(by=['line_x'])

x_shift = 0.0
y_shift = 0.0

l = ax1.plot(line['line_x'], line['line_y'], color='black', linestyle='--', linewidth=1.5, zorder=0)

p2 = ax1.scatter(par2['X']-x_shift, par2['Y']-y_shift, c=par2['par2'], cmap=c_par2, s=5, zorder=5)

plt.annotate('(a)',
ha = 'center', va = 'bottom',
xytext = (0.075, 11.5),
xy = (0.11, 10.7),
fontsize=22)

plt.annotate('(b)',
ha = 'center', va = 'bottom',
xytext = (0.21, 5.9),
xy = (0.11, 10.7),
fontsize=22)

plt.annotate('(c)',
ha = 'center', va = 'bottom',
xytext = (0.39, 2.85),
xy = (0.11, 10.7),
fontsize=22)

pos2=fig.add_axes([0.75, 0.15, 0.03, 0.7])

cb2 = plt.colorbar(p2, shrink=0.2, cax=pos2)

cb2.set_label(r'$\theta$', size=16)

plt.show()
