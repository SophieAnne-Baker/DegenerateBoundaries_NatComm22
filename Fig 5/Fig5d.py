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
c_par1 = matplotlib.colors.ListedColormap(['limegreen', 'deeppink', 'blue'])


# Plotting
fig = plt.figure()
ax1 = fig.add_axes([0.15,0.15,0.55,0.75])

plt.ylim(None, 17)
plt.yticks(np.arange(0, 17, step=5))
plt.xlim(None, 0.5)
ax1.tick_params(labelsize=18)
plt.xlabel('mean error', fontsize=26)
plt.ylabel('mean decision time', fontsize=26)
plt.title(r'c/W', fontsize=28)
#plt.xticks(np.arange(0, 0., step=0.1))

line = line.sort_values(by=['line_x'])

x_shift = 0.05
y_shift = 2

l = ax1.plot(line['line_x'], line['line_y'], color='black', linestyle='--', linewidth=1.5, zorder=0)

p1 = ax1.scatter(par1['X'], par1['Y'], c=par1['par1'], cmap=c_par1, s=5, zorder=1, edgecolors='face')

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

pos1=fig.add_axes([0.75, 0.15, 0.03, 0.7])

cb1 = plt.colorbar(p1, cax=pos1)

cb1.set_label(r'c/W', size=16)
cb1.set_ticks([0.025,0.05,0.08])
cb1.set_ticklabels([0.001,0.04,0.1])
cb1.ax.tick_params(labelsize=16)

fig = plt.gcf()
plt.show()
