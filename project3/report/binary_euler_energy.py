import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

filenames = [("benchmark/binary_euler"+str(i)+".dat") for i in range(2,7)]
labels = ['$N = 10^2$', '$N = 10^3$', '$N = 10^4$', '$N = 10^5$', '$N = 10^6$']
colors = ['indianred', 'yellowgreen', 'darkturquoise', 'royalblue','black']

plt.figure(figsize=(8,8))
fig = plt.figure(1)
axes = plt.gca()
axes.set_xlim([0,10.0])
axes.tick_params(labelsize=14)
#axes.annotate('* = no Coulomb interaction',xy=(4.8,0.051),xytext=(4.8,0.051))

for i in range(0,5):
	file = np.loadtxt(filenames[i],unpack=True)
	plt.plot(file[0],file[5],linestyle='-',linewidth=2,label=labels[i],color=colors[i])

plt.legend(loc=1, shadow=True)
plt.xlabel(r'$x$', fontsize=16, weight='normal', family='serif')
plt.ylabel(r'$y$', fontsize=16, weight='normal', family='serif')
plt.title(r'Energy of Earth (Forward Euler Method)', fontsize=16, weight='normal', family='serif')
plt.grid()

figname = 'binary_euler_energy.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)