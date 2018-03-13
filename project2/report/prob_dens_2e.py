import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

filenames = [("benchmark/dots"+str(i)+".dat") for i in range(1,5)]
labels = ['$\omega_r = 0.01$', '$\omega_r = 0.5$', '$\omega_r = 1$', '$\omega_r = 5$']
colors = ['skyblue', 'mediumseagreen', 'darksalmon','plum']

plt.figure(figsize=(10,8))
fig = plt.figure(1)
axes = plt.gca()
axes.set_xlim([0,7])
axes.tick_params(labelsize=14)
axes.annotate('* = no Coulomb interaction',xy=(4.8,0.051),xytext=(4.8,0.051))

for i in range(0,4):
	file = np.loadtxt(filenames[i],unpack=True)
	plt.plot(file[0],file[1],linewidth=2,label=labels[i],color=colors[i])
	plt.plot(file[0],file[2],linewidth=2,linestyle='--',label=labels[i]+'*',color=colors[i])


plt.legend(loc=1, shadow=True)
plt.xlabel(r'$\rho$', fontsize=16, weight='normal', family='serif')
plt.ylabel(r'$|\psi(\rho)|^2$', fontsize=16, weight='normal', family='serif')
plt.title(r'Ground State Probability Densities for Various $\omega_r$', fontsize=16, weight='normal', family='serif')
plt.grid()

figname = 'prob_dens_2e.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)