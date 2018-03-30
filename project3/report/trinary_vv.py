import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

Xfiles = [("benchmark/trinary_euler_"+str(i)+"X.dat") for i in range(2,5)]
Yfiles = [("benchmark/trinary_euler_"+str(i)+"Y.dat") for i in range(2,5)]
Efiles = [("benchmark/trinary_euler_"+str(i)+"E.dat") for i in range(2,5)]

labels = ['$N = 10^2$', '$N = 10^3$', '$N = 10^4$']
colors = ['darkturquoise', 'indianred', 'black']

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.set_xlim([-10,10])
axes.set_ylim([-10,10])
axes.tick_params(labelsize=12)

for i in range(0,3):
	Xfile = np.loadtxt(Xfiles[i],unpack=True)
	Yfile = np.loadtxt(Yfiles[i],unpack=True)
	plt.plot(Xfile[2],Yfile[2],linewidth=1,linestyle='-',label=labels[i]+'(Earth)',color=colors[i])
	plt.plot(Xfile[3],Yfile[3],linewidth=1,linestyle='--',label=labels[i]+'(Jupiter)',color=colors[i])

plt.xlabel(r'$x$', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'$y$', fontsize=12, weight='normal', family='serif')
plt.title(r'Orbits of Earth and Jupiter (Velocity Verlet)', fontsize=12, weight='normal', family='serif')
plt.grid()
plt.tight_layout()

figname = 'trinary_vv_orbit.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)
plt.clf()

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.set_xlim([0,14])
axes.tick_params(labelsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

for i in range(0,3):
	Efile = np.loadtxt(Efiles[i],unpack=True)
	plt.plot(Efile[0],Efile[2],linewidth=1,linestyle='-',label=labels[i],color=colors[i])
	plt.plot(Efile[0],Efile[3],linewidth=1,linestyle='--',label=labels[i],color=colors[i])


plt.xlabel(r'$t$ (year)', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'$E$', fontsize=12, weight='normal', family='serif')
plt.title(r'Energies of Earth and Jupiter (Velocity Verlet)', fontsize=12, weight='normal', family='serif')
plt.grid()
plt.tight_layout()

figname = 'trinary_vv_energy.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)