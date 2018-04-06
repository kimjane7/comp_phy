import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

##############################################
#            ORBIT WITH FIXED SUN            #
##############################################

Xfiles = [("benchmark/binary_fixed_vv_"+str(i)+"_X.dat") for i in range(2,7)]
Yfiles = [("benchmark/binary_fixed_vv_"+str(i)+"_Y.dat") for i in range(2,7)]
Efiles = [("benchmark/binary_fixed_vv_"+str(i)+"_E.dat") for i in range(2,7)]

labels = ['$\Delta t = 5$ year','$\Delta t = 5E-1$ year','$\Delta t = 5E-2$ year','$\Delta t = 5E-3$ year','$\Delta t = 5E-4$ year']
colors = ['indianred', 'yellowgreen', 'darkturquoise', 'royalblue', 'black',]

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.set_xlim([-1.4,1.4])
axes.set_ylim([-1.4,1.4])
axes.tick_params(labelsize=12)

for i in range(0,5):
	Xfile = np.loadtxt(Xfiles[i],unpack=True)
	Yfile = np.loadtxt(Yfiles[i],unpack=True)
	plt.plot(Xfile[2],Yfile[2],linewidth=1,label=labels[i],color=colors[i])

plt.legend(loc=1, shadow=True)
plt.xlabel(r'$x$', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'$y$', fontsize=12, weight='normal', family='serif')
plt.title(r'Orbit of Earth (Velocity Verlet)', fontsize=12, weight='normal', family='serif')
plt.grid()
plt.tight_layout()

figname = 'binary_fixed_vv_orbit.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)
plt.clf()

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.tick_params(labelsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

for i in range(0,6):
	Efile = np.loadtxt(Efiles[i],unpack=True)
	plt.plot(Efile[0],Efile[2],linewidth=1,label=labels[i],color=colors[i])

plt.legend(loc=1, shadow=True)
plt.xlabel(r'$t$ (year)', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'$E$', fontsize=12, weight='normal', family='serif')
plt.title(r'Energy of Earth (Velocity Verlet)', fontsize=12, weight='normal', family='serif')
plt.grid()
plt.tight_layout()

figname = 'binary_fixed_vv_energy.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)