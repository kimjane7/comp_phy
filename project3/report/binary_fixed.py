import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter


#######################################
###########  FORWARD EULER  ###########
#######################################

'''
Xfiles = [("benchmark/binary_fixed_euler"+str(i)+"_X.dat") for i in range(3,8)]
Yfiles = [("benchmark/binary_fixed_euler"+str(i)+"_Y.dat") for i in range(3,8)]
Efiles = [("benchmark/binary_fixed_euler"+str(i)+"_E.dat") for i in range(3,8)]

labels = ['$\Delta t = 2.5\cdot 10^{-2}$ yr', '$\Delta t = 2.5\cdot 10^{-3}$ yr', '$\Delta t = 2.5\cdot 10^{-4}$ yr', '$\Delta t = 2.5\cdot 10^{-5}$ yr', '$\Delta t = 2.5\cdot 10^{-6}$ yr']
colors = ['indianred','yellowgreen','cyan','dodgerblue','purple']

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.set_xlim([-3,3])
axes.set_ylim([-3,3])
axes.tick_params(labelsize=12)

for i in range(0,5):
	Xfile = np.loadtxt(Xfiles[i],unpack=True)
	Yfile = np.loadtxt(Yfiles[i],unpack=True)
	plt.plot(Xfile[2],Yfile[2],linewidth=1,label=labels[i],color=colors[i])

plt.legend(loc=1, shadow=True)
plt.xlabel(r'$x$', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'$y$', fontsize=12, weight='normal', family='serif')
plt.title(r'Orbit of Earth (Forward Euler)', fontsize=12, weight='normal', family='serif')
plt.grid()
plt.tight_layout()

figname = 'binary_fixed_euler_orbit.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)
plt.clf()

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.tick_params(labelsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

for i in range(0,5):
	Efile = np.loadtxt(Efiles[i],unpack=True)
	plt.plot(Efile[0],Efile[2],linewidth=1,label=labels[i],color=colors[i])

plt.legend(loc=1, shadow=True)
plt.xlabel(r'$t$ (yr)', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'$E$', fontsize=12, weight='normal', family='serif')
plt.title(r'Energy of Earth (Forward Euler)', fontsize=12, weight='normal', family='serif')
plt.grid()
plt.tight_layout()

figname = 'binary_fixed_euler_energy.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)
plt.clf()

'''
#######################################
##########  VELOCITY VERLET  ##########
#######################################


Xfiles = [("benchmark/binary_fixed_vv"+str(i)+"_X.dat") for i in range(1,8)]
Yfiles = [("benchmark/binary_fixed_vv"+str(i)+"_Y.dat") for i in range(1,8)]
Efiles = [("benchmark/binary_fixed_vv"+str(i)+"_E.dat") for i in range(1,8)]

labels = ['$\Delta t = 2.5$ yr','$\Delta t = 0.25$ yr','$\Delta t = 2.5\cdot 10^{-2}$ yr', '$\Delta t = 2.5\cdot 10^{-3}$ yr', '$\Delta t = 2.5\cdot 10^{-4}$ yr', '$\Delta t = 2.5\cdot 10^{-5}$ yr', '$\Delta t = 2.5\cdot 10^{-6}$ yr']
colors = ['indianred','salmon','orange','yellowgreen','cyan','dodgerblue','purple']

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.set_xlim([-3,1.5])
axes.set_ylim([-1.5,3])
axes.tick_params(labelsize=12)

for i in range(1,7):
	Xfile = np.loadtxt(Xfiles[i],unpack=True)
	Yfile = np.loadtxt(Yfiles[i],unpack=True)
	plt.plot(Xfile[2],Yfile[2],linewidth=1,label=labels[i],color=colors[i])

plt.legend(loc=2, shadow=True)
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
axes.set_xlim([0,3])
axes.set_ylim([-5.95e-05,-5.9e-05])
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

for i in range(1,7):
	Efile = np.loadtxt(Efiles[i],unpack=True)
	plt.plot(Efile[0],Efile[2],linewidth=2,label=labels[i],color=colors[i])

plt.legend(loc=4, shadow=True)
plt.xlabel(r'$t$ (yr)', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'$E$', fontsize=12, weight='normal', family='serif')
plt.title(r'Energy of Earth (Velocity Verlet)', fontsize=12, weight='normal', family='serif')
plt.grid()
plt.tight_layout()

figname = 'binary_fixed_vv_energy.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)
plt.clf()