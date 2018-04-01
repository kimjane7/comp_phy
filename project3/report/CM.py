import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

colors = ['indianred', 'yellowgreen', 'darkturquoise', 'royalblue','black']

Xfile = np.loadtxt("benchmark/CM_euler_X.dat",unpack=True)
Yfile = np.loadtxt("benchmark/CM_euler_Y.dat",unpack=True)

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.tick_params(labelsize=12)

for i in range(2,11):
	plt.plot(Xfile[i],Yfile[i],linewidth=1,color=colors[3])

plt.xlabel(r'$x$', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'$y$', fontsize=12, weight='normal', family='serif')
plt.title(r'Orbits of Planets in CM frame (Forward Euler)', fontsize=12, weight='normal', family='serif')
plt.grid()
plt.tight_layout()

figname = 'CM_euler_orbit.pdf'
plt.savefig(figname, format='pdf')
plt.clf()

Xfile = np.loadtxt("benchmark/CM_vv_X.dat",unpack=True)
Yfile = np.loadtxt("benchmark/CM_vv_Y.dat",unpack=True)

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.tick_params(labelsize=12)

for i in range(2,11):
	plt.plot(Xfile[i],Yfile[i],linewidth=1,color=colors[3])

plt.xlabel(r'$x$', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'$y$', fontsize=12, weight='normal', family='serif')
plt.title(r'Orbits of Planets in CM frame (Velocity Verlet)', fontsize=12, weight='normal', family='serif')
plt.grid()
plt.tight_layout()

figname = 'CM_vv_orbit.pdf'
plt.savefig(figname, format='pdf')