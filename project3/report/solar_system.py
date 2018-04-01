import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

Xfile = np.loadtxt("benchmark/test_euler_X.dat",unpack=True)
Yfile = np.loadtxt("benchmark/test_euler_Y.dat",unpack=True)

colors = ['indianred', 'coral', 'orange', 'gold', 'yellowgreen', 'darkturquoise', 'royalblue','darkorchid', 'magenta', 'pink']

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.tick_params(labelsize=12)

for i in range(2,11):
	plt.plot(Xfile[i],Yfile[i],linewidth=1,color=colors[3])

plt.xlabel(r'$x$', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'$y$', fontsize=12, weight='normal', family='serif')
plt.title(r'Orbits of Planets (Forward Euler)', fontsize=12, weight='normal', family='serif')
plt.grid()
plt.tight_layout()

figname = 'solar_system_euler_orbit.pdf'
plt.savefig(figname, format='pdf')
plt.clf()

Xfile = np.loadtxt("benchmark/test_vv_X.dat",unpack=True)
Yfile = np.loadtxt("benchmark/test_vv_Y.dat",unpack=True)

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.tick_params(labelsize=12)

for i in range(2,11):
	plt.plot(Xfile[i],Yfile[i],linewidth=1,color=colors[i-2])

plt.xlabel(r'$x$', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'$y$', fontsize=12, weight='normal', family='serif')
plt.title(r'Orbits of Planets (Velocity Verlet)', fontsize=12, weight='normal', family='serif')
plt.grid()
plt.tight_layout()

figname = 'solar_system_vv_orbit.pdf'
plt.savefig(figname, format='pdf')