import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

Xfile = np.loadtxt("benchmark/solarsystem_CM_euler_X.dat",unpack=True)
Yfile = np.loadtxt("benchmark/solarsystem_CM_euler_Y.dat",unpack=True)
Zfile = np.loadtxt("benchmark/solarsystem_CM_euler_Z.dat",unpack=True)

colors = ['indianred','salmon','orange','gold','yellowgreen','cyan','dodgerblue','slateblue','purple','magenta']
labels = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']

plt.figure(figsize=(9,6))
fig = plt.figure(1)
axes = fig.gca(projection='3d')
axes.tick_params(labelsize=12)

for i in range(10,0,-1):
	axes.plot(Xfile[i],Yfile[i],Zfile[i],linewidth=1,color=colors[i-1],label=labels[i-1])

axes.set_xlim3d(-30, 30)
axes.set_ylim3d(-30, 30)
axes.set_zlim3d(-30, 30)
axes.set_xlabel('x')
axes.set_ylabel('y')
axes.set_zlabel('z')
plt.title(r'Orbits of Planets (Forward Euler)', fontsize=12, weight='normal', family='serif')
plt.legend(loc=2, shadow=True)
#plt.tight_layout()

figname = 'solar_system_euler_orbit.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)
plt.clf()

Xfile = np.loadtxt("benchmark/solarsystem_CM_vv_X.dat",unpack=True)
Yfile = np.loadtxt("benchmark/solarsystem_CM_vv_Y.dat",unpack=True)
Zfile = np.loadtxt("benchmark/solarsystem_CM_vv_Z.dat",unpack=True)

plt.figure(figsize=(9,6))
fig = plt.figure(1)
axes = fig.gca(projection='3d')
axes.tick_params(labelsize=12)

for i in range(10,0,-1):
	axes.plot(Xfile[i],Yfile[i],Zfile[i],linewidth=1,color=colors[i-1],label=labels[i-1])

axes.set_xlim3d(-30, 30)
axes.set_ylim3d(-30, 30)
axes.set_zlim3d(-30, 30)
axes.set_xlabel('x')
axes.set_ylabel('y')
axes.set_zlabel('z')
plt.title(r'Orbits of Planets (Velocity Verlet)', fontsize=12, weight='normal', family='serif')
plt.legend(loc=2, shadow=True)
#plt.tight_layout()

figname = 'solar_system_vv_orbit.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)