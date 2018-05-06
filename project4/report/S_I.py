import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.tick_params(labelsize=12)

det = np.loadtxt("benchmark/deterministic_A.dat",unpack=True)
sto = np.loadtxt("benchmark/montecarlo_A0.dat",unpack=True)


plt.plot(sto[1],sto[2],linewidth=1,linestyle='-',color='darkcyan')
plt.plot(det[1],det[2],linewidth=1,linestyle='-',color='k')

plt.legend(loc=1, shadow=True)
plt.xlabel(r'Susceptible People', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'Infected People', fontsize=12, weight='normal', family='serif')
plt.title(r'Population A', fontsize=12, weight='normal', family='serif')
#plt.grid()
plt.tight_layout()

figname = 'S_I.png'
plt.savefig(figname, format='png')
os.system('okular '+figname)
plt.clf()
