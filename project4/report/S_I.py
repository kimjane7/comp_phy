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
#axes.set_xlim([0,500])
#axes.set_ylim([0,500])
axes.tick_params(labelsize=12)

det = np.loadtxt("benchmark/deterministic_test.dat",unpack=True)
sto = np.loadtxt("benchmark/stochastic_test1.dat",unpack=True)


plt.plot(det[1],det[2],linewidth=1,linestyle='-',color='darkgray')
plt.plot(sto[1],sto[2],linewidth=1,linestyle='-',color='darkcyan')

plt.legend(loc=1, shadow=True)
plt.xlabel(r'Susceptible People', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'Infected People', fontsize=12, weight='normal', family='serif')
plt.title(r'SIRS Model', fontsize=12, weight='normal', family='serif')
#plt.grid()
plt.tight_layout()

figname = 'S_I.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)
plt.clf()
