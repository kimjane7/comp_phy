import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
 

labels = ['Susceptible','Infected','Resistant']
colors = ['yellowgreen','indianred','dodgerblue']

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.set_xlim([0,15])
axes.set_ylim([0,400])
axes.tick_params(labelsize=12)

detA = np.loadtxt("benchmark/deterministic_A.dat",unpack=True)
detB = np.loadtxt("benchmark/deterministic_B.dat",unpack=True)
detC = np.loadtxt("benchmark/deterministic_C.dat",unpack=True)
detD = np.loadtxt("benchmark/deterministic_D.dat",unpack=True)
sto = np.loadtxt("benchmark/montecarlo_A0.dat",unpack=True)



for i in range(1,4):
	plt.plot(detA[0],detA[i],linewidth=2,linestyle='-',label=labels[i-1],color=colors[i-1])
	plt.plot(detB[0],detB[i],linewidth=2,linestyle='--',color=colors[i-1])
	plt.plot(detC[0],detC[i],linewidth=2,linestyle='-.',color=colors[i-1])
	plt.plot(detD[0],detD[i],linewidth=2,linestyle=':',color=colors[i-1])
	#plt.plot(sto[0],sto[i],linewidth=1,linestyle=(0,(5,1)),label=labels[i-1],color=colors[i-1])

plt.legend(loc=1, shadow=True)
plt.xlabel(r'Time', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'Number of People', fontsize=12, weight='normal', family='serif')
plt.title(r'SIRS Model', fontsize=12, weight='normal', family='serif')
#plt.grid()
plt.tight_layout()

figname = 'SIR_t.png'
plt.savefig(figname, format='png')
os.system('okular '+figname)
plt.clf()
