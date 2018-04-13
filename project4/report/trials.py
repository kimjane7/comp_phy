import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
 
files = ['benchmark/stochastic_test'+str(i)+'.dat' for i in range(0,100)]
labels = ['Susceptible','Infected','Resistant']
colors = ['yellowgreen','indianred','dodgerblue']

plt.figure(figsize=(6,6))
fig = plt.figure(1)
axes = plt.gca()
axes.set_xlim([0,15])
#axes.set_ylim([-6,8])
axes.tick_params(labelsize=12)


for j in range(0,100):
	sto = np.loadtxt(files[j],unpack=True)
	for i in range(1,4):
		plt.plot(sto[0],sto[i],linewidth=1,linestyle='-',color=colors[i-1],alpha=0.3)

det = np.loadtxt("benchmark/deterministic_test.dat",unpack=True)

for i in range(1,3):
	plt.plot(det[0],det[i],linewidth=1,linestyle='-',label=labels[i-1],color=colors[i-1])
	plt.plot(det[0],det[i],linewidth=1,linestyle='-',color='k')
plt.plot(det[0],det[3],linewidth=1,linestyle='-',label=labels[2],color=colors[2])
plt.plot(det[0],det[3],linewidth=1,linestyle='-',label='Deterministic Solutions',color='k')

plt.legend(loc=1, shadow=True, fontsize=12)
plt.xlabel(r'Time', fontsize=12, weight='normal', family='serif')
plt.ylabel(r'Number of People', fontsize=12, weight='normal', family='serif')
plt.title(r'SIRS Model (100 Samples)', fontsize=12, weight='normal', family='serif')
#plt.grid()
plt.tight_layout()

figname = 'trials.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)
plt.clf()


