import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

pops = ['A','B','C','D']
a = [r'a = 4',r'a = 4',r'a = 4',r'a = 4']
b = [r'b = 1',r'b = 2',r'b = 3',r'b = 4']
c = [r'c = 0.5',r'c = 0.5',r'c = 0.5',r'c = 0.5']

labels = ['Susceptible','Infected','Resistant']
colors = ['yellowgreen','indianred','dodgerblue']


for i in range(0,4):

	files = ['benchmark/montecarlo_'+pops[i]+str(j)+'.dat' for j in range(0,100)]

	plt.figure(figsize=(6,6))
	fig = plt.figure(1)
	axes = plt.gca()
	axes.set_xlim([0,20])
	axes.set_ylim([0,400])
	axes.tick_params(labelsize=12)

	for j in range(0,100):
		sto = np.loadtxt(files[j],unpack=True)
		for k in range(1,4):
			plt.plot(sto[0],sto[k],linewidth=1,linestyle='-',color=colors[k-1],alpha=0.3)

	stats = np.loadtxt("benchmark/montecarlo_"+pops[i]+"100_stats.dat",unpack=True)

	plt.plot(stats[0],stats[1],linewidth=1,linestyle='-',color=colors[0],label=labels[0])
	plt.plot(stats[0],stats[3],linewidth=1,linestyle='-',color=colors[1],label=labels[1])
	plt.plot(stats[0],stats[5],linewidth=1,linestyle='-',color=colors[2],label=labels[2])

	plt.plot(stats[0],stats[1],linewidth=1,linestyle='-',color='k',label='Average')
	plt.plot(stats[0],stats[3],linewidth=1,linestyle='-',color='k')
	plt.plot(stats[0],stats[5],linewidth=1,linestyle='-',color='k')

	'''
	plt.plot(stats[0],stats[1]+stats[2],linewidth=1,linestyle=':',label=r'1 $\sigma$',color='k')
	plt.plot(stats[0],stats[1]-stats[2],linewidth=1,linestyle=':',color='k')
	plt.plot(stats[0],stats[3]+stats[4],linewidth=1,linestyle=':',color='k')
	plt.plot(stats[0],stats[3]-stats[4],linewidth=1,linestyle=':',color='k')
	plt.plot(stats[0],stats[5]+stats[6],linewidth=1,linestyle=':',color='k')
	plt.plot(stats[0],stats[5]-stats[6],linewidth=1,linestyle=':',color='k')
	'''

	plt.legend(loc=1, shadow=True, fontsize=12)
	plt.xlabel(r'Time', fontsize=12, weight='normal', family='serif')
	plt.ylabel(r'Number of People', fontsize=12, weight='normal', family='serif')
	plt.title(r'SIRS Model for Population '+pops[i]+' (100 Samples)', fontsize=12, weight='normal', family='serif')
	#plt.grid()
	plt.tight_layout()

	figname = 'trials_'+pops[i]+'100.png'
	plt.savefig(figname, format='png')
	os.system('okular '+figname)
	plt.clf()


