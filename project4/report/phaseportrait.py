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

	files = ['benchmark/phaseportrait_'+pops[i]+str(j)+'.dat' for j in range(0,12)]
	plt.figure(figsize=(6,6))
	fig = plt.figure(1)
	axes = plt.gca()
	axes.set_xlim([0,1])
	axes.set_ylim([0,1])
	axes.tick_params(labelsize=12)

	# curves
	for j in range(0,12):
		file = np.loadtxt(files[j],unpack=True)
		plt.plot(file[1]/400,file[2]/400,linewidth=1,linestyle='-',color='darkcyan')

	# triangle boundary
	x = np.arange(0.0,1.0,0.01)
	y = 1.0-x
	plt.plot(x,y,linewidth=1,linestyle='--',color='k')

	# details
	plt.text(0.70, 0.8, a[i],fontsize=14)
	plt.text(0.70, 0.75, b[i],fontsize=14)
	plt.text(0.70, 0.7, c[i],fontsize=14)

	plt.legend(loc=1, shadow=True, fontsize=12)
	plt.xlabel(r'$S/N$', fontsize=12, weight='normal', family='serif')
	plt.ylabel(r'$I/N$', fontsize=12, weight='normal', family='serif')
	plt.title(r'Phase Portrait for Population '+pops[i], fontsize=12, weight='normal', family='serif')
	#plt.grid()
	plt.tight_layout()

	figname = 'phaseportrait_'+pops[i]+'.png'
	plt.savefig(figname, format='png')
	os.system('okular '+figname)
	plt.clf()


