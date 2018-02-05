import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

'''
os.system('g++ LU_decomp.cpp -o LU -std=c++11 -larmadillo')
os.system('g++ gen_tridiag.cpp -o gen -std=c++11')
os.system('g++ symm_tridiag.cpp -o symm -std=c++11')
os.system('g++ simple_tridiag.cpp -o simple -std=c++11')
'''
os.system('LU 3 LU')
os.system('gen 6 gen')
os.system('symm 7 symm')
os.system('simple 8 simple')


files = ['LU','gen','symm', 'simple']
labels = ['n = 1000', 'n = 100', 'n = 10']
colors = ['skyblue', 'mediumseagreen', 'darksalmon']
borders = ['darkcyan', 'darkgreen', 'maroon']


i = 3

plt.figure(figsize=(10,8))
fig = plt.figure(1)
axes = plt.gca()
axes.set_xlim([0,1])
axes.tick_params(labelsize=14)

for j in range(0,3):
	file = np.loadtxt(files[i]+str(3-j)+'.dat', unpack=True)
	plt.plot(file[0],file[1],markerfacecolor=colors[j],markeredgecolor=borders[j],label=labels[j],markersize=10,marker='o',markeredgewidth=1)

exact = np.loadtxt('simple6.dat', unpack=True)
plt.plot(exact[0],exact[2],color='k',linewidth=2,linestyle='-',label='exact')

legend = plt.legend(loc='upper right', shadow=True)
plt.xlabel(r'$x$', fontsize=16, weight='normal', family='serif')
plt.ylabel(r'$v(x)$', fontsize=16, weight='normal', family='serif')
plt.title(r'Exact and Approximate Solutions of 1-D Poisson Equation', fontsize=16, weight='normal', family='serif')

figname = 'compare_'+files[i]+'.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)




