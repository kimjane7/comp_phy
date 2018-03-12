import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
from pylab import *
matplotlib.rcParams['font.family'] = "serif"
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

filenames = [("benchmark/oscillator"+str(i)+".dat") for i in range(1,4)]
labels = ["ground state","first excited state","second excited state"]
colors = ['skyblue', 'mediumseagreen', 'darksalmon']

plt.figure(figsize=(10,8))
fig = plt.figure(1)
axes = plt.gca()
axes.set_xlim([-0.01,5])
axes.tick_params(labelsize=14)

for i in range(0,3):
	file = np.loadtxt(filenames[i],unpack=True)
	plt.plot(file[0],file[1],linewidth=2,label=labels[i],color=colors[i])


legend = plt.legend(loc='upper right', shadow=True)
plt.xlabel(r'$\rho$', fontsize=16, weight='normal', family='serif')
plt.ylabel(r'$u(\rho)$', fontsize=16, weight='normal', family='serif')
plt.title(r'Wavefunctions for First Three Energy States', fontsize=16, weight='normal', family='serif')
plt.grid()

figname = 'wvfunc_1e.pdf'
plt.savefig(figname, format='pdf')
os.system('okular '+figname)