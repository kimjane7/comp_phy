import matplotlib.pyplot as plt 
from matplotlib import rc
import numpy as np
import os

filenames = [("benchmark/binary_vv"+str(i)+".dat") for i in range(2,7)]
labels = ['$N = 10^2$', '$N = 10^3$', '$N = 10^4$', '$N = 10^5$', '$N = 10^6$']
colors = ['indianred', 'yellowgreen', 'darkturquoise', 'royalblue','black']


plt.figure(1)
plt.figure(figsize=(16,8))
plt.rc('text', usetex=True)


orbit = plt.subplot(221)

plt.xlabel(r'Trajectory of Earth', fontsize=16, weight='normal')
plt.setp(orbit.get_xticklabels(),fontsize=16)
plt.setp(orbit.get_yticklabels(),fontsize=16)
for i in range(0,3):
	file = np.loadtxt(filenames[i],unpack=True)
	plt.plot(file[1],file[2],linewidth=1,label=labels[i],color=colors[i])

energy = plt.subplot(222)

plt.xlabel(r'Energy of Earth', fontsize=16, weight='normal')
plt.setp(energy.get_xticklabels(),fontsize=16)
plt.setp(energy.get_yticklabels(),fontsize=16)
for i in range(0,3):
	file = np.loadtxt(filenames[i],unpack=True)
	plt.plot(file[0],file[5],linestyle='-',linewidth=1,label=labels[i],color=colors[i])


plt.savefig('binary_vv.pdf',format='pdf')
os.system('okular binary_vv.pdf')
