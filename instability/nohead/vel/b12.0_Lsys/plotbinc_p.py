from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,minimize
import pdb
from numpy.linalg import norm
from numpy.linalg import eig
import pdb
from matplotlib.lines import Line2D
import matplotlib
import non_dimensionalize as ndm
matplotlib.rcParams.update({'font.size': 13})
font = {'size':13}
legend_elements1 = []
legend_elements2 = []
colors = ['b','r','g','y','c','m','k']
linestyles = ['.-','*-']
pointstyles = ['.','x','+','v','<','s','*']

fig1,ax1 = plt.subplots(1,1)
#fig2,ax2 = plt.subplots(1,1)
Dirs = ['single_w0.4','pair_w0.4','single_w0.8','pair_w0.8','single_w1.0','pair_w1.0']
labels = ['single','pair']
Omega = [0.4,0.8,1.0]
B_list = []
for cn,direct in enumerate(Dirs):
    Tau = []
    B = []
    with open(direct+'/B.data','r') as fb:
        for line in fb:
            tau,b1m,b2m = [float(item) for item in line.strip().split()]
            Tau.append(tau)
            B.append((b1m+b2m)/2)
    Tau = np.array(Tau)
    ax1.plot(Tau*Omega[int(cn/2)],B,linestyles[cn%2],color=colors[int(cn/2)],label=r'%s $\omega=%.1f$'%(labels[cn%2],ndm.ndm_omega(Omega[int(cn/2)])),markersize=8,linewidth=1)
    B_list.append(np.array(B))
    cl = Line2D([0],[0],color=colors[int(cn/2)],label=r'$\omega/\omega_\mathrm{ref}=%.1f$'%(myround(ndm.ndm_omega(Omega[int(cn/2)]))))
    if cn%2==0:
        legend_elements1.append(cl)
cl1 = Line2D([0],[0],color='w',marker='.',markerfacecolor='k',label=r'single',markersize=12)
cl2 = Line2D([0],[0],color='w',marker='*',markerfacecolor='k',label=r'pair',markersize=12)
legend_elements2.append(cl1)
legend_elements2.append(cl2)
Tau = np.array(Tau)
#for i,omega in enumerate(Omega):
#    ax2.plot(Tau*Omega[i],B_list[i*2+1]-B_list[i*2],'.-',color=colors[i],label=r'$\omega=%.1f$'%(Omega[int(i/2)],))
lg1 = ax1.legend(handles=legend_elements1,loc='lower right',ncol=1)
lg2 = ax1.legend(handles=legend_elements2,loc='lower left',ncol=1)
ax1.add_artist(lg1)
ax1.set_xlabel('De')
ax1.set_ylabel('B')
#ax1.legend(loc='best',ncol=4)
ax1.set_ylim((1,6))
#ax2.set_xlabel('De')
#ax2.set_ylabel(r'$\Delta $B')
#ax2.legend(loc='best',ncol=2)


plt.show()
