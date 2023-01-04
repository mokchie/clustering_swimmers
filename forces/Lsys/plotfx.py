from variable3 import *
from readvel3 import *
import pdb
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def coarse(lst):
    return np.array(coarseave(lst,50))
fig,(ax1,ax2) = plt.subplots(2,1)
D = [4.4,4.8,5.2,5.6,6.0]
dirs = ['sdist'+str(d) for d in D]
colors = ['b','r','g','y','m','c','k']*10
Sm = [0,5,7,10,15,25]
Fxc1 = [[] for i in range(len(Sm))]
Fxc2 = [[] for i in range(len(Sm))]
for di,direct in enumerate(dirs):
    infiles = [direct+'/'+i for i in ['in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run','in-7.run','in-8.run','in-9.run']]
    Tau = []
    De = []
    F = [[],[]]
    for infile in infiles:
        Dt = readvar(infile,'Dt')
        taus = readvar(infile,'tau')
        omega = -readvar(infile,'omega')
        for tau in taus:
            Tau.append(tau)
            De.append(tau*omega)
            linestyles = ['-','--']
            file_list = zip([direct+'/data/fs1-%s.data'%myround(tau), direct+'/data/fs2-%s.data'%myround(tau)],['c_cfs1','c_cfs2'])
            for cn,(f_fs,c_cf) in enumerate(file_list):
                dst,ds = readrow(f_fs)

                t = coarse(np.array(dst['TimeStep'])*Dt)

                fs = np.array(ds[c_cf])
                fx = fs[:,0]
                fy = fs[:,1]
                F[cn].append(np.average(fx[int(len(fx)/3.0):]))
    Fx1 = F[0]
    Fx2 = F[1]
    De,Fx1,Fx2,Tau = np.transpose(sorted(zip(De,Fx1,Fx2,Tau)))
    ax1.plot(De,Fx1,'o-',color=colors[di],label='d=%s'%D[di])
    ax1.plot(De,Fx2,'o--',color=colors[di])
    for si,sm in enumerate(Sm):
        Fxc1[si].append(Fx1[sm])
        Fxc2[si].append(Fx2[sm])

ax1.set_xlabel('De')
ax1.set_ylabel('f')
ax1.legend(loc='best',ncol=3)
for si,sm in enumerate(Sm):
    ax2.plot(D,Fxc1[si],'o-',color=colors[si],label=r'De$=%s$'%myround(De[sm]))
    ax2.plot(D,Fxc2[si],'o--',color=colors[si])

ax2.set_xlabel('d')
ax2.set_ylabel('f')
ax2.legend(loc='best',ncol=3)
plt.show()

