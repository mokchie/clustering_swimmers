from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def coarse(lst):
    return np.array(coarseave(lst,50))
fig,ax1 = plt.subplots(1,1)
infiles = ['in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run','in-7.run','in-8.run','in-9.run']
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
        for cn,(f_fs,c_cf) in enumerate(zip(['data/fs1-%s.data'%myround(tau), 'data/fs2-%s.data'%myround(tau)],['c_cfs1','c_cfs2'])):
            dst,ds = readrow(f_fs)

            t = coarse(np.array(dst['TimeStep'])*Dt)

            fs = np.array(ds[c_cf])
            fx = fs[:,0]
            fy = fs[:,1]
            F[cn].append(np.average(fy[int(len(fy)/3.0):]))

Fx = F[0]
Fy = F[1]
De,Fx,Fy,Tau = np.transpose(sorted(zip(De,Fx,Fy,Tau)))
ax1.plot(De,Fx,'bo-')
ax1.plot(De,Fy,'ro-')
ax1.set_xlabel('De')
ax1.set_ylabel('$F_y$')
plt.show()

