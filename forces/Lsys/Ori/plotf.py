from variable3 import *
from readvel3 import *
import pdb
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
import non_dimensionalize as ndm
matplotlib.rcParams.update({'font.size':13})
def coarse(lst):
    return np.array(coarseave(lst,50))
fig1,ax1 = plt.subplots(1,1)
fig2,ax2 = plt.subplots(1,1)
D = np.array([-0.02,0.0,0.02,0.04])
dirs = ['phi'+str(d) for d in D]
colors = ['b','r','g','y','m','c','k','C0','C1','C2','C3','C4','C5','C6','C7']*10
Sm = [0,5,7,10,15,25]
Fyc1 = [[] for i in range(len(Sm))]
Fyc2 = [[] for i in range(len(Sm))]
for di,direct in enumerate(dirs):
    infiles = [direct+'/'+i for i in ['in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run','in-7.run','in-8.run','in-9.run']]
    Tau = []
    De = []
    F = [[],[]]
    for infile in infiles:
        Dt = readvar(infile,'Dt')
        taus = readvar(infile,'tau')
        omega = -readvar(infile,'omega')
        k = readvar(infile,'k')
        eta = readvar(infile,'mu')
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
                F[cn].append(np.average(fy[int(len(fy)/3.0):]))
    Fy1 = F[0]
    Fy2 = F[1]
    De,Fy1,Fy2,Tau = np.transpose(sorted(zip(De,Fy1,Fy2,Tau)))
    ax1.plot(De,ndm.ndm_force(Fy1)*ndm.ndm_k(k)/ndm.ndm_omega(omega)/ndm.ndm_eta(eta),'.-',color=colors[di],label=r'$\phi/\pi=%.3f$'%(D[di]/np.pi),markersize=8,linewidth=1)
    ax1.plot([0,32],[0,0],'k--')
    #ax1.plot(De,Fy2,'+--',color=colors[di])
    for si,sm in enumerate(Sm):
        Fyc1[si].append(Fy1[sm])
        Fyc2[si].append(Fy2[sm])

ax1.set_xlabel('De')
ax1.set_ylabel('$F_vk/\omega/\eta$')
ax1.legend(loc='best',ncol=1)
for si,sm in enumerate(Sm):
    ax2.plot(D/np.pi,Fyc1[si],'.-',color=colors[si],label=r'De$=%s$'%myround(De[sm]))
    ax2.plot(D/np.pi,Fyc2[si],'+--',color=colors[si])
ax1.set_ylim((-0.4,1.0))
ax2.set_xlabel('$\phi/\pi$')
ax2.set_ylabel('$F_vk/\omega/\eta$')
ax2.legend(loc='best',ncol=2)
plt.show()

