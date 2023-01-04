from variable3 import *
from readvel3 import *
import pdb
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib
import non_dimensionalize as ndm
matplotlib.rcParams.update({'font.size':13})
def coarse(lst):
    return np.array(coarseave(lst,50))
fig1,ax1 = plt.subplots(1,1)
fig2,ax2 = plt.subplots(1,1)
fig3,ax3 = plt.subplots(1,1)
D = np.array([2.4,2.8,3.2,3.6,4.0,4.4])
dirs = ['sdist'+str(d) for d in D]
colors = ['b','r','g','y','m','c','k','C0','C1','C2','C3','C4','C5']*10
Sm = [0,5,10,15,20,25]
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
        k = readvar(infile,'k')
        omega = -readvar(infile,'omega')
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
    if di not in [2,4]:
        ax1.plot(De,ndm.ndm_force(np.array(Fy1))*ndm.ndm_k(k)/ndm.ndm_omega(np.abs(omega))/ndm.ndm_eta(eta),'.-',color=colors[di],label='$kr_v$=%.2f'%(D[di]*k),markersize=8,linewidth=1)
        #ax1.plot(De,Fy2,'.--',color=colors[di])
    for si,sm in enumerate(Sm):
        Fyc1[si].append(Fy1[sm])
        Fyc2[si].append(Fy2[sm])
Fyc1 = np.array(Fyc1)
Fyc2 = np.array(Fyc2)
ax1.set_xlabel('De')
ax1.set_ylabel('$F_vk/\omega/\eta$')
ax1.legend(loc='upper right',ncol=2)
D = np.array(D)*k
def ipow(x,p,k,b):
    return k/x**p+b
for si,sm in enumerate(Sm):
    ax2.plot(D,ndm.ndm_force(np.array(Fyc1[si]))*ndm.ndm_k(k)/ndm.ndm_omega(np.abs(omega))/ndm.ndm_eta(eta),'.-',color=colors[si],label=r'De$=%s$'%myround(De[sm]),markersize=8,linewidth=1)
    popt,pconv = curve_fit(ipow,D,np.abs(Fyc1[si]),p0=[4,1,50])
    ax3.plot(D,np.abs(Fyc1[si]),'.',color=colors[si],label=r'De$=%s$'%myround(De[sm]),markersize=8,linewidth=1)
    ax3.plot(D,ipow(D,popt[0],popt[1],popt[2]))
    print(popt[0])
ax1.plot(De,np.zeros_like(De),'k--')    
ax2.plot(D,np.zeros_like(D),'k--')
ax2.set_xlabel('$kr_v$')
ax2.set_ylabel('$F_vk/\omega/\eta$')
ax2.legend(loc='best',ncol=2)
ax1.set_ylim((-2.7,1))
ax2.set_ylim((-2.7,1))
plt.show()

