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
fig1,(ax1,ax3) = plt.subplots(2,1)
fig2,ax2 = plt.subplots(1,1)
fig3,ax4 = plt.subplots(1,1)
Dxs = ['0_16','1_16','2_16','3_16','4_16','5_16','6_16','7_16','8_16']
Dxsa = np.array([0,1/16,2/16,3/16,4/16,5/16,6/16,7/16,8/16])
dirs = ['dxs'+d for d in Dxs]
colors = ['b','r','g','y','m','c','k','C0','C1','C2','C3']*10
Fdx = []
Sm = [0,1,5,10,15,25]
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
                F[cn].append(np.average(fx[int(len(fx)/3.0):]))
    Fx1 = np.array(F[0])
    Fx2 = np.array(F[1])
    De,Fx1,Fx2,Tau = np.transpose(sorted(zip(De,Fx1,Fx2,Tau)))
    ax1.plot(De,ndm.ndm_force(Fx1)*ndm.ndm_k(k)/ndm.ndm_omega(omega)/ndm.ndm_eta(eta),'.-',color=colors[di],label=r'$\Delta x_h=%s\lambda$'%Dxsa[di],markersize=8)
    ax1.plot(De,ndm.ndm_force(Fx2)*ndm.ndm_k(k)/ndm.ndm_omega(omega)/ndm.ndm_eta(eta),'.--',color=colors[di],markersize=8)
    ax3.plot(De,Fx1,'.-',color=colors[di],label=r'$\Delta x_h=%s\lambda$'%Dxsa[di])
    Fdx.append(Fx1)
    for si,sm in enumerate(Sm):
        Fxc1[si].append(Fx1[sm])
        Fxc2[si].append(Fx2[sm])

ax1.set_xlabel('De')
ax1.set_ylabel(r'$F_h$')
ax1.legend(loc='best',ncol=3)
Fxc1 = np.array(Fxc1)
Fxc2 = np.array(Fxc2)
for si,sm in enumerate(Sm):
    ax2.plot(Dxsa*2*np.pi,ndm.ndm_force(Fxc1[si])*ndm.ndm_k(k)/ndm.ndm_omega(omega)/ndm.ndm_eta(eta),'.-',color=colors[si],label=r'De$=%s$'%myround(De[sm]),markersize=8,linewidth=1)
    #ax2.plot(Dxsa,Fxc2[si],'o--',color=colors[si])
Fdx = np.transpose(Fdx)
Grad = []
for de,fdx in zip(De,Fdx):
    found = False
    for i,f in enumerate(fdx[0:-1]):
        if f*fdx[i+1]<=0:
          Grad.append((fdx[i+1]-fdx[i])/(Dxsa[i+1]-Dxsa[i])/(2*np.pi/k))
          ax4.plot([de,],[ndm.ndm_grad_force(Grad[-1]),],'b.',markersize=8)
          found = True
          break
    if not found:
        Grad.append(0)
        ax4.plot([de,],[0,],'rx')

ax2.plot([np.min(Dxsa)*2*np.pi,np.max(Dxsa)*2*np.pi],[0,0],'k--')
ax2.set_xlabel(r'$kr_h$')
ax2.set_ylabel(r'$F_hk/\omega/\eta$')
ax3.set_xlabel(r'De')
ax3.set_ylabel(r'$F_h$')
ax4.set_xlabel(r'De')
ax4.set_ylabel(r'$\nabla F_h r_c/F_\mathrm{ref}$')
#ax1.set_ylim((-60,80))
#ax3.set_ylim((-20,110))
ax2.legend(loc='best',ncol=2)
ax3.legend(loc='best',ncol=3)
#ax2.set_ylim((ndm.ndm_force(-15),ndm.ndm_force(75)))
plt.show()

