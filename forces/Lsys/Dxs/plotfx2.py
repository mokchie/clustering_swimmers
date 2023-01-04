from variable3 import *
from readvel3 import *
import pdb
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def coarse(lst):
    return np.array(coarseave(lst,50))
fig1,ax1 = plt.subplots(1,1)
fig2,ax3 = plt.subplots(1,1)
fig3,(ax2,ax4) = plt.subplots(2,1)
Dxs = ['-4_16','-3_16','-2_16','-1_16','0_16','1_16','2_16','3_16','4_16']
Dxsa = [-4/16,-3/16,-2/16,-1/16,0/16,1/16,2/16,3/16,4/16,]
dirs = ['dxs'+d for d in Dxs]
colors = ['b','r','g','y','m','c','k','C0','C1','C2','C3','C4']*10
Sm = [0,1,2,3,5,7,10,15,25]
Fd = []
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
    ax1.plot(De,Fx1,'.-',color=colors[di],label=r'$\Delta x_h=%s\lambda$'%Dxsa[di])
    ax1.plot(De,Fx2,'.--',color=colors[di])
    ax3.plot(De,Fx1-Fx2,'.-',color=colors[di],label=r'$\Delta x_h=%s\lambda$'%Dxsa[di])
    Fd.append(Fx1-Fx2)
    for si,sm in enumerate(Sm):
        Fxc1[si].append(Fx1[sm])
        Fxc2[si].append(Fx2[sm])

Fxc1 = np.array(Fxc1)
Fxc2 = np.array(Fxc2)
ax1.set_xlabel('De')
ax1.set_ylabel(r'$F_x$')
ax1.legend(loc='best',ncol=3)
for si,sm in enumerate(Sm):
    ax2.plot(Dxsa,Fxc1[si]-Fxc2[si],'.-',color=colors[si],label=r'De$=%s$'%myround(De[sm]))
    #ax2.plot(Dxsa,Fxc2[si],'o--',color=colors[si])

Fd = np.transpose(Fd)
Grad = []
def flin(x,k,b):
    return k*x+b
for de,fdx in zip(De,Fd):
    #popt,pconv = curve_fit(flin,Dxsa[4:7],fdx[4:7])
    if fdx[4]*fdx[5]>0:
        ax4.plot([de,],[0,],'rx')
        Grad.append(0)
    else:
        grad = (fdx[5]-fdx[4])/(Dxsa[5]-Dxsa[4])/(2*np.pi/k)
        ax4.plot([de,],[grad,],'b.')
        Grad.append(grad)
#ax4.plot(De,Grad,'b.')
ax2.plot([Dxsa[0],Dxsa[-1]],[0,0],'b--')

ax2.set_xlabel(r'$\Delta x_h/\lambda$')
ax2.set_ylabel(r'$\Delta F_x$')
ax4.set_xlabel(r'De')
ax4.set_ylabel(r'$\nabla F_x$')
ax3.set_xlabel(r'De')
ax3.set_ylabel(r'$\Delta F_x$')
ax2.set_ylim((-120,300))
ax2.legend(loc='best',ncol=3)
ax3.legend(loc='best',ncol=3)
plt.show()

