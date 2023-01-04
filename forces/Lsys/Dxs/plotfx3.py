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
#fig1,ax1 = plt.subplots(1,1)
#fig2,ax3 = plt.subplots(1,1)
fig3,ax2 = plt.subplots(1,1)
fig4,ax4 = plt.subplots(1,1)
Dxs = ['-4_16','-3_16','-2_16','-1_16','0_16','1_16','2_16','3_16','4_16']
Dxsa = np.array([-4/16,-3/16,-2/16,-1/16,0/16,1/16,2/16,3/16,4/16,])
dirs = ['dxs'+d for d in Dxs]
colors = ['b','r','g','y','m','c','k','C0','C1','C2','C3','C4']*10
#Sm = [0,1,2,3,5,7,10,15,25]
Sm = [0,2,5,10,18,25]
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
    #ax1.plot(De,Fx1,'.-',color=colors[di],label=r'$\Delta x_h=%s\lambda$'%Dxsa[di])
    #ax1.plot(De,Fx2,'.--',color=colors[di])
    #ax3.plot(De,Fx1-Fx2,'.-',color=colors[di],label=r'$\Delta x_h=%s\lambda$'%Dxsa[di])
    Fd.append(Fx1-Fx2)
    for si,sm in enumerate(Sm):
        Fxc1[si].append(Fx1[sm])
        Fxc2[si].append(Fx2[sm])

Fxc1 = np.array(Fxc1)
Fxc2 = np.array(Fxc2)

for si,sm in enumerate(Sm):
    ax2.plot(Dxsa*2*np.pi,ndm.ndm_force(Fxc1[si])*ndm.ndm_k(k)/ndm.ndm_omega(omega)/ndm.ndm_eta(eta),'.-',color=colors[si],label=r'De$=%s$'%myround(De[sm]),markersize=8,linewidth=1)
    #ax2.plot(Dxsa,Fxc2[si],'o--',color=colors[si])

Fd = np.transpose(Fd)
Grad = []
def flin(x,k,b):
    return k*x+b
for de,fdx in zip(De,Fd):
    #popt,pconv = curve_fit(flin,Dxsa[4:7],fdx[4:7])
    if fdx[4]*fdx[5]>0:
        ax4.plot([de,],[0,],'rx',markersize=8)
        Grad.append(0)
    else:
        grad = (fdx[5]-fdx[4])/(Dxsa[5]-Dxsa[4])/(2*np.pi/k)
        ax4.plot([de,],[ndm.ndm_grad_force(grad),],'b.',markersize=8)
        Grad.append(grad)
#ax4.plot(De,Grad,'b.')
ax2.plot([2*np.pi*Dxsa[0],2*np.pi*Dxsa[-1]],[0,0],'k--')




# Dxs = ['-4_16','-3_16','-2_16','-1_16','-0.5_16','0_16','0.5_16','1_16','2_16','3_16','4_16']
# labels = ['-8/32','-6/32','-4/32','-2/32','-1/32','0','1/32','2/32','4/32','6/32','8/32']
# Dxsa = [-4/16,-3/16,-2/16,-1/16,-0.5/16,0/16,0.5/16,1/16,2/16,3/16,4/16,]
# dirs = ['../Dxs_db0.0/dxs'+d for d in Dxs]
# colors = ['b','r','g','y','m','c','k','C0','C1','C2','C3','C4']*10
# Sm = [0,1,2,3,5,7,10,15,25]
# Fd = []
# Fxc1 = [[] for i in range(len(Sm))]
# Fxc2 = [[] for i in range(len(Sm))]
# for di,direct in enumerate(dirs):
#     infiles = [direct+'/'+i for i in ['in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run','in-7.run','in-8.run','in-9.run']]
#     Tau = []
#     De = []
#     F = [[],[]]
#     for infile in infiles:
#         Dt = readvar(infile,'Dt')
#         taus = readvar(infile,'tau')
#         k = readvar(infile,'k')
#         omega = -readvar(infile,'omega')
#         for tau in taus:
#             Tau.append(tau)
#             De.append(tau*omega)
#             linestyles = ['-','--']
#             file_list = zip([direct+'/data/fs1-%s.data'%myround(tau), direct+'/data/fs2-%s.data'%myround(tau)],['c_cfs1','c_cfs2'])
#             for cn,(f_fs,c_cf) in enumerate(file_list):
#                 dst,ds = readrow(f_fs)

#                 t = coarse(np.array(dst['TimeStep'])*Dt)

#                 fs = np.array(ds[c_cf])
#                 fx = fs[:,0]
#                 fy = fs[:,1]
#                 F[cn].append(np.average(fx[int(len(fx)/3.0):]))
#     Fx1 = np.array(F[0])
#     Fx2 = np.array(F[1])
#     De,Fx1,Fx2,Tau = np.transpose(sorted(zip(De,Fx1,Fx2,Tau)))
#     Fd.append(Fx1-Fx2)
#     for si,sm in enumerate(Sm):
#         Fxc1[si].append(Fx1[sm])
#         Fxc2[si].append(Fx2[sm])

# Fxc1 = np.array(Fxc1)
# Fxc2 = np.array(Fxc2)

# Fd = np.transpose(Fd)
# Grad = []
# def flin(x,k,b):
#     return k*x+b
# for de,fdx in zip(De,Fd):
#     popt,pconv = curve_fit(flin,Dxsa[5:8],fdx[5:8])
#     Grad.append(popt[0])
# Grad = np.array(Grad)/(2*np.pi/k)
# ax4.plot(De,Grad,'r.')

ax2.set_xlabel(r'$kr_h$')
ax2.set_ylabel(r'$F_hk/\omega/\eta$')
ax4.set_xlabel(r'De')
ax4.set_ylabel(r'$\nabla F_h r_c/F_\mathrm{ref}$')
ax2.legend(loc='best',ncol=2)
ax2.set_ylim((-0.3,0.8))
plt.show()

