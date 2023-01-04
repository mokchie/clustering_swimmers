from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib
import non_dimensionalize as ndm
matplotlib.rcParams.update({'font.size': 13})
infiles = ['in1.run','in2.run','in3.run','in4.run','in5.run']
colors = ['b-','k-','r-','g-','y-','c-','m-','b-','k-.','r-.','g-.','y-.','c-.','m-.']*3
linestyles = ['b.-','k.-','r.-','g.-','y.-','c.-','rx-','gx-','yx-','mx-','bx-','cx-']*2
def flin(x,k,b):
    return k*x+b
fig,ax1 = plt.subplots(1,1)
fig2,ax2 = plt.subplots(1,1)

Dirs = ['omega0.4','omega0.6','omega0.8','omega1.2','omega1.6']
De_all = []
for cn,dir in enumerate(Dirs):
    De = []
    Vx = []
    TAU = []
    for infile in [dir+'/'+i for i in infiles]:
        dt = readvar(infile,'Dt')
        b = readvar(infile,'b')
        k = readvar(infile,'k')
        omega = readvar(infile,'omega')
        period = -2.0*np.pi/omega
        fn = readvar(infile,'fn')
        Ly = readvar(infile,'Ly')
        tau = readvar(infile,'tau')
        etas = readvar(infile,'mu')
        Np = readvar(infile,'Np')
        kBT = readvar(infile,'T')
        rho0 = readvar(infile,'rho0')
        etap = tau*rho0*Np*kBT
        eta = etap+etas
        print(eta,etap,tau,Np)
        files = []
        vxfiles = []
        for ta in tau:
            if ta<1:
                files.append(dir+'/swimmer1-pos-%.1f.data'%ta)
                vxfiles.append(dir+'/vxf-%.1f.data'%ta)
            else:
                files.append(dir+'/swimmer1-pos-%.0f.data'%ta)
                vxfiles.append(dir+'/vxf-%.0f.data'%ta)
        for i,(filename,vxfile) in enumerate(zip(files,vxfiles)):
            da = readave(vxfile)
            Time = np.array(da['TimeStep'])*dt
            vxt = da['c_cvxf']
            ax2.plot(Time,vxt)
            vxt = vxt[int(len(vxt)/2):]
            Vx.append(np.average(vxt))
            de = -tau[i]*omega
            De.append(de)
            TAU.append(tau[i])
    De,Vx = np.transpose(sorted(zip(De,Vx),key=lambda tup:tup[0]))
    mask = De<50
    De = De[mask]
    Vx = Vx[mask]
    De_all.append(De)
    ax1.plot(De,ndm.ndm_v(Vx)*ndm.ndm_k(k)/ndm.ndm_omega(np.abs(omega)),linestyles[cn],markersize=8,linewidth=1,label=r'$\omega/\omega_\mathrm{ref}=%.1f$'%(-ndm.ndm_omega(omega),))
    #ax1.plot(De,Vx*k/np.abs(omega)/(1+etas/eta*De**2),linestyles[cn],markersize=8,linewidth=1,label=r'$\omega/\omega_\mathrm{ref}=%.1f$'%(-ndm.ndm_omega(omega),))
def Unn(de,omega,b,k,etas,eta):
    return 1/2.0*b**2*k**2*(1+etas/eta*de**2)/(1+de**2)
    #return 1/2.0*b**2*k**2*(1)/(1+de**2)
    
X = np.linspace(np.min(list(map(np.min,De_all))),np.max(list(map(np.max,De_all))),200)
Y = Unn(X,1,b,k,etas,eta)
ax1.plot(X,Y,'k--',label="theory Eq.(14)",linewidth=1.2)
ax1.legend(loc='best',ncol=2)
ax1.set_xlabel('De')
ax1.set_ylabel(r'$Vk/\omega$')
ax1.set_xlim((0,35))
plt.show()
