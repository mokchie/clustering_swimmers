from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib
import non_dimensionalize as ndm
matplotlib.rcParams.update({'font.size': 13})
colors = ['b','k','r','g','y','c','m']*5
linestyles = ['r.','g.','y.','m.','b.','b+','g+','y+','m+','k+']
def flin(x,k,b):
    return k*x+b
fig,ax1 = plt.subplots(1,1)
b = 1.8 #approximately the beating amplitude is 1.8, but it is actually dependent to the beating particle mass and De number
mu = 100
mass = [1.0,0.5,0.25,0.125,0.0625,0.03125]
directs = ['mass'+str(i) for i in mass]
for cn,dir in enumerate(directs):
    infiles = ['in1.run','in2.run','in3.run','in4.run','in5.run']
    Infiles = [dir+'/'+i for i in infiles]
    De = []
    Vx = []
    Vf = []
    TAU = []

    for infile in Infiles:
        rho = readvar(infile,'rho0')
        dt = readvar(infile,'Dt')
        omega = readvar(infile,'omega')
        k = readvar(infile,'k')
        period = -2.0*np.pi/omega
        fn = readvar(infile,'fn')
        Ly = readvar(infile,'Ly')
        tau = readvar(infile,'tau')
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
            dc,lc = readaverows(filename)
            da = readave(vxfile)
            Time = np.array(dc['TimeStep'])*dt
            dstep = dc['TimeStep'][1]-dc['TimeStep'][0]
            Xc = np.array([d['c_mc1'][0] for d in lc])
            Yc = np.array([d['c_mc1'][1] for d in lc])
            Tp = round(period/dt/dstep)
            Tp = 1
            Time = np.array(coarseave(Time,Tp))
            Xc = np.array(coarseave(Xc,Tp))
            Yc = np.array(coarseave(Yc,Tp))
            Xc -= Xc[0]
            Yc -= Yc[0]
            Dc = np.sqrt(Xc**2+Yc**2)
            de = -tau[i]*omega
            De.append(de)
            st = int(len(Time)/3*2)
            popt,pconv = curve_fit(flin,Time[st:],Dc[st:])
            Vx.append(popt[0]+np.average(da['c_cvxf']))
            Vf.append(popt[0])
            TAU.append(tau[i])
    De,Vx = np.transpose(sorted(zip(De,Vx),key=lambda tup: tup[0]))
    mask = De<35
    De = De[mask]
    Vx = Vx[mask]
    Re = np.abs(omega)/2/np.pi*b**2*rho*mass[cn]/mu
    ax1.plot(De,ndm.ndm_v(Vx)*ndm.ndm_k(k)/ndm.ndm_omega(np.abs(omega)),'.-',color=colors[cn],label=r'$\rho/\rho_\mathrm{ref}=%s,\mathrm{Re}\approx%.3f$'%(ndm.ndm_rho(mass[cn]*rho),Re),markersize=8,linewidth=1)
#    ax1.plot(De,Vf,'--',color=colors[cn])
ax1.legend(loc='best',ncol=1)
ax1.set_xlabel('De')
ax1.set_ylabel('$Vk/\omega$')
ax1.set_ylim((0.018,0.03))
plt.show()
