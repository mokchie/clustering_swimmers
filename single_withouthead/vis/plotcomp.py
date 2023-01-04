from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
import matplotlib
import non_dimensionalize as ndm
font = {'size':11}
matplotlib.rcParams.update({'font.size': 13})
colors = ['b','r','g','y','c','m']*5
linestyles1 = ['r.','g.','y.','m.','b.','b+','g+','y+','m+','k+']
kappa = 2811.0
mu = 100.0
def flin(x,k,b):
    return k*x+b
fig,ax1 = plt.subplots(1,1)
Omega = [0.2,0.3,0.4,0.6,1.0]
legend_elements1 = []
legend_elements2 = []
####################### without head ################################

directs = ['tau_omega'+str(i) for i in Omega]
for cn,dir in enumerate(directs):
    infiles = ['in1.run','in2.run','in3.run','in4.run','in5.run']
    Infiles = [dir+'/'+i for i in infiles]
    De = []
    Vx = []
    Vf = []
    TAU = []
    for infile in Infiles:
        dt = readvar(infile,'Dt')
        omega = readvar(infile,'omega')
        k = readvar(infile,'k')
        period = -2.0*np.pi/omega
        fn = readvar(infile,'fn')
        Ly = readvar(infile,'Ly')
        Ls = readvar(infile,'Ls')
        tau = readvar(infile,'tau')
        files = []
        vxfiles = []
        for ta in tau:
            files.append(dir+'/swimmer1-pos-%s.data'%myround(ta))
            vxfiles.append(dir+'/vxf-%s.data'%myround(ta))

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
            st = int(len(Time)/2)
            popt,pconv = curve_fit(flin,Time[st:],Dc[st:])
            vxt = da['c_cvxf']
            Vx.append(popt[0]+np.average(vxt))
            Vf.append(popt[0])
            TAU.append(tau[i])
    De,Vx,Vf = np.transpose(sorted(zip(De,Vx,Vf), key=lambda tup: tup[0]))
    Sp = (mu*Omega[cn]/kappa)**(1.0/3)*Ls/3.0/np.pi
    ax1.plot(De,ndm.ndm_v(Vx)/ndm.ndm_omega(Omega[cn])*ndm.ndm_k(k),'.-',color=colors[cn],markersize=8,linewidth=1)
    cl = Line2D([0],[0],color=colors[cn],label=r'$\omega/\omega_\mathrm{ref}=%.1f, \mathrm{Sp}=%.2f$'%(ndm.ndm_omega(Omega[cn]),Sp))
    legend_elements1.append(cl)
#    ax1.plot(De,Vf,'--',color=colors[cn])

###################### with head ########################
#colors = ['r','c','g','c','m']*5

Omega = [0.2,0.3,0.4,0.6,1.0]#[0.3,1.0]
directs = ['../../single_withhead/var/vis/m0.125/tau_omega'+str(i) for i in Omega]
for cn,dir in enumerate(directs):
    infiles = ['in1.run','in2.run','in5.run','in3.run','in4.run']
    Infiles = [dir+'/'+i for i in infiles]
    De = []
    Vx = []
    Vf = []
    TAU = []

    for infile in Infiles:
        dt = readvar(infile,'Dt')
        k = readvar(infile,'k')
        omega = readvar(infile,'omega')
        period = -2.0*np.pi/omega
        fn = readvar(infile,'fn')
        Ly = readvar(infile,'Ly')
        Ls = readvar(infile,'Ls')
        tau = readvar(infile,'tau')
        files = []
        vxfiles = []
        for ta in tau:
            files.append(dir+'/swimmer1-pos-%s.data'%myround(ta))
            vxfiles.append(dir+'/vxf-%s.data'%myround(ta))
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
            st = int(len(Time)/2)
            popt,pconv = curve_fit(flin,Time[st:],Dc[st:])
            Vx.append(popt[0]+np.average(da['c_cvxf']))
            Vf.append(popt[0])
            TAU.append(tau[i])
    De,Vx = np.transpose(sorted(zip(De,Vx),key=lambda tup: tup[0]))
    Sp = (mu*Omega[cn]/kappa)**(1.0/3)*Ls/3.0/np.pi
    ax1.plot(De,ndm.ndm_v(Vx)/ndm.ndm_omega(Omega[cn])*ndm.ndm_k(k),'*-',color=colors[cn],markersize=8,linewidth=1)
cl1 = Line2D([0],[0],color='w',marker='.',markerfacecolor='k',label=r'without a head',markersize=12)
cl2 = Line2D([0],[0],color='w',marker='*',markerfacecolor='k',label=r'with a head',markersize=12)
legend_elements2.append(cl1)
legend_elements2.append(cl2)
lg1 = ax1.legend(handles=legend_elements1,loc='upper right',ncol=1)
lg2 = ax1.legend(handles=legend_elements2,loc='lower right',ncol=1)
ax1.add_artist(lg1)
ax1.set_xlabel('De')
ax1.set_ylabel(r'$Vk/\omega$')
ax1.set_ylim((0,0.08))
plt.show()
