from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
infiles = ['in-1.run','in-2.run','in-3.run','in-4.run']
colors = ['b-','k-','r-','g-','y-','c-','m-','b-','k-.','r-.','g-.','y-.','c-.','m-.']*3
linestyles = ['r.','g.','y.','m.','b.','b+','g+','y+','m+','k+']
def flin(x,k,b):
    return k*x+b
fig,ax1 = plt.subplots(1,1)

Dirs = ['single','pair']
for direct in Dirs:
    De = []
    Vx = []
    Vf = []
    Tau = []

    for infile in [direct+'/'+i for i in infiles]:
        try:
            dt = readvar(infile,'Dt')
            omega = readvar(infile,'omega')
            period = -2.0*np.pi/omega
            fn = readvar(infile,'fn')
            Ly = readvar(infile,'Ly')
            tau = readvar(infile,'tau')
            files = []
            vxfiles = []
            for ta in tau:
                if ta<1:
                    files.append(direct+'/data/swimmer1-pos-%.1f.data'%ta)
                    vxfiles.append(direct+'/data/vxf-%.1f.data'%ta)
                else:
                    files.append(direct+'/data/swimmer1-pos-%.0f.data'%ta)
                    vxfiles.append(direct+'/data/vxf-%.0f.data'%ta)
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
                popt,pconv = curve_fit(flin,Time,Dc)
                Vx.append(popt[0]+np.average(da['c_cvxf']))
                Vf.append(popt[0])
                Tau.append(tau[i])
        except FileNotFoundError as err:
            print(err)
    ax1.plot(Tau,Vx,'.-',label=direct)
    #ax2.plot(Tau,Vf,'--')


ax1.legend(loc='best',ncol=4)
ax1.set_xlabel(r'$\tau$')
ax1.set_ylabel('V')
plt.show()
