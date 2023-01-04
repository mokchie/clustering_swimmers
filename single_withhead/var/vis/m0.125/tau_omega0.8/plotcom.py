from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
infiles = ['in1.run','in2.run','in3.run','in4.run']
colors = ['b-','k-','r-','g-','y-','c-','m-','b-','k-.','r-.','g-.','y-.','c-.','m-.']*3
linestyles = ['r.','g.','y.','m.','b.','b+','g+','y+','m+','k+']
def flin(x,k,b):
    return k*x+b
fig,(ax1,ax2) = plt.subplots(2,1)
cn = 0
De = []
Vx = []
Vf = []
TAU = []

for infile in infiles:
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
            files.append('swimmer1-pos-%.1f.data'%ta)
            vxfiles.append('vxf-%.1f.data'%ta)
        else:
            files.append('swimmer1-pos-%.0f.data'%ta)
            vxfiles.append('vxf-%.0f.data'%ta)
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
        de = -tau[i]*omega/2/np.pi
        De.append(de)
        ax1.plot(Time,Dc,colors[cn],label=r'$\mathrm{De}=%.2f$'%de)
        st = 600
        popt,pconv = curve_fit(flin,Time[st:],Dc[st:])
        ax1.plot(Time[st:],flin(Time[st:],popt[0],popt[1]),'--',color=colors[cn][0])
        Vx.append(popt[0]+np.average(da['c_cvxf']))
        Vf.append(popt[0])
        TAU.append(tau[i])
        cn+=1

ax2.plot(De,Vx,'.-')
ax2.plot(De,Vf,'--')
ax1.legend(loc='best',ncol=4)
ax1.set_xlabel('time')
ax1.set_ylabel('displacement')
plt.show()
