from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
infiles = ['in.run',]
colors = ['b','r','g','y','c','m','k']*10
linestyles = ['-','--','-.']
pointstyles = ['.','x','+']
def flin(x,k,b):
    return k*x+b
fig1,ax1 = plt.subplots(1,1)
fig2,ax2 = plt.subplots(1,1)
cn = 0
De = []
TAU = []

for infile in infiles:
    dt = readvar(infile,'Dt')
    omega = readvar(infile,'omega')
    period = -2.0*np.pi/omega
    fn = readvar(infile,'fn')
    Ly = readvar(infile,'Ly')
    tau = readvar(infile,'tau')
    files1 = []
    files2 = []
    vxfiles = []
    for fni in fn:
        files1.append('data/swimmer1-pos-%s.data'%fni)
        files2.append('data/swimmer2-pos-%s.data'%fni)
        vxfiles.append('data/vxf-%s.data'%fni)
    for i,(filename1,filename2,vxfile) in enumerate(zip(files1,files2,vxfiles)):
        da = readave(vxfile)
        XC = []
        YC = []
        for j,filename in enumerate((filename1,filename2)):
            dc,lc = readaverows(filename)
            Time = np.array(dc['TimeStep'])*dt
            dstep = dc['TimeStep'][1]-dc['TimeStep'][0]
            Xc = np.array([d['c_mc%d'%(j+1)][0] for d in lc])
            Yc = np.array([d['c_mc%d'%(j+1)][1] for d in lc])
            Tp = round(period/dt/dstep)
            Tp = 1
            Time = np.array(coarseave(Time,Tp))
            Xc = np.array(coarseave(Xc,Tp))
            Yc = np.array(coarseave(Yc,Tp))
            XC.append(Xc)
            YC.append(Yc)
            Xc -= Xc[0]
            Yc -= Yc[0]
            Dc = np.sqrt(Xc**2+Yc**2)
            de = -tau[i]*omega/2/np.pi
            De.append(de)
            st = int(len(Time)/2)
            ax1.plot(Xc,Yc,linestyles[j],color=colors[cn],label=r'$\mathrm{De}=%.2f$'%de)
            #popt,pconv = curve_fit(flin,Time[st:],Dc[st:])
            #ax1.plot(Time[st:],flin(Time[st:],popt[0],popt[1]),'-.',color=colors[cn])
            Time2 = (Time[0:-1]+Time[1:])/2
            V = np.sqrt((Xc[1:]-Xc[0:-1])**2+(Yc[1:]-Yc[0:-1])**2)/(Time[1:]-Time[0:-1])
            Time2 = coarseave(Time2,4000)
            V = coarseave(V,4000)
            ax2.plot(Time2,V,linestyles[j],color=colors[cn],label=r'$\mathrm{De}=%.2f$'%de)
        cn+=1

ax1.legend(loc='best',ncol=4)
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax2.set_xlabel('time')
ax2.set_ylabel('vel')
plt.show()
