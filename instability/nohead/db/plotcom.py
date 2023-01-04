from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
from scipy.optimize import curve_fit
import matplotlib
import non_dimensionalize as ndm
matplotlib.rcParams.update({'font.size':13})

colors = ['b','r','g','y','c','m','k']*10
linestyles = ['.-','.--','.-.']*10
pointstyles = ['.','x','+','v','<','s','*']*10
def cstyle(i,colors=colors,linestyles=linestyles):
    c = i%len(colors)
    p = int(i/len(colors))
    p = p%len(linestyles)
    return colors[c]+linestyles[p]

def flin(x,k,b):
    return k*x+b

infiles = ['in-%d.run'%i for i in range(1,5)]
IND = [8.0,8.25,8.5,8.75,9.0,9.25,9.5,9.75,10.0,10.25,10.5,10.75,11.0,11.125,11.25,11.5,11.75,12.0]+[10.875,11.0625]
Dirs = ['db'+str(i) for i in IND]

Db = []
TAU = []
Syn = []
for cn,dir in enumerate(Dirs):
    for infile in infiles:
        infile = dir+'/'+infile
        dt = readvar(infile,'Dt')
        omega = readvar(infile,'omega')
        b1 = readvar(infile,'b1')
        b2 = readvar(infile,'b2')

        print('dir:',dir)
        k = readvar(infile,'k')
        period = -2.0*np.pi/omega
        fn = readvar(infile,'fn')
        Ly = readvar(infile,'Ly')
        Tau = readvar(infile,'tau')
        for fni,tau in zip(fn,Tau):
            Db.append((b2-b1))
            TAU.append(tau)
            if tau<5.0: 
                filename1 = dir+'/data/swimmer1-pos-%.1f.data'%tau
                filename2 = dir+'/data/swimmer2-pos-%.1f.data'%tau
                vxfile = dir+'/data/vxf-%.1f.data'%tau
            else:
                filename1 = dir+'/data/swimmer1-pos-%d.data'%tau
                filename2 = dir+'/data/swimmer2-pos-%d.data'%tau
                vxfile = dir+'/data/vxf-%d.data'%tau
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
                Dc = np.sqrt((Xc-Xc[0])**2+(Yc-Yc[0])**2)
                de = tau*omega
                st = int(len(Time)/2)
                #ax2.plot(Xc,Yc,linestyles[j],color=colors[cn])
            ds = np.sqrt((XC[0]-XC[1])**2+(YC[0]-YC[1])**2)
            #ax1.plot(Time,ds,color=colors[cn],label=r'$\Delta b=%.2f$'%(np.abs(b1-b2)/(b1+b2)*2*100.0,)+'%')
            if ds[-1]>Ly/2.0:
                Syn.append(0)
            else:
                Syn.append(1)
            
fig,ax = plt.subplots(1,1)
De = -np.array(TAU)*omega
for x,y,z in zip(De,Db,Syn):
    if z>0:
        ax.plot(x,ndm.ndm_angle(np.array(y)),'b.')
    else:
        ax.plot(x,ndm.ndm_angle(np.array(y)),'rx')

ax.set_xlabel(r'De')
ax.set_ylabel(r'$\Delta \theta_b/\pi$')
plt.savefig('figure.png')
plt.show()
