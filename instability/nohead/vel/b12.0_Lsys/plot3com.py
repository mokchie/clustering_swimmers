from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
infiles = ['in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run','in-7.run','in-8.run',]
colors = ['b','r','g','y','c','m','k']*10
linestyles = ['o-','s--']*10
labels = ['single','pair']
pointstyles = ['.','x']*10
def flin(x,k,b):
    return k*x+b
fig,(ax1,ax2) = plt.subplots(2,1)
Ns = [11,20]
Vs = []
Vp = []
Om = []
Dirs = ['single_w0.4','pair_w0.4','single_w0.6','pair_w0.6','single_w0.8','pair_w0.8','single_w1.0','pair_w1.0','single_w1.2','pair_w1.2',]#'single_w1.0','pair_w1.0','single_w1.2','pair_w1.2',]
V_list = []
for cn,direct in enumerate(Dirs):
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
                    files.append(direct+'/data/swimmer1-pos-%s.data'%myround(ta))
                    vxfiles.append(direct+'/data/vxf-%s.data'%myround(ta))
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
                Tau.append(tau[i])
        except FileNotFoundError as err:
            print(err)
    De,Vx,Tau = np.transpose(sorted(zip(De,Vx,Tau),key=lambda tup: tup[0]))
    V_list.append(Vx)
    if cn%2==0:
        Om.append(-omega)
        Vs.append(Vx)
    else:
        Vp.append(Vx)
    if cn%2==1:
        ax1.plot(De,Vx-V_list[-2],linestyles[cn%2],color=colors[int(cn/2)],label=r'$\omega=%s$'%(myround(-omega),))
    #ax1.plot(De,Vx,linestyles[cn%2],color=colors[int(cn/2)],label=r'%s $\omega=%s$'%(labels[cn%2],myround(-omega)))
Vs = np.array(Vs)
Vp = np.array(Vp)
for i,ns in enumerate(Ns):
    ax2.plot(Om,Vs[:,ns],'o-',color=colors[i],label=r'single $\tau=%s$'%myround(Tau[ns]))
    ax2.plot(Om,Vp[:,ns],'s-',color=colors[i],label=r'pair $\tau=%s$'%myround(Tau[ns]))
ax1.legend(loc='best',ncol=3)
ax1.set_xlabel(r'De')
ax1.set_ylabel(r'$\Delta$V')
ax2.legend(loc='best',ncol=2)
ax2.set_xlabel(r'$\omega$')
ax2.set_ylabel('V')
print('tau = ',Tau[Ns])
plt.show()
