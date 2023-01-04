from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib
import non_dimensionalize as ndm
matplotlib.rcParams.update({'font.size': 13})
colors = ['b','r','g','y','c','m','k']*10
linestyles = ['r.','g.','y.','m.','b.','b+','g+','y+','m+','k+']
def flin(x,k,b):
    return k*x+b
fig,ax1 = plt.subplots(1,1)

directs = ['omega_Ls'+str(i) for i in (20,22,24,26,28,30,32)]
for cn,direct in enumerate(directs):
    infiles = [direct+'/'+i for i in ['in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run']]
    De = []
    Vx = []
    Vf = []
    Omega = []
    for infile in infiles:
        dt = readvar(infile,'Dt')
        omega = -readvar(infile,'omega')
        b = readvar(infile,'b')
        period = 2.0*np.pi/omega
        fn = readvar(infile,'fn')
        Ly = readvar(infile,'Ly')
        Ls = readvar(infile,'Ls')
        tau = readvar(infile,'tau')
        files = []
        vxfiles = []
        for n in fn:
            files.append(direct+'/data/swimmer1-pos-%d.data'%n)
            vxfiles.append(direct+'/data/vxf-%d.data'%n)
        for i,(filename,vxfilename) in enumerate(zip(files,vxfiles)):
            Omega.append(omega[i])
            da = readave(vxfilename)
            dc,lc = readaverows(filename)
            Time = np.array(dc['TimeStep'])*dt
            dstep = dc['TimeStep'][1]-dc['TimeStep'][0]
            Xc = np.array([d['c_mc1'][0] for d in lc])
            Yc = np.array([d['c_mc1'][1] for d in lc])
            Tp = round(period[i]/dt/dstep)
            Tp = 1
            Time = np.array(coarseave(Time,Tp))
            Xc = np.array(coarseave(Xc,Tp))
            Yc = np.array(coarseave(Yc,Tp))
            Xc -= Xc[0]
            Yc -= Yc[0]
            Dc = np.sqrt(Xc**2+Yc**2)
            de = tau*omega[i]
            De.append(de)
            st = int(len(Time)/3)
            popt,pconv = curve_fit(flin,Time[st:],Dc[st:])
            if Xc[-1]<Xc[0]:
                Vx.append(popt[0]+np.average(da['c_cvxf']))
                Vf.append(popt[0])
            else:
                Vx.append(-(popt[0]-np.average(da['c_cvxf'])))
                Vf.append(-popt[0])

    ax1.plot(ndm.ndm_omega(np.array(Omega)),ndm.ndm_v(np.array(Vx)),'.-',color=colors[cn],label=r'$L/r_c=%.1f$'%ndm.ndm_length(Ls),markersize=8,linewidth=1)

ax1.legend(loc='best',ncol=2)
ax1.set_xlabel(r'$\omega/\omega_\mathrm{ref}$')
ax1.set_ylabel(r'$V/v_\mathrm{ref}$')
plt.show()
