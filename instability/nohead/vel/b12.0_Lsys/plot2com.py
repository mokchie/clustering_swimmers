from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
import matplotlib
import non_dimensionalize as ndm
matplotlib.rcParams.update({'font.size': 13})
font = {'size':13}
infiles = ['in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run','in-7.run','in-8.run',]
colors = ['b','k','r','g','y','c','m']*10
colors2 = ['b','r','g','y','c','m']*10
linestyles = ['.-','*-']*10
labels = ['single','pair']
pointstyles = ['.','*']*10
def flin(x,k,b):
    return k*x+b
mu = 100
kappa = 2188
fig1,ax1 = plt.subplots(1,1)
fig2,ax2 = plt.subplots(1,1)
legend_elements1 = []
legend_elements2 = []
legend_elements3 = []
legend_elements4 = []
Ns = [11,20]
Vs = []
Vp = []
Om = []
Dirs = ['single_w0.4','pair_w0.4','single_w0.6','pair_w0.6','single_w0.8','pair_w0.8','single_w1.0','pair_w1.0','single_w1.2','pair_w1.2',]#'single_w1.0','pair_w1.0','single_w1.2','pair_w1.2',]
#Dirs = ['single_w0.4','single_w0.6','single_w0.8','single_w1.0','single_w1.2',]
for cn,direct in enumerate(Dirs):
    De = []
    Vx = []
    Vf = []
    Tau = []
    for infile in [direct+'/'+i for i in infiles]:
        try:
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
    if cn%2==0:
        Om.append(-omega)
        Vs.append(Vx)
    else:
        Vp.append(Vx)
    if cn in [2,3,8,9]:
        continue
    ax1.plot(De,ndm.ndm_v(Vx)*ndm.ndm_k(k)/ndm.ndm_omega(np.abs(omega)),linestyles[cn%2],color=colors[int(cn/2)],label=r'$\omega/\omega_\mathrm{ref}=%.1f$'%(myround(ndm.ndm_omega(-omega)),),markersize=8,linewidth=1)
    cl = Line2D([0],[0],color=colors[int(cn/2)],label=r'$\omega/\omega_\mathrm{ref}=%.1f$'%(myround(ndm.ndm_omega(-omega)),))
    if cn%2==0:
        legend_elements1.append(cl)
cl1 = Line2D([0],[0],color='w',marker='.',markerfacecolor='k',label=r'single',markersize=12)
cl2 = Line2D([0],[0],color='w',marker='*',markerfacecolor='k',label=r'pair',markersize=12)
legend_elements2.append(cl1)
legend_elements2.append(cl2)
        
Vs = np.array(Vs)
Vp = np.array(Vp)
for i,ns in enumerate(Ns):
    ax2.plot(ndm.ndm_omega(np.array(Om)),ndm.ndm_v(np.array(Vs[:,ns])),'.-',color=colors2[i],label=r'single $\tau/t_\mathrm{ref}=%.1f$'%myround(ndm.ndm_time(Tau[ns])),markersize=8,linewidth=1)
    ax2.plot(ndm.ndm_omega(np.array(Om)),ndm.ndm_v(np.array(Vp[:,ns])),'*-',color=colors2[i],label=r'pair $\tau/\t_mathrm{ref}=%.1f$'%myround(ndm.ndm_time(Tau[ns])),markersize=8,linewidth=1)
    cl = Line2D([0],[0],color=colors2[i],label=r'$\tau/t_\mathrm{ref}=%.1f$'%myround(ndm.ndm_time(Tau[ns])))
    legend_elements3.append(cl)

legend_elements4.append(cl1)
legend_elements4.append(cl2)
    
lg1 = ax1.legend(handles=legend_elements1,loc='upper left',ncol=2)
lg2 = ax1.legend(handles=legend_elements2,loc='lower right',ncol=1)
lg3 = ax2.legend(handles=legend_elements3,loc='lower right',ncol=1)
lg4 = ax2.legend(handles=legend_elements4,loc='upper left',ncol=1)
ax1.add_artist(lg1)
ax2.add_artist(lg3)
#ax1.legend(loc='upper left',ncol=2)
ax1.set_xlabel(r'De')
ax1.set_ylabel(r'$Vk/\omega$')
#ax2.legend(loc='best',ncol=2)
ax2.set_xlabel(r'$\omega/\omega_\mathrm{ref}$')
ax2.set_ylabel(r'$V/v_\mathrm{ref}$')
ax1.set_ylim((0.01,0.07))
#ax2.set_ylim((0.05,0.12))
#ax1.text(45,0.09,'(a)',fontdict=font)
#ax2.text(1.1,0.09,'(b)',fontdict=font)
print('tau = ',Tau[Ns])
plt.show()
