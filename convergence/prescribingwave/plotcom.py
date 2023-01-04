from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib
matplotlib.rcParams.update({'font.size':14,'font.family':'sans-serif'})

colors = ['b','r','g','m','c']
linestyles = ['r.','g.','y.','m.','b.','b+','g+','y+','m+','k+']
def flin(x,k,b):
    return k*x+b
fig,ax1 = plt.subplots(1,1)
cmap = plt.get_cmap('jet')
Rho = [4.17,6.25,8.33,10.42,12.50]
nlay = [2,3,4,5,6]
dirs = ['rho'+str(i) for i in Rho]
for di,direct in enumerate(dirs):
    infiles = [direct+'/'+i for i in ('in0.run','in1.run','in2.run','in3.run','in4.run')]
    De = []
    Vx = []
    for infile in infiles:
        dt = readvar(infile,'Dt')
        omega = -readvar(infile,'omega')
        Ly = readvar(infile,'Ly')
        tau = readvar(infile,'tau')
        k = readvar(infile,'k')
        files = []
        vxfiles = []
        for ta in tau:
            files.append(direct+'/swimmer1-pos-%.1f.data'%ta)
            vxfiles.append(direct+'/vxf-%.1f.data'%ta)
        for i,(filename,vxfile) in enumerate(zip(files,vxfiles)):
            da = readave(vxfile)
            Time = np.array(da['TimeStep'])*dt
            vxt = da['c_cvxf']
            vxt = vxt[int(len(vxt)*2/3):]
            Vx.append(np.average(vxt)*k/omega)
            de = tau[i]*omega
            De.append(de)
    ax1.plot(De,Vx,'.-',color=colors[di],label=r'$N_\mathrm{lay}=%s,d_0=%s$'%(nlay[di],Rho[di]))
ax1.legend(loc='best',ncol=1)
ax1.set_xlabel('De')
ax1.set_ylabel(r'$Vk/\omega$')
ax1.set_ylim((0.018,0.035))
plt.show()
