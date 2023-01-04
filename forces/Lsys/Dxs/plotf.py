from variable3 import *
from readvel3 import *
import pdb
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
def coarse(lst):
    return np.array(coarseave(lst,50))
fig,(ax1,ax2) = plt.subplots(2,1)
Dxs = ['-4_16','-3_16','-2_16','-1_16','0_16','1_16','2_16','3_16','4_16']
Dxsa = [-4/16,-3/16,-2/16,-1/16,0/16,1/16,2/16,3/16,4/16,]
dirs = ['dxs'+d for d in Dxs]
colors = ['b','r','g','y','m','c','k','C0','C1','C2','C3']*10
Sm = [0,1,3,5,10,15,25]
Fyc1 = [[] for i in range(len(Sm))]
Fyc2 = [[] for i in range(len(Sm))]
for di,direct in enumerate(dirs):
    infiles = [direct+'/'+i for i in ['in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run','in-7.run','in-8.run','in-9.run']]
    Tau = []
    De = []
    F = [[],[]]
    for infile in infiles:
        Dt = readvar(infile,'Dt')
        taus = readvar(infile,'tau')
        omega = -readvar(infile,'omega')
        for tau in taus:
            Tau.append(tau)
            De.append(tau*omega)
            linestyles = ['-','--']
            file_list = zip([direct+'/data/fs1-%s.data'%myround(tau), direct+'/data/fs2-%s.data'%myround(tau)],['c_cfs1','c_cfs2'])
            for cn,(f_fs,c_cf) in enumerate(file_list):
                dst,ds = readrow(f_fs)

                t = coarse(np.array(dst['TimeStep'])*Dt)

                fs = np.array(ds[c_cf])
                fx = fs[:,0]
                fy = fs[:,1]
                F[cn].append(np.average(fy[int(len(fy)/3.0):]))
    Fy1 = F[0]
    Fy2 = F[1]
    De,Fy1,Fy2,Tau = np.transpose(sorted(zip(De,Fy1,Fy2,Tau)))
    ax1.plot(De,Fy1,'o-',color=colors[di],label='dl=%s'%Dxs[di])
    ax1.plot(De,Fy2,'o--',color=colors[di])
    for si,sm in enumerate(Sm):
        Fyc1[si].append(Fy1[sm])
        Fyc2[si].append(Fy2[sm])
Fyc1 = np.array(Fyc1)
Fyc2 = np.array(Fyc2)
ax1.set_xlabel('De')
ax1.set_ylabel('$F_y$')
ax1.legend(loc='best',ncol=3)
for si,sm in enumerate(Sm):
    ax2.plot(Dxsa,Fyc1[si]-Fyc2[si],'o-',color=colors[si],label=r'De$=%s$'%myround(De[sm]))
    #ax2.plot(Dxsa,Fyc2[si],'o--',color=colors[si])

ax2.set_xlabel(r'$d/\lambda$')
ax2.set_ylabel('$F_y$')
ax2.legend(loc='best',ncol=3)
plt.show()

