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
fig,ax1 = plt.subplots(1,1)
cn = 0
De = []
Vx = []
TAU = []

for infile in infiles:
    dt = readvar(infile,'Dt')
    omega = readvar(infile,'omega')
    period = -2.0*np.pi/omega
    Ly = readvar(infile,'Ly')
    tau = readvar(infile,'tau')
    files = []
    vxfiles = []
    for ta in tau:
        if ta<1:
            files.append('swimmer1-pos-%.1f.data'%ta)
            vxfiles.append('vxf-%.1f.data'%ta)
        else:
            files.append('swimmer1-pos-%.1f.data'%ta)
            vxfiles.append('vxf-%.1f.data'%ta)
    for i,(filename,vxfile) in enumerate(zip(files,vxfiles)):
        da = readave(vxfile)
        Time = np.array(da['TimeStep'])*dt
        Vx.append(np.average(da['c_cvxf']))
        de = -tau[i]*omega
        De.append(de)
        TAU.append(tau[i])
        cn+=1

ax1.plot(De,Vx,'.-')
ax1.legend(loc='best',ncol=4)
ax1.set_xlabel('De')
ax1.set_ylabel('Vx')
plt.show()
