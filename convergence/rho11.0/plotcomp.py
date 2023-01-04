from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.integrate import trapezoid

infile = 'in.run'
fig1,ax1 = plt.subplots(1,1)
fig2,ax2 = plt.subplots(1,1)
omega = readvar(infile,'omega')
y01 = readvar(infile,'y01')
y02 = readvar(infile,'y02')
pdiff = readvar(infile,'pdiff')
dt = readvar(infile,'Dt')
kt = readvar(infile,'k')
Ls = readvar(infile,'Ls')
Ly = readvar(infile,'Ly')
Lx = readvar(infile,'Lx')
dx = readvar(infile,'dx')
tot = readvar(infile,'tot')
Rc = readvar(infile,'rc')
files = ['lammpstrj-%.1f.data'%rc for rc in Rc]
linestyles = ['b.','g.','y.','m.','k.','b+','g+','y+','m+','k+']
def fsin(x,b,k,theta0):
    return b*np.sin(k*x+theta0)
L1error = []
for j,filename in enumerate(files):
    print( 'processing %s ... '%filename )
    Time = []
    Timea = []
    Dpx = []
    Hpx1 = []
    Hpx2 = []
    data = readxyz(filename,[2,3])
    count = 0
    for d in data:
        try:
            timestep = d['timestep']
            X = d['x']
            Y = d['y']
            Typ = d['type']
            ID = d['id']
            X1,Y1,Typ1,ID1 = np.transpose(select(zip(X,Y,Typ,ID), lambda tup: tup[2]==2))
            X2,Y2,Typ2,ID2 = np.transpose(select(zip(X,Y,Typ,ID), lambda tup: tup[2]==3))
            X1,Y1,ID1 = np.transpose(sorted(zip(X1,Y1,ID1),key=lambda tup: tup[0]))
            headx1,heady1,headid1 = sorted(zip(X1,Y1,ID1),key=lambda tup: tup[2])[1]
            X2,Y2,ID2 = np.transpose(sorted(zip(X2,Y2,ID2),key=lambda tup: tup[0]))
            headx2,heady2,headid2 = sorted(zip(X2,Y2,ID2),key=lambda tup: tup[2])[1]
            ave1 = np.average(Y1)
            count += 1
            Time.append(timestep*dt)            
            Hpx1.append(headx1+Ls/2)
            Hpx2.append(headx2+Ls/2)
            dhx = headx1-headx2
            dhx -= np.round(dhx/Lx)*Lx
            Dpx.append(dhx)
        except RuntimeError:
            pass
    ax1.plot(Time,Dpx,'b-')



    FdT = []
    Fddy = []
    with open('fd.data','r') as fp:
        for line in fp:
            fdt,fddy = [float(a) for a in line.strip().split()]
            FdT.append(fdt)
            Fddy.append(fddy)
    f = interpolate.interp1d(Time,Dpx)
    ax1.plot(FdT,f(FdT),'-')
    L1error.append(trapezoid(np.abs(Fddy-f(FdT)),FdT))
ax1.plot(FdT,Fddy,'k--')
ax1.set_xlabel('time')
ax1.set_ylabel('delta y')
ax2.plot(Rc,L1error,'ro-')
ax2.set_xlabel(r'$r_c$')
ax2.set_ylabel('L1 error')
plt.show()



        
        

