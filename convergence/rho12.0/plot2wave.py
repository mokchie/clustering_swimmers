from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.integrate import trapezoid

infile = 'in.run'
fig,(ax1,ax2,ax3) = plt.subplots(3,1)
fn = readvar(infile,'fn')
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
files = ['lammpstrj-%d.data'%n for n in fn]
linestyles = ['b.','g.','y.','m.','k.','b+','g+','y+','m+','k+']
def fsin(x,b,k,theta0):
    return b*np.sin(k*x+theta0)

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
    ax1.plot(X1,Y1,'b.')

    ax1.plot(X2,Y2,'r.')

    ax1.plot([headx1,],[heady1,],'go')
    ax1.plot([headx2,],[heady2,],'go')

    ax2.plot(Time,Dpx,'b-')

    ax3.plot(Time,Hpx1,'b-')
    ax3.plot(Time,Hpx2,'r-')

    ax1.set_ylim((y02-2,y01+3))
    ax1.set_xlim((-Ls/2-2,Ls/2+2))
    ax1.text(-Ls/2+2,y01+1,'T = %.0f\n(totT = %.0f)'%(timestep*dt,tot*dt))
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax2.set_xlabel('time')
    ax2.set_ylabel('delta y')
    ax3.set_ylabel('y')
    FdT = []
    Fddy = []
    with open('fd.data','r') as fp:
        for line in fp:
            fdt,fddy = [float(a) for a in line.strip().split()]
            FdT.append(fdt)
            Fddy.append(fddy)
    ax2.plot(FdT,Fddy,'k--')
    f = interpolate.interp1d(Time,Dpx)
    ax2.plot(FdT,f(FdT),'b-')
    print(trapezoid(np.abs(Fddy-f(FdT)),FdT))
    plt.show()



        
        

