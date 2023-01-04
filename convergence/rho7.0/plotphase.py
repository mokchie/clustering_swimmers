from variable import *
from readvel import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

infile = 'in.run'
fig,ax = plt.subplots(1,1)
kt = readvar(infile,'k')
fn = readvar(infile,'fn')
omega = readvar(infile,'omega')
dt = readvar(infile,'Dt')
sdist = readvar(infile,'sdist')
rho = readvar(infile,'rho0')
mu = readvar(infile,'mu')
files = ['lammpstrj-%d.data'%n for n in fn]
linestyles = ['b.','g.','y.','m.','r.','b+','g+','y+','m+','k+']
def fcos(x,b,k,theta0):
    return b*np.cos(k*x+theta0)

for j,filename in enumerate(files):
    print 'processing %s ... '%filename
    Time = []
    Dtheta = []
    data = readxyz(filename,[2,3])
    B1 = []
    B2 = []
    for d in data:
        timestep = d['timestep']
        X = d['x']
        Y = d['y']
        Typ = d['type']
        X1,Y1,Typ1 = np.transpose(select(zip(X,Y,Typ), lambda tup: tup[2]==2))
        X2,Y2,Typ2 = np.transpose(select(zip(X,Y,Typ), lambda tup: tup[2]==3))
        X1,Y1 = np.transpose(sorted(zip(X1,Y1)))
        X2,Y2 = np.transpose(sorted(zip(X2,Y2)))
        popt,pcov = curve_fit(fcos,X1,Y1-np.average(Y1))
        b1,k1,theta01 = popt[0:3]
        if b1<0:
            b1=-b1
            theta01 -= np.pi

        B1.append(b1)
        popt,pcov = curve_fit(fcos,X2,Y2-np.average(Y2))
        b2,k2,theta02 = popt[0:3]
        if b2<0:
            b2=-b2
            theta02 -= np.pi
        B2.append(b2)
        Time.append(timestep*dt)
        dtheta = theta01-theta02

        dtheta -= round(dtheta/2.0/np.pi)*2.0*np.pi

        Dtheta.append(dtheta)
    Tp = -2*np.pi/omega[j]
    ax.plot(Time,-np.abs(Dtheta),linestyles[j],label=r'$\omega=%.2f$'%(-omega[j],))
    bave = np.average((np.average(B1),np.average(B2)))
    print 'b/lambda = ',bave/(2*np.pi/kt)
    print 'h/lambda = ',sdist/(2*np.pi/kt)
    print 'Re = ',-bave**2*omega[j]/2/np.pi*rho/mu
ax.set_xlabel('time')
ax.set_ylabel(r'phase diff $\theta_0$')
plt.legend(loc='best',ncol=2)
plt.show()

        
        

