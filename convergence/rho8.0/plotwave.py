from variable import *
from readvel import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

infile = 'in.run'
fig,(ax1,ax2,ax3) = plt.subplots(3,1)
fn = readvar(infile,'fn')
fn = [1,1,1,1,1,1,1,1,1,1]
omega = readvar(infile,'omega')
y01 = readvar(infile,'y01')
y02 = readvar(infile,'y02')
pdiff = readvar(infile,'pdiff')
dt = readvar(infile,'Dt')
kt = readvar(infile,'k')
Ls = readvar(infile,'Ls')
Ly = readvar(infile,'Ly')
dx = readvar(infile,'dx')
tot = readvar(infile,'tot')
files = ['lammpstrj-%d.data'%n for n in fn]
linestyles = ['b.','g.','y.','m.','k.','b+','g+','y+','m+','k+']
def fsin(x,b,k,theta0):
    return b*np.sin(k*x+theta0)

for j,filename in enumerate(files):
    print 'processing %s ... '%filename
    Time = []
    Timea = []
    Dtheta1 = []
    Dtheta2 = []
    Dpx = []
    Dpxa = []
    Hpx1 = []
    Hpx1a = []
    Hpx2 = []
    Hpx2a = []
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
            popt,pcov = curve_fit(fsin,X1,Y1-ave1,p0=[1.0,kt,0.0])
            b1,k1,theta01 = popt[0:3]
            if b1<0:
                b1=-b1
                theta01 -= np.pi

            ave2 = np.average(Y2)
            popt,pcov = curve_fit(fsin,X2,Y2-ave2,p0=[1.0,kt,0.0])
            b2,k2,theta02 = popt[0:3]
            if b2<0:
                b2=-b2
                theta02 -= np.pi

            count += 1

            ax1.clear()
            ax2.clear()
            ax3.clear()

            ax1.plot(X1,Y1,'b.')
            ax1.plot(X1,fsin(X1,b1,k1,theta01)+ave1,'b-')

            ax1.plot(X2,Y2,'r.')
            ax1.plot(X2,fsin(X2,b2,k2,theta02)+ave2,'r-')
            if timestep<=tot/2:
                Time.append(timestep*dt)            
                Hpx1.append(headx1+Ls/2)
                Hpx2.append(headx2+Ls/2)
                Dpx.append(headx1-headx2)
            else:
                Timea.append((tot-timestep)*dt)
                Hpx1a.append(headx1+Ls/2)
                Hpx2a.append(headx2+Ls/2)
                Dpxa.append(headx1-headx2)
            ax1.plot([headx1,],[heady1,],'go')
            ax1.plot([headx2,],[heady2,],'go')

            ax2.plot(Time,Dpx,'b-')
            ax2.plot(Timea,Dpxa,'b--')
            ax2.plot((Time[0],Time[-1]),(1,1),'b--')
            ax2.plot((Time[0],Time[-1]),(-1,-1),'b--')

            ax3.plot(Time,Hpx1,'b-')
            ax3.plot(Timea,Hpx1a,'b--')
            
            ax3.plot(Time,Hpx2,'r-')
            ax3.plot(Timea,Hpx2a,'r--')

            ax1.set_ylim((y02-2,y01+3))
            ax1.set_xlim((-Ls/2-2,Ls/2+2))
#            ax2.set_ylim((-1.5,1.5))
#            ax3.set_ylim((-1,2.5))
            ax1.text(-Ls/2+2,y01+1,'T = %.0f\n(totT = %.0f)'%(timestep*dt,tot*dt))
            ax1.set_xlabel('x')
            ax1.set_ylabel('y')
            ax2.set_xlabel('time')
            ax2.set_ylabel('phase lag')
            fig.canvas.draw()
            plt.pause(.001)
        except RuntimeError:
            pass



        
        

