from variable import *
from readvel import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

infile = 'in.run'
fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
fn = readvar(infile,'fn')
fn = [1,]
kt = readvar(infile,'k')
omega = readvar(infile,'omega')
dt = readvar(infile,'Dt')
Ls = readvar(infile,'Ls')
files = ['lammpstrj-%d.data'%n for n in fn]
linestyles = ['b.','g.','y.','m.','k.','b+','g+','y+','m+','k+']
def fcos(x,b,k,theta0):
    return b*np.cos(k*x+theta0)

for j,filename in enumerate(files):
    B = []
    K = []
    Yp1 = []
    Yp2 = []
    avey1 = []
    avey2 = []
    print 'processing %s ... '%filename
    Time = []
    Dtheta = []
    data = readxyz(filename,[2,3])
    for d in data:
        timestep = d['timestep']
        X = d['x']
        Y = d['y']
        Typ = d['type']
        ID = d['id']
        X1,Y1,ID1,Typ1 = np.transpose(select(zip(X,Y,ID,Typ), lambda tup: tup[3]==2))
        X2,Y2,ID2,Typ2 = np.transpose(select(zip(X,Y,ID,Typ), lambda tup: tup[3]==3))

        Yp1.append(sorted(zip(ID1,Y1))[0][1])
        Yp2.append(sorted(zip(ID2,Y2))[0][1])
        X1,Y1 = np.transpose(sorted(zip(X1,Y1)))
        X2,Y2 = np.transpose(sorted(zip(X2,Y2)))
        ave1 = np.average(Y1)
        popt,pcov = curve_fit(fcos,X1,Y1-ave1,p0=[1.0,kt,0.0])

        b1,k1,theta01 = popt[0:3]
        if b1<0:
            b1=-b1
            theta01 -= np.pi

        ave2 = np.average(Y2)
        popt,pcov = curve_fit(fcos,X2,Y2-ave2,p0=[1.0,kt,0.0])
        b2,k2,theta02 = popt[0:3]
        if b2<0:
            b2=-b2
            theta02 -= np.pi

        Time.append(timestep*dt)
        avey1.append(np.average(Y1))
        avey2.append(np.average(Y2))
        B.append(b1)
        K.append(k1)
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    ax1.plot(Time,B,'.-')
    ax2.plot(Time,K,'.-')
    ax2.plot((Time[0],Time[-1]),(kt,kt),'--')
    ax3.plot(Time,Yp1,'b-')
    ax3.plot(Time,Yp2,'r-')
    ax4.plot(Time,avey1,'-')
    ax4.plot(Time,avey2,'-')
#        ax.set_ylim((-1.5,1.5))
#        ax.set_xlim((-Ls/2-2,Ls/2+2))
    ax1.set_xlabel('Time')
    ax1.set_ylabel('b')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('k')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('y')

#        fig.canvas.draw()
#        plt.pause(.001)
plt.show()


        
        

