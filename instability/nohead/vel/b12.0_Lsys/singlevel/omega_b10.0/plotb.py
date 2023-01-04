from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pdb

infile = 'in.run'
fig,ax1 = plt.subplots(1,1)
fn = readvar(infile,'fn')
B = readvar(infile,'b')
omega = readvar(infile,'omega')
pdiff = readvar(infile,'pdiff')
dt = readvar(infile,'Dt')
kt = readvar(infile,'k')
Ls = readvar(infile,'Ls')
Lx = readvar(infile,'Lx')
Ly = readvar(infile,'Ly')
dx = readvar(infile,'dx')
tot = readvar(infile,'tot')
Ns = int(round(Ls/dx))
files = ['data/sheetstrj-%d.data'%n for n in fn]
linestyles = ['b.','g.','y.','m.','k.','b+','g+','y+','m+','k+']
def fsin(x,b,k,theta0):
    return b*np.sin(k*x+theta0)
def flin(x,k):
    return k*x
def recover(X,L):
    for i,x in enumerate(X[1:]):
        if np.abs(x-X[i])>L/2:
            if x>X[i]:
                X[i+1]-=L
            else:
                X[i+1]+=L
def getsign(x):
    if x!=0:
        return x/np.abs(x)
    else:
        return 1.0
def rotate(x,y,theta):
    t0 = np.arctan2(y,x)
    r = np.sqrt(x**2+y**2)
    t1 = t0+theta
    x,y = np.cos(t1)*r,np.sin(t1)*r
    return (x,y)
for j,filename in enumerate(files):
    if B[j]!=12.0:
        continue
    print( 'processing %s ... '%filename )
    Time = []
    dth10 = 0
    dth20 = 0
    dth1 = 0
    dth2 = 0
    data = readxyz(filename,[2,3])
    count = 0
    for d in data:
        try:
            timestep = d['timestep']
            X = d['x']
            Y = d['y']
            Typ = d['type']
            ID = d['id']
            X1,Y1,ID1,Typ1 = np.transpose(select(zip(X,Y,ID,Typ), lambda tup: tup[3]==2))
            X1,Y1,ID1 = np.transpose(sorted(zip(X1,Y1,ID1),key=lambda tup: tup[2]))
            l1 = Ns
            X1 = X1[l1:2*l1]
            Y1 = Y1[l1:2*l1]
            ID1 = ID1[l1:2*l1]

            headx1,heady1,headid1 = sorted(zip(X1,Y1,ID1),key=lambda tup: tup[2])[0]
            recover(X1,Lx)
            recover(Y1,Ly)
            X1-=X1[0]
            Y1-=Y1[0]

            pr1 = np.arctan2(Y1[-1],X1[-1])
            X1,Y1 = rotate(X1,Y1,-pr1)

            popt1,pcon1 = curve_fit(flin,X1,Y1)
            if np.abs(popt1[0])>1:
                xa1 = getsign(np.average(Y1)*popt1[0])
            else:
                xa1 = getsign(np.average(X1))
            dth1 = np.arctan2(flin(xa1,popt1[0]),xa1)

            X1,Y1 = rotate(X1,Y1,-dth1)
            p1 = pr1+dth1

            ax1.clear()
            ax1.plot(X1,Y1,'b.')

            Time.append(timestep*dt)            

            ax1.set_xlim((-2,Lx/2))
            ax1.set_ylim((-5,5))
            ax1.text(-Lx/2+4,Ly/2-4,'T = %.0f\n(totT = %.0f)'%(timestep*dt,tot*dt))
            ax1.set_xlabel('x')
            ax1.set_ylabel('y')

            fig.canvas.draw()
            plt.pause(.01)
        except RuntimeError:
            pass

