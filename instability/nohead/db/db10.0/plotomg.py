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
Tau = readvar(infile,'tau')
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
colors = ['b','r','g','y','m','c','k']*10
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
    print( 'processing %s ... '%filename )
    Time = []
    dth1 = 0
    dth2 = 0
    data = readxyz(filename,[2,3])
    count = 0
    Or1 = []
    Or2 = []
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
            X2,Y2,ID2,Typ2 = np.transpose(select(zip(X,Y,ID,Typ), lambda tup: tup[3]==3))
            X2,Y2,ID2 = np.transpose(sorted(zip(X2,Y2,ID2),key=lambda tup: tup[2]))
            l2 = Ns
            X2 = X2[l2:2*l2]
            Y2 = Y2[l2:2*l2]
            ID2 = ID2[l2:2*l2]

            headx1,heady1,headid1 = sorted(zip(X1,Y1,ID1),key=lambda tup: tup[2])[0]
            headx2,heady2,headid2 = sorted(zip(X2,Y2,ID2),key=lambda tup: tup[2])[0]
            recover(X1,Lx)
            recover(Y1,Ly)
            recover(X2,Lx)
            recover(Y2,Ly)
            X1-=X1[0]
            Y1-=Y1[0]
            X2-=X2[0]
            Y2-=Y2[0]
            X10,Y10,X20,Y20 = X1,Y1,X2,Y2
            
            pr1 = np.arctan2(Y1[-1],X1[-1])
            pr2 = np.arctan2(Y2[-1],X2[-1])
            X1,Y1 = rotate(X1,Y1,-pr1)
            X2,Y2 = rotate(X2,Y2,-pr2)

            popt1,pcon1 = curve_fit(flin,X1,Y1)
            popt2,pcon2 = curve_fit(flin,X2,Y2)
            if np.abs(popt1[0])>1:
                xa1 = getsign(np.average(Y1)*popt1[0])
            else:
                xa1 = getsign(np.average(X1))
            if np.abs(popt2[0])>1:
                xa2 = getsign(np.average(Y2)*popt2[0])
            else:
                xa2 = getsign(np.average(X2))
            dth1 = np.arctan2(flin(xa1,popt1[0]),xa1)
            dth2 = np.arctan2(flin(xa2,popt2[0]),xa2)

            X1,Y1 = rotate(X1,Y1,-dth1)
            X2,Y2 = rotate(X2,Y2,-dth2)
            p1 = pr1+dth1
            p2 = pr2+dth2
            Or1.append(p1)
            Or2.append(p2)
            Time.append(timestep*dt)
        except IOError as err:
            print(err)
    Time = np.array(Time)
    Or1 = np.array(Or1)
    Or2 = np.array(Or2)
    Om1 = (Or1[1:]-Or1[0:-1])
    Om2 = (Or2[1:]-Or2[0:-1])
    for i,om in enumerate(Om1):
        if om>np.pi:
            Om1[i]-=2*np.pi
        elif om<-np.pi:
            Om1[i]+=2*np.pi
    for i,om in enumerate(Om2):
        if om>np.pi:
            Om2[i]-=2*np.pi
        elif om<-np.pi:
            Om2[i]+=2*np.pi
            
    Om1/=(Time[1:]-Time[0:-1])
    Om2/=(Time[1:]-Time[0:-1])
    Timem = (Time[1:]+Time[0:-1])/2
    aven = 20
    ax1.plot(coarseave(Timem,aven),coarseave(Om1,aven),'b-',color=colors[j],label=r'$\tau=%s$'%Tau[j]) 
    ax1.plot((Timem[0],Timem[-1]),(np.average(Om1),np.average(Om1)),'-.',color='b')
    ax1.plot(coarseave(Timem,aven),coarseave(Om2,aven),'r--',)
    ax1.plot((Timem[0],Timem[-1]),(np.average(Om2),np.average(Om2)),':',color='r')
    print('ave omg1: ',np.average(Om1))
    print('ave omg2: ',np.average(Om2))
ax1.set_xlabel('Time')
ax1.set_ylabel(r'$\omega$')
ax1.legend(loc='best',ncol=2)
plt.show()

