from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,minimize
import pdb
from numpy.linalg import norm
from numpy.linalg import eig
import pdb
colors = ['b','r','g','y','c','m','k']
linestyles = ['-','--','-.']
pointstyles = ['.','x','+','v','<','s','*']
def cstyle(i,colors=colors,linestyles=pointstyles):
    c = i%len(colors)
    p = int(i/len(colors))
    p = p%len(linestyles)
    return colors[c]+linestyles[p]
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
def get_angles(X,Y):
    X12 = zip(X[0:-1],X[1:])
    Y12 = zip(Y[0:-1],Y[1:])
    Alpha = []
    for x12,y12 in zip(X12,Y12):
        x1,x2 = x12
        y1,y2 = y12
        alpha = np.arctan2(y2-y1,x2-x1)
        Alpha.append(alpha)
    return Alpha
def rebuild(Xs,Alpha):
    X = [0,]
    Y = [0,]
    Ds = Xs[1:]-Xs[0:-1]
    for ds,alpha in zip(Ds,Alpha):
        X.append(X[-1]+ds*np.cos(alpha))
        Y.append(Y[-1]+ds*np.sin(alpha))
    return np.array(X),np.array(Y)
        
fig,ax1 = plt.subplots(1,1)

DS = []
TAU = []
direct='.'
infiles = [direct+'/'+i for i in ('in-1.run','in-2.run','in-3.run','in-4.run')]
for infile in infiles:
    try:
        fn = readvar(infile,'fn')
        Tau = readvar(infile,'tau')
        omega = readvar(infile,'omega')
        pdiff = readvar(infile,'pdiff')
        dt = readvar(infile,'Dt')
        Tp = -2*np.pi/omega
        kt = readvar(infile,'k')
        Ls = readvar(infile,'Ls')
        Lx = readvar(infile,'Lx')
        Ly = readvar(infile,'Ly')
        dx = readvar(infile,'dx')
        tot = readvar(infile,'tot')
        tottime = tot*dt
    except FileNotFoundError as err:
        print(err)
        continue
    Ns = int(round(Ls/dx))
    files = []
    for n in Tau:
        if n<1:
            files.append(direct+'/data/sheetstrj-%.1f.data'%n)
        else:
            files.append(direct+'/data/sheetstrj-%d.data'%n)
    for j,filename in enumerate(files):
        TAU.append(Tau[j])
        print( 'processing %s ... '%filename )
        Time = []
        Tim = []
        Alpha = []
        Yaw = []
        Hp = []
        Ds = []
        data = readxyz(filename,[2,3])
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


                X2,Y2,ID2,Typ2 = np.transpose(select(zip(X,Y,ID,Typ), lambda tup: tup[3]==3))
                X2,Y2,ID2 = np.transpose(sorted(zip(X2,Y2,ID2),key=lambda tup: tup[2]))
                l2 = Ns
                X2 = X2[l2:2*l2]
                Y2 = Y2[l2:2*l2]
                ID2 = ID2[l2:2*l2]
                headx2,heady2,headid2 = sorted(zip(X2,Y2,ID2),key=lambda tup: tup[2])[0]

                recover(X2,Lx)
                recover(Y2,Ly)
                d1 = np.average(X1)-np.average(X2)
                d1 = d1-np.round(d1/Lx)*Lx
                d2 = np.average(Y1)-np.average(Y2)
                d2 = d2-np.round(d2/Ly)*Ly
                ds = [d1,d2]
                
                yaw = np.arctan2(Y1[0]-Y1[1],X1[0]-X1[1])
                hpx = X1[0]-np.average(X1)
                hpy = Y1[0]-np.average(Y1)
                x0 = X1[0]
                y0 = Y1[0]
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
                Time.append(timestep*dt)

                if Time[-1]>tottime-Tp:
                    Tim.append(Time[-1])
                    Alpha.append(get_angles(X1,Y1))
                    Yaw.append(yaw)
                    Hp.append([hpx,hpy])
                    Ds.append(ds)
            except IOError:
                pass

        Hp = np.array(Hp)
        Ds = np.array(Ds)
        DS.append(np.average(np.sqrt(Ds[:,0]**2+Ds[:,1]**2)))

ax1.plot(TAU,DS)
ax1.set_xlabel('t')
ax1.set_ylabel('d')
plt.show()
