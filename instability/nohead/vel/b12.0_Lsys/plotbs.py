from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pdb
from numpy.linalg import norm

stau = 10.0
colors = ['b','r','g','y','c','m','k']
linestyles = ['.-','.--','.-.']
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
    X123 = zip(X[0:-2],X[1:-1],X[2:])
    Y123 = zip(Y[0:-2],Y[1:-1],Y[2:])
    Alpha = []
    for x123,y123 in zip(X123,Y123):
        x1,x2,x3 = x123
        y1,y2,y3 = y123
        v1 = np.array([x2-x1,y2-y1])
        v2 = np.array([x3-x2,y3-y2])
        alpha = np.arccos(v1.dot(v2)/norm(v1)/norm(v2))
        Alpha.append(alpha)
    return Alpha


        
fig,ax1 = plt.subplots(1,1)
cn = 0

Dirs = ['single_w0.8','pair_w0.8']
for direct in Dirs:
    infiles = [direct+'/'+i for i in ('in-1.run','in-2.run','in-3.run','in-4.run')]
    for infile in infiles:
        try:
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
            Tp = -2*np.pi/omega
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
            if Tau[j]!=stau:
                continue
            print( 'processing %s ... '%filename )
            Time = []
            data = readxyz(filename,[2,])
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
                    x0 = X1[0]
                    y0 = Y1[0]
                    X1-=x0
                    Y1-=y0

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
                        ax1.plot(X1+x0-np.average(X1),Y1+y0-np.average(Y1),cstyle(cn))
                except IOError:
                    pass
            cn+=1
#ax1.set_xlim((-2,Lx/2))
#ax1.set_ylim((-3,3))
ax1.set_xlabel('x')
ax1.set_ylabel('y')
plt.show()
