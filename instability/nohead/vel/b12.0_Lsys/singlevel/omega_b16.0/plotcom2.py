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

direct='.'

infiles = [direct+'/'+i for i in ('in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run')]
for infile in infiles:
    try:
        fn = readvar(infile,'fn')
        Tau = readvar(infile,'tau')
        Omega = -readvar(infile,'omega')
        pdiff = readvar(infile,'pdiff')
        dt = readvar(infile,'Dt')
        kt = readvar(infile,'k')
        Ls = readvar(infile,'Ls')
        Lx = readvar(infile,'Lx')
        Ly = readvar(infile,'Ly')
        dx = readvar(infile,'dx')
        Tp = 2*np.pi/Omega
        tot = readvar(infile,'tot')
    except FileNotFoundError as err:
        print(err)
        continue
    Ns = int(round(Ls/dx))
    files1 = []
    files2 = []
    for n in fn:
        files1.append(direct+'/data/swimmer1-pos-%d.data'%n)
        files2.append(direct+'/data/swimmer2-pos-%d.data'%n)
    for i,filename1 in enumerate(files1):
        print( 'processing %s ... '%(filename1,) )
        XC = []
        YC = []
        for j,filename in enumerate((filename1,)):
            dc,lc = readaverows(filename)
            Time = np.array(dc['TimeStep'])*dt
            dstep = Time[1]-Time[0]
            Xc = np.array([d['c_mc%d'%(j+1)][0] for d in lc])
            Yc = np.array([d['c_mc%d'%(j+1)][1] for d in lc])
            Xc -= Xc[0]
            Yc -= Yc[0]
            st = 0#len(Time)-int(Tp[j]/dstep)*3
            XC.append(Xc)
            YC.append(Yc)
            Dc = np.sqrt(Xc**2+Yc**2)
            ax1.plot(Xc[st:],Yc[st:],linestyles[i],label=r'$\tau=%.1f$'%Omega[i])

ax1.legend(loc='best',ncol=3)
plt.show()
