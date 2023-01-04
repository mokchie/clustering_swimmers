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

def plot_shape(ax,direct,tau,N,ln,cn,lab):
    infile = direct+'/in-1.run'
    Lx = readvar(infile,'Lx')
    Ly = readvar(infile,'Ly')
    filename = direct+'/data/sheetstrj-%s.data'%myround(tau)
    data = readxyz(filename,[2,3])
    found = False
    for d in data:
        timestep = d['timestep']
        if timestep != N:
            continue
        else:
            found = True
            X = d['x']
            Y = d['y']
            Typ = d['type']
            ID = d['id']
            X1,Y1,ID1,Typ1 = np.transpose(select(zip(X,Y,ID,Typ), lambda tup: tup[3]==2))
            X1,Y1,ID1 = np.transpose(sorted(zip(X1,Y1,ID1),key=lambda tup: tup[2]))
            l1 = int(len(X1)/3)
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
            ax.plot(X1,Y1,linestyles[ln],color=colors[cn],label=lab)
    if not found:
        print('specified timestep can not be found')

Direct1 = ['single_w0.4','single_w0.6','single_w0.8','single_w1.0','single_w1.2']
Direct2 = ['pair_w0.4','pair_w0.6','pair_w0.8','pair_w1.0','pair_w1.2']
Omega = [0.4,0.6,0.8,1.0,1.2]
N = 400000
tau = 15
fig,ax1 = plt.subplots(1,1)
for cn,(direct1,direct2) in enumerate(zip(Direct1,Direct2)):
    plot_shape(ax1,direct1,tau,N,0,cn,r'single $\omega=%.1f$'%Omega[cn])
    plot_shape(ax1,direct2,tau,N,2,cn,r'pair $\omega=%.1f$'%Omega[cn])
ax1.legend(loc='best',ncol=4)
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
plt.show()
