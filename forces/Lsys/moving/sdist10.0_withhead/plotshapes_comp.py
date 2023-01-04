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
import matplotlib
matplotlib.rcParams.update({'font.size': 13})
colors = ['b','r','g','y','c','m','k']*50
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
    l = len(X)
    for i,x in enumerate(X[1:]):
        if np.abs(x-X[i])>L/2:
            if x>X[i]:
                for j in range(i+1,l):
                    X[j]-=L
            else:
                for j in range(i+1,l):
                    X[j]+=L
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

def read_shape(direct,tau):
    infile = direct+'/in-1.run'
    Lx = readvar(infile,'Lx')
    Ly = readvar(infile,'Ly')
    Ls = readvar(infile,'Ls')
    dx = readvar(infile,'dx')
    dt = readvar(infile,'Dt')
    Ns = int(round(Ls/dx))
    filename = direct+'/data/sheetstrj-%s.data'%myround(tau)
    print(filename)
    data = readxyz(filename,[2,3])
    Time = []
    P1 = []
    P2 = []

    for i,d in enumerate(data):
        timestep = d['timestep']
        Time.append(timestep*dt)
        X = d['x']
        Y = d['y']
        Typ = d['type']
        ID = d['id']

        X1,Y1,ID1,Typ1 = np.transpose(select(zip(X,Y,ID,Typ), lambda tup: tup[3]==2))
        X1,Y1,ID1 = np.transpose(sorted(zip(X1,Y1,ID1),key=lambda tup: tup[2]))

        X2,Y2,ID2,Typ2 = np.transpose(select(zip(X,Y,ID,Typ), lambda tup: tup[3]==3))
        X2,Y2,ID2 = np.transpose(sorted(zip(X2,Y2,ID2),key=lambda tup: tup[2]))

        l1 = Ns
        X1 = X1[l1:2*l1]
        Y1 = Y1[l1:2*l1]
        ID1 = ID1[l1:2*l1]

        l2 = Ns
        X2 = X2[l2:2*l2]
        Y2 = Y2[l2:2*l2]
        ID2 = ID2[l2:2*l2]

        recover(X1,Lx)
        recover(X2,Lx)
        recover(Y1,Ly)
        recover(Y2,Ly)

        P1.append([X1,Y1])
        P2.append([X2,Y2])
    return (Time,P1,P2,Lx,Ly)
Stau = [0.1,2.55,25]
Shifty = [0,16,32]

fig1,ax1 = plt.subplots(1,1)
for cn,tau in enumerate(Stau):
    Time,P1,P2,Lx,Ly = read_shape('.',tau)
    P1 = np.array(P1)
    P2 = np.array(P2)
    for i in range(np.shape(P1)[2]):
        recover(P1[:,0,i],Lx)
        recover(P1[:,1,i],Ly)
    for i in range(np.shape(P2)[2]):
        recover(P2[:,0,i],Lx)
        recover(P2[:,1,i],Ly)

    t = Time[-1]
    p1 = P1[-1]
    p2 = P2[-1]
    pcx = (np.average(p1[0])+np.average(p2[0]))/2
    pcy = (np.average(p1[1])+np.average(p2[1]))/2
    p1[0]-=pcx
    p2[0]-=pcx
    p1[1]-=pcy-Shifty[cn]
    p2[1]-=pcy-Shifty[cn]
    ax1.plot(p1[0],p1[1],'-',color=colors[cn])
    ax1.plot((p1[0][0],),(p1[1][0],),'o',color=colors[cn])
    ax1.plot(p2[0],p2[1],'-',color=colors[cn])
    ax1.plot((p2[0][0],),(p2[1][0],),'o',color=colors[cn])
    ax1.text(0,Shifty[cn]+4,r'$\mathrm{De}=%s$'%(myround(tau*np.abs(readvar('in-1.run','omega')))))
#ax1.legend(loc='best',ncol=4)
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
#ax1.set_ylim((-8,32))
ax1.set_aspect('equal')
ax1.set_axis_off()
plt.show()
