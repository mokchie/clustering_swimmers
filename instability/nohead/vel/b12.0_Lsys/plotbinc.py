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

colors = ['b','r','g','y','c','m','k']*10
linestyles = ['o-','s--']
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

Dirs = ['single_w0.4','pair_w0.4','single_w0.6','pair_w0.6','single_w0.8','pair_w0.8','single_w1.0','pair_w1.0','single_w1.2','pair_w1.2',]
Dirs = ['single_w1.2','pair_w1.2',]
for cn,direct in enumerate(Dirs):
    infiles = [direct+'/'+i for i in ('in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run','in-7.run','in-8.run')]
    Tau = []
    B1m = []
    B2m = []
    for infile in infiles:
        try:
            fn = readvar(infile,'fn')
            tau = readvar(infile,'tau')
            omega = readvar(infile,'omega')
            pdiff = readvar(infile,'pdiff')
            dt = readvar(infile,'Dt')
            kt = readvar(infile,'k')
            Ls = readvar(infile,'Ls')
            Lx = readvar(infile,'Lx')
            Ly = readvar(infile,'Ly')
            dx = readvar(infile,'dx')
            tot = readvar(infile,'tot')
            tottime = tot*dt
            Tp = -2*np.pi/omega
        except FileNotFoundError as err:
            print(err)
            continue
        Ns = int(round(Ls/dx))
        files = []
        for n in tau:
            files.append(direct+'/data/sheetstrj-%s.data'%myround(n))
            Tau.append(n)
        for j,filename in enumerate(files):
            print( 'processing %s ... '%filename )
            Tim = []
            Alpha = []
            Time = []
            cr = 0
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
                        if cr%1==0:
                            Alpha.append(get_angles(X1,Y1))
                            Tim.append(timestep*dt)
                        cr+=1
                except IOError:
                    pass
            Alpha = np.array(Alpha)
            alpha0 = np.average(Alpha,axis=0)
            A = Alpha - np.tile(alpha0,(Alpha.shape[0],1))
            AT = np.transpose(A)
            eigw,eigv = eig(AT.dot(A))
            xs = np.array(range(0,Alpha.shape[1]+1))*dx
            Xa,Ya = rebuild(xs,alpha0)

            Ph = []
            Phx = []
            Phy = []
            for alpha in Alpha:
                def f12(b):
                    return norm(b[0]*eigv[:,0]+b[1]*eigv[:,1]+alpha0-alpha)
                res = minimize(f12,[1,1])
                if res.success == False:
                    print('unable to converge')
                else:
                    Phx.append(res.x[0])
                    Phy.append(res.x[1])
                    Ph.append(np.arctan2(res.x[1],res.x[0]))
            #ax0.plot(Phx,Phy,pointstyles[cn],color=colors[cn],label=direct)
            B1m.append(np.max(Phx))
            B2m.append(np.max(Phy))
    B1m = np.array(B1m)
    B2m = np.array(B2m)
    Tau,B1m,B2m = np.transpose(sorted(zip(Tau,B1m,B2m)))
    with open(direct+'/B.data','w') as fb:
        for ta,b1m,b2m in zip(Tau,B1m,B2m):
            fb.write('%f %f %f\n'%(ta,b1m,b2m))
    ax1.plot(Tau,B1m+B2m/2,linestyles[cn%2],color=colors[int(cn/2)])

ax1.set_xlabel('tau')
ax1.set_ylabel('B')
#ax1.legend(loc='best')


plt.show()
