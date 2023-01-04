from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,minimize
import pdb
from numpy.linalg import norm
from numpy.linalg import eig
import matplotlib
import non_dimensionalize as ndm
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams.update({'font.size': 13})
Stau = [0.1,1.3,3,40]
colors = ['b','r','g','y','c','m','k','C0','C1','C2','C3','C4','C5']*10
linestyles = ['-','--','-.']
pointstyles = ['.','x','+','v','<','s','*']
Time_list = []
Ds_list = []
Tau_list = []
De_list = []
def tanfun(t,k,t0,phi0,b):
    return k*np.arctan(np.tan(phi0/2)*np.exp(-t/t0))+b
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
        
fig1,ax1 = plt.subplots(1,1)
fig2,ax2 = plt.subplots(1,1)

direct='.'
infiles = [direct+'/'+i for i in ('in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run','in-7.run','in-8.run','in-9.run')]
cn = -1
for infile in infiles:
    try:
        fn = readvar(infile,'fn')
        Tau = readvar(infile,'tau')
        Tau = np.array([myround(ii) for ii in Tau])
#        pdb.set_trace()
        omega = -readvar(infile,'omega')
        pdiff = readvar(infile,'pdiff')
        dt = readvar(infile,'Dt')
        kt = readvar(infile,'k')
        Ls = readvar(infile,'Ls')
        Lx = readvar(infile,'Lx')
        Ly = readvar(infile,'Ly')
        dx = readvar(infile,'dx')
        Tp = -2*np.pi/omega
        tot = readvar(infile,'tot')
    except FileNotFoundError as err:
        print(err)
        continue
    Ns = int(round(Ls/dx))
    filenames = []
    for n in Tau:
        filenames.append(direct+'/data/sheetstrj-%s.data'%myround(n))
    for i,filename in enumerate(filenames):
#         if Tau[i] not in Stau:
#             continue
        cn += 1
        print( 'processing %s ... '%(filename,) )
        Time = []
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
                del_x = np.average(X1)-np.average(X2)
                del_y = np.average(Y1)-np.average(Y2)
                del_x -= np.round(del_x/Lx)*Lx
                del_y -= np.round(del_y/Ly)*Ly
                ds = np.abs(del_y)#np.sqrt(del_x**2+del_y**2)
                Time.append(timestep*dt)
                Ds.append(ds)
            except IOError:
                pass
        Time = np.array(Time)
        Ds = np.array(Ds)
#        pdb.set_trace()
#        ax1.plot(Time,Ds,'-',label=r'$De=%.1f$'%(Tau[i]*omega))
        Time_list.append(Time)
        Ds_list.append(Ds)
        De_list.append(Tau[i]*omega)
        Tau_list.append(Tau[i])

DeTimeDs = sorted(zip(De_list,Time_list,Ds_list,Tau_list))
T0 = []
DE = []
TAU = []
cn = 0
for i,(De,Time,Ds,tau) in enumerate(DeTimeDs):
    Time = Time[80:]
    Ds = Ds[80:]
    popt,conv = curve_fit(tanfun,Time,Ds,p0=[(np.max(Ds)-np.min(Ds))/np.pi*2,(Time[-1]-Time[0])/2,np.pi-0.01,np.min(Ds)])
    if tau in Stau:
        ax1.plot(np.array(coarseave(Time,10))*omega,ndm.ndm_length(np.array(coarseave(Ds,10))),'.',color = colors[cn],label=r'$\mathrm{De}=%s$'%(myround(De),),markersize=8,linewidth=1)
        ax1.plot(np.array(Time)*omega,ndm.ndm_length(np.array(tanfun(Time,popt[0],popt[1],popt[2],popt[3]))),'-',color=colors[cn],markersize=8,linewidth=1)
        cn+=1
    T0.append(popt[1])
    DE.append(De)
ax2.plot(DE,omega*np.array(T0),'b.-',label='without a head',markersize=8,linewidth=1)
print(omega)
Time_list = []
Ds_list = []
Tau_list = []
De_list = []
direct='../moving_withhead'
infiles = [direct+'/'+i for i in ('in-1.run','in-2.run','in-3.run','in-4.run','in-5.run','in-6.run','in-7.run','in-8.run','in-9.run')]
cn = -1
for infile in infiles:
    try:
        fn = readvar(infile,'fn')
        Tau = readvar(infile,'tau')
        Tau = np.array([myround(ii) for ii in Tau])
#        pdb.set_trace()
        omega = -readvar(infile,'omega')
        pdiff = readvar(infile,'pdiff')
        dt = readvar(infile,'Dt')
        kt = readvar(infile,'k')
        Ls = readvar(infile,'Ls')
        Lx = readvar(infile,'Lx')
        Ly = readvar(infile,'Ly')
        dx = readvar(infile,'dx')
        Tp = -2*np.pi/omega
        tot = readvar(infile,'tot')
    except FileNotFoundError as err:
        print(err)
        continue
    Ns = int(round(Ls/dx))
    filenames = []
    for n in Tau:
        filenames.append(direct+'/data/sheetstrj-%s.data'%myround(n))
    for i,filename in enumerate(filenames):
#         if Tau[i] not in Stau:
#             continue
        cn += 1
        print( 'processing %s ... '%(filename,) )
        Time = []
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
                del_x = np.average(X1)-np.average(X2)
                del_y = np.average(Y1)-np.average(Y2)
                del_x -= np.round(del_x/Lx)*Lx
                del_y -= np.round(del_y/Ly)*Ly
                ds = np.sqrt(del_x**2+del_y**2)
                Time.append(timestep*dt)
                Ds.append(ds)
            except IOError:
                pass
        Time = np.array(Time)
        Ds = np.array(Ds)
#        pdb.set_trace()
#        ax1.plot(Time,Ds,'-',label=r'$De=%.1f$'%(Tau[i]*omega))
        Time_list.append(Time)
        Ds_list.append(Ds)
        De_list.append(Tau[i]*omega)
        Tau_list.append(Tau[i])

DeTimeDs = sorted(zip(De_list,Time_list,Ds_list,Tau_list))
T0 = []
DE = []
TAU = []
for i,(De,Time,Ds,tau) in enumerate(DeTimeDs):
    Time = Time[80:]
    Ds = Ds[80:]
    popt,conv = curve_fit(tanfun,Time,Ds,p0=[(np.max(Ds)-np.min(Ds))/np.pi*2,(Time[-1]-Time[0])/2,np.pi-0.01,np.min(Ds)])
    T0.append(popt[1])
    DE.append(De)
ax2.plot(DE,omega*np.array(T0),'r.-',label='with a head',markersize=8,linewidth=1)
print(omega)

ax2.set_xlabel('$\mathrm{De}$')
ax2.set_ylabel(r'$\omega\tau_s$')
ax1.set_xlabel('$\omega t$')
ax1.set_ylabel(r'$\Delta D/r_c$')
ax1.legend(loc='best',ncol=2)
ax2.legend(loc='best',ncol=1)
#ax1.text(150,3,'(a)',fontdict={'size':18})
#ax2.text(10,105,'(b)',fontdict={'size':18})
plt.show()
