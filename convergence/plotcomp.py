from variable3 import *
from readvel3 import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.integrate import trapezoid
import matplotlib
matplotlib.rcParams.update({'font.size':14,'font.family':'sans-serif'})
Rho = [6.0,7.0,8.0,9.0,10.0,11.0,12.0]
Rhodirs = ['rho'+str(i) for i in Rho]
fig1,ax1 = plt.subplots(1,1)
fig2,ax2 = plt.subplots(1,1)
linestyles = ['b.','g.','y.','m.','k.','b+','g+','y+','m+','k+']
colors = ['C'+str(i) for i in range(10)]
cmap = plt.get_cmap('jet')
L1errors = []
def fsin(x,b,k,theta0):
    return b*np.sin(k*x+theta0)
for rhodir in Rhodirs:
    infile = rhodir+'/in.run'
    omega = readvar(infile,'omega')
    y01 = readvar(infile,'y01')
    y02 = readvar(infile,'y02')
    pdiff = readvar(infile,'pdiff')
    dt = readvar(infile,'Dt')
    kt = readvar(infile,'k')
    Ls = readvar(infile,'Ls')
    Ly = readvar(infile,'Ly')
    Lx = readvar(infile,'Lx')
    Tp = readvar(infile,'Tp')
    dx = readvar(infile,'dx')
    tot = readvar(infile,'tot')
    Rc = readvar(infile,'rc')
    files = [rhodir+'/lammpstrj-%.1f.data'%rc for rc in Rc]
    l1error = []
    for j,filename in enumerate(files):
        print( 'processing %s ... '%filename )
        Time = []
        Timea = []
        Dpx = []
        Hpx1 = []
        Hpx2 = []
        data = readxyz(filename,[2,3])
        count = 0
        for d in data:
            try:
                timestep = d['timestep']
                X = d['x']
                Y = d['y']
                Typ = d['type']
                ID = d['id']
                X1,Y1,Typ1,ID1 = np.transpose(select(zip(X,Y,Typ,ID), lambda tup: tup[2]==2))
                X2,Y2,Typ2,ID2 = np.transpose(select(zip(X,Y,Typ,ID), lambda tup: tup[2]==3))
                X1,Y1,ID1 = np.transpose(sorted(zip(X1,Y1,ID1),key=lambda tup: tup[0]))
                headx1,heady1,headid1 = sorted(zip(X1,Y1,ID1),key=lambda tup: tup[2])[1]
                X2,Y2,ID2 = np.transpose(sorted(zip(X2,Y2,ID2),key=lambda tup: tup[0]))
                headx2,heady2,headid2 = sorted(zip(X2,Y2,ID2),key=lambda tup: tup[2])[1]
                ave1 = np.average(Y1)
                count += 1
                Time.append(timestep*dt)            
                Hpx1.append(headx1+Ls/2)
                Hpx2.append(headx2+Ls/2)
                dhx = headx1-headx2
                dhx -= np.round(dhx/Lx)*Lx
                Dpx.append(dhx)
            except RuntimeError:
                pass

        FdT = []
        Fddy = []
        with open(rhodir+'/fd.data','r') as fp:
            for line in fp:
                fdt,fddy = [float(a) for a in line.strip().split()]
                FdT.append(fdt)
                Fddy.append(fddy)
        f = interpolate.interp1d(Time,Dpx)
        Fddy = np.array(Fddy)
        FdT = np.array(FdT)
        l1error.append(np.average(np.abs(Fddy/(Tp/2/np.pi)-f(FdT)/(Tp/2/np.pi))))

    L1errors.append(l1error)
L1errors = np.array(L1errors)
for i,rc in enumerate(Rc):
    ax1.plot(Rho,L1errors[:,i],'o-',color=cmap(i/len(Rc)),label=r'$r_c=%.1f$'%rc)
Nng_all = []
l1_all = []
for i,rho in enumerate(Rho):
    print(rho)
    print(Rc)
    Nng = Rc**2*np.pi*rho
    print(Nng)    
    mask = Nng<100
    Nng_all += list(Nng[mask])
    l1_all += list(L1errors[i,:][mask])
    ax2.plot(Nng[mask],L1errors[i,:][mask],'o',color=cmap(i/len(Rho)),label=r'$d_0=%.1f$'%rho)
def fpow(x,k,n):
    return k*x**n
popt,pconv = curve_fit(fpow,Nng_all,l1_all)
print(popt)
Nng_all.sort()
x0 = np.array(Nng_all)
y0 = fpow(x0,popt[0],popt[1])
ax2.plot(x0,y0,'k--')
ax1.set_xlabel(r'$d_0$')
ax1.set_ylabel('L1 error')
ax1.legend(loc='best',ncol=2)
ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.set_xlabel(r'$N_\mathrm{neig}$')
ax2.set_ylabel(r'$L_1\omega/(2V)$')
ax2.legend(loc='best',ncol=2)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim((0.004,0.4))
plt.show()



        
        

